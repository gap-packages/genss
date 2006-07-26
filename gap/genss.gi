# A little try to implement a generic Schreier-Sims using the orb package

# Configuration variables:

GENSS.InitialHashSize := NextPrimeInt(20000);
GENSS.Report := 10000;
GENSS.SchreierPatience := 200;
GENSS.SchreierPatienceBaseImage := infinity;
GENSS.VerificationElements := 10;  # this makes for error prob. 1/1024
GENSS.ShortOrbitsNrRandoms := 10;
GENSS.OrbitLengthLimit := 10000000;
GENSS.FindVectorsWithShortOrbitFunc := GENSS_FindVectorsWithShortOrbit;
GENSS.ShortOrbitsOrbLimit := 20000;
GENSS.ShortOrbitsStartLimit := 100;
GENSS.NumberPrevOrbitPoints := 30;
GENSS.NumberSchreierGensWait := 5;
GENSS.RandomStabGens := 5;

# A stabiliser chain is a component object with components:
# Its type is IsStabilizerChain and IsStabilizerChainByOrb
#
#  gens    - the group generators
#  gensi   - the inverses of the group generators
#  size    - if not equal to fail, this is the known size of the group
#  orb     - the orbit
#  cand    - a record for the candidates, consisting of
#            .list : a list of candidates for base points
#            .ops  : a corresponding list of operations
#            .used : number of already used candidates
#  stab    - information for the stabiliser, can be:
#               false                     if the stabiliser is trivial
#               a record like this here   if there is a nontrivial stabiliser
#  base    - a list shared by all records in the chain with the base images
#  depth   - the depth of the record in the chain, 1 is outermost
#  ntries  - number of unsuccessful tries to increase the stabilizer
#  opt     - an options record shared by all records in the chain, comps:
#            SchreierPatience : number of useless Schreier gens before
#                               giving up sifting Schreier gens
#            SchreierPatienceBaseImage : number of useless Schreier gens
#                               before giving up sifting base images
#
# The orbits have a Schreier vector but do *not* compute the stabiliser

InstallGlobalFunction( GENSS_NextBasePoint, 
  function( gens, cand, S )
    local NotFixedUnderAllGens,i;

    NotFixedUnderAllGens := function( gens, x, op )
      return ForAny( gens, g->op(x,g) <> x );
    end;

    # S can be false or a stabilizer chain record
    if S <> false then  # try points in previous orbit
        for i in [2..Minimum(Length(S!.orb),GENSS.NumberPrevOrbitPoints)] do
            if NotFixedUnderAllGens(gens,S!.orb[i],S!.orb!.op) then
                return rec( point := S!.orb[i], op := S!.orb!.op, 
                            cand := cand );
            fi;
        od;
    fi;

    repeat
        if cand.used >= Length(cand.points) then
            cand := FindBasePointCandidates(Group(gens),1);
        fi;
        cand.used := cand.used + 1;
    until NotFixedUnderAllGens(gens,cand.points[cand.used],cand.ops[cand.used]);

    return rec( point := cand.points[cand.used], op := cand.ops[cand.used],
                cand := cand );
  end );

InstallGlobalFunction( GENSS_CreateStabChainRecord,
  function( gens, size, depth, nextpoint, nextop, base, cand, opt )
    local S,hashsize;

    Info( InfoGenSS, 1, "Creating new stab chain record..." );

    gens := ShallowCopy(gens);
    # Note that we do ShallowCopy twice such that this list and the
    # one in the orbit record are different from each other and from the
    # list given as an argument.
    S := rec( gens := ShallowCopy(gens), gensi := List(gens,x->x^-1),
              stab := false, cand := cand, size := size, base := base,
              nrtries := 0, opt := opt, stabgens := [], stabwords := [],
              depth := depth );

    if IsInt(size) then
        hashsize := NextPrimeInt(Minimum(size,GENSS.InitialHashSize));
        S.orb := Orb( gens, nextpoint, nextop, hashsize, 
                      rec( schreier := true, storenumbers := true,
                           schreiergenaction := GENSS_CollectSchreierGens,
                           report := GENSS.Report, grpsizebound := size ) );
    else
        hashsize := GENSS.InitialHashSize;
        S.orb := Orb( gens, nextpoint, nextop, hashsize, 
                      rec( schreier := true, storenumbers := true,
                           schreiergenaction := GENSS_CollectSchreierGens,
                           report := GENSS.Report ) );
    fi;
    Add(base,nextpoint);
    S.orb!.stabilizerchain := S;
    Objectify( StabChainByOrbType, S );

    return S;
  end );

InstallGlobalFunction( GENSS_StartNewOrbitEnumeration,
  function( S )
    # Really starts the orbit enumeration for the layer below
    # thus filling S!.stab with a new layer:
    local next;
    Info( InfoGenSS, 1, "Creating new stabilizer layer..." );
    next := GENSS_NextBasePoint( S!.stabgens, S!.cand, S );
    S!.stab := GENSS_CreateStabChainRecord(S!.stabgens,false,
               S!.depth+1,next.point,next.op,S!.base,next.cand,S!.opt);
    Info( InfoGenSS, 1, "Entering orbit enumeration depth ",
          S!.stab!.depth, "..." );
    repeat
        Enumerate(S!.stab!.orb,GENSS.OrbitLengthLimit);
        if not(IsClosed(S!.stab!.orb)) then
            Error("Orbit too long, increase GENSS.OrbitLengthLimit!");
        fi;
    until IsClosed(S!.stab!.orb);
    Info( InfoGenSS, 1, "Done orbit enumeration depth ",S!.stab!.depth,
          "..." );
    S!.orb!.schreiergenaction := GENSS_ConsiderNewSchreierGenerator;
    S!.orb!.stabsize := Size(S!.stab);
  end );

InstallGlobalFunction( GENSS_CollectSchreierGens,
  function( o, i, j, pos )
    # Called for the first few Schreier generators
    # After a certain number of generators we switch to 
    # GENSS_CreateNewStabilizer
    local S,sgen,wordb,wordf;

    S := o!.stabilizerchain;

    # For randomisation we stop looking at Schreier generators eventually:
    if S!.nrtries > S!.opt.SchreierPatience then return false; fi;
    
    # First create the new stabilizer element:
    wordf := TraceSchreierTreeForward(o,i);
    wordb := TraceSchreierTreeForward(o,pos);
    sgen := EvaluateWord(o!.gens,wordf) * o!.gens[j] * 
            EvaluateWord(o!.gens,wordb)^-1;
    S!.nrtries := S!.nrtries + 1;
    if not(IsOne(sgen)) then  # found a stabilizer element, begin new orbit
        Add( S!.stabwords, Concatenation( wordf, [j], ORB_InvWord(wordb) ) );
        Add( S!.stabgens, sgen );
        if Length(S!.stabgens) >= GENSS.NumberSchreierGensWait then

            GENSS_StartNewOrbitEnumeration(S);
            return true;

        fi;
    fi;
    return false;
  end );

InstallGlobalFunction( GENSS_ConsiderNewSchreierGenerator,
  function( o, i, j, pos )
    # Called for all further Schreier generators
    # A new generator is created and sifted
    local S,new,sgen,wordb,wordf;

    S := o!.stabilizerchain;

    # For randomisation we stop looking at Schreier generators eventually:
    if S!.nrtries > S!.opt.SchreierPatience then return false; fi;
    
    # First create the new stabilizer element:
    wordf := TraceSchreierTreeForward(o,i);
    wordb := TraceSchreierTreeForward(o,pos);
    sgen := EvaluateWord(o!.gens,wordf) * o!.gens[j] * 
            EvaluateWord(o!.gens,wordb)^-1;
    S!.nrtries := S!.nrtries + 1;
    if not(IsOne(sgen)) then   # found a stabilizer element, add generator:
        Info( InfoGenSS, 2, "Considering new Schreier generator in depth ",
              S!.depth,"..." );
        new := AddGeneratorToStabilizerChain(S!.stab,sgen,false);
        if new then 
            Add(S!.stabgens,sgen);
            Add(S!.stabwords,Concatenation(wordf,[j],ORB_InvWord(wordb)));
            S!.nrtries := 0; 
            o!.stabsize := Size(S!.stab);
        fi;
        return new;
    fi;
    return false;  # nothing new
  end );

InstallMethod( AddGeneratorToStabilizerChain, 
  "for a stabchain, a group element, and a boolean value",
  [ IsStabilizerChain and IsStabilizerChainByOrb, IsObject, IsBool ],
  function( S, el, complete )
    # Increases the set represented by S by the generator el.
    local r;
    r := SiftGroupElement( S, el );
    # if this is one then sgen is already contained in the stabilizer
    if r.isone then     # already in the group!
        return false;
    fi;
    # Now there remain two cases:
    #  (1) the sift stopped somewhere and we have to add a generator there
    #  (2) the sift ran all through the chain and the element still was not
    #      the identity, then we have to prolong the chain
    if r.S <> false then   # case (1)
        S := r.S;
        Info( InfoGenSS, 1, "Adding new generator to stabilizer chain..." );
        Add(S!.gens,r.rem);
        Add(S!.gensi,r.rem^-1);
        AddGeneratorsToOrbit(S!.orb,[r.rem]);
        S!.nrtries := 0;   # again look at a few Schreier generators!
        Info( InfoGenSS, 1, "Entering orbit enumeration depth ",S!.depth,
              "..." );
        repeat
            Enumerate(S!.orb,GENSS.OrbitLengthLimit);
            if not(IsClosed(S!.orb)) then
                Error("Orbit too long, increase GENSS.OrbitLengthLimit!");
            fi;
        until IsClosed(S!.orb);
        Info( InfoGenSS, 1, "Done orbit enumeration depth ",S!.depth,
              "..." );
    else   # case (2)
        S := r.preS;   # this is the last stabilizer in the chain
        Add( S!.stabgens, r.rem );
        Add( S!.stabwords, "new" );
        if Length(S!.stabgens) >= GENSS.NumberSchreierGensWait then
            GENSS_StartNewOrbitEnumeration(S);
        fi;
    fi;
    if complete then
        GENSS_CompleteStabilizerChain(S);
    fi;
    return true;
  end );

InstallMethod( SiftGroupElement, "for a stabilizer chain and a group element",
  [ IsStabilizerChain and IsStabilizerChainByOrb, IsObject ],
  function( S, x )
    local o,p,po,preS;
    #Info( InfoGenSS, 1, "Sifting group element..." );
    preS := false;
    while S <> false do
        o := S!.orb;
        p := o!.op(o[1],x);
        po := Position(o,p);
        if po = fail then   # not in current stabilizer
            return rec( isone := false, rem := x, S := S, preS := preS );
        fi;
        # Now sift through Schreier tree:
        while po > 1 do
            x := x * S!.gensi[o!.schreiergen[po]];
            po := o!.schreierpos[po];
        od;
        preS := S;
        S := S!.stab;
    od;
    return rec( isone := IsOne(x), rem := x, S := false, preS := preS );
  end );

InstallMethod( SiftGroupElement2, "for a stabilizer chain and a group element",
  [ IsStabilizerChain and IsStabilizerChainByOrb, IsObject ],
  function( S, x )
    local o,p,po,preS;
    #Info( InfoGenSS, 1, "Sifting group element..." );
    preS := false;
    while S <> false do
        o := S!.orb;
        p := o!.op(o[1],x);
        po := Position(o,p);
        if po = fail then   # not in current stabilizer
            return rec( isone := false, rem := x, S := S, preS := preS );
        fi;
        # Now sift through Schreier tree:
        while po > 1 do
            x := x * S!.orb!.gensi[o!.schreiergen[po]];
            po := o!.schreierpos[po];
        od;
        preS := S;
        S := S!.stab;
    od;
    return rec( isone := IsOne(x), rem := x, S := false, preS := preS );
  end );

InstallGlobalFunction( GENSS_CompleteStabilizerChain,
  function( S )
    while S!.stab <> false do S := S!.stab; od;
    if Length(S!.stabgens) > 0 then
        GENSS_StartNewOrbitEnumeration(S);
    fi;
  end );

InstallMethod( Size, "for a stabilizer chain",
  [ IsStabilizerChain and IsStabilizerChainByOrb ],
  function( S )
    local size;
    size := 1;
    while S <> false do
        size := size * Length(S!.orb);
        S := S!.stab;
    od;
    return size;
  end );

InstallMethod( FindBasePointCandidates, "for a matrix group over a FF",
  [ IsGroup and IsMatrixGroup, IsInt ],
  function( grp, i )
    local bv,cand,d,F,v,w;
    F := DefaultFieldOfMatrixGroup(grp);
    if not(IsFinite(F)) then
        TryNextMethod();
    fi;
    d := DimensionOfMatrixGroup(grp);
    Info( InfoGenSS, 2, "Finding nice base points..." );
    if IsObjWithMemory(GeneratorsOfGroup(grp)[1]) then
        grp := Group(StripMemory(GeneratorsOfGroup(grp)));
    fi;
    if i = 0 and Size(F)^d > 100000 then
        bv := GENSS.FindVectorsWithShortOrbitFunc(grp);
    else
        bv := One(grp);
    fi;
    cand := rec( points := [], ops := [], used := 0 );
    for v in bv do
        w := ORB_NormalizeVector(ShallowCopy(v));
        if Size(F) <> 2 then
            Add(cand.points,w);
            Add(cand.ops,OnLines);
            Add(cand.points,v);
            Add(cand.ops,OnRight);
        else
            Add(cand.points,w);
            Add(cand.ops,OnRight);
        fi;
    od;
    return cand;
  end );

InstallMethod( FindBasePointCandidates, "for a permutation group",
  [ IsGroup and IsPermGroup, IsInt ],
  function( grp, i )
    local ops,points;
    if i = 0 then
        points := [1..Minimum(20,LargestMovedPoint(grp))];
    else
        points := [1..LargestMovedPoint(grp)];
    fi;
    ops := List([1..Length(points)],x->OnPoints);
    return rec( points := points, ops := ops, used := 0 );
  end );
    
InstallMethod( StabilizerChain, "for a group object", [ IsGroup, IsRecord ],
  function( grp, opt )
    # Computes a stabilizer chain for the group grp
    local S,cand,i,next,res,x;
    # Add some default options:
    if not(IsBound(opt.SchreierPatience)) then
        opt.SchreierPatience := GENSS.SchreierPatience;
    fi;
    if not(IsBound(opt.VerificationElements)) then
        opt.VerificationElements := GENSS.VerificationElements;
    fi;
    if not(IsBound(opt.SchreierPatienceBaseImage)) then
        opt.SchreierPatienceBaseImage := GENSS.SchreierPatienceBaseImage;
    fi;
    if IsBound(opt.cand) then
        cand := opt.cand;
    else
        cand := FindBasePointCandidates( grp, 0 );
        if IsMatrixGroup(grp) and IsBound(opt.ShortOrbit) then
            i := opt.ShortOrbit;
            repeat
                i := i - 1;
                res := GENSS_FindShortOrbit(grp);
            until res <> fail or i = 0;
            if res <> fail then
                if Size(DefaultFieldOfMatrixGroup(grp)) > 2 then
                    Add(cand.points,res[1],1);
                    Add(cand.ops,OnLines,1);
                    Add(cand.points,res[1],2);
                    Add(cand.ops,OnPoints,2);
                else
                    Add(cand.points,res[1],1);
                    Add(cand.ops,OnPoints,1);
                fi;
            fi;
        fi;
    fi;
    next := GENSS_NextBasePoint(GeneratorsOfGroup(grp),cand,false);
    if HasSize(grp) then
        S := GENSS_CreateStabChainRecord(GeneratorsOfGroup(grp),Size(grp),1,
                 next.point,next.op,[],cand,opt);
    else
        S := GENSS_CreateStabChainRecord(GeneratorsOfGroup(grp),false,1,
                 next.point,next.op,[],cand,opt);
    fi;
    Info( InfoGenSS, 1, "Entering orbit enumeration depth 1..." );
    repeat
        Enumerate(S!.orb,GENSS.OrbitLengthLimit);
        if not(IsClosed(S!.orb)) then
            Error("Orbit too long, increase GENSS.OrbitLengthLimit!");
        fi;
    until IsClosed(S!.orb);
    Info( InfoGenSS, 1, "Done orbit enumeration depth 1..." );
    GENSS_CompleteStabilizerChain(S);
    if IsInt(opt.SchreierPatience) then
        if HasSize(grp) then
            while Size(S) < Size(grp) do
                x := PseudoRandom(grp);
                if AddGeneratorToStabilizerChain(S,x,true) then
                    Info( InfoGenSS, 1, "Verification found error ... ",
                          "stabilizer chain was improved." );
                fi;
            od;
        else
            # Do some verification here:
            i := 0; 
            while i < opt.VerificationElements do
                i := i + 1;
                #Info(InfoGenSS,1,"i=",i," size=",Size(S));
                x := PseudoRandom(grp);
                if AddGeneratorToStabilizerChain(S,x,true) then
                    Info( InfoGenSS, 1, "Verification found error ... ",
                          "stabilizer chain was improved." );
                    i := 0;
                fi;
            od;
        fi;
    fi;
    return S;
  end );

InstallMethod( StabilizerChain, "for a group object", [ IsGroup ],
  function( grp )
    return StabilizerChain( grp, rec() );
  end );

InstallMethod( ViewObj, "for a stabilizer chain",
  [ IsStabilizerChain and IsStabilizerChainByOrb ],
  function( S )
    local i;
    if not(IsBound(GENSS.VIEWDEPTH)) then
        GENSS.VIEWDEPTH := 0;
    else
        GENSS.VIEWDEPTH := GENSS.VIEWDEPTH + 1;
    fi;
    for i in [1..GENSS.VIEWDEPTH] do Print(" "); od;
    Print("<stabilizer chain size=",Size(S)," orblen=",Length(S!.orb),
          " layer=",S!.layer," Schreier depth=",S!.orb!.depth,">");
    if S!.stab <> false then
        Print("\n");
        ViewObj(S!.stab);
    fi;
    GENSS.VIEWDEPTH := GENSS.VIEWDEPTH - 1;
    if GENSS.VIEWDEPTH < 0 then
        Unbind(GENSS.VIEWDEPTH);
    fi;
  end );

InstallGlobalFunction( GENSS_MapBaseImage,
  function( bi, el, S )
    # bi must be a base image belonging to S, el a group element
    local i,l,res;
    l := Length(bi);
    res := 0*[1..l];
    i := 1;
    while true do
        res[i] := S!.orb!.op(bi[i],el);
        i := i + 1;
        if i = l+1 then break; fi;
        S := S!.stab;
    od;
    return res;
  end );

InstallGlobalFunction( GENSS_FindVectorsWithShortOrbit,
  function(g)
    local c,f,i,inters,j,l,nw,sortfun,v,vv,w,wb,ww;
    l := ShallowCopy(GeneratorsOfGroup(g));
    f := DefaultFieldOfMatrixGroup(g);
    for i in [1..GENSS.ShortOrbitsNrRandoms] do
        Add(l,PseudoRandom(g));
    od;
    c := List(l,x->Set(Factors(CharacteristicPolynomial(x,1))));
    v := [];
    for i in [1..Length(l)] do
        for j in [1..Length(c[i])] do
            vv := [];
            Add(vv,[VectorSpace(f,NullspaceMat(Value(c[i][j],l[i]))),
                    Degree(c[i][j]),
                    WeightVecFFE(CoefficientsOfLaurentPolynomial(c[i][j])[1]),
                    1]);
        od;
        Add(v,vv);
    od;
    Info(InfoGenSS,2,"Have eigenspaces.");
    # Now collect a list of all those spaces together with all
    # possible intersects:
    w := [];
    for i in [1..Length(l)] do
        nw := [];
        for j in [1..Length(v[i])] do
            for ww in w do
                inters := Intersection(ww[1],v[i][j][1]);
                if Dimension(inters) > 0 then
                    Add(nw,[inters,Minimum(ww[2],v[i][j][2]),
                            Minimum(ww[3],v[i][j][3]),ww[4]+v[i][j][4]]);
                fi;
            od;
            Add(nw,v[i][j]);
        od;
        Append(w,nw);
    od;
    sortfun := function(a,b)
        if a[2] < b[2] then return true;
        elif a[2] > b[2] then return false;
        elif a[3] < b[3] then return true;
        elif a[3] > b[3] then return false;
        elif a[4] < b[4] then return true;
        elif a[4] > b[4] then return false;
        elif Dimension(a[1]) < Dimension(b[1]) then return true;
        else return false;
        fi;
    end;
    Sort(w,sortfun);
    wb := List(w,ww->Basis(ww[1])[1]);
    Info(InfoGenSS,2,"Have ",Length(wb)," vectors for possibly short orbits.");
    for ww in [1..Length(wb)] do
        if not(IsMutable(wb[ww])) then
            wb[ww] := ShallowCopy(wb[ww]);
        fi;
        ORB_NormalizeVector(wb[ww]);
    od;
    return wb;
end );

InstallGlobalFunction( GENSS_FindShortOrbit,
  function( g )
    local ThrowAwayOrbit,found,gens,hashlen,i,j,limit,newnrorbs,nrorbs,o,wb;

    wb := GENSS_FindVectorsWithShortOrbit(g);

    # Now we have a list of vectors with (hopefully) short orbits.
    # We start enumerating all those orbits, but first only 50 elements:
    nrorbs := Minimum(Length(wb),64);  # take only the 64 first
    gens := GeneratorsOfGroup(g);
    o := [];
    hashlen := NextPrimeInt(QuoInt(GENSS.ShortOrbitsOrbLimit * 3,2));
    for i in [1..nrorbs] do
        Add(o,Orb(gens,ShallowCopy(wb[i]),OnLines,hashlen));
    od;
    limit := GENSS.ShortOrbitsStartLimit;
    i := 1;               # we start to work on the first one

    ThrowAwayOrbit := function(i)
        # This removes orbit number i from o, thereby handling nrorbs and
        # Length(o) correctly. If you want to use o[i] further, please
        # make a copy (of the reference) before calling this function.
        if Length(o) > nrorbs then
            o[i] := o[nrorbs+1];
            o{[nrorbs+1..Length(o)-1]} := o{[nrorbs+2..Length(o)]};
            Unbind(o[Length(o)]);
        else
            o{[i..nrorbs-1]} := o{[i+1..nrorbs]};
            Unbind(o[nrorbs]);
            nrorbs := nrorbs-1;
        fi;
    end;

    repeat
        Enumerate(o[i],limit);
        found := IsClosed(o[i]);
        if Length(o[i]) = 1 then
            Info(InfoGenSS,2,"Orbit Number ",i," has length 1.");
            found := false;
            # Now throw away this orbit:
            ThrowAwayOrbit(i);
            # we intentionally do not increase i here!
        elif not(found) then
            i := i + 1;
        fi;
        if i > nrorbs then
          Info(InfoGenSS,2,"Done ",nrorbs,
               " orbit(s) to limit ",limit,".");
          limit := limit * 2;
          if limit > GENSS.ShortOrbitsOrbLimit then
              Info(InfoGenSS,2,"Limit reached, giving up.");
              return fail;
          fi;
          i := 1;
          if nrorbs < i then
              Info(InfoGenSS,2,"No orbits left, giving up.");
              return fail;
          fi;
          if nrorbs > 1 then
              newnrorbs := QuoInt((nrorbs+1),2);
              for j in [newnrorbs+1..nrorbs] do
                  Unbind(o[j]);
              od;
              nrorbs := newnrorbs;
          fi;
        fi;
    until found;
    Info(InfoGenSS,2,"Found orbit of length ",Length(o[i])," (#",i,").");
    return o[i];
  end );   

# Second try:

# Some notes:
# - get rid of nrtries entry?
# - do we need stabgens and stabwords
# - call depth "layer" ?

InstallGlobalFunction( GENSS_CreateStabChainRecord2,
  function( gens, size, layer, nextpoint, nextop, base, cand, opt )
    local S,hashsize;

    Info( InfoGenSS, 2, "Creating new stab chain record..." );

    gens := ShallowCopy(gens);
    # Note that we do ShallowCopy twice such that this list and the
    # one in the orbit record are different from each other and from the
    # list given as an argument.
    S := rec( stab := false, cand := cand, size := size, base := base,
              nrtries := 0, opt := opt, stabgens := [], stabwords := [],
              layer := layer );

    if IsInt(size) then
        hashsize := NextPrimeInt(Minimum(size,GENSS.InitialHashSize));
    else
        hashsize := GENSS.InitialHashSize;
    fi;
    S.orb := Orb( gens, nextpoint, nextop,
                  rec( hashlen := hashsize, schreier := true, log := true,
                       report := GENSS.Report ) );
    Add(base,nextpoint);
    S.orb!.stabilizerchain := S;
    Objectify( StabChainByOrbType, S );

    return S;
  end );

InstallMethod( StabilizerChain2, "for a group object", [ IsGroup ],
  function( grp )
    return StabilizerChain2( grp, rec() );
  end );

InstallGlobalFunction( GENSS_StabilizerChain2Inner,
  function( gens, size, layer, cand, opt, parentS )
    # Computes a stabilizer chain for the group generated by gens
    # with known size size (can be false if not known). This will be
    # layer layer in the final stabilizer chain. cand is a (shared)
    # record for base point candidates and opt the (shared) option
    # record. This is called in StabilizerChain2 and calls itself.
    # It also can be called if a new layer is needed.
    local S,i,next,pr,r,stabgens,x;
    Info(InfoGenSS,2,"Entering GENSS_StabilizerChain2Inner layer=",layer);
    next := GENSS_NextBasePoint(gens,cand,parentS);
    if size <> false then
        S := GENSS_CreateStabChainRecord2(gens,size,layer,
                 next.point,next.op,[],cand,opt);
    else
        S := GENSS_CreateStabChainRecord2(gens,false,layer,
                 next.point,next.op,[],cand,opt);
    fi;
    Info( InfoGenSS, 2, "Entering orbit enumeration layer ",layer,"..." );
    repeat
        Enumerate(S!.orb,GENSS.OrbitLengthLimit);
        if not(IsClosed(S!.orb)) then
            Error("Orbit too long, increase GENSS.OrbitLengthLimit!");
        fi;
    until IsClosed(S!.orb);
    Info(InfoGenSS, 1, "Layer ", layer, ": Orbit length is ", Length(S!.orb), 
         ", creating shallow Schreier tree...");
    MakeSchreierTreeShallow(S!.orb);
    S!.orb!.gensi := List(S!.orb!.gens,x->x^-1);
    Info( InfoGenSS, 2, "Done orbit enumeration layer ",layer );

    # Are we done?
    if size <> false and Length(S!.orb) = size then
        S!.proof := true;
        Info(InfoGenSS,2,"Leaving GENSS_StabilizerChain2Inner layer=",layer);
        return S;
    fi;

    # Now create a few random stabilizer elements:
    Info(InfoGenSS,2,"Creating ",GENSS.RandomStabGens,
         " random elements of the point stabilizer...");
    pr := ProductReplacer( gens );
    stabgens := [];
    for i in [1..GENSS.RandomStabGens] do
        x := Next(pr);
        r := SiftGroupElement2( S, x );   # this goes down only one step!
        if not(r.isone) then
            # Now r.rem is the remainder, a stabilizer element
            Add(stabgens,r.rem);
        fi;
    od;
    if Length(stabgens) > 0 then   # there is a non-trivial stabiliser
        Info(InfoGenSS,2,"Found ",Length(stabgens)," non-trivial ones.");
        if size <> false then
            S!.stab := GENSS_StabilizerChain2Inner(stabgens,size/Length(S!.orb),
                                                   layer+1,cand,opt,S);
        else
            S!.stab := GENSS_StabilizerChain2Inner(stabgens,false,
                                                   layer+1,cand,opt,S);
        fi;
        S!.proof := S!.stab!.proof;   # hand up information
    else
        # We are not sure that the next stabiliser is trivial, but we believe!
        Info(InfoGenSS,2,"Found no non-trivial ones.");
        S!.proof := false;
    fi;

    Info(InfoGenSS,2,"Leaving GENSS_StabilizerChain2Inner layer=",layer);
    return S;
  end );

InstallMethod( AddGeneratorToStabilizerChain2,
  "for a stabilizer chain and a new generator",
  [ IsStabilizerChain and IsStabilizerChainByOrb, IsObject ],
  function( S, el )
    # Increases the set represented by S by the generator el.
    local r,subsize;
    r := SiftGroupElement2( S, el );
    # if this is one then sgen is already contained in the stabilizer
    if r.isone then     # already in the group!
        return false;
    fi;
    # Now there remain two cases:
    #  (1) the sift stopped somewhere and we have to add a generator there
    #  (2) the sift ran all through the chain and the element still was not
    #      the identity, then we have to prolong the chain
    if r.S <> false then   # case (1)
        S := r.S;
        Info( InfoGenSS, 1, "Adding new generator to stabilizer chain..." );
        AddGeneratorsToOrbit(S!.orb,[r.rem]);
        Add(S!.orb!.gensi,r.rem^-1);
        Info( InfoGenSS, 2, "Entering orbit enumeration layer ",S!.layer,
              "..." );
        repeat
            Enumerate(S!.orb,GENSS.OrbitLengthLimit);
            if not(IsClosed(S!.orb)) then
                Error("Orbit too long, increase GENSS.OrbitLengthLimit!");
            fi;
        until IsClosed(S!.orb);
        MakeSchreierTreeShallow(S!.orb);
        Info( InfoGenSS, 2, "Done orbit enumeration layer ",S!.layer,
              "..." );
        S!.proof := false;
    else   # case (2)
        S := r.preS;   # this is the last stabilizer in the chain
        subsize := Order(r.rem);
        S!.stab := GENSS_StabilizerChain2Inner( [r.rem], subsize, S!.layer+1,
                                                S!.cand, S!.opt, S );
        S!.proof := false;
    fi;
    return true;   # something new
  end );

InstallMethod( StabilizerChain2, "for a group object", [ IsGroup, IsRecord ],
  function( grp, opt )
    # Computes a stabilizer chain for the group grp
    local S,cand,i,pr,res,x;

    # First a few preparations, then we delegate to GENSS_StabilizerChain2Inner:

    # Add some default options:
    if not(IsBound(opt.VerificationElements)) then
        opt.VerificationElements := GENSS.VerificationElements;
    fi;
    # Find base point candidates:
    if IsBound(opt.cand) then
        cand := opt.cand;
    else
        cand := FindBasePointCandidates( grp, 0 );
        if IsMatrixGroup(grp) and IsBound(opt.ShortOrbit) then
            i := opt.ShortOrbit;
            repeat
                i := i - 1;
                res := GENSS_FindShortOrbit(grp);
            until res <> fail or i = 0;
            if res <> fail then
                if Size(DefaultFieldOfMatrixGroup(grp)) > 2 then
                    Add(cand.points,res[1],1);
                    Add(cand.ops,OnLines,1);
                    Add(cand.points,res[1],2);
                    Add(cand.ops,OnPoints,2);
                else
                    Add(cand.points,res[1],1);
                    Add(cand.ops,OnPoints,1);
                fi;
            fi;
        fi;
    fi;

    if HasSize(grp) then
        S := GENSS_StabilizerChain2Inner(GeneratorsOfGroup(grp), Size(grp), 1,
                                         cand, opt, false);
    else
        S := GENSS_StabilizerChain2Inner(GeneratorsOfGroup(grp), false, 1, 
                                         cand, opt, false);
    fi;

    # Do we already have a proof?
    if S!.proof then return S; fi;

    # Now a possible verification phase:
    if S!.size <> false then
        pr := ProductReplacer(GeneratorsOfGroup(grp));
        while Size(S) < Size(grp) do
            Info(InfoGenSS,1,"Known size not reached, throwing in a random ",
                 "element...");
            x := Next(pr);
            if AddGeneratorToStabilizerChain2(S,x) then
                Info( InfoGenSS, 1, "Increased size to ",Size(S) );
            fi;
        od;
    else
        # Do some verification here:
        Info(InfoGenSS,1,"Doing randomized verification...");
        pr := ProductReplacer(GeneratorsOfGroup(grp));
        i := 0; 
        while i < opt.VerificationElements do
            i := i + 1;
            #Info(InfoGenSS,1,"i=",i," size=",Size(S));
            x := Next(pr);
            if AddGeneratorToStabilizerChain2(S,x) then
                Info( InfoGenSS, 1, "Verification found error ... ",
                      "new size ", Size(S) );
                i := 0;
            fi;
        od;
    fi;
    return S;
  end );

