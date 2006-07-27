# A little try to implement a generic Schreier-Sims using the orb package

# Configuration variables, these are taken as defaults for the options:

GENSS.InitialHashSize := NextPrimeInt(20000);
GENSS.Report := 10000;
GENSS.ShortOrbitsNrRandoms := 10;
GENSS.OrbitLengthLimit := 10000000;
GENSS.ShortOrbitsOrbLimit := 20000;
GENSS.ShortOrbitsStartLimit := 100;
GENSS.NumberPrevOrbitPoints := 30;
GENSS.RandomStabGens := 5;
GENSS.StabGenScramble := 30;
GENSS.StabGenScrambleFactor := 6;
GENSS.StabGenAddSlots := 3;
GENSS.StabGenMaxDepth := 200;
GENSS.VerifyScramble := 100;
GENSS.VerifyScrambleFactor := 10;
GENSS.VerifyAddSlots := 10;
GENSS.VerifyMaxDepth := 400;
GENSS.VerifyElements := 10;   # this amounts to an error bound of 1/1024
GENSS.DeterministicVerification := false;

InstallGlobalFunction( GENSS_CopyDefaultOptions,
  function( defopt, opt )
    local n;
    for n in RecFields(defopt) do
        if not(IsBound(opt.(n))) then
            opt.(n) := defopt.(n);
        fi;
    od;
  end );


# A stabiliser chain is a component object with components:
# Its type is IsStabilizerChain and IsStabilizerChainByOrb
#
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
#  layer   - the layer of the record in the chain, 1 is outermost
#  opt     - an options record shared by all records in the chain, comps:
#            defaults are in the global record GENSS
#
# The orbits have a Schreier vector but do *not* compute the stabiliser
# they write a log.

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
            x := x * S!.orb!.gensi[o!.schreiergen[po]];
            po := o!.schreierpos[po];
        od;
        preS := S;
        S := S!.stab;
    od;
    return rec( isone := IsOne(x), rem := x, S := false, preS := preS );
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
        bv := GENSS_FindVectorsWithShortOrbit(grp);
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
    Print("<stabchain size=",Size(S)," orblen=",Length(S!.orb),
          " layer=",S!.layer," SchreierDepth=",S!.orb!.depth,">");
    if S!.stab <> false then
        Print("\n");
        ViewObj(S!.stab);
    fi;
    GENSS.VIEWDEPTH := GENSS.VIEWDEPTH - 1;
    if GENSS.VIEWDEPTH < 0 then
        Unbind(GENSS.VIEWDEPTH);
    fi;
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


InstallGlobalFunction( GENSS_CreateStabChainRecord,
  function( gens, size, layer, nextpoint, nextop, base, cand, opt )
    local S,hashsize;

    Info( InfoGenSS, 2, "Creating new stab chain record..." );

    gens := ShallowCopy(gens);
    # Note that we do ShallowCopy twice such that this list and the
    # one in the orbit record are different from each other and from the
    # list given as an argument.
    S := rec( stab := false, cand := cand, size := size, base := base,
              nrtries := 0, opt := opt, layer := layer );

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

InstallMethod( StabilizerChain, "for a group object", [ IsGroup ],
  function( grp )
    return StabilizerChain( grp, rec() );
  end );

InstallGlobalFunction( GENSS_StabilizerChainInner,
  function( gens, size, layer, cand, opt, parentS )
    # Computes a stabilizer chain for the group generated by gens
    # with known size size (can be false if not known). This will be
    # layer layer in the final stabilizer chain. cand is a (shared)
    # record for base point candidates and opt the (shared) option
    # record. This is called in StabilizerChain and calls itself.
    # It also can be called if a new layer is needed.
    local gen,S,i,next,pr,r,stabgens,x;
    Info(InfoGenSS,2,"Entering GENSS_StabilizerChainInner layer=",layer);
    next := GENSS_NextBasePoint(gens,cand,parentS);
    if size <> false then
        S := GENSS_CreateStabChainRecord(gens,size,layer,
                 next.point,next.op,[],cand,opt);
    else
        S := GENSS_CreateStabChainRecord(gens,false,layer,
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
        Info(InfoGenSS,2,"Leaving GENSS_StabilizerChainInner layer=",layer);
        return S;
    fi;

    # Now create a few random stabilizer elements:
    Info(InfoGenSS,2,"Creating ",GENSS.RandomStabGens,
         " random elements of the point stabilizer...");
    pr := ProductReplacer( gens,
                      rec( scramble := S!.opt.StabGenScramble,
                           scramblefactor := S!.opt.StabGenScrambleFactor,
                           addslots := S!.opt.StabGenAddSlots,
                           maxdepth := S!.opt.StabGenMaxDepth ));
    S!.pr := pr;   # for later use
    stabgens := [];
    for i in [1..GENSS.RandomStabGens] do
        x := Next(pr);
        r := SiftGroupElement( S, x );   # this goes down only one step!
        if not(r.isone) then
            # Now r.rem is the remainder, a stabilizer element
            Add(stabgens,r.rem);
        fi;
    od;
    if Length(stabgens) > 0 then   # there is a non-trivial stabiliser
        Info(InfoGenSS,2,"Found ",Length(stabgens)," non-trivial ones.");
        if size <> false then
            S!.stab := GENSS_StabilizerChainInner(stabgens,size/Length(S!.orb),
                                                   layer+1,cand,opt,S);
        else
            S!.stab := GENSS_StabilizerChainInner(stabgens,false,
                                                   layer+1,cand,opt,S);
        fi;
        S!.proof := S!.stab!.proof;   # hand up information
    else
        # We are not sure that the next stabiliser is trivial, but we believe!
        Info(InfoGenSS,2,"Found no non-trivial ones.");
        S!.proof := false;
    fi;

    # Now we sift our generators through to ensure that no generator
    # fixes the whole basis:
    #Info(InfoGenSS,2,"Sifting the generators in layer ",S!.layer);
    #for gen in S!.orb!.gens do
    #    AddGeneratorToStabilizerChain(S,gen);
    #od;

    Info(InfoGenSS,2,"Leaving GENSS_StabilizerChainInner layer=",layer);
    return S;
  end );

InstallMethod( AddGeneratorToStabilizerChain,
  "for a stabilizer chain and a new generator",
  [ IsStabilizerChain and IsStabilizerChainByOrb, IsObject ],
  function( S, el )
    # Increases the set represented by S by the generator el.
    local i,r,SS,subsize;
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
        SS := r.S;
        Info( InfoGenSS, 1, "Adding new generator to stabilizer chain ",
              "in layer ", SS!.layer, "..." );
        AddGeneratorsToOrbit(SS!.orb,[r.rem]);
        Add(SS!.orb!.gensi,r.rem^-1);
        Info( InfoGenSS, 2, "Entering orbit enumeration layer ",SS!.layer,
              "..." );
        repeat
            Enumerate(SS!.orb,GENSS.OrbitLengthLimit);
            if not(IsClosed(SS!.orb)) then
                Error("Orbit too long, increase GENSS.OrbitLengthLimit!");
            fi;
        until IsClosed(SS!.orb);
        #MakeSchreierTreeShallow(S!.orb);
        #for i in [Length(S!.orb!.gensi)+1..Length(S!.orb!.gens)] do
        #    S!.orb!.gensi[i] := S!.orb!.gens[i]^-1;
        #od;
        Info( InfoGenSS, 2, "Done orbit enumeration layer ",SS!.layer,
              "..." );
        SS!.proof := false;
    else   # case (2)
        # Note that we do not create a pr instance here for one
        # generator, this will be done later on as needed...
        SS := r.preS;
        subsize := Order(r.rem);
        SS!.stab := GENSS_StabilizerChainInner([r.rem],subsize,
                           SS!.layer+1,SS!.cand, SS!.opt, SS );
        SS!.proof := false;
        SS := SS!.stab;
    fi;
    # Now we have added a new generator (or a new layer) at layer SS,
    # the new gen came from layer S (we were called here, after all),
    # thus we have to check, whether all the orbits between S (inclusively)
    # and SS (exclusively) are also closed under the new generator r.rem,
    # we add it to all these orbits, thereby also making the Schreier trees
    # shallower:
    while S!.layer <> SS!.layer do
        Info(InfoGenSS,2,"Adding new generator to orbit at layer ",S!.layer);
        AddGeneratorsToOrbit(S!.orb,[r.rem]);
        Add(S!.orb!.gensi,r.rem^-1);
        S := S!.stab;
    od;
    return true;
  end );

InstallMethod( StabilizerChain, "for a group object", [ IsGroup, IsRecord ],
  function( grp, opt )
    # Computes a stabilizer chain for the group grp
    local S,cand,i,pr,prob,res,x;

    # First a few preparations, then we delegate to GENSS_StabilizerChain2Inner:

    # Add some default options:
    GENSS_CopyDefaultOptions(GENSS,opt);

    # Old style error probability for compatibility:
    if IsBound(opt.random) then
        if opt.random = 0 then
            opt.VerifyElements := 0;
        elif opt.random = 1000 then
            opt.DeterministicVerification := true;
        else
            prob := 1/2;
            opt.VerifyElements := 1;
            while prob * (1000-opt.random) >= 1 do
                prob := prob / 2;
                opt.VerifyElements := opt.VerifyElements + 1;
            od;
        fi;
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
        S := GENSS_StabilizerChainInner(GeneratorsOfGroup(grp), Size(grp), 1,
                                         cand, opt, false);
    else
        S := GENSS_StabilizerChainInner(GeneratorsOfGroup(grp), false, 1, 
                                         cand, opt, false);
    fi;

    # Do we already have a proof?
    if S!.proof then return S; fi;

    Info(InfoGenSS,1,"Current size found: ",Size(S));
    # Now a possible verification phase:
    if S!.size <> false then   # we knew the size in advance
        Info(InfoGenSS,1,"Doing verification via known size...");
        pr := ProductReplacer(GeneratorsOfGroup(grp),
                      rec( scramble := opt.StabGenScramble,
                           scramblefactor := opt.StabGenScrambleFactor,
                           addslots := opt.StabGenAddSlots,
                           maxdepth := opt.StabGenMaxDepth ));
        while Size(S) < Size(grp) do
            Info(InfoGenSS,1,"Known size not reached, throwing in a random ",
                 "element...");
            x := Next(pr);
            if AddGeneratorToStabilizerChain(S,x) then
                Info( InfoGenSS, 1, "Increased size to ",Size(S) );
            fi;
        od;
    else
        # Do some verification here:
        Info(InfoGenSS,1,"Doing randomized verification...");
        pr := ProductReplacer(GeneratorsOfGroup(grp),
                      rec( scramble := opt.VerifyScramble,
                           scramblefactor := opt.VerifyScrambleFactor,
                           addslots := opt.VerifyAddSlots,
                           maxdepth := opt.VerifyMaxDepth ));
        i := 0; 
        while i < opt.VerifyElements do
            i := i + 1;
            x := Next(pr);
            if AddGeneratorToStabilizerChain(S,x) then
                Info( InfoGenSS, 1, "Verification found error ... ",
                      "new size ", Size(S) );
                i := 0;
            fi;
        od;
    fi;
    return S;
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

