#############################################################################
##
##  genss.gi              genss package           
##                                                           Max Neunhoeffer
##                                                              Felix Noeske
##
##  Copyright 2006 Lehrstuhl D für Mathematik, RWTH Aachen
##
##  Implementation stuff for generic Schreier-Sims
##
#############################################################################

#############################################################################
# The following global record contains default values for options for the 
# main function "StabilizerChain":
#############################################################################

# Initial hash size for orbits:
GENSS.InitialHashSize := NextPrimeInt(20000);
# Number of points to process before reporting:
GENSS.Report := 30000;
# Number of random elements to consider for the determination of short orbits:
GENSS.ShortOrbitsNrRandoms := 10;
# Absolute limit for orbit length in search for short orbits:
GENSS.ShortOrbitsOrbLimit := 80000;
# First limit for the parallel enumeration of orbits looking for short orbits:
GENSS.ShortOrbitsStartLimit := 400;
# Absolute limit for single orbit length:
GENSS.OrbitLengthLimit := 10000000;
# Number of points in the previous orbit to consider for the next base point:
GENSS.NumberPrevOrbitPoints := 30;
# Number of (evenly distributed) random generators for the stabilizer:
GENSS.RandomStabGens := 4;
# Product replacement parameters for the creation of random elements
# for short orbit search:
GENSS.ShortOrbScramble := 0;
GENSS.ShortOrbScrambleFactor := 0;
GENSS.ShortOrbAddSlots := 0;
GENSS.ShortOrbMaxDepth := 100;
# Product replacement parameters for the stabilizer element generation:
GENSS.StabGenScramble := 30;
GENSS.StabGenScrambleFactor := 6;
GENSS.StabGenAddSlots := 3;
GENSS.StabGenMaxDepth := 200;
# Product replacement parameters for the verification phase:
GENSS.VerifyScramble := 100;
GENSS.VerifyScrambleFactor := 10;
GENSS.VerifyAddSlots := 10;
GENSS.VerifyMaxDepth := 400;
# Number of random elements used for verification,
# note that this is changed by the "random" and "ErrorBound" options!
GENSS.VerifyElements := 10;   # this amounts to an error bound of 1/1024
GENSS.DeterministicVerification := false;
# Number of random elements used to do immediate verification:
GENSS.ImmediateVerificationElements := 3;
# Are we working projectively?
GENSS.Projective := false;
# To find a very short orbit we try a random vector:
GENSS.VeryShortOrbLimit := 1000;
# Never consider more than this number of candidates for short orbits:
GENSS.LimitShortOrbCandidates := 100;

#############################################################################
# A few helper functions needed elsewhere:
#############################################################################

InstallGlobalFunction( GENSS_CopyDefaultOptions,
  function( defopt, opt )
    local n;
    for n in RecFields(defopt) do
        if not(IsBound(opt.(n))) then
            opt.(n) := defopt.(n);
        fi;
    od;
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
  function(g,opt)
    # Needs opt components "ShortOrbitsNrRandoms"
    local c,f,i,inters,j,l,nw,pr,sortfun,v,vv,w,wb,ww;
    Info(InfoGenSS,2,"Trying Murray/O'Brien heuristics...");
    l := ShallowCopy(GeneratorsOfGroup(g));
    f := DefaultFieldOfMatrixGroup(g);
    if HasPseudoRandomSeed(g) then
        for i in [1..opt.ShortOrbitsNrRandoms] do
            Add(l,PseudoRandom(g));
        od;
    else
        pr := ProductReplacer(GeneratorsOfGroup(g),
                      rec( scramble := opt.ShortOrbScramble,
                           scramblefactor := opt.ShortOrbScrambleFactor,
                           addslots := opt.ShortOrbAddSlots,
                           maxdepth := opt.ShortOrbMaxDepth ));
        for i in [1..opt.ShortOrbitsNrRandoms] do
            Add(l,Next(pr));
        od;
    fi;
    if IsObjWithMemory(l[1]) then
        ForgetMemory(l);
    fi;
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
    i := 1;
    while i <= Length(l) and Length(w) < GENSS.LimitShortOrbCandidates do
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
        i := i + 1;
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
  function( g,opt )
    # Needs opt components:
    #  "ShortOrbitsNrRandoms"  (because it uses GENSS_FindVectorsWithShortOrbit)
    #  "ShortOrbitsOrbLimit"
    #  "ShortOrbitsStartLimit"
    local ThrowAwayOrbit,found,gens,hashlen,i,j,limit,newnrorbs,nrorbs,o,wb;

    wb := GENSS_FindVectorsWithShortOrbit(g,opt);

    # Now we have a list of vectors with (hopefully) short orbits.
    # We start enumerating all those orbits, but first only 50 elements:
    nrorbs := Minimum(Length(wb),64);  # take only the 64 first
    gens := GeneratorsOfGroup(g);
    o := [];
    hashlen := NextPrimeInt(QuoInt(opt.ShortOrbitsOrbLimit * 3,2));
    for i in [1..nrorbs] do
        Add(o,Orb(gens,ShallowCopy(wb[i]),OnLines,hashlen));
    od;
    limit := opt.ShortOrbitsStartLimit;
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
          if limit > opt.ShortOrbitsOrbLimit then
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
    Info(InfoGenSS,1,"Found short orbit of length ",Length(o[i])," (#",i,").");
    return o[i];
  end );   

InstallGlobalFunction( GENSS_IsOneProjective,
  function(el)
    local s;
    s := el[1][1];
    if IsZero(s) then return false; fi;
    s := s^-1;
    return IsOne( s*el );
  end );


#############################################################################
# Now to the heart of the method, the Schreier-Sims:
#############################################################################

InstallMethod( FindBasePointCandidates, "for a matrix group over a FF",
  [ IsGroup and IsMatrixGroup and IsFinite, IsRecord, IsInt ],
  function( grp, opt, i )
    local bv,cand,d,F,o,gens,v,w,res;
    F := DefaultFieldOfMatrixGroup(grp);
    d := DimensionOfMatrixGroup(grp);
    Info( InfoGenSS, 2, "Finding nice base points..." );
    if IsObjWithMemory(GeneratorsOfGroup(grp)[1]) then
        grp := Group(StripMemory(GeneratorsOfGroup(grp)));
    fi;

    # Try a random vector to find a very short orbit:
    Info( InfoGenSS, 3, "Trying random vector for very short orbit..." );
    gens := GeneratorsOfGroup(grp);
    v := ShallowCopy( gens[1][1] );
    Randomize(v);
    o := Orb(gens,v,OnRight);
    Enumerate(o,opt.VeryShortOrbLimit);
    if Length(o) > 1 and Length(o) < opt.VeryShortOrbLimit then
        Info( InfoGenSS, 3, "Found orbit of length ",Length(o) );
        if Size(F) = 2 then
            cand := rec( points := [v], ops := [OnPoints], used := 0 );
        else 
            v := ORB_NormalizeVector(v);
            cand := rec( points := [v], ops := [OnLines],
                         used := 0 );
            # Note that if we work non-projectively, then the same
            # point will be taken next with OnRight automatically!
        fi;
        return cand;
    fi;
    Info( InfoGenSS, 3, "Found no very short orbit up to limit ",
          opt.VeryShortOrbLimit );

    # Next possibly "TryShortOrbit":
    cand := rec( points := [], used := 0 );
    if IsBound(opt.TryShortOrbit) and opt.TryShortOrbit > 0 then
        repeat
            opt.TryShortOrbit := opt.TryShortOrbit - 1;
            Info(InfoGenSS,1,"Looking for short orbit (",opt.TryShortOrbit,
                 ")...");
            res := GENSS_FindShortOrbit(grp,opt);
        until res <> fail or opt.TryShortOrbit = 0;
        if res <> fail then
            if Size(F) > 2 then
                cand.points := [res[1]];
                cand.ops := [OnLines];
                # Note that if we work non-projectively, then the same
                # point will be taken next with OnRight automatically!
            else
                cand.points := [res[1]];
                cand.ops := [OnRight];
            fi;
            return cand;
        fi;
    fi;

    # Otherwise standard Murry/O'Brien heuristics:
    if i = 0 and 
       ((opt!.Projective = false and Size(F)^d > 300000) or
        (opt!.Projective = true and Size(F)^(d-1) > 300000)) then
        bv := GENSS_FindVectorsWithShortOrbit(grp,opt);
    else
        bv := One(grp);
    fi;
    cand := rec( points := [], ops := [], used := 0 );
    for v in bv do
        w := ORB_NormalizeVector(ShallowCopy(v));
        if Size(F) = 2 then
            Add(cand.points,w);
            Add(cand.ops,OnRight);
        else
            Add(cand.points,w);
            Add(cand.ops,OnLines);
            # Note that if we work non-projectively, then the same
            # point will be taken next with OnRight automatically!
        fi;
    od;
    return cand;
  end );

InstallMethod( FindBasePointCandidates, "for a permutation group",
  [ IsGroup and IsPermGroup, IsRecord, IsInt ],
  function( grp, opt, i )
    local ops,points;
    if i = 0 then
        points := [1..Minimum(20,LargestMovedPoint(grp))];
    else
        points := [1..LargestMovedPoint(grp)];
    fi;
    ops := List([1..Length(points)],x->OnPoints);
    return rec( points := points, ops := ops, used := 0 );
  end );
    
InstallGlobalFunction( GENSS_NextBasePoint, 
  function( gens, cand, opt, S )
    local NotFixedUnderAllGens,i;

    NotFixedUnderAllGens := function( gens, x, op )
      return ForAny( gens, g->op(x,g) <> x );
    end;

    # S can be false or a stabilizer chain record
    if S <> false and not(opt.StrictlyUseCandidates) then  
        # try points in previous orbit
        # Maybe we can take the last base point now acting non-projectively?
        if opt.Projective = false and IsIdenticalObj(S!.orb!.op,OnLines) and 
           NotFixedUnderAllGens(gens,S!.orb[1],OnRight) then
            return rec( point := S!.orb[1], op := OnRight, cand := cand );
        fi;
        for i in [2..Minimum(Length(S!.orb),opt.NumberPrevOrbitPoints)] do
            if NotFixedUnderAllGens(gens,S!.orb[i],S!.orb!.op) then
                return rec( point := S!.orb[i], op := S!.orb!.op, 
                            cand := cand );
            fi;
        od;
    fi;

    repeat
        if cand.used >= Length(cand.points) then
            cand := FindBasePointCandidates(Group(gens),opt,1);
            opt.StrictlyUseCandidates := false;
        fi;
        cand.used := cand.used + 1;
    until NotFixedUnderAllGens(gens,cand.points[cand.used],cand.ops[cand.used]);

    return rec( point := cand.points[cand.used], op := cand.ops[cand.used],
                cand := cand );
  end );



InstallGlobalFunction( GENSS_CreateStabChainRecord,
  function( gens, size, layer, nextpoint, nextop, base, cand, opt )
    local S,hashsize;

    Info( InfoGenSS, 2, "Creating new stab chain record..." );

    gens := ShallowCopy(gens);
    # Note that we do ShallowCopy such that the original list and the
    # one in the orbit record are different from each other.
    S := rec( stab := false, cand := cand, size := size, base := base,
              opt := opt, layer := layer );

    if IsInt(size) then
        hashsize := NextPrimeInt(Minimum(size,opt.InitialHashSize));
    else
        hashsize := opt.InitialHashSize;
    fi;
    S.orb := Orb( gens, nextpoint, nextop,
                  rec( hashlen := hashsize, schreier := true, log := true,
                       report := opt.Report ) );
    Add(base,nextpoint);
    S.orb!.stabilizerchain := S;
    Objectify( StabChainByOrbType, S );

    return S;
  end );

InstallGlobalFunction( GENSS_StabilizerChainInner,
  function( gens, size, layer, cand, opt, parentS )
    # Computes a stabilizer chain for the group generated by gens
    # with known size size (can be false if not known). This will be
    # layer layer in the final stabilizer chain. cand is a (shared)
    # record for base point candidates and opt the (shared) option
    # record. This is called in StabilizerChain and calls itself.
    # It also can be called if a new layer is needed.
    local base,gen,S,i,next,pr,r,stabgens,x;

    Info(InfoGenSS,2,"Entering GENSS_StabilizerChainInner layer=",layer);
    next := GENSS_NextBasePoint(gens,cand,opt,parentS);
    cand := next.cand;   # This could have changed
    if parentS <> false then
        base := parentS!.base;
    else
        base := [];
    fi;
    S := GENSS_CreateStabChainRecord(gens,size,layer,
             next.point,next.op,base,next.cand,opt);

    Info( InfoGenSS, 3, "Entering orbit enumeration layer ",layer,"..." );
    repeat
        Enumerate(S!.orb,opt.OrbitLengthLimit);
        if not(IsClosed(S!.orb)) then
            Error("Orbit too long, increase opt.OrbitLengthLimit");
        fi;
    until IsClosed(S!.orb);
    Info(InfoGenSS, 1, "Layer ", layer, ": Orbit length is ", Length(S!.orb)); 
    if Length(S!.orb) > 50 or S!.orb!.depth > 5 then
        Info(InfoGenSS, 2, "Trying to make Schreier tree shallower...");
        MakeSchreierTreeShallow(S!.orb);
        Info(InfoGenSS, 2, "Depth is now ",S!.orb!.depth);
    fi;
    Info( InfoGenSS, 3, "Done orbit enumeration layer ",layer );
    S!.orb!.gensi := List(S!.orb!.gens,x->x^-1);

    # Are we done?
    if size <> false and Length(S!.orb) = size then
        S!.proof := true;
        Info(InfoGenSS,2,"Leaving GENSS_StabilizerChainInner layer=",layer);
        return S;
    fi;

    # Now create a few random stabilizer elements:
    Info(InfoGenSS,2,"Creating ",opt.RandomStabGens,
         " random elements of the point stabilizer...");
    pr := ProductReplacer( gens,
                      rec( scramble := S!.opt.StabGenScramble,
                           scramblefactor := S!.opt.StabGenScrambleFactor,
                           addslots := S!.opt.StabGenAddSlots,
                           maxdepth := S!.opt.StabGenMaxDepth ));
    S!.pr := pr;   # for later use
    stabgens := [];
    for i in [1..opt.RandomStabGens] do
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
        if opt.ImmediateVerificationElements > 0 then
            Info(InfoGenSS,1,"Doing immediate verification in layer ",
                 S!.layer," (",opt.ImmediateVerificationElements,
                 " elements)...");
            i := 0;
            while i < opt.ImmediateVerificationElements do
                i := i + 1;
                x := Next(pr);
                if AddGeneratorToStabilizerChain(S,x) then
                    Info( InfoGenSS, 1, "Immediate verification found error ",
                          "(layer ",S!.layer,")..." );
                    i := 0;
                fi;
            od;
        fi;

        S!.proof := S!.stab!.proof;   # hand up information
    else
        # We are not sure that the next stabiliser is trivial, but we believe!
        Info(InfoGenSS,2,"Found no non-trivial ones.");
        S!.proof := false;
    fi;

    Info(InfoGenSS,2,"Leaving GENSS_StabilizerChainInner layer=",layer);
    return S;
  end );

InstallGlobalFunction( GENSS_DeriveCandidatesFromStabChain,
  function( S )
    local cand;
    cand := rec( points := [], ops := [], used := 0 );
    while S <> false do
        Add( cand.points, S!.orb[1] );
        Add( cand.ops,    S!.orb!.op );
        S := S!.stab;
    od;
    return cand;
  end );

InstallMethod( StabilizerChain, "for a group object", [ IsGroup ],
  function( grp )
    return StabilizerChain( grp, rec() );
  end );

InstallMethod( StabilizerChain, "for a group object", [ IsGroup, IsRecord ],
  function( grp, opt )
    # Computes a stabilizer chain for the group grp
    # Possible options:
    #   random:     compatibility to StabChain, works as there
    #   ErrorBound: rational giving the error bound, probability of an
    #               error will be proven to be lower than this
    #               if given, this sets VerifyElements and
    #                                   DeterministicVerification
    #   VerifyElements: number of random elements for verification
    #   DeterministicVerification: flag, whether to do or not
    #                              (not yet implemented)
    #   Base:       specify a known base, this can either be the 
    #               StabilizerChain of a known supergroup or a list of
    #               points, if "BaseOps" is not given, it is either
    #               set to OnPoints if Projective is not set or false,
    #               or is set to OnLines if Projective is set to true,
    #               if Base is set
    #   BaseOps:    a list of operation functions for the base points
    #   Projective: if set to true indicates that the matrices given
    #               as generators are to be thought as projective
    #               elements, that is, the whole StabilizerChain will
    #               be one for a projective group!
    #   StrictlyUseCandidates: if set to true exactly the points in
    #               opt.cand will be used as base points and no other
    #               points from previous orbits are used first
    #               this is set when Base was specified
    #   Cand:       initializer for cand, which are base points
    #               candidates
    #               must be a record with components "points", "ops",
    #               "used", the latter indicates the largest index
    #               that has already been used.
    #   TryShortOrbit: Number of tries for the short orbit finding alg.
    #   StabGenScramble,
    #   StabGenScrambleFactor,
    #   StabGenAddSlots,
    #   StabGenMaxDepth:   parameters for product replacer for generating
    #                      stabilizer elements when size is known
    #   VerifyScramble,
    #   VerifyScrambleFactor,
    #   VerifyAddSlots,
    #   VerifyMaxDepth:    parameters for product replacer for verification
    #                      phase
    #   VeryShortOrbLimit: when looking for short orbits try a random
    #                      vector and enumerate its orbit until this limit
    #   
    #   ... to be continued
    #
    local S,cand,i,pr,prob,x;

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
    # The new style error bounds:
    if IsBound(opt.ErrorBound) then
        prob := 1/2;
        opt.VerifyElements := 1;
        while prob >= opt.ErrorBound do
            prob := prob / 2;
            opt.VerifyElements := opt.VerifyElements + 1;
        od;
    fi;
                
    # Find base point candidates:
    if IsBound(opt.Base) then
        if IsStabilizerChain(opt.Base) then
            cand := GENSS_DeriveCandidatesFromStabChain(opt.Base);
        else  # directly take the base points:
            cand := rec( points := opt.Base, used := 0 );
            if not(IsBound(opt.BaseOps)) then
                # Let's guess the ops:
                if IsBound(opt.Projective) and opt.Projective = true then
                    cand.ops := ListWithIdenticalEntries(Length(opt.Base),
                                                         OnLines);
                else
                    cand.ops := ListWithIdenticalEntries(Length(opt.Base),
                                                         OnPoints);
                fi;
            else
                cand.ops := opt.BaseOps;
            fi;
        fi;
        opt.StrictlyUseCandidates := true;
    elif IsBound(opt.Cand) then
        cand := opt.Cand;
    else
        # Otherwise try Murray/O'Brien:
        cand := FindBasePointCandidates( grp, opt, 0 );
    fi;
    if not(IsBound(opt.StrictlyUseCandidates)) then
        opt.StrictlyUseCandidates := false;
    fi;
    if HasSize(grp) and not(opt.Projective) then
        S := GENSS_StabilizerChainInner(GeneratorsOfGroup(grp), Size(grp), 1,
                                         cand, opt, false);
    elif IsBound(opt.Size) then
        S := GENSS_StabilizerChainInner(GeneratorsOfGroup(grp), opt.Size, 1,
                                         cand, opt, false);
    else
        S := GENSS_StabilizerChainInner(GeneratorsOfGroup(grp), false, 1, 
                                         cand, opt, false);
    fi;

    # Do we already have a proof?
    if S!.proof or opt.VerifyElements = 0 then return S; fi;

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
        S!.proof := true;
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
        Info( InfoGenSS, 3, "Entering orbit enumeration layer ",SS!.layer,
              "..." );
        repeat
            Enumerate(SS!.orb,S!.opt.OrbitLengthLimit);
            if not(IsClosed(SS!.orb)) then
                Error("Orbit too long, increase S!.opt.OrbitLengthLimit!");
            fi;
        until IsClosed(SS!.orb);
        Info( InfoGenSS, 3, "Done orbit enumeration layer ",SS!.layer,"..." );
        SS!.proof := false;
    else   # case (2)
        # Note that we do not create a pr instance here for one
        # generator, this will be done later on as needed...
        SS := r.preS;
        subsize := Order(r.rem);
        SS!.stab := GENSS_StabilizerChainInner([r.rem],subsize,
                           SS!.layer+1,SS!.cand, SS!.opt, SS );
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

InstallMethod( SiftGroupElement, "for a stabilizer chain and a group element",
  [ IsStabilizerChain and IsStabilizerChainByOrb, IsObject ],
  function( S, x )
    local o,p,po,preS,r;
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
    r := rec( rem := x, S := false, preS := preS );
    if preS!.opt.Projective then
        r.isone := GENSS_IsOneProjective(x);
    else
        r.isone := IsOne(x);
    fi;
    return r;
  end );

InstallMethod( SiftGroupElementSLP, 
  "for a stabilizer chain and a group element",
  [ IsStabilizerChain and IsStabilizerChainByOrb, IsObject ],
  function( S, x )
    local nrstrong,o,p,po,preS,r,slp;
    preS := false;
    slp := [];     # will be reversed in the end
    nrstrong := 0; # we count here the number of generators
    while S <> false do
        o := S!.orb;
        p := o!.op(o[1],x);
        po := Position(o,p);
        if po = fail then   # not in current stabilizer
            return rec( isone := false, rem := x, S := S, preS := preS,
                        slp := fail );
        fi;
        # Now sift through Schreier tree:
        while po > 1 do
            x := x * S!.orb!.gensi[o!.schreiergen[po]];
            Add(slp,1);
            Add(slp,nrstrong + o!.schreiergen[po]);
            po := o!.schreierpos[po];
        od;
        nrstrong := nrstrong + Length(S!.orb!.gens);
        preS := S;
        S := S!.stab;
    od;
    r := rec( rem := x, S := false, preS := preS );
    if preS!.opt.Projective then
        r.isone := GENSS_IsOneProjective(x);
    else
        r.isone := IsOne(x);
    fi;
    if r.isone then
        if Length(slp) = 0 then   # element was the identity!
            r.slp := StraightLineProgramNC([[1,0]],nrstrong);
        else
            r.slp := StraightLineProgramNC(
                           [slp{[Length(slp),Length(slp)-1..1]}],nrstrong);
        fi;
    else
        r.slp := fail;
    fi;
    return r;
  end );

InstallMethod( StrongGenerators, "for a stabilizer chain",
  [ IsStabilizerChain and IsStabilizerChainByOrb ],
  function( S )
    local gens;
    gens := [];
    while S <> false do
        Append(gens,S!.orb!.gens);
        S := S!.stab;
    od;
    return gens;
  end );

InstallMethod( NrStrongGenerators, "for a stabilizer chain",
  [ IsStabilizerChain and IsStabilizerChainByOrb ],
  function( S )
    local nrgens;
    nrgens := 0;
    while S <> false do
        nrgens := nrgens + Length(S!.orb!.gens);
        S := S!.stab;
    od;
    return nrgens;
  end );

InstallMethod( ForgetMemory, "for a stabilizer chain",
  [ IsStabilizerChain and IsStabilizerChainByOrb ],
  function( S )
    while S <> false do
        ForgetMemory(S!.orb);
        S := S!.stab;
    od;
  end );

#############################################################################
# The following operations apply to stabilizer chains:
#############################################################################

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

InstallMethod( \in, "for a group element and a stabilizer chain",
  [ IsObject, IsStabilizerChain and IsStabilizerChainByOrb ],
  function( el, S )
    local r;
    r := SiftGroupElement(S,el);
    if r.isone then
        return true;
    else
        return false;
    fi;
  end );
    
GENSS_VIEWDEPTH := 0;  # to please the interpreter
InstallMethod( ViewObj, "for a stabilizer chain",
  [ IsStabilizerChain and IsStabilizerChainByOrb ],
  function( S )
    local i;
    if not(IsBound(GENSS_VIEWDEPTH)) then
        GENSS_VIEWDEPTH := 0;
    else
        GENSS_VIEWDEPTH := GENSS_VIEWDEPTH + 1;
    fi;
    for i in [1..GENSS_VIEWDEPTH] do Print(" "); od;
    Print("<stabchain size=",Size(S)," orblen=",Length(S!.orb),
          " layer=",S!.layer," SchreierDepth=",S!.orb!.depth,">");
    if S!.stab <> false then
        Print("\n");
        ViewObj(S!.stab);
    fi;
    GENSS_VIEWDEPTH := GENSS_VIEWDEPTH - 1;
    if GENSS_VIEWDEPTH < 0 then
        Unbind(GENSS_VIEWDEPTH);
    fi;
  end );
Unbind(GENSS_VIEWDEPTH);


#############################################################################
# We can store a stabilizer chain in the (mutable) attribute
# StabChainMutable, then the following methods for matrix group apply:
#############################################################################

InstallMethod( Size, "for a finite matrix group with a stabilizer chain",
  [ IsMatrixGroup and IsFinite and HasStabChainMutable ],
  function( g )
    return Size(StabChainMutable(g));
  end );

InstallMethod( \in, "for a matrix and a finite matrix group with a stabchain",
  [ IsMatrix and IsFFECollColl, 
    IsMatrixGroup and IsFinite and HasStabChainMutable ],
  function( el, g )
    local S,r;
    S := StabChainMutable(g);
    r := SiftGroupElement(S,el);
    if r.isone then
        return true;
    else
        return false;
    fi;
  end );

InstallMethod( IsProved, "for a stabilizer chain",
  [ IsStabilizerChain and IsStabilizerChainByOrb ],
  function( S )
    return S!.proof;
  end );

InstallMethod( MakeGAPStabChain, "for a stabilizer chain",
  [ IsStabilizerChain and IsStabilizerChainByOrb ],
  function( S )
    local first,last,s,ss;
    if not(IsPermOnIntOrbitRep(S!.orb)) then
        Error("Can only work with permutations acting on integers!");
        return fail;
    fi;
    s := rec();
    if S!.stab <> false then
        ss := MakeGAPStabChain(S!.stab);
        s.labels := ss.labels;
        s.stabilizer := ss;
    else
        s.labels := [];
        s.stabilizer := rec( labels := s.labels, genlabels := [],
                             generators := [], identity := S!.orb!.gens[1]^0 );
    fi;
    s.orbit := List(S!.orb);
    first := Length(s.labels)+1;
    Append(s.labels,List(S!.orb!.gens,x->x^-1));
    last := Length(s.labels);
    s.genlabels := [first..last];
    s.generators := s.labels{s.genlabels};
    s.identity := s.generators[1]^0;
    s.translabels := [1];
    Append(s.translabels,
           List(S!.orb!.schreiergen{[2..Length(S!.orb!.schreiergen)]},
                x->x+first-1));
    s.transversal := s.labels{s.translabels};
    return s;
  end );
 
