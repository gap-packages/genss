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
# Do not throw Errors but return fail:
GENSS.FailInsteadOfError := false;
# Number of Schreier generators to create in TC verification:
GENSS.NumberSchreierGens := 20;
# Maximal number of Schreier generators to create in TC verification:
GENSS.MaxNumberSchreierGens := 100;
# By default do 10 short orbit tries:
GENSS.TryShortOrbit := 10;

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
    Info(InfoGenSS,3,"Trying Murray/O'Brien heuristics...");
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
            Add(vv,[NullspaceMat(Value(c[i][j],l[i])),
                    Degree(c[i][j]),
                    WeightVecFFE(CoefficientsOfLaurentPolynomial(c[i][j])[1]),
                    1]);
        od;
        Add(v,vv);
    od;
    Info(InfoGenSS,3,"Have eigenspaces.");
    # Now collect a list of all those spaces together with all
    # possible intersects:
    w := [];
    i := 1;
    while i <= Length(l) and Length(w) < GENSS.LimitShortOrbCandidates do
        nw := [];
        for j in [1..Length(v[i])] do
            for ww in w do
                inters := SumIntersectionMat(ww[1],v[i][j][1])[2];
                if Length(inters) > 0 then
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
        elif Length(a[1]) < Length(b[1]) then return true;
        else return false;
        fi;
    end;
    Sort(w,sortfun);
    wb := List(w,ww->ww[1][1]);
    Info(InfoGenSS,3,"Have ",Length(wb)," vectors for possibly short orbits.");
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
            Info(InfoGenSS,3,"Orbit Number ",i," has length 1.");
            found := false;
            # Now throw away this orbit:
            ThrowAwayOrbit(i);
            # we intentionally do not increase i here!
        elif not(found) then
            i := i + 1;
        fi;
        if i > nrorbs then
          Info(InfoGenSS,3,"Done ",nrorbs,
               " orbit(s) to limit ",limit,".");
          limit := limit * 2;
          if limit > opt.ShortOrbitsOrbLimit then
              Info(InfoGenSS,3,"Limit reached, giving up.");
              return fail;
          fi;
          i := 1;
          if nrorbs < i then
              Info(InfoGenSS,3,"No orbits left, giving up.");
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
    Info(InfoGenSS,2,"Found short orbit of length ",Length(o[i])," (#",i,").");
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
    Info( InfoGenSS, 3, "Finding nice base points..." );
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

    Info( InfoGenSS, 3, "Creating new stab chain record..." );

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

InstallGlobalFunction( GENSS_ComputeStrongBelowNumbers,
  function( S )
    if S!.stab = false then
        S!.strongbelow := 0;
        S!.nrstrong := Length(S!.orb!.gens);
        return S!.nrstrong;
    else
        S!.strongbelow := GENSS_ComputeStrongBelowNumbers(S!.stab);
        S!.nrstrong := S!.strongbelow + Length(S!.orb!.gens);
        return S!.nrstrong;
    fi;
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

    Info(InfoGenSS,4,"Entering GENSS_StabilizerChainInner layer=",layer);
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
            if opt.FailInsteadOfError then
                return "Orbit too long, increase opt.OrbitLengthLimit";
            else
                Error("Orbit too long, increase opt.OrbitLengthLimit");
            fi;
        fi;
    until IsClosed(S!.orb);
    Info(InfoGenSS, 1, "Layer ", layer, ": Orbit length is ", Length(S!.orb)); 
    if Length(S!.orb) > 50 or S!.orb!.depth > 5 then
        Info(InfoGenSS, 3, "Trying to make Schreier tree shallower...");
        MakeSchreierTreeShallow(S!.orb);
        Info(InfoGenSS, 3, "Depth is now ",S!.orb!.depth);
    fi;
    Info( InfoGenSS, 3, "Done orbit enumeration layer ",layer );
    S!.orb!.gensi := List(S!.orb!.gens,x->x^-1);

    # Are we done?
    if size <> false and Length(S!.orb) = size then
        S!.proof := true;
        Info(InfoGenSS,4,"Leaving GENSS_StabilizerChainInner layer=",layer);
        return S;
    fi;

    # Now create a few random stabilizer elements:
    Info(InfoGenSS,3,"Creating ",opt.RandomStabGens,
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
        Info(InfoGenSS,3,"Found ",Length(stabgens)," non-trivial ones.");
        if size <> false then
            S!.stab := GENSS_StabilizerChainInner(stabgens,size/Length(S!.orb),
                                                   layer+1,cand,opt,S);
        else
            S!.stab := GENSS_StabilizerChainInner(stabgens,false,
                                                   layer+1,cand,opt,S);
        fi;
        if IsString(S!.stab) then return S!.stab; fi; 
        if opt.ImmediateVerificationElements > 0 then
            Info(InfoGenSS,2,"Doing immediate verification in layer ",
                 S!.layer," (",opt.ImmediateVerificationElements,
                 " elements)...");
            i := 0;
            while i < opt.ImmediateVerificationElements do
                i := i + 1;
                x := Next(pr);
                if AddGeneratorToStabilizerChain(S,x) then
                    Info( InfoGenSS, 2, "Immediate verification found error ",
                          "(layer ",S!.layer,")..." );
                    i := 0;
                fi;
            od;
        fi;

        S!.proof := S!.stab!.proof;   # hand up information
    else
        # We are not sure that the next stabiliser is trivial, but we believe!
        Info(InfoGenSS,3,"Found no non-trivial ones.");
        S!.proof := false;
    fi;

    Info(InfoGenSS,4,"Leaving GENSS_StabilizerChainInner layer=",layer);
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

InstallGlobalFunction( GENSS_TrivialOp,
  function( p, x )
    return p;
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
    local S,cand,i,pr,prob,x,gens;

    # First a few preparations, then we delegate to GENSS_StabilizerChainInner:

    # Add some default options:
    GENSS_CopyDefaultOptions(GENSS,opt);

    # Check for the identity group:
    gens := GeneratorsOfGroup(grp);
    if Length(gens) = 0 or ForAll(gens,IsOne) then
        # Set up a trivial stabilizer chain record:
        S := GENSS_CreateStabChainRecord(gens,1,1,1,GENSS_TrivialOp,
                                         [],false,opt);
        Enumerate(S!.orb);
        S!.proof := true;
        S!.trivialgroup := true;
        GENSS_ComputeStrongBelowNumbers(S);
        return S;
    fi;
    
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
    if IsString(S) then return S; fi;
    GENSS_ComputeStrongBelowNumbers(S);

    # Do we already have a proof?
    if S!.proof or opt.VerifyElements = 0 then return S; fi;

    Info(InfoGenSS,2,"Current size found: ",Size(S));
    # Now a possible verification phase:
    if S!.size <> false then   # we knew the size in advance
        Info(InfoGenSS,2,"Doing verification via known size...");
        pr := ProductReplacer(GeneratorsOfGroup(grp),
                      rec( scramble := opt.StabGenScramble,
                           scramblefactor := opt.StabGenScrambleFactor,
                           addslots := opt.StabGenAddSlots,
                           maxdepth := opt.StabGenMaxDepth ));
        while Size(S) < Size(grp) do
            Info(InfoGenSS,2,"Known size not reached, throwing in a random ",
                 "element...");
            x := Next(pr);
            if AddGeneratorToStabilizerChain(S,x) then
                Info( InfoGenSS, 2, "Increased size to ",Size(S) );
            fi;
        od;
        S!.proof := true;
    else
        # Do some verification here:
        Info(InfoGenSS,2,"Doing randomized verification...");
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
                Info( InfoGenSS, 2, "Verification found error ... ",
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
    local i,r,SS,subsize,origS;
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
        Info( InfoGenSS, 2, "Adding new generator to stabilizer chain ",
              "in layer ", SS!.layer, "..." );
        AddGeneratorsToOrbit(SS!.orb,[r.rem]);
        Add(SS!.orb!.gensi,r.rem^-1);
        Info( InfoGenSS, 4, "Entering orbit enumeration layer ",SS!.layer,
              "..." );
        repeat
            Enumerate(SS!.orb,S!.opt.OrbitLengthLimit);
            if not(IsClosed(SS!.orb)) then
                if S!.opt.FailInsteadOfError then
                    return "Orbit too long, increase S!.opt.OrbitLengthLimit";
                else
                    Error("Orbit too long, increase S!.opt.OrbitLengthLimit!");
                fi;
            fi;
        until IsClosed(SS!.orb);
        Info( InfoGenSS, 4, "Done orbit enumeration layer ",SS!.layer,"..." );
        SS!.proof := false;
    else   # case (2)
        # Note that we do not create a pr instance here for one
        # generator, this will be done later on as needed...
        SS := r.preS;
        subsize := Order(r.rem);
        SS!.stab := GENSS_StabilizerChainInner([r.rem],false,
                           SS!.layer+1,SS!.cand, SS!.opt, SS );
        if IsString(SS!.stab) then return SS!.stab; fi; 
        SS := SS!.stab;
    fi;
    # Now we have added a new generator (or a new layer) at layer SS,
    # the new gen came from layer S (we were called here, after all),
    # thus we have to check, whether all the orbits between S (inclusively)
    # and SS (exclusively) are also closed under the new generator r.rem,
    # we add it to all these orbits, thereby also making the Schreier trees
    # shallower:
    origS := S;
    while S!.layer <> SS!.layer do
        Info(InfoGenSS,3,"Adding new generator to orbit at layer ",S!.layer);
        AddGeneratorsToOrbit(S!.orb,[r.rem]);
        Add(S!.orb!.gensi,r.rem^-1);
        S := S!.stab;
    od;
    GENSS_ComputeStrongBelowNumbers(origS);
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
    local o,p,po,preS,r,slp,nrstrong;
    nrstrong := S!.nrstrong;
    preS := false;
    slp := [];     # will be reversed in the end
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
            Add(slp,S!.strongbelow + o!.schreiergen[po]);
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
        Add(gens,S!.orb!.gens);
        S := S!.stab;
    od;
    return Concatenation(Reversed(gens));
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
        ForgetMemory(S!.orb!.gensi);
        S := S!.stab;
    od;
  end );

InstallMethod( GENSS_CreateSchreierGenerator, 
  "for a stabilizer chain, a position in the first orbit and a gen number",
  [ IsStabilizerChain and IsStabilizerChainByOrb, IsPosInt, IsPosInt ],
  function( S, i, j )
    local o,p,w1,w2;
    o := S!.orb;
    p := Position(o,o!.op(o[i],o!.gens[j]));
    if o!.schreierpos[p] = i and o!.schreiergen[p] = j then
        return fail;
    fi;
    w1 := TraceSchreierTreeForward(o,i);
    w2 := TraceSchreierTreeBack(o,p);
    return [w1,j,w2];
  end );

InstallGlobalFunction( GENSS_Prod,
  function( gens, word )
    if Length(word) = 0 then 
        return gens[1]^0; 
    else
        return Product(gens{word});
    fi;
  end );
  
InstallGlobalFunction( VerifyStabilizerChainTC,
  function( S )
    local Grels,Hrels,Prels,MakeSchreierGens,ct,f,gens,gensi,i,j,k,l,li,max,
          newpres,nrgens,nrschr,o,ords,pres,sb,sgs,slp,subgens,v,w,x,
          cosetnrlimitfactor;
    if S!.stab <> false then
        pres := VerifyStabilizerChainTC(S!.stab);
        if IsList(pres) then return pres; fi;
    else
        pres := StraightLineProgram([[]],0);
    fi;
    Info(InfoGenSS,1,"Verifying stabilizer chain in layer ",S!.layer);
    # First create a few Schreier generators:
    sgs := [];
    i := 1;
    j := 1;
    o := S!.orb;
    nrgens := Length(o!.gens);
    sb := S!.strongbelow;
    MakeSchreierGens := function(n)
        local sg;
        Info(InfoGenSS,3,"Creating ",n," Schreier generators...");
        while Length(sgs) < n and
              i <= Length(o) do
            sg := GENSS_CreateSchreierGenerator(S,i,j);
            j := j + 1;
            if j > nrgens then
                j := 1;
                i := i + 1;
            fi;
            if sg <> fail then
                Add(sgs,sg);
            fi;
        od;
    end;

    nrschr := S!.opt.NumberSchreierGens;
    MakeSchreierGens(nrschr);
    f := FreeGroup(S!.nrstrong);
    gens := GeneratorsOfGroup(f);
    gensi := List(gens,x->x^-1);
    subgens := gens{[1..S!.strongbelow]};
    Hrels := ResultOfStraightLineProgram(pres,subgens);
    if S!.opt.Projective then
        ords := List([1..nrgens],i->ProjectiveOrder(o!.gens[i]));
    else
        ords := List([1..nrgens],i->Order(o!.gens[i]));
    fi;
    Prels := List([1..nrgens],i->gens[i+sb]^ords[i]);
    Grels := [];
    cosetnrlimitfactor := 4;
    while true do   # will be left by return eventually
        for k in [Length(Grels)+1..Length(sgs)] do
            Grels[k] := GENSS_Prod(gens,sgs[k][1]+sb) * gens[sgs[k][2]+sb] * 
                        GENSS_Prod(gensi,sgs[k][3]+sb);
            x := GENSS_Prod(o!.gens,sgs[k][1]) * o!.gens[sgs[k][2]] * 
                 GENSS_Prod(o!.gensi,sgs[k][3]);
            if S!.stab <> false then
                slp := SiftGroupElementSLP(S!.stab,x);
                if not(slp.isone) then
                    return [fail,S!.layer];
                fi;
                Grels[k] := Grels[k] / ResultOfStraightLineProgram(slp.slp,
                                                                   subgens);
                sgs[k][4] := slp.slp;
            else
                if not(IsOne(x)) then
                    return [fail,S!.layer];
                fi;
                sgs[k][4] := false;
            fi;
        od;
        Info(InfoGenSS,2,"Doing coset enumeration with limit ",
             cosetnrlimitfactor*Length(o));
        ct := CosetTableFromGensAndRels(gens,Concatenation(Prels,Hrels,Grels),
                   subgens:max := cosetnrlimitfactor*Length(o),silent);
        if ct = fail then   # did not close!
            cosetnrlimitfactor := QuoInt(cosetnrlimitfactor*3,2);
            Info(InfoGenSS,2,"Coset enumeration did not finish!");
        #Error(1);
            if nrschr > Length(sgs) # or
               # nrschr > S!.opt.MaxNumberSchreierGens 
               then   # we are done!
                # Something is wrong!
                return [fail, S!.layer];
            fi;
        else
            Info(InfoGenSS,2,"Coset enumeration found ",Length(ct[1]),
                 " cosets.");
        #Error(2);
            if Length(ct[1]) = Length(o) then
                # Verification is OK, now build a presentation:
                l := GeneratorsWithMemory(
                       ListWithIdenticalEntries(S!.nrstrong,()));
                li := List(l,x->x^-1);
                newpres := ResultOfStraightLineProgram(pres,
                                   l{[1..S!.strongbelow]});
                for k in [1..nrgens] do
                    Add(newpres,l[k+sb]^ords[k]);
                od;
                for k in [1..Length(sgs)] do
                    if sgs[k][4] <> false then
                        Add(newpres,
                            GENSS_Prod(l,sgs[k][1]+sb)*l[sgs[k][2]+sb]*
                            GENSS_Prod(li,sgs[k][3]+sb)*
                            ResultOfStraightLineProgram(sgs[k][4],l)^-1);
                    else
                        Add(newpres,GENSS_Prod(l,sgs[k][1]+sb)*l[sgs[k][2]+sb]*
                                    GENSS_Prod(li,sgs[k][3]+sb));
                    fi;
                od;
                Info(InfoGenSS,2,"Found presentation for layer ",S!.layer,
                     " using ",Length(newpres)," relators.");
                return SLPOfElms(newpres);
            fi;
        fi;
        nrschr := QuoInt(nrschr*4,2);
        MakeSchreierGens(nrschr);
    od;
  end);



InstallGlobalFunction( VerifyStabilizerChainTC2,
  function( S )
    local Grels,Hrels,MakeSchreierGens,Prels,ace,acenrGrels,
          cosetnrlimitfactor,f,gens,gensi,i,j,k,l,li,max,newpres,
          nrcosets,nrgens,nrschr,o,ords,pres,sb,sgs,slp,subgens,x;
    if S!.stab <> false then
        pres := VerifyStabilizerChainTC2(S!.stab);
        if IsList(pres) then return pres; fi;
    else
        pres := StraightLineProgram([[]],0);
    fi;
    Info(InfoGenSS,1,"Verifying stabilizer chain in layer ",S!.layer);
    # First create a few Schreier generators:
    sgs := [];
    i := 1;
    j := 1;
    o := S!.orb;
    nrgens := Length(o!.gens);
    sb := S!.strongbelow;
    MakeSchreierGens := function(n)
        local sg;
        Info(InfoGenSS,3,"Creating ",n," Schreier generators...");
        while Length(sgs) < n and
              i <= Length(o) do
            sg := GENSS_CreateSchreierGenerator(S,i,j);
            j := j + 1;
            if j > nrgens then
                j := 1;
                i := i + 1;
            fi;
            if sg <> fail then
                Add(sgs,sg);
            fi;
        od;
    end;

    nrschr := S!.opt.NumberSchreierGens;
    MakeSchreierGens(nrschr);
    f := FreeGroup(S!.nrstrong);
    gens := GeneratorsOfGroup(f);
    gensi := List(gens,x->x^-1);
    subgens := gens{[1..S!.strongbelow]};
    Hrels := ResultOfStraightLineProgram(pres,subgens);
    if S!.opt.Projective then
        ords := List([1..nrgens],i->ProjectiveOrder(o!.gens[i]));
    else
        ords := List([1..nrgens],i->Order(o!.gens[i]));
    fi;
    Prels := List([1..nrgens],i->gens[i+sb]^ords[i]);
    Grels := [];
    cosetnrlimitfactor := 100;
    ace := false;
    while true do   # will be left by return eventually
        for k in [Length(Grels)+1..Length(sgs)] do
            Grels[k] := GENSS_Prod(gens,sgs[k][1]+sb) * gens[sgs[k][2]+sb] * 
                        GENSS_Prod(gensi,sgs[k][3]+sb);
            x := GENSS_Prod(o!.gens,sgs[k][1]) * o!.gens[sgs[k][2]] * 
                 GENSS_Prod(o!.gensi,sgs[k][3]);
            if S!.stab <> false then
                slp := SiftGroupElementSLP(S!.stab,x);
                if not(slp.isone) then
                    if ace <> false then
                        ACEQuit(ace);
                    fi;
                    return [fail,S!.layer];
                fi;
                Grels[k] := Grels[k] / ResultOfStraightLineProgram(slp.slp,
                                                                   subgens);
                sgs[k][4] := slp.slp;
            else
                if not(IsOne(x)) then
                    return [fail,S!.layer];
                fi;
                sgs[k][4] := false;
            fi;
        od;
        Info(InfoGenSS,2,"Doing ACE coset enumeration with limit ",
             cosetnrlimitfactor*Length(o)," and ",Length(Hrels),
             "+",Length(Prels),"+",Length(Grels)," relations...");
        if ace = false then
            ace := ACEStart(gens,Concatenation(Hrels,Prels,Grels),subgens:
                            max := cosetnrlimitfactor * Length(o) );
        else
            ACEAddRelators(ace,Grels{[acenrGrels+1..Length(Grels)]});
        fi;
        acenrGrels := Length(Grels);
        if not(IsCompleteACECosetTable(ace))  then   # did not close!
            #cosetnrlimitfactor := QuoInt(cosetnrlimitfactor*3,2);
            Info(InfoGenSS,2,"Coset enumeration did not finish!");
            if nrschr > Length(sgs) # or
               # nrschr > S!.opt.MaxNumberSchreierGens 
               then   # we are done!
                # Something is wrong!
                Error();
                ACEQuit(ace);
                return [fail, S!.layer];
            fi;
        else
            nrcosets := ACEStats(ace).index;
            Info(InfoGenSS,2,"Coset enumeration found ",nrcosets," cosets.");
            if nrcosets = Length(o) then
                ACEQuit(ace);
                # Verification is OK, now build a presentation:
                l := GeneratorsWithMemory(
                       ListWithIdenticalEntries(S!.nrstrong,()));
                li := List(l,x->x^-1);
                newpres := ResultOfStraightLineProgram(pres,
                                   l{[1..S!.strongbelow]});
                for k in [1..nrgens] do
                    Add(newpres,l[k+sb]^ords[k]);
                od;
                for k in [1..Length(sgs)] do
                    if sgs[k][4] <> false then
                        Add(newpres,
                            GENSS_Prod(l,sgs[k][1]+sb)*l[sgs[k][2]+sb]*
                            GENSS_Prod(li,sgs[k][3]+sb)*
                            ResultOfStraightLineProgram(sgs[k][4],l)^-1);
                    else
                        Add(newpres,GENSS_Prod(l,sgs[k][1]+sb)*l[sgs[k][2]+sb]*
                                    GENSS_Prod(li,sgs[k][3]+sb));
                    fi;
                od;
                Info(InfoGenSS,2,"Found presentation for layer ",S!.layer,
                     " using ",Length(newpres)," relators.");
                return SLPOfElms(newpres);
            fi;
        fi;
        # nrschr := nrschr + S!.opt.NumberSchreierGens;
        nrschr := QuoInt(nrschr*4,3);
        MakeSchreierGens(nrschr);
    od;
  end);

VerifyStabilizerChainTC3 := 
  function( S )
    local Grels,Hrels,MakeSchreierGens,Prels,TestGenerationSGens,ace,
          acenrGrels,cosetnrlimitfactor,dr,f,gens,gensi,i,j,k,l,li,max,
          newpres,nrcosets,nrgens,nrschr,o,ok,ords,pres,sb,sgs,slp,subgens,x;
    if S!.stab <> false then
        pres := VerifyStabilizerChainTC3(S!.stab);
        if IsList(pres) then return pres; fi;
    else
        pres := StraightLineProgram([[]],0);
    fi;
    Info(InfoGenSS,1,"Verifying stabilizer chain in layer ",S!.layer);
    # First create a few Schreier generators:
    sgs := [];
    i := 1;
    j := 1;
    o := S!.orb;
    nrgens := Length(o!.gens);
    sb := S!.strongbelow;
    MakeSchreierGens := function(n)
        local sg;
        Info(InfoGenSS,3,"Creating ",n," Schreier generators...");
        while Length(sgs) < n and
              i <= Length(o) do
            sg := GENSS_CreateSchreierGenerator(S,i,j);
            j := j + 1;
            if j > nrgens then
                j := 1;
                i := i + 1;
            fi;
            if sg <> fail then
                Add(sgs,sg);
            fi;
        od;
    end;
    TestGenerationSGens := function()
        local SS,g,k,l;
        l := [];
        for k in [1..Length(sgs)] do
            Add(l,GENSS_Prod(o!.gens,sgs[k][1]) * o!.gens[sgs[k][2]] * 
                  GENSS_Prod(o!.gensi,sgs[k][3]));
        od;
        g := Group(l);
        SS := StabilizerChain(g,rec( Base := S!.stab, 
                                     ErrorBound := 1/3,
                                     Projective := S!.opt.Projective ));
        Info( InfoGenSS,3,"Schreier gen test gave ",Size(SS)," of ",
              Size(S!.stab) );
        return Size(SS) = Size(S!.stab);
    end;        
            
    nrschr := S!.opt.NumberSchreierGens;
    while true do;
        MakeSchreierGens(nrschr);
        ok := S!.stab = false or TestGenerationSGens();
        if ok then break; fi;
        if i > Length(o) then
            Error( "this cannot have happenend" );
        fi;
        nrschr := QuoInt(nrschr*4,3);
    od;
    # We are now pretty sure that the Schreier generators generate the
    # stabilizer, just make some more and then start coset enumeration:
    nrschr := nrschr * 2;
    MakeSchreierGens(nrschr);

    f := FreeGroup(S!.nrstrong);
    gens := GeneratorsOfGroup(f);
    gensi := List(gens,x->x^-1);
    subgens := gens{[1..S!.strongbelow]};
    Hrels := ResultOfStraightLineProgram(pres,subgens);
    if S!.opt.Projective then
        ords := List([1..nrgens],i->ProjectiveOrder(o!.gens[i]));
    else
        ords := List([1..nrgens],i->Order(o!.gens[i]));
    fi;
    Prels := List([1..nrgens],i->gens[i+sb]^ords[i]);
    Grels := [];
    cosetnrlimitfactor := 4;
    ace := false;
    # Xrels := [];
    while true do   # will be left by return eventually
        for k in [Length(Grels)+1..Length(sgs)] do
            Grels[k] := GENSS_Prod(gens,sgs[k][1]+sb) * gens[sgs[k][2]+sb] * 
                        GENSS_Prod(gensi,sgs[k][3]+sb);
            x := GENSS_Prod(o!.gens,sgs[k][1]) * o!.gens[sgs[k][2]] * 
                 GENSS_Prod(o!.gensi,sgs[k][3]);
            if S!.stab <> false then
                slp := SiftGroupElementSLP(S!.stab,x);
                if not(slp.isone) then
                    if ace <> false then
                        ACEQuit(ace);
                    fi;
                    return [fail,S!.layer];
                fi;
                Grels[k] := Grels[k] / ResultOfStraightLineProgram(slp.slp,
                                                                   subgens);
                sgs[k][4] := slp.slp;
            else
                if not(IsOne(x)) then
                    return [fail,S!.layer];
                fi;
                sgs[k][4] := false;
            fi;
        od;
        Info(InfoGenSS,2,"Doing ACE coset enumeration with limit ",
             cosetnrlimitfactor*Length(o)," and ",Length(Hrels),
             "+",Length(Prels),"+",Length(Grels)," relations...");
        if ace = false then
            ace := ACEStart(gens,Concatenation(Hrels,Prels,Grels),subgens:
                            max := cosetnrlimitfactor * Length(o), hlt );
            dr := ACEDataRecord(ace);
            nrcosets := dr.stats.index;
        else
            # k := Random([Length(o)+1..nrcosts]).
            if Length(Grels) > acenrGrels then
                ACEAddRelators(ace,Grels{[acenrGrels+1..Length(Grels)]}:
                               max := cosetnrlimitfactor * Length(o), hlt );
            else
                ACEContinue(ace : max := cosetnrlimitfactor * Length(o));
            fi;
            dr := ACEDataRecord(ace);
            nrcosets := dr.stats.index;
        fi;
        acenrGrels := Length(Grels);
        if not(IsCompleteACECosetTable(ace))  then   # did not close!
            cosetnrlimitfactor := QuoInt(cosetnrlimitfactor*3,2);
            Info(InfoGenSS,2,"Coset enumeration did not finish!");
            #if nrschr > Length(sgs) # or
            #   # nrschr > S!.opt.MaxNumberSchreierGens 
            #   then   # we are done!
            #    # Something is wrong!
            #    Error("wrong1");
            #    ACEQuit(ace);
            #    return [fail, S!.layer];
            #fi;
        else
            Info(InfoGenSS,2,"Coset enumeration found ",nrcosets," cosets.");
            if nrcosets = Length(o) then
                ACEQuit(ace);
                # Verification is OK, now build a presentation:
                l := GeneratorsWithMemory(
                       ListWithIdenticalEntries(S!.nrstrong,()));
                li := List(l,x->x^-1);
                newpres := ResultOfStraightLineProgram(pres,
                                   l{[1..S!.strongbelow]});
                for k in [1..nrgens] do
                    Add(newpres,l[k+sb]^ords[k]);
                od;
                for k in [1..Length(sgs)] do
                    if sgs[k][4] <> false then
                        Add(newpres,
                            GENSS_Prod(l,sgs[k][1]+sb)*l[sgs[k][2]+sb]*
                            GENSS_Prod(li,sgs[k][3]+sb)*
                            ResultOfStraightLineProgram(sgs[k][4],l)^-1);
                    else
                        Add(newpres,GENSS_Prod(l,sgs[k][1]+sb)*l[sgs[k][2]+sb]*
                                    GENSS_Prod(li,sgs[k][3]+sb));
                    fi;
                od;
                Info(InfoGenSS,2,"Found presentation for layer ",S!.layer,
                     " using ",Length(newpres)," relators.");
                return SLPOfElms(newpres);
            fi;
        fi;
        # nrschr := nrschr + S!.opt.NumberSchreierGens;
        nrschr := QuoInt(nrschr*4,3);
        MakeSchreierGens(nrschr);
    od;
  end;


VerifyStabilizerChainTC4 := 
  function( S )
    local Grels,Hrels,MakeSchreierGen,Prels,ace,cosetlimit,done,dr,f,
          gens,gensi,guck1,guck2,hlt,i,j,k,l,li,max,newpres,nrcosets,
          nrgens,o,ords,pres,sb,sg,sgs,slp,st,subgens,x,y,y1,y2;

    if S!.stab <> false then
        pres := VerifyStabilizerChainTC4(S!.stab);
        if IsList(pres) then return pres; fi;
    else
        pres := StraightLineProgram([[]],0);
    fi;
    Info(InfoGenSS,1,"Verifying stabilizer chain in layer ",S!.layer);

    # The following are global to "MakeSchreierGen":
    i := 1;
    j := 1;
    MakeSchreierGen := function()
        local sg;
        Info(InfoGenSS,4,"Creating Schreier generator... i=",i," j=",j);
        while i <= Length(o) do
            sg := GENSS_CreateSchreierGenerator(S,i,j);
            j := j + 1;
            if j > nrgens then
                j := 1;
                i := i + 1;
            fi;
            if sg <> fail then return sg; fi;
        od;
        return fail;
    end;

    o := S!.orb;
    nrgens := Length(o!.gens);
    sb := S!.strongbelow;
            
    if S!.nrstrong > 26 then
        f := FreeGroup(List([1..S!.nrstrong],String));
    else
        gens := [];
        st := "a";
        for k in [1..S!.nrstrong] do
            Add(gens,ShallowCopy(st));
            st[1] := CHAR_INT(INT_CHAR(st[1])+1);
        od;
        f := FreeGroup(gens);
    fi;
    gens := GeneratorsOfGroup(f);
    gensi := List(gens,x->x^-1);
    subgens := gens{[1..S!.strongbelow]};
    Hrels := ResultOfStraightLineProgram(pres,subgens);
    if S!.opt.Projective then
        ords := List([1..nrgens],i->ProjectiveOrder(o!.gens[i]));
    else
        ords := List([1..nrgens],i->Order(o!.gens[i]));
    fi;
    Prels := List([1..nrgens],i->gens[i+sb]^ords[i]);
    sgs := [];
    Grels := [];

    # Now start up a coset enumeration:
    cosetlimit := QuoInt(5 * Length(o),4);
    Info(InfoGenSS,2,"Starting ACE coset enumeration with limit ",
         cosetlimit," and ",Length(Hrels),
         "+",Length(Prels),"+",Length(Grels)," relations...");
    ace := ACEStart(gens,Concatenation(Hrels,Prels,Grels),subgens:
                    max := cosetlimit, hlt := true );
    done := IsCompleteACECosetTable(ace);
    if done then
        dr := ACEDataRecord(ace);
        nrcosets := dr.stats.index;
    fi;
        
    while true do   # will be left by return eventually
        if done then
            Info(InfoGenSS,2,"Coset enumeration found ",nrcosets," cosets.");
            if nrcosets = Length(o) then
                ACEQuit(ace);
                # Verification is OK, now build a presentation:
                l := GeneratorsWithMemory(
                       ListWithIdenticalEntries(S!.nrstrong,()));
                li := List(l,x->x^-1);
                newpres := ResultOfStraightLineProgram(pres,
                                   l{[1..S!.strongbelow]});
                for k in [1..nrgens] do
                    Add(newpres,l[k+sb]^ords[k]);
                od;
                for k in [1..Length(sgs)] do
                    if sgs[k][4] <> false then
                        Add(newpres,
                            GENSS_Prod(l,sgs[k][1]+sb)*l[sgs[k][2]+sb]*
                            GENSS_Prod(li,sgs[k][3]+sb)*
                            ResultOfStraightLineProgram(sgs[k][4],l)^-1);
                    else
                        Add(newpres,GENSS_Prod(l,sgs[k][1]+sb)*l[sgs[k][2]+sb]*
                                    GENSS_Prod(li,sgs[k][3]+sb));
                    fi;
                od;
                Info(InfoGenSS,2,"Found presentation for layer ",S!.layer,
                     " using ",Length(newpres)," relators.");
                return SLPOfElms(newpres);
            elif nrcosets < Length(o) then
                Error("This cannot possible have happened!");
            else   # nrcosets > Length(o)
                Info(InfoGenSS,2,"Too many cosets, we must have forgotten ",
                     "another relation!");
                done := false;
            fi;
        fi;
        if not(done) then
            while true do    # will be left by break
                sg := MakeSchreierGen();
                if sg = fail then
                    Error("Something wrong, have processed all Schreier gens");
                fi;
                # Sift residue:
                x := GENSS_Prod(o!.gens,sg[1]) * o!.gens[sg[2]] * 
                     GENSS_Prod(o!.gensi,sg[3]);
                if not(IsOne(x)) then
                    if S!.stab <> false then
                        slp := SiftGroupElementSLP(S!.stab,x);
                        if not(slp.isone) then
                            ACEQuit(ace);
                            return [fail,S!.layer,x];
                        fi;
                    fi;
                    sg[4] := slp.slp;
                else
                    sg[4] := false;
                fi;
                y1 := GENSS_Prod(gens,sg[1]+sb) * gens[sg[2]+sb];
                y2 := GENSS_Prod(gens,Reversed(sg[3]+sb));
                # Now check with ACE:
                guck1 := ACETraceWord(ace,1,y1);
                guck2 := ACETraceWord(ace,1,y2);
                if guck1 = fail or guck2 = fail or guck1 <> guck2 then
                    y := y1/y2;
                    if sg[4] <> false then
                        y := y / ResultOfStraightLineProgram(sg[4],subgens);
                    fi;
                    Add(sgs,sg);
                    Add(Grels,y);
                    break;
                fi;
                # Otherwise go to next Schreier generator.
            od;
            Info(InfoGenSS,2,"Redoing ACE coset enumeration with limit ",
                 cosetlimit," and ",Length(Hrels),
                 "+",Length(Prels),"+",Length(Grels)," relations...");
            ACEAddRelators(ace,[y]);
            done := IsCompleteACECosetTable(ace);
            if done then
                dr := ACEDataRecord(ace);
                nrcosets := dr.stats.index;
            fi;
        fi;
    od;
  end;

GENSS_CosetTableFromGensAndRelsInit :=
function(fgens,grels,fsgens,limit)
    local   next,  prev,            # next and previous coset on lists
            firstFree,  lastFree,   # first and last free coset
            firstDef,   lastDef,    # first and last defined coset
            table,                  # columns in the table for gens
            rels,                   # representatives of the relators
            relsGen,                # relators sorted by start generator
            subgroup,               # rows for the subgroup gens
            i, gen, inv,            # loop variables for generator
            g,                      # loop variable for generator col
            rel,                    # loop variables for relation
            p, p1, p2,              # generator position numbers
            app,                    # arguments list for 'MakeConsequences'
            j,                      # integer variable
            length, length2,        # length of relator (times 2)
            cols,
            nums,
            l,
            nrdef,                  # number of defined cosets
            nrmax,                  # maximal value of the above
            nrdel,                  # number of deleted cosets
            nrinf;                  # number for next information message

    # give some information
    Info( InfoFpGroup, 3, "    defined deleted alive   maximal");
    nrdef := 1;
    nrmax := 1;
    nrdel := 0;
    nrinf := 1000;

    # define one coset (1)
    firstDef  := 1;  lastDef  := 1;
    firstFree := 2;  lastFree := limit;

    # make the lists that link together all the cosets
    next := [ 2 .. limit + 1 ];  next[1] := 0;  next[limit] := 0;
    prev := [ 0 .. limit - 1 ];  prev[2] := 0;

    # compute the representatives for the relators
    rels := RelatorRepresentatives( grels );

    # make the columns for the generators
    table := [];
    for gen  in fgens  do
        g := ListWithIdenticalEntries( limit, 0 );
        Add( table, g );
        if not ( gen^2 in rels or gen^-2 in rels ) then
            g := ListWithIdenticalEntries( limit, 0 );
        fi;
        Add( table, g );
    od;

    # make the rows for the relators and distribute over relsGen
    relsGen := RelsSortedByStartGen( fgens, rels, table, true );

    # make the rows for the subgroup generators
    subgroup := [];
    for rel  in fsgens  do
      #T this code should use ExtRepOfObj -- its faster
      # cope with SLP elms
      if IsStraightLineProgElm(rel) then
        rel:=EvalStraightLineProgElm(rel);
      fi;
        length := Length( rel );
        length2 := 2 * length;
        nums := [ ]; nums[length2] := 0;
        cols := [ ]; cols[length2] := 0;

        # compute the lists.
        i := 0;  j := 0;
        while i < length do
            i := i + 1;  j := j + 2;
            gen := Subword( rel, i, i );
            p := Position( fgens, gen );
            if p = fail then
                p := Position( fgens, gen^-1 );
                p1 := 2 * p;
                p2 := 2 * p - 1;
            else
                p1 := 2 * p - 1;
                p2 := 2 * p;
            fi;
            nums[j]   := p1;  cols[j]   := table[p1];
            nums[j-1] := p2;  cols[j-1] := table[p2];
        od;
        Add( subgroup, [ nums, cols ] );
    od;

    # make the structure that is passed to 'MakeConsequences'
    app := [ table, next, prev, relsGen, subgroup, 
             firstFree, lastFree, firstDef, lastDef, 0, 0 ];

    # we do not want minimal gaps to be marked in the coset table
    app[12] := 0;

    return rec( nrdef := nrdef, nrmax := nrmax, nrdel := nrdel, nrinf := nrinf,
                app := app, closed := false, table := table, limit := limit,
                fgens := fgens );
end;

GENSS_DoToddCoxeter := function( r )
    local app,fgens,firstDef,firstFree,gen,i,inv,lastDef,lastFree,next,nrdef,
          nrdel,nrinf,nrmax,prev,relsGen,subgroup,table;

    fgens := r.fgens;
    nrdef := r.nrdef;
    nrmax := r.nrmax;
    nrdel := r.nrdel;
    nrinf := r.nrinf;
    app := r.app;
    table := app[1];
    next := app[2];
    prev := app[3];
    relsGen := app[4];
    subgroup := app[5];
    firstFree := app[6];
    lastFree := app[7];
    firstDef := app[8];
    lastDef := app[9];

    # run over all the cosets:
    while firstDef <> 0  do

        # run through all the rows and look for undefined entries
        for i  in [ 1 .. Length( table ) ]  do
            gen := table[i];

            if gen[firstDef] <= 0  then

                inv := table[i + 2*(i mod 2) - 1];

                # if necessary expand the table
                if firstFree = 0  then
                    app[6] := firstFree;
                    app[7] := lastFree;
                    app[8] := firstDef;
                    app[9] := lastDef;
                    r.nrdef := nrdef;
                    r.nrmax := nrmax;
                    r.nrdel := nrdel;
                    r.nrinf := nrinf;
                    return fail;
                fi;

                # update the debugging information
                nrdef := nrdef + 1;
                if nrmax <= firstFree  then
                    nrmax := firstFree;
                fi;

                # define a new coset
                gen[firstDef]   := firstFree;
                inv[firstFree]  := firstDef;
                next[lastDef]   := firstFree;
                prev[firstFree] := lastDef;
                lastDef         := firstFree;
                firstFree       := next[firstFree];
                next[lastDef]   := 0;

                # set up the deduction queue and run over it until it's empty
                app[6] := firstFree;
                app[7] := lastFree;
                app[8] := firstDef;
                app[9] := lastDef;
                app[10] := i;
                app[11] := firstDef;
                nrdel := nrdel + MakeConsequences( app );
                firstFree := app[6];
                lastFree := app[7];
                firstDef := app[8];
                lastDef  := app[9];

                # give some information
                if nrinf <= nrdef+nrdel then
                    Info( InfoFpGroup, 3, "\t", nrdef, "\t", nrinf-nrdef,
                          "\t", 2*nrdef-nrinf, "\t", nrmax );
                    nrinf := ( Int(nrdef+nrdel)/1000 + 1 ) * 1000;
                fi;

            fi;
        od;

        firstDef := next[firstDef];
    od;

    Info( InfoFpGroup, 2, "\t", nrdef, "\t", nrdel, "\t", nrdef-nrdel, "\t",
          nrmax );

    # separate pairs of identical table columns.                  
    for i in [ 1 .. Length( fgens ) ] do                             
        if IsIdenticalObj( table[2*i-1], table[2*i] ) then
            table[2*i] := StructuralCopy( table[2*i-1] );       
        fi;                                                                
    od;

    # standardize the table
    StandardizeTable( table );

    # return the success result
    r.closed := true;
    return true;
end;

GENSS_ToddCoxeterExpandLimit := function(r,newlimit)
    local app,g,l,limit,next,prev,table;
    if newlimit <= r.limit then return r; fi;
    app := r.app;
    table := app[1];
    next := app[2];
    prev := app[3];
    limit := r.limit;
    next[newlimit] := 0;
    prev[newlimit] := newlimit-1;
    for g  in table  do g[newlimit] := 0;  od;
    for l  in [ limit+2 .. newlimit-1 ]  do
        next[l] := l+1;
        prev[l] := l-1;
        for g  in table  do g[l] := 0;  od;
    od;
    next[limit+1] := limit+2;
    prev[limit+1] := 0;
    for g  in table  do g[limit+1] := 0;  od;
    app[6] := limit+1;   # this is firstFree
    r.limit := newlimit;
    app[7] := newlimit;     # this is lastFree
end;

GENSS_TCAddRelators := function(r,newrels)
  local i,relsGen;
  newrels := RelatorRepresentatives(newrels);
  newrels := RelsSortedByStartGen( r.fgens, newrels, r.table, true );
  relsGen := r.app[4];
  for i in [1..Length(relsGen)] do
      Append(relsGen[i],newrels[i]);
  od;
end;


VerifyStabilizerChainTC5 := 
  function( S )
    local FindTwoWords,Grels,Hrels,Prels,allgens,cosetlimit,done,el,f,
          gens,gensi,hom,k,l,li,newpres,newrel,nrcosets,nrgens,o,opgens,
          ords,pres,r,rels,sb,stronggens,strongi,subgens,tc,words;

    if S!.stab <> false then
        pres := VerifyStabilizerChainTC5(S!.stab);
        if IsList(pres) then return pres; fi;
    else
        pres := StraightLineProgram([[]],0);
    fi;
    Info(InfoGenSS,1,"Verifying stabilizer chain in layer ",S!.layer);

    o := S!.orb;
    nrgens := Length(o!.gens);
    sb := S!.strongbelow;
            
    f := FreeGroup(sb+nrgens);
    gens := GeneratorsOfGroup(f);
    gensi := List(gens,x->x^-1);
    allgens := 0*[1..2*Length(gens)];
    allgens{[1,3..2*Length(gens)-1]} := gens;
    allgens{[2,4..2*Length(gens)]} := gensi;
    subgens := gens{[1..S!.strongbelow]};
    Hrels := ResultOfStraightLineProgram(pres,subgens);
    if S!.opt.Projective then
        ords := List([1..nrgens],i->ProjectiveOrder(o!.gens[i]));
    else
        ords := List([1..nrgens],i->Order(o!.gens[i]));
    fi;
    Prels := List([1..nrgens],i->gens[i+sb]^ords[i]);
    Grels := [];
    stronggens := StrongGenerators(S);
    strongi := List(stronggens,x->x^-1);
    opgens := 0*[1..2*Length(stronggens)];
    opgens{[1,3..2*Length(stronggens)-1]} := stronggens;
    opgens{[2,4..2*Length(stronggens)]} := strongi;

    # Now start up a coset enumeration:
    cosetlimit := Maximum(QuoInt(7 * Length(o),6),Length(o)+5);
    Info(InfoGenSS,2,"Starting coset enumeration with limit ",
         cosetlimit," and ",Length(Hrels),
         "+",Length(Prels),"+",Length(Grels)," relations...");
    rels := Concatenation(Hrels,Prels);
    tc := GENSS_CosetTableFromGensAndRelsInit(gens,rels,subgens,cosetlimit);
    done := GENSS_DoToddCoxeter(tc);
        
    FindTwoWords := function(o,opgens,table)
        local TraceWord,cosets,cosetsrevtab,i,j,new,nrcosets,pt,pts,
              ptsrev,schgen,schpt,w1,w2,x,y;
        nrcosets := Length(table[1]);
        cosets := [1];
        cosetsrevtab := 0*[1..nrcosets];
        cosetsrevtab[1] := 1;
        pts := [o[1]];
        ptsrev := 0*[1..Length(o)];
        ptsrev[1] := 1;
        schpt := [fail];    # the Schreier tree
        schgen := [fail];
        i := 1;
        TraceWord := function(pos)
            local w;
            w := [];
            while pos > 1 do
                Add(w,schgen[pos]);
                pos := schpt[pos];
            od;
            return Reversed(w);
        end;
        while i <= Length(cosets) do
            for j in [1..Length(opgens)] do
                x := table[j][cosets[i]];
                if x <> 0 then   # image is defined:
                    if cosetsrevtab[x] = 0 then   # not visited
                        Add(cosets,x);
                        new := Length(cosets);
                        cosetsrevtab[x] := new;
                        schpt[new] := i;
                        schgen[new] := j;
                        pt := o!.op(pts[i],opgens[j]);
                        y := Position(o,pt);
                        if ptsrev[y] = 0 then
                            Add(pts,pt);
                            ptsrev[y] := new;
                        else
                            # We have reached a new coset by a word that
                            # maps the starting point of the orbit to the 
                            # same point as the one of another coset!
                            w1 := TraceWord(ptsrev[y]);
                            w2 := TraceWord(new);
                            return [w1,w2];
                        fi;
                    fi;
                fi;
            od;
            i := i + 1;
        od;
        Error("Bad, this should never have been reached!");
        return fail;
    end;

    while true do   # will be left by return eventually
        if done = true then
            nrcosets := Length(tc.table[1]);
            Info(InfoGenSS,2,"Coset enumeration found ",nrcosets," cosets.");
            if nrcosets = Length(o) then
                # Verification is OK, now build a presentation:
                l := GeneratorsWithMemory(
                       ListWithIdenticalEntries(S!.nrstrong,()));
                newpres := ResultOfStraightLineProgram(pres,
                                   l{[1..S!.strongbelow]});
                for k in [1..nrgens] do
                    Add(newpres,l[k+sb]^ords[k]);
                od;
                hom := GroupHomomorphismByImagesNC(f,Group(l),gens,l);
                for k in Grels do
                    Add(newpres,ImageElm(hom,k));
                od;
                Info(InfoGenSS,2,"Found presentation for layer ",S!.layer,
                     " using ",Length(newpres)," relators.");
                return SLPOfElms(newpres);
            elif nrcosets < Length(o) then
                Error("This cannot possibly have happened!");
                return [fail,S!.layer];
            else   # nrcosets > Length(o)
                Info(InfoGenSS,2,"Too many cosets, we must have forgotten ",
                     "another relation!");
                Error("This cannot possibly have happened2!");
                return [fail,S!.layer];
            fi;
        fi;
        # Now we have to find another relation, we do a breadth-first
        # search through the already defined cosets to find two cosets
        # that are still different but ought to be equal because the
        # corresponding orbit points are equal:
        words := FindTwoWords(o,opgens,tc.table);
        el := EvaluateWord(opgens,words[1])/EvaluateWord(opgens,words[2]);
        r := SiftGroupElementSLP(S,el);
        if not(r.isone) then
            # Error, we found a new stabilizer element!
            return [fail,S!.layer,el];
        fi;
        newrel := [(EvaluateWord(allgens,words[1])
                    /EvaluateWord(allgens,words[2]))
                   / ResultOfStraightLineProgram(r.slp,gens)];
        GENSS_TCAddRelators(tc,newrel);
        Add(Grels,newrel[1]);
        if Length(words[1]) > 0 then
            tc.app[10] := words[1][1];
        else
            # More difficult:
            words := ExtRepOfObj(newrel[1]);
            if words[2] > 0 then
                tc.app[10] := 2*words[1]-1;
            else
                tc.app[10] := 2*words[1];
            fi;
        fi;
        #Print("<\c");
        #for i in [1..tc.limit] do
        #    for j in [1..Length(allgens)] do
        #        if tc.table[j][i] <> 0 then
        #            tc.app[11] := i;
        #            tc.app[10] := j;
        #            tc.nrdel := tc.nrdel + MakeConsequences( tc.app );
        #        fi;
        #    od;
        #od;
        tc.app[11] := 1;
        tc.nrdel := tc.nrdel + MakeConsequences( tc.app );
        #Print("-\c");
        done := GENSS_DoToddCoxeter(tc);
        #Print(">\c");
        if Length(Grels) mod 100 = 0 then
            #Print("\n");
            Info(InfoGenSS,2,"Currently using ",Length(Hrels),"+",
                 Length(Prels),"+",Length(Grels)," relations.");
        fi;
    od;
  end;







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
 
InstallGlobalFunction( GENSS_FindGensStabilizer,
  function( S, pt, opt )
    # Assumes that S is a correct stabilizer chain for the group generated by
    # its gens (in S!.orb!.gens) and that pt lies in the first orbit S!.orb.
    # Finds an SLP to express generators of the stabilizer of pt in
    # the group in terms of the gens.
    # Assumes that the group and the first stabilizer in the chain are 
    # non-trivial.
    local gens,g,pos,mgens,word,cosetrep,invgens,pr,stabgens,size,stabsize,
          randel,newpt,ptpos,ptcosetrep,el,SS,stab,i;
    gens := S!.orb!.gens;
    g := GroupWithGenerators(gens);
    SetSize(g,Size(S));
    pos := Position(S!.orb,pt);
    if pos = fail then
        Error("point pt must lie in first orbit");
        return fail;
    fi;
    if not(IsBound(opt.scramble)) then opt.scramble := 30; fi;
    if not(IsBound(opt.scramblefactor)) then opt.scramblefactors := 0; fi;
    if not(IsBound(opt.maxdepth)) then opt.maxdepth := 200; fi;
    if not(IsBound(opt.addslots)) then 
        opt.addslots := Maximum(0,5-Length(gens)); 
    fi;
    # Give the gens a (new) memory:
    if IsObjWithMemory(gens[1]) then
        mgens := GeneratorsWithMemory(List(gens,x->x!.el));
    else
        mgens := GeneratorsWithMemory(gens);
    fi;
    # FIXME: use GENSS_Prod
    if pos = 1 then
        cosetrep := mgens[1]^0;
    else
        word := TraceSchreierTreeForward(S!.orb,pos);
        cosetrep := Product(mgens{word});
    fi;
    invgens := List(mgens,x->x^-1);

    # Set up the product replacer:
    pr := ProductReplacer( mgens, opt );
    stabgens := [];   # here we collect the stabilizer generators
    size := 1;
    stabsize := Size(S!.stab);   # this has to exist!
    repeat
        randel := Next(pr);
        newpt := S!.orb!.op(pt,randel!.el);
        ptpos := Position(S!.orb,newpt);
        word := TraceSchreierTreeBack(S!.orb,ptpos);
        if Length(word) > 0 then
            ptcosetrep := Product(invgens{word});
            el := randel * ptcosetrep;
        else
            el := randel;
        fi;
        if pos <> 1 then
            el := el * cosetrep;
        fi;
        # now el is an element in the stabilizer of pt 
        if IsOne(el) then
            Info(InfoGenSS,3,"Found trivial stabilizer element.");
        else
            Add(stabgens,el);
            stab := GroupWithGenerators(stabgens);
            SS := StabilizerChain(stab,rec( Base := S ));
            Info(InfoGenSS,1,"Have group size ",Size(SS)," (of ",
                 stabsize,")");
            if Size(SS) = size then
                Remove(stabgens);
            else
                size := Size(SS);
            fi;
        fi;
    until size = stabsize;
    # Try leaving out the first generators found:
    i := 1;
    repeat
        Info(InfoGenSS,1,"Need ",Length(stabgens)+1-i," generators.");
        i := i + 1;
        if i < Length(stabgens) then
            stab := Group(stabgens{[i..Length(stabgens)]});
            SS := StabilizerChain(stab,rec(Base := S));
            size := Size(SS);
        else
            size := 0;
        fi;
    until size < stabsize;
    return SLPOfElms(stabgens{[i-1..Length(stabgens)]});
  end );

InstallGlobalFunction( GENSS_FindShortGensStabilizerOld,
  function( S, pt )
    # Assumes that S is a correct stabilizer chain for the group generated by
    # its gens (in S!.orb!.gens) and that pt lies in the first orbit S!.orb.
    # Finds an SLP to express generators of the stabilizer of pt in
    # the group in terms of the gens.
    # Assumes that the group and the first stabilizer in the chain are 
    # non-trivial.
    local gens,g,pos,mgens,word,cosetrep,invgens,pr,stabgens,size,stabsize,
          randel,newpt,ptpos,ptcosetrep,el,SS,stab,i,j;
    gens := S!.orb!.gens;
    g := GroupWithGenerators(gens);
    SetSize(g,Size(S));
    pos := Position(S!.orb,pt);
    if pos = fail then
        Error("point pt must lie in first orbit");
        return fail;
    fi;
    # Give the gens a (new) memory:
    if IsObjWithMemory(gens[1]) then
        mgens := GeneratorsWithMemory(List(gens,x->x!.el));
    else
        mgens := GeneratorsWithMemory(gens);
    fi;
    if pos = 1 then
        cosetrep := mgens[1]^0;
    else
        word := TraceSchreierTreeForward(S!.orb,pos);
        cosetrep := Product(mgens{word});
    fi;
    invgens := List(mgens,x->x^-1);

    # Set up the search:
    pr := ShallowCopy(mgens);
    i := 1;   # counts points
    j := 1;   # counts gens
    stabgens := [];   # here we collect the stabilizer generators
    size := 1;
    stabsize := Size(S!.stab);   # this has to exist!
    repeat
        Add(pr,pr[i]*mgens[j]);
        j := j + 1; 
        if j > Length(mgens) then
            j := 1;
            i := i + 1;
        fi;
        randel := pr[Length(pr)];
        newpt := S!.orb!.op(pt,randel!.el);
        ptpos := Position(S!.orb,newpt);
        word := TraceSchreierTreeBack(S!.orb,ptpos);
        if Length(word) > 0 then
            ptcosetrep := Product(invgens{word});
            el := randel * ptcosetrep;
        else
            el := randel;
        fi;
        if pos <> 1 then
            el := el * cosetrep;
        fi;
        # now el is an element in the stabilizer of pt 
        if IsOne(el) then
            Info(InfoGenSS,2,"Found trivial stabilizer element.");
        else
            Add(stabgens,el);
            stab := GroupWithGenerators(stabgens);
            SS := StabilizerChain(stab,rec( Base := S ));
            Info(InfoGenSS,1,"Have group size ",Size(SS)," (of ",
                 stabsize,")");
            if Size(SS) = size then
                Remove(stabgens);
            else
                size := Size(SS);
            fi;
        fi;
    until size = stabsize;
    # Try leaving out the first generators found:
    i := 1;
    repeat
        Info(InfoGenSS,1,"Need ",Length(stabgens)+1-i," generators.");
        i := i + 1;
        if i < Length(stabgens) then
            stab := Group(stabgens{[i..Length(stabgens)]});
            SS := StabilizerChain(stab,rec(Base := S));
            size := Size(SS);
        else
            size := 0;
        fi;
    until size < stabsize;
    return SLPOfElms(stabgens{[i-1..Length(stabgens)]});
  end );

InstallGlobalFunction( GENSS_FindShortGensStabilizer,
  function( gens, pt, op, grpsize, stabsize, S )
    # Assumes that S is a correct stabilizer chain a supergroup of the
    # group g generated by gens. pt and op is an action for g.
    # grpsize is the Size of g and stabsize is the size of the stabilizer.
    # Finds an SLP to express generators of the stabilizer of pt in
    # g in terms of the gens.
    # Assumes that g and the orbit of pt under g are non-trivial.
    local SS,el,g,i,invgens,j,mgens,newpt,orb,pr,ptcosetrep,ptpos,randel,size,
          stab,stabgens,word;
    g := GroupWithGenerators(gens);
    Info(InfoGenSS,2,"Enumerating orbit of length ",grpsize/stabsize);
    orb := Orb(gens,pt,op,rec(schreier := true,stabsizebound := stabsize));
    Enumerate(orb);
    # Give the gens a (new) memory:
    if IsObjWithMemory(gens[1]) then
        mgens := GeneratorsWithMemory(List(gens,x->x!.el));
    else
        mgens := GeneratorsWithMemory(gens);
    fi;
    invgens := List(mgens,x->x^-1);

    # Set up the search:
    pr := ShallowCopy(mgens);
    i := 1;   # counts points
    j := 1;   # counts gens
    stabgens := [];   # here we collect the stabilizer generators
    size := 1;
    repeat
        Add(pr,pr[i]*mgens[j]);
        j := j + 1; 
        if j > Length(mgens) then
            j := 1;
            i := i + 1;
        fi;
        randel := pr[Length(pr)];
        newpt := op(pt,randel!.el);
        ptpos := Position(orb,newpt);
        word := TraceSchreierTreeBack(orb,ptpos);
        if Length(word) > 0 then
            ptcosetrep := Product(invgens{word});
            el := randel * ptcosetrep;
        else
            el := randel;
        fi;
        # now el is an element in the stabilizer of pt 
        if IsOne(el) then
            Info(InfoGenSS,3,"Found trivial stabilizer element.");
        else
            Add(stabgens,el);
            stab := GroupWithGenerators(stabgens);
            SS := StabilizerChain(stab,rec( Base := S ));
            Info(InfoGenSS,1,"Have group size ",Size(SS)," (of ",
                 stabsize,")");
            if Size(SS) = size then
                Remove(stabgens);
            else
                size := Size(SS);
            fi;
        fi;
    until size >= stabsize;
    if size > stabsize then
        Error("Something is wrong, stabsize limit was too small");
        return fail;
    fi;
    # Try leaving out the first generators found:
    i := 1;
    Info(InfoGenSS,1,"Need ",Length(stabgens)," generators.");
    while i < Length(stabgens) do
       stab := Group(stabgens{Concatenation([1..i-1],[i+1..Length(stabgens)])});
       SS := StabilizerChain(stab,rec(Base := S));
       if Size(SS) < stabsize then
           i := i + 1;   # keep generator
       else
           Remove(stabgens,i);
           Info(InfoGenSS,1,"Need ",Length(stabgens)," generators.");
       fi;
    od;
    return SLPOfElms(stabgens);
  end );

InstallGlobalFunction( SLPChainStabilizerChain,
  function( S, gens )
    # S a stabilizer chain assumed to be correct for the group g generated
    # by generators gens.
    # Returns a list of straight line programs expressing successively
    # the stabilizers in the chain, each in terms of the generators of the
    # previous, the first in terms of gens.
    local l,ll,slp;
    l := [];
    ll := [Size(S)];
    while S!.stab <> false do
        Info(InfoGenSS,1,"Working on group with size ",Size(S),"...");
        slp := GENSS_FindShortGensStabilizer( 
             gens, S!.orb[1], S!.orb!.op, Size(S), Size(S!.stab), S );
        Add(l,slp);
        Add(ll,Size(S!.stab));
        gens := ResultOfStraightLineProgram(slp,gens);
        Info(InfoGenSS,1,"Found SLP with ",
             Length(LinesOfStraightLineProgram(slp))," lines to compute ",
             Length(gens)," generators.");
        S := S!.stab;
    od;
    return rec(slps := l,sizes := ll);
  end );

InstallGlobalFunction( GENSS_MakeIterRecord,
  function( S )
    local SS,Slist,it;
    Slist := [S];
    SS := S!.stab;
    while SS <> false do
        Add(Slist,SS);
        SS := SS!.stab;
    od;
    it := rec( NextIterator := GENSS_GroupNextIterator,
               IsDoneIterator := GENSS_GroupIsDoneIterator,
               ShallowCopy := GENSS_GroupShallowCopy,
               pos := ListWithIdenticalEntries( Length( S!.base ), 1 ),
               Slist := Slist,
               cache := ListWithIdenticalEntries( Length(S!.base), 
                                                  One(S!.orb!.gens[1]) ),
               one := One(S!.orb!.gens[1]),
               lens := List(Slist,x->Length(x!.orb)) );
    return it;
  end );

InstallMethod( GroupIteratorByStabilizerChain, "for a stabilizer chain",
  [ IsStabilizerChain ],
  function( S )
    return IteratorByFunctions(GENSS_MakeIterRecord(S));
  end );

InstallGlobalFunction( GENSS_GroupNextIterator,
  function( iter )
    local done,i,j,res,word;
    if iter!.pos[1] > iter!.lens[1] then return fail; fi;
    res := iter!.cache[Length(iter!.cache)];
    # Now step forwards:
    i := Length(iter!.Slist);
    repeat
        iter!.pos[i] := iter!.pos[i] + 1;
        if iter!.pos[i] > iter!.lens[i] then
            iter!.pos[i] := 1;
            i := i - 1;
            done := false;
        else
            done := true;
        fi;
    until done or i = 0;
    if i = 0 then 
        iter!.pos[1] := iter!.lens[1]+1; 
        return res; 
    fi;   # we are done
    word := TraceSchreierTreeForward(iter!.Slist[i]!.orb,iter!.pos[i]);
    iter!.cache[i] := Product(iter!.Slist[i]!.orb!.gens{word});
    if i > 1 then iter!.cache[i] := iter!.cache[i] * iter!.cache[i-1]; fi;
    for j in [i+1..Length(iter!.Slist)] do
        iter!.cache[j] := iter!.cache[j-1];
    od;
    return res;
  end );

InstallGlobalFunction( GENSS_GroupIsDoneIterator,
  function( iter )
    return iter!.pos[1] > iter!.lens[1];
  end );

InstallGlobalFunction( GENSS_GroupShallowCopy,
  function( iter )
    return GENSS_MakeIterRecord( iter!.Slist[1] );
  end );

InstallGlobalFunction( GENSS_ImageElm,
  function( data, x )
    local r;
    r := SiftGroupElementSLP(data!.Sg,x);
    if not(r.isone) then return fail; fi;
    return ResultOfStraightLineProgram(r.slp,data!.stronggims);
  end );

InstallGlobalFunction( GENSS_PreImagesRepresentative,
  function( data, x )
    local r;
    r := SiftGroupElementSLP(data!.Si,x);
    if not(r.isone) then return fail; fi;
    return ResultOfStraightLineProgram(r.slp,data!.strongipre);
  end );

InstallGlobalFunction( GroupHomomorphismByImagesNCStabilizerChain,
  function( g, h, images, opt )
    local Sg,Si,data,gm,im,slpstrongg,slpstrongi,strongg,stronggims,
          strongi,strongipre;
    gm := GroupWithMemory(g);
    Sg := StabilizerChain(gm,opt);
    strongg := StrongGenerators(Sg);
    ForgetMemory(Sg);
    slpstrongg := SLPOfElms(strongg);
    stronggims := ResultOfStraightLineProgram(slpstrongg,images);
    im := GroupWithMemory(images);
    Si := StabilizerChain(im,opt);
    strongi := StrongGenerators(Si);
    ForgetMemory(Si);
    slpstrongi := SLPOfElms(strongi);
    strongipre := ResultOfStraightLineProgram(slpstrongi,GeneratorsOfGroup(g));
    data := rec( Sg := Sg, stronggims := stronggims,
                 Si := Si, strongipre := strongipre,
                 slpstrongg := slpstrongg, slpstrongi := slpstrongi );
    return GroupHomByFuncWithData( g, h, GENSS_ImageElm, false,
                                   GENSS_PreImagesRepresentative, data );
  end );


