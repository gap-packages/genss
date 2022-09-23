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
GENSS.InitialHashSize := NextPrimeInt(1000);
# Number of points to process before reporting:
GENSS.Report := 30000;
# Number of random elements to consider for the determination of short orbits:
GENSS.ShortOrbitsNrRandoms := 10;
# Absolute limit for orbit length in search for short orbits:
GENSS.ShortOrbitsOrbLimit := 80000;
# First limit for the parallel enumeration of orbits looking for short orbits:
GENSS.ShortOrbitsInitialLimit := 400;
# Absolute limit for single orbit length:
GENSS.OrbitLengthLimit := 10000000;
# Number of points in the previous orbit to consider for the next base point:
GENSS.NumberPrevOrbitPoints := 10;
# Number of (evenly distributed) random generators for the stabilizer:
GENSS.RandomStabGens := 3;
# Product replacement parameters for the stabilizer element generation:
# Now actually used for the generation of all random elements on top
# level. Random elements further down are created on top and then
# sifted.
GENSS.StabGenScramble := 30;
GENSS.StabGenScrambleFactor := 6;
GENSS.StabGenAddSlots := 3;
GENSS.StabGenMaxDepth := 400;
# Number of random elements used for verification,
# note that this is changed by the "random" and "ErrorBound" options!
GENSS.VerifyElements := 10;   # this amounts to an error bound of 1/1024
GENSS.DeterministicVerification := false;
# Number of random elements used to do immediate verification:
GENSS.ImmediateVerificationElements := 3;
# Are we working projectively?
GENSS.Projective := false;
# To find a very short orbit we try two basis vectors and a random vector:
GENSS.VeryShortOrbLimit := 500;
# Never consider more than this number of candidates for short orbits:
GENSS.LimitShortOrbCandidates := 50;
# Do not throw Errors but return fail:
GENSS.FailInsteadOfError := false;
# Number of Schreier generators to create in TC verification:
GENSS.NumberSchreierGens := 20;
# Maximal number of Schreier generators to create in TC verification:
GENSS.MaxNumberSchreierGens := 100;
# By default do 3 rounds of birthday paradox method:
GENSS.TryBirthdayParadox := 3;
# By default do 1 short orbit tries:
GENSS.TryShortOrbit := 1;
# Limit for orbit length during orbit estimation using birthday paradox:
GENSS.OrbitLimitBirthdayParadox := 1000000;
# We immediately take an orbit if its estimate is lower than this limit:
GENSS.OrbitLimitImmediatelyTake := 10000;
# Limit for number of random elements during orbit estimation:
GENSS.NrRandElsBirthdayParadox := 6000;
# Set this to false to allow orbits of length 1:
GENSS.Reduced := true;

# Product replacement parameters for Stab:
GENSS.StabScramble := 10;
GENSS.StabScrambleFactor := 1;
GENSS.StabAddSlots := 1;
GENSS.StabMaxDepth := 400;
# Initial limit for orbit length for Stab computation:
GENSS.StabInitialLimit := 1000;
# Patience to find random elements mapping the orbit into itself:
GENSS.StabInitialPatience := 10;
# Maximal amount of memory used for the orbit:
GENSS.StabOrbitLimit := 1000000;
# If the probability of a wrong stabiliser is smaller than this, do no
# longer try to create stabiliser elements:
GENSS.StabAssumeCompleteLimit := 1/(10^7);
# The following switches off the log used in all orbits for stabilizer
# chains if set to false. Do not do this unless you know what you are
# doing since since could prevent Schreier trees from being shallow.
GENSS.OrbitsWithLog := true;

#############################################################################
# A few helper functions needed elsewhere:
#############################################################################

InstallGlobalFunction( GENSS_CopyDefaultOptions,
  function( defopt, opt )
    local n;
    for n in RecNames(defopt) do
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
  # This implements Murray/O'Brien-like heuristics.
  # It produces new random elements using GENSS_RandomElementFromAbove
  # and stores them in opt.FindBasePointCandidatesData.randpool for
  # later usage.
  function(g,opt,parentS)
    # Needs opt components "ShortOrbitsNrRandoms"
    local l, f, data, x, c, onlydegs, v, vv, w, i, nw, inters, sortfun, 
          wb, j, ww;
    Info(InfoGenSS,4,"Trying Murray/O'Brien heuristics...");
    l := ShallowCopy(GeneratorsOfGroup(g));
    f := DefaultFieldOfMatrixGroup(g);
    data := opt.FindBasePointCandidatesData;
    for i in [1..opt.ShortOrbitsNrRandoms] do
        if parentS = false then
            x := GENSS_RandomElementFromAbove(opt,0);
        else
            x := GENSS_RandomElementFromAbove(parentS,parentS!.layer);
        fi;
        Add(l,x);
        Add(data.randpool,x);
    od;
    ForgetMemory(l);
    c := List(l,x->Set(Factors(CharacteristicPolynomial(x,1):
                               onlydegs := [1..3])));
    v := [];
    for i in [1..Length(l)] do
        for j in [1..Length(c[i])] do
            vv := EmptyPlist(Length(c[i]));
            Add(vv,[NullspaceMat(Value(c[i][j],l[i])),
                    Degree(c[i][j]),
                    WeightVecFFE(CoefficientsOfLaurentPolynomial(c[i][j])[1]),
                    1]);
        od;
        Add(v,vv);
    od;
    Info(InfoGenSS,4,"Have eigenspaces.");
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
  function( g, opt, parentS )
    # Needs opt components:
    #  "ShortOrbitsNrRandoms"  (because it uses GENSS_FindVectorsWithShortOrbit)
    #  "ShortOrbitsOrbLimit"
    #  "ShortOrbitsInitialartLimit"
    local ThrowAwayOrbit,found,gens,hashlen,i,j,limit,newnrorbs,nrorbs,o,wb;

    wb := GENSS_FindVectorsWithShortOrbit(g,opt,parentS);
    if Length(wb) = 0 then return fail; fi;

    # Now we have a list of vectors with (hopefully) short orbits.
    # We start enumerating all those orbits, but first only 50 elements:
    nrorbs := Minimum(Length(wb),64);  # take only the 64 first
    gens := GeneratorsOfGroup(g);
    o := [];
    hashlen := NextPrimeInt(QuoInt(opt.ShortOrbitsOrbLimit,2));
    for i in [1..nrorbs] do
        Add(o,Orb(gens,ShallowCopy(wb[i]),OnLines,
                  rec( treehashsize := hashlen )));
    od;
    limit := opt.ShortOrbitsInitialLimit;
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
    if not(IsOne(s)) then
        s := s^-1;
        el := s * el;
    fi;
    return IsOne( el );
  end );


#############################################################################
# Now to the heart of the method, the Schreier-Sims:
#############################################################################

InstallMethod( FindBasePointCandidates, 
  "for a scalar matrix group over a FF",
  [ IsGroup and IsMatrixGroup and IsFinite, IsRecord, IsInt, IsObject ],
  50,  # highest weight,
  function( grp, opt, i, parentS )
    local F, q, d, gens, v;
    Info( InfoGenSS, 3, "Finding nice base points (scalar)..." );
    F := DefaultFieldOfMatrixGroup(grp);
    q := Size(F);
    d := DimensionOfMatrixGroup(grp);
    gens := GeneratorsOfGroup(grp);
    if IsObjWithMemory(gens[1]) then
        gens := StripMemory(gens);
    fi;
    if ForAny(gens,x->not(GENSS_IsOneProjective(x))) then
        opt.FindBasePointCandidatesData := rec( randpool := [], vecs := [] );
        # This is needed for communication between different methods
        TryNextMethod();
    fi;
    v := ZeroMutable(gens[1][1]);
    v[1] := PrimitiveRoot(F);   # to indicate the field!
    return rec( points := [v], ops := [OnRight], used := 0 );
  end );

InstallMethod( FindBasePointCandidates, 
  "for a matrix group over a FF, very short orbit",
  [ IsGroup and IsMatrixGroup and IsFinite, IsRecord, IsInt, IsObject ],
  40,  # highest weight,
  function( grp, opt, i, parentS )
    local F, q, d, gens, op, v, vv, k, kk, o, cand, j;
    Info( InfoGenSS, 3, "Finding nice base points (very short)..." );
    F := DefaultFieldOfMatrixGroup(grp);
    q := Size(F);
    d := DimensionOfMatrixGroup(grp);
    gens := GeneratorsOfGroup(grp);
    if IsObjWithMemory(gens[1]) then
        gens := StripMemory(gens);
    fi;

    # Try two standard basis vectors and a random vector to find a very 
    # short orbit:
    if q = 2 then
        op := OnPoints;
    else
        op := OnLines;
    fi;
    v := [];
    # Find first standard basis vector that is moved:
    vv := ZeroMutable( gens[1][1] );
    k := 1;
    while k <= d do
        vv[k] := One(F);
        if ForAny(gens,x->op(vv,x) <> vv) then
            Add(v,vv);
            break;
        fi;
        vv[k] := Zero(F);
        k := k + 1;
    od;
    # Find last standard basis vector that is moved:
    vv := ZeroMutable( gens[1][1] );
    kk := d;
    while kk >= k do
        vv[kk] := One(F);
        if ForAny(gens,x->op(vv,x) <> vv) then
            Add(v,vv);
            break;
        fi;
        vv[kk] := Zero(F);
        kk := kk - 1;
    od;
    # Pick a random vector:
    vv := ZeroMutable( gens[1][1] );
    if IsPlistRep(vv) then
        for j in [1..Length(vv)] do
            vv[j] := Random(F);
        od;
    else
        Randomize(vv);
    fi;
    ORB_NormalizeVector(vv);
    Add(v,vv);
    # Now investigate these up to a certain limit:
    for j in [1..Length(v)] do
        o := Orb(gens,v[j],op, 
                 rec(treehashsize := QuoInt(opt.VeryShortOrbLimit,2)+1));
        Enumerate(o,opt.VeryShortOrbLimit);
        if Length(o) > 1 and Length(o) < opt.VeryShortOrbLimit then
            Info( InfoGenSS, 3, "Found orbit of length ",Length(o) );
            cand := rec( points := [v[j]], ops := [op], used := 0 );
            # Note that if we work non-projectively, then the same
            # point will be taken next with OnRight automatically!
            return cand;
        fi;
    od;
    Info( InfoGenSS, 3, "Found no very short orbit up to limit ",
          opt.VeryShortOrbLimit );
    Append(opt.FindBasePointCandidatesData.vecs,v);  # hand on vectors
    TryNextMethod();
  end );

# Method with rank 30 taking some commutators and invariant spaces

InstallMethod( FindBasePointCandidates,
  "for a matrix group over a FF, using birthday paradox method",
  [ IsGroup and IsMatrixGroup and IsFinite, IsRecord, IsInt, IsObject ], 20,
  function( grp, opt, mode, parentS )
    local F, q, d, randels, immorblimit, orblimit, data, op, v, l, c, e, ht, 
          val, x, w, cand, minest, minpos, round, i, j, gens;
    F := DefaultFieldOfMatrixGroup(grp);
    q := Size(F);
    d := DimensionOfMatrixGroup(grp);
    randels := opt.NrRandElsBirthdayParadox;
    immorblimit := opt.OrbitLimitImmediatelyTake;
    orblimit := opt.OrbitLimitBirthdayParadox;

    Info( InfoGenSS, 3, "Finding base points (birthday paradox, limit=",
                        orblimit,", randels=",randels,")..." );
    data := opt.FindBasePointCandidatesData; # this we get from earlier methods
    if q = 2 then
        op := OnPoints;
    else
        op := OnLines;
    fi;
    gens := GeneratorsOfGroup(grp);
    if IsObjWithMemory(gens[1]) then
        gens := StripMemory(gens);
    fi;
    for round in [1..opt.TryBirthdayParadox] do
        v := Set(GENSS_FindVectorsWithShortOrbit(grp,opt,parentS));
        if round = 1 then
            Append(v,data.vecs);   # take previously tried ones as well
        fi;
        v := Filtered(v,vv->ForAny(gens,x-> vv <> op(vv,x)));
        l := Length(v);
        if l = 0 then
            # Find a vector on which grp acts non-trivially.
            # At least one basis vector is guaranteed to be moved,
            # unless grp is diagonal and op is OnLines, in which case
            # they might all be fixed. However, in that case, the
            # vector [1,1,...,1] is moved.
            v := OneMutable(gens[1]); # List of basis vectors
            Add(v, Sum(v)); # add vector [1,1,...,1]
            v := Filtered(v,vv->ForAny(gens,x-> vv <> op(vv,x)));
            l := Length(v);
        fi;
        c := 0*[1..l];    # the number of coincidences
        e := ListWithIdenticalEntries(l,infinity);   # the current estimates
        ht := HTCreate(v[1]*PrimitiveRoot(F),
                       rec(hashlen := NextPrimeInt(l * randels * 4)));
        for i in [1..l] do
            val := HTValue(ht,v[i]);
            if val = fail then
                HTAdd(ht,v[i],[i]);
            else
                AddSet(val,i);
            fi;
        od;
        for i in [1..randels] do
            if parentS = false then
                x := GENSS_RandomElementFromAbove(opt,0);
            else
                x := GENSS_RandomElementFromAbove(parentS,parentS!.layer);
            fi;
            Add(data.randpool,x);
            for j in [1..l] do
                if IsObjWithMemory(x) then
                    w := op(v[j],x!.el);
                else
                    w := op(v[j],x);
                fi;
                val := HTValue(ht,w);
                if val <> fail then   # we know this point!
                    if j in val then    # a coincidence!
                        c[j] := c[j] + 1;
                        e[j] := QuoInt(i^2,2*c[j]);
                        if (c[j] >= 3 and e[j] <= immorblimit) or
                           (c[j] >= 15 and e[j] <= orblimit) then
                             Info( InfoGenSS, 2, "Found orbit with estimated ",
                                   "length ",e[j]," (coinc=",c[j],")" );
                             cand := rec(points := [v[j]], ops := [op], 
                                         used := 0);
                             for i in [1..l] do
                                 if i <> j and c[i] >= 10 and
                                    e[i] <= orblimit then
                                     Add(cand.points,v[i]);
                                     Add(cand.ops,op);
                                 fi;
                             od;
                             if Length(cand.points) > 1 then
                                 Info( InfoGenSS, 2, "Adding ", 
                                       Length(cand.points)-1, " more vectors",
                                       " to candidates.");
                             fi;
                             return cand;
                        fi;
                    else
                        AddSet(val,j);
                    fi;
                else
                    HTAdd(ht,w,[j]);
                fi;
            od;
        od;
        minest := Minimum(e);
        minpos := Position(e,minest);
        Info( InfoGenSS,2,"Birthday #", round, ": no small enough estimate. ",
              "MinEst=",minest," Coinc=",c[j] );
        randels := randels * 2;
        orblimit := orblimit * 4;
    od;
    TryNextMethod();
  end );


InstallMethod( FindBasePointCandidates, 
  "for a matrix group over a FF, original try short orbit method",
  [ IsGroup and IsMatrixGroup and IsFinite, IsRecord, IsInt, IsObject ], 10,
  function( grp, opt, i, parentS )
    local F, d, data, cand, res;
    F := DefaultFieldOfMatrixGroup(grp);
    d := DimensionOfMatrixGroup(grp);
    data := opt.FindBasePointCandidatesData;
    Info( InfoGenSS, 3, "Finding nice base points (TryShortOrbit)..." );

    # Next possibly "TryShortOrbit":
    cand := rec( points := [], used := 0 );
    if IsBound(opt.TryShortOrbit) and opt.TryShortOrbit > 0 then
        repeat
            opt.TryShortOrbit := opt.TryShortOrbit - 1;
            Info(InfoGenSS,1,"Looking for short orbit (",opt.TryShortOrbit,
                 ")...");
            res := GENSS_FindShortOrbit(grp,opt,parentS);
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
    TryNextMethod();
  end );

InstallMethod( FindBasePointCandidates, 
  "for a matrix group over a FF, traditional Murray/O'Brien", 
  [ IsGroup and IsMatrixGroup and IsFinite, IsRecord, IsInt, IsObject ],
  function( grp, opt, i, parentS )
    local F, d, bv, cand, w, v;
    F := DefaultFieldOfMatrixGroup(grp);
    d := DimensionOfMatrixGroup(grp);
    Info( InfoGenSS, 3, "Finding nice base points (Murray/O'Brien)..." );

    # Standard Murray/O'Brien heuristics:
    if i = 0 and 
       ((opt!.Projective = false and Size(F)^d > 300000) or
        (opt!.Projective = true and Size(F)^(d-1) > 300000)) then
        bv := GENSS_FindVectorsWithShortOrbit(grp,opt,parentS);
        if Length(bv) < 3 then
            bv := MutableCopyMat(One(grp));
        fi;
        bv := bv{[1..3]};   # just take 3 of them
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
  [ IsGroup and IsPermGroup, IsRecord, IsInt, IsObject ],
  function( grp, opt, i, parentS )
    local ops,points;
    if i = 0 then
        points := [1..Minimum(20,LargestMovedPoint(grp))];
    else
        points := [1..LargestMovedPoint(grp)];
    fi;
    ops := List([1..Length(points)],x->OnPoints);
    return rec( points := points, ops := ops, used := 0 );
  end );
    
GENSS_HACK := OnPoints;

InstallGlobalFunction( GENSS_OpFunctionMaker, function(op,index)
  local name,s,f;
  name := NAME_FUNC(op);
  if name = "unknown" or name = "GENSS_HACK" then return fail; fi;
  s := Concatenation( "GENSS_HACK := function(x,el) return ",
                      name, "(x,el[",String(index),"]); end;" );
  f := InputTextString(s);
  Read(f);
  return GENSS_HACK;
end);

InstallMethod( FindBasePointCandidates, "for a direct product",
  [ IsGroup, IsRecord, IsInt, IsObject ],
  function( grp, opt, i, parentS )
    local gens,l,cand,j,factgens,fac,cand2,k,op,S2,op2;
    gens := GeneratorsOfGroup(grp);
    if not(ForAll(gens,IsDirectProductElement)) or Length(gens) = 0 then
        TryNextMethod();
    fi;
    l := Length(gens[1]);
    cand := rec( points := [], ops := [] );
    for j in [1..l] do
        factgens := List(gens,x->x[j]);
        fac := Group(factgens);
        S2 := StabilizerChain(fac);
        cand2 := BaseStabilizerChain(S2);
        for k in [1..Length(cand2.ops)] do
            op := cand2.ops[k];
            op2 := GENSS_OpFunctionMaker(op,j);
            if op2 <> fail then
                Add(cand.points,cand2.points[k]);
                Add(cand.ops,op2);
            fi;
        od;
    od;
    cand.used := 0;
    return cand;
  end );


InstallGlobalFunction( GENSS_NextBasePoint, 
  function( gens, cand, opt, S )
    local NotFixedUnderAllGens,i,notfixed;

    if IsBound(opt.FindBasePointCandidatesData) then
        Unbind(opt.FindBasePointCandidatesData);
    fi;

    NotFixedUnderAllGens := function( gens, x, op )
      if IsObjWithMemory(gens[1]) then
        return ForAny( gens, g->op(x,g!.el) <> x );
      else
        return ForAny( gens, g->op(x,g) <> x );
      fi;
    end;

    # S can be false or a stabilizer chain record
    if S <> false and not(opt.StrictlyUseCandidates) then  
        # try points in previous orbit
        for i in [2..Minimum(Length(S!.orb),opt.NumberPrevOrbitPoints)] do
            if NotFixedUnderAllGens(gens,S!.orb[i],S!.orb!.op) then
                Info(InfoGenSS,3,"Taking another point in same orbit.");
                return rec( point := S!.orb[i], op := S!.orb!.op, 
                            cand := cand );
            fi;
        od;
        # Maybe we can take the last base point now acting non-projectively?
        if opt.Projective = false and IsIdenticalObj(S!.orb!.op,OnLines) and 
           NotFixedUnderAllGens(gens,S!.orb[1],OnRight) then
            Info(InfoGenSS,3,"Taking same point in non-projective action.");
            return rec( point := S!.orb[1], op := OnRight, cand := cand );
        fi;
    fi;

    # Now use the candidates we already have:
    repeat
        if cand.used >= Length(cand.points) then
            if IsBound(cand.points2) and Length(cand.points2) > 0 then
                cand := rec( points := cand.points2, ops := cand.ops2,
                             used := 0 );
            else
                cand := FindBasePointCandidates(Group(gens),opt,1,S);
                opt.StrictlyUseCandidates := false;
                opt.Reduced := true;
            fi;
        fi;
        cand.used := cand.used + 1;
        notfixed := NotFixedUnderAllGens(gens,
                          cand.points[cand.used],cand.ops[cand.used]);
        if notfixed = false and opt.Reduced = true then # everything is fixed!
            if IsBound(cand.points2) then
                # Maybe our current stabilizer is too small, we still
                # might want to use these points again:
                Add(cand.points2,cand.points[cand.used]);
                Add(cand.ops2,cand.ops[cand.used]);
            fi;
        fi;
    until opt.Reduced = false or notfixed;
          
    return rec( point := cand.points[cand.used], op := cand.ops[cand.used],
                cand := cand );
  end );


InstallGlobalFunction( GENSS_CreateStabChainRecord,
  function( parentS, gens, size, nextpoint, nextop, cand, opt )
    # parentS can be false or the parent in the stabiliser chain
    local base, layer, stronggens, layergens, nr, orb, S, hashsize;

    Info( InfoGenSS, 4, "Creating new stab chain record..." );

    if parentS = false then
        base := [];
        layer := 1;
        stronggens := ShallowCopy(gens);
        layergens := [1..Length(gens)];   # indices in stronggens
    else
        base := parentS!.base;
        layer := parentS!.layer + 1;
        stronggens := parentS!.stronggens;
        nr := Length(stronggens);
        Append(stronggens,gens);
        layergens := [nr+1..Length(stronggens)];
    fi;

    gens := ShallowCopy(gens);
    # Note that we do ShallowCopy such that the original list and the
    # one in the orbit record are different from each other.
    if IsInt(size) then
        hashsize := NextPrimeInt(Minimum(size,opt.InitialHashSize));
    else
        hashsize := opt.InitialHashSize;
    fi;
    orb := Orb( gens, nextpoint, nextop,
                rec( treehashsize := hashsize, schreier := true, 
                     log := opt.OrbitsWithLog,
                     report := opt.Report ) );
    S := rec( stab := false, orb := orb, cand := cand, base := base,
              opt := opt, layer := layer, parentS := parentS,
              stronggens := stronggens, layergens := layergens,
              size := size, randpool := [], IsOne := opt.IsOne );
    if parentS = false then
        S!.topS := S;
    else
        S!.topS := parentS!.topS;
    fi;
    if IsBound(opt.FindBasePointCandidatesData) then
        S.randpool := opt.FindBasePointCandidatesData.randpool;
    fi;

    Add(base,nextpoint);
    S.orb!.stabilizerchain := S;
    Objectify( StabChainByOrbType, S );

    return S;
  end );

InstallGlobalFunction( GENSS_RandomElementFromAbove,
  function( S, i )
    # This function provides a random element for the i-th stabiliser
    # in the chain "from above". These elements eventually come from
    # the one product replacer object in the opt record for the whole
    # group. They are sifted through i layers of the stabiliser chain
    # infrastructure until they stabilise the first i points (i=0 simply
    # means a random element in the whole group). If we find out underways
    # that an orbit is too small (that is, a stabiliser was not complete),
    # we fix the stabiliser chain as we go. Random elements that have
    # been generated in a layer j < i previously and are not yet used
    # as strong generators can be taken from a pool that was kept in
    # the stabiliser chain. S must either be layer i of the stabiliser 
    # or (for i=0) it can be equal to the options record opt chain.
    local SS, x, topS, j, o, p, po;
    if i = 0 then
        return Next(S.pr);    # in this case we got the options record
    fi;
    if S!.layer <> i then
        Error("i must be equal to the layer");
    fi;
    SS := S;
    while true do   # will be left by break
        if Length(SS!.randpool) > 0 then
            x := Remove(SS!.randpool,Length(SS!.randpool));
            topS := SS!.topS;
            break;
        fi;
        if SS!.parentS = false then   # we have reached the top
            if IsBound(SS!.opt.RandomElmFunc) then
                x := SS!.opt.RandomElmFunc();
            else
                x := Next(SS!.opt.pr);
            fi;
            topS := SS;   # remember top
            break;
        fi;
        SS := SS!.parentS;
    od;
    # We now have a random element x in layer SS!.layer, that is,
    # it is contained in the SS!.layer-1-th stabiliser, sift it down:
    j := SS!.layer-1;
    while j < i do
        o := SS!.orb;
        if IsObjWithMemory(x) then
          p := o!.op(o[1],x!.el);
        else
          p := o!.op(o[1],x);
        fi;
        po := Position(o,p);
        if po = fail then   # not in current stabilizer
            Info(InfoGenSS,3,"Random element from top found error in layer ",
                             SS!.layer);
            AddGeneratorToStabilizerChain(topS,x);
            # We add it at the top to add it in every orbit above except
            # the first one!
            return GENSS_RandomElementFromAbove(S,i);
        fi;
        # Now sift through Schreier tree:
        while po > 1 do
            x := x * SS!.orb!.gensi[o!.schreiergen[po]];
            po := o!.schreierpos[po];
        od;
        SS := SS!.stab;
        j := j + 1;
    od;
    # After this, we have successfully reached i-th stabilizer
    return x;
  end );
   
InstallGlobalFunction( GENSS_StabilizerChainInner,
  function( gens, size, layer, cand, opt, parentS )
    # Computes a stabilizer chain for the group generated by gens
    # with known size size (can be false if not known). This will be
    # layer layer in the final stabilizer chain. cand is a (shared)
    # record for base point candidates and opt the (shared) option
    # record. This is called in StabilizerChain and calls itself.
    # It also can be called if a new layer is needed.
    local base,gen,S,i,merk,merk2,next,pr,r,stabgens,x;

    Info(InfoGenSS,4,"Entering GENSS_StabilizerChainInner layer=",layer);
    next := GENSS_NextBasePoint(gens,cand,opt,parentS);
    cand := next.cand;   # This could have changed
    S := GENSS_CreateStabChainRecord(parentS,gens,size,
                                     next.point,next.op,next.cand,opt);
    base := S!.base;

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
    Info(InfoGenSS, 2, "Layer ", layer, ": Orbit length is ", Length(S!.orb)); 

    if layer > 1 then
        parentS!.stab := S;   # such that from now on random element
                              # generation works!
    else
        if (Length(S!.orb) > 50 or S!.orb!.depth > 5) and
           S!.opt.OrbitsWithLog then
            Info(InfoGenSS, 3, "Trying to make Schreier tree shallower (depth=",
                 S!.orb!.depth,")...");
            merk := Length(S!.orb!.gens);
            merk2 := Length(S!.stronggens);
            MakeSchreierTreeShallow(S!.orb);
            Append(S!.stronggens,S!.orb!.gens{[merk+1..Length(S!.orb!.gens)]});
            Append(S!.layergens,[merk2+1..Length(S!.stronggens)]);
            Info(InfoGenSS, 3, "Depth is now ",S!.orb!.depth);
        fi;
    fi;
    S!.orb!.gensi := List(S!.orb!.gens,x->x^-1);
 

    # Are we done?
    if size <> false and Length(S!.orb) = size then
        S!.proof := true;
        Info(InfoGenSS,4,"Leaving GENSS_StabilizerChainInner layer=",layer);
        return S;
    fi;

    # Now create a few random stabilizer elements:
    stabgens := EmptyPlist(opt.RandomStabGens);
    for i in [1..opt.RandomStabGens] do
        x := GENSS_RandomElementFromAbove(S,layer);
        if not(S!.IsOne(x)) then
            Add(stabgens,x);
        fi;
    od;
    Info(InfoGenSS,3,"Created ",opt.RandomStabGens,
         " random stab els, ",
         Length(stabgens)," non-trivial.");
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
#Error("foobar");
            while i < opt.ImmediateVerificationElements do
                i := i + 1;
                x := GENSS_RandomElementFromAbove(S,layer);
                if AddGeneratorToStabilizerChain(S!.topS,x) then
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

InstallMethod(VerifyStabilizerChainMC, 
  "for a stabilizer chain and a positive integer",
  [ IsStabilizerChainByOrb, IsInt ],
  function( S, nrels )
    local i,x;

    Info(InfoGenSS,2,"Doing randomized verification...");
    i := 0; 
    while i < nrels do
        i := i + 1;
        if IsBound(S!.opt.RandomElmFunc) then
            x := S!.opt.RandomElmFunc();
        else
            x := Next(S!.opt.pr);
        fi;
        if AddGeneratorToStabilizerChain(S,x) then
            Info( InfoGenSS, 2, "Verification found error ... ",
                  "new size ", Size(S) );
            i := 0;
        fi;
    od;
  end );

InstallMethod( StabilizerChain, "for a group object and a record", 
  [ IsGroup, IsRecord ],
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
    #   Reduced:    if set to true no orbits of length 1 will occur
    #               (except for the trivial group), is true by default
    #   Cand:       initializer for cand, which are base points
    #               candidates
    #               must be a record with components "points", "ops",
    #               "used", the latter indicates the largest index
    #               that has already been used. "used" is automatically
    #               set to 0 if not set.
    #   TryShortOrbit: Number of tries for the short orbit finding alg.
    #   StabGenScramble,
    #   StabGenScrambleFactor,
    #   StabGenAddSlots,
    #   StabGenMaxDepth:   parameters for product replacer for generating
    #                      random elements
    #   VeryShortOrbLimit: when looking for short orbits try a random
    #                      vector and enumerate its orbit until this limit
    #   
    #   ... to be continued
    #
    local S,cand,i,pr,prob,x,gens,SS;

    # First a few preparations, then we delegate to GENSS_StabilizerChainInner:

    # Add some default options:
    if (HasSize(grp) and not(IsBound(opt.Projective) and opt.Projective))
       or IsBound(opt.Size) then
        if not(IsBound(opt.ImmediateVerificationElements)) then
            opt.ImmediateVerificationElements := 0;
        fi;
    fi;
    if IsBound(opt.Projective) and opt.Projective then
        opt.IsOne := GENSS_IsOneProjective;
    elif not(IsBound(opt.IsOne)) then
        opt.IsOne := IsOne;
    fi;
    # Now opt.IsOne is set to a function to check whether or not a group
    # element is equal to the identity.
    GENSS_CopyDefaultOptions(GENSS,opt);

    # Check for the identity group:
    gens := GeneratorsOfGroup(grp);
    if Length(gens) = 0 or 
       ForAll(gens,opt.IsOne) then
        # Set up a trivial stabilizer chain record:
        S := GENSS_CreateStabChainRecord(false,gens,1,1,GENSS_TrivialOp,
                                         rec( points := [], ops := [],
                                              used := 0),opt);
        Enumerate(S!.orb);
        S!.orb!.gensi := List(S!.orb!.gens,x->x^-1);
        S!.proof := true;
        S!.trivialgroup := true;
        if (not(IsBound(opt.Projective)) or
            opt.Projective = false) and
           IsIdenticalObj(opt.IsOne,IsOne) then
            SetStoredStabilizerChain(grp,S);
        fi;
        return S;
    fi;
    
    # Setup a random element generator for the whole group and store
    # it in the opt record:
    opt.pr := ProductReplacer( gens,
                      rec( scramble := opt.StabGenScramble,
                           scramblefactor := opt.StabGenScrambleFactor,
                           addslots := opt.StabGenAddSlots,
                           maxdepth := opt.StabGenMaxDepth ));

    # Old style error probability for compatibility:
    if IsBound(opt.random) then
        if opt.random = 0 then
            opt.VerifyElements := 0;
        elif opt.random = 1000 then
            opt.DeterministicVerification := true;
            Info(InfoGenSS,1,"Warning: Deterministic verification not yet ",
                 "implemented!");
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
            cand := rec( points := opt.Base, used := 0, points2 := [],
                         ops2 := [] );
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
        if not(IsBound(cand.used)) then
            cand.used := 0;
        fi;
        cand.points2 := [];
        cand.ops2 := [];
    else
        # Otherwise try different things later using generic methods:
        cand := rec( points := [], ops := [], used := 0 );
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

    # Do we already have a proof?
    if S!.proof then
        if not(IsBound(opt.Projective)) or opt.Projective = false then
            SetStoredStabilizerChain(grp,S);
        fi;
        return S;
    fi;
    if opt.VerifyElements = 0 then return S; fi;

    Info(InfoGenSS,2,"Current size found: ",Size(S));
    # Now a possible verification phase:
    if S!.size <> false then   # we knew the size in advance
        Info(InfoGenSS,2,"Doing verification via known size...");
        while Size(S) < S!.size do
            Info(InfoGenSS,2,"Known size not reached, throwing in a random ",
                 "element...");
            if IsBound(opt.RandomElmFunc) then
                x := opt.RandomElmFunc();
            else
                x := Next(opt.pr);
            fi;
            if AddGeneratorToStabilizerChain(S,x) then
                Info( InfoGenSS, 2, "Increased size to ",Size(S) );
            fi;
        od;
        S!.proof := true;
        if (not(IsBound(opt.Projective)) or opt.Projective = false) and
           IsIdenticalObj(opt.IsOne,IsOne) then
            SetStoredStabilizerChain(grp,S);
        fi;
    else
        # Do some verification here:
        VerifyStabilizerChainMC(S,opt.VerifyElements);
    fi;
    # Now clean up the random element pools to save memory:
    SS := S;
    while SS <> false do
        if IsBound(SS!.randpool) and Length(SS!.randpool) > 0 then
            SS!.randpool := [];
        fi;
        SS := SS!.stab;
    od;
    return S;
  end );

InstallMethod( AddGeneratorToStabilizerChain,
  "for a stabilizer chain and a new generator",
  [ IsStabilizerChain and IsStabilizerChainByOrb, IsObject ],
  function( S, el )
    # Increases the set represented by S by the generator el.
    local SS, r, n, pr, i, newstrongnr;
    if IsBound(S!.trivialgroup) and S!.trivialgroup then
        if S!.IsOne(el) then
            return false;
        fi;
        SS := StabilizerChain(Group(el),S!.opt);
        if IsString(SS) then return SS; fi;
        for n in NamesOfComponents(SS) do
            S!.(n) := SS!.(n);
        od;
        Unbind(S!.trivialgroup);
        return true;
    fi;

    r := SiftGroupElement( S, el );
    # if this is one then el is already contained in the stabilizer chain
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
        Add(SS!.stronggens,r.rem);
        Add(SS!.layergens,Length(SS!.stronggens));
        AddGeneratorsToOrbit(SS!.orb,[r.rem]);
        Add(SS!.orb!.gensi,r.rem^-1);
        newstrongnr := Length(SS!.stronggens);
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
        newstrongnr := Length(SS!.stronggens)+1;  # r.rem will end up there !
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
    while S!.layer < SS!.layer do
        Info(InfoGenSS,2,"Adding new generator to orbit in layer ",S!.layer);
        Add(S!.layergens,newstrongnr);
        AddGeneratorsToOrbit(S!.orb,[r.rem]);
        Add(S!.orb!.gensi,r.rem^-1);
        S := S!.stab;
    od;
    # Finally, we have to add it to the product replacer!
    AddGeneratorToProductReplacer(S!.opt!.pr,r.rem);
    return true;
  end );

InstallMethod( SiftGroupElement, "for a stabilizer chain and a group element",
  [ IsStabilizerChain and IsStabilizerChainByOrb, IsObject ],
  function( S, x )
    local o,p,po,preS,r,isone;
    isone := S!.IsOne;
    preS := false;
    while S <> false do
        o := S!.orb;
        if IsObjWithMemory(x) then
          p := o!.op(o[1],x!.el);
        else
          p := o!.op(o[1],x);
        fi;
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
    r := rec( rem := x, S := false, preS := preS, isone := isone(x) );
    return r;
  end );

InstallMethod( SiftGroupElementSLP, 
  "for a stabilizer chain and a group element",
  [ IsStabilizerChain and IsStabilizerChainByOrb, IsObject ],
  function( S, x )
    local preS, nrstrong, slp, o, p, po, r, isone;
    preS := false;
    isone := S!.IsOne;
    nrstrong := Length(S!.stronggens);
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
            Add(slp,S!.layergens[o!.schreiergen[po]]);
            po := o!.schreierpos[po];
        od;
        preS := S;
        S := S!.stab;
    od;
    r := rec( rem := x, S := false, preS := preS, isone := isone(x) );
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
    return S!.stronggens;
  end );

InstallMethod( NrStrongGenerators, "for a stabilizer chain",
  [ IsStabilizerChain and IsStabilizerChainByOrb ],
  function( S )
    return Length(S!.stronggens);
  end );

InstallMethod( ForgetMemory, "for a stabilizer chain",
  [ IsStabilizerChain and IsStabilizerChainByOrb ],
  function( S )
    ForgetMemory(S!.stronggens);
    while S <> false do
        ForgetMemory(S!.orb);
        ForgetMemory(S!.orb!.gensi);
        S := S!.stab;
    od;
  end );

InstallMethod( AddNormalizingGenToLayer,
  "for a stabilizer chain, a group element, and a prime number",
  [ IsStabilizerChain and IsStabilizerChainByOrb, IsObject, IsPosInt ],
  function( S, x, p )
    # This is a very specialised function which is not officially documented.
    # It is used for the matrix PCGS code where a stabiliser chain is
    # built up using a known base in a generator by generator fashion
    # where the generators are the generators of the PCGS. Thus a lot
    # is known about them.
    # This function assumes that x is a new generator normalising the
    # group generated by the previous generators in the layer described
    # by S. So in particular, S is not necessarily the top of the 
    # stabiliser chain but the layer where actually something happens.
    # Namely, since it is known that the new generator generates
    # a group which is p times as big as the previous one (p a prime)
    # it is known that the orbit in the given layer will become
    # p times as big and we can enumerate it relatively cheaply.
    local lmp,o,oldnrgens;
    o := S!.orb;
    oldnrgens := Length(o!.gens);
    Add(o!.gens,x);
    Add(o!.gensi,x^-1);
    if IsPermOnIntOrbitRep(o) then
        lmp := LargestMovedPoint(o!.gens);
        if lmp > Length(o!.tab) then
            Append(o!.tab,ListWithIdenticalEntries(lmp-Length(o!.tab),0));
        fi;
    fi;
    ResetFilterObj(o,IsClosed);
    o!.pos := 1;
    o!.genstoapply := [oldnrgens+1..Length(o!.gens)];
    repeat
        Enumerate(o,S!.opt.OrbitLengthLimit);
        if not(IsClosed(o)) then
            if S!.opt.FailInsteadOfError then
                return "Orbit too long, increase S!.opt.OrbitLengthLimit";
            else
                Error("Orbit too long, increase S!.opt.OrbitLengthLimit!");
            fi;
        fi;
    until IsClosed(S!.orb);
    o!.genstoapply := [1..Length(o!.gens)];

    # Now fix up the stabilizer chain record:
    if not x in S!.stronggens then
      Add(S!.stronggens,x);
      Add(S!.layergens,Length(S!.stronggens));
    fi;
    Unbind(S!.size);   # whatever we knew before, we know no longer
 
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
  function( S ) Error("Currently not functional!"); return fail; end );
##    function( S )
##      local Grels,Hrels,Prels,MakeSchreierGens,ct,f,gens,gensi,i,j,k,l,li,max,
##            newpres,nrgens,nrschr,o,ords,pres,sb,sgs,slp,subgens,v,w,x,
##            cosetnrlimitfactor;
##  
##      allstrong := function(S)
##        st := Set(S!.layergens);
##        while S!.stab <> false do
##            S := S!.stab;
##            UniteSet(st,S!.layergens);
##        od;
##        return st;
##      end;
##      if S!.stab <> false then
##          pres := VerifyStabilizerChainTC(S!.stab);
##          if IsList(pres) then return pres; fi;
##          strongbelow := allstrong(S!.stab);
##      else
##          pres := StraightLineProgram([[]],0);
##          strongbelow := [];
##      fi;
##      strong := allstrong(S);
##      Info(InfoGenSS,1,"Verifying stabilizer chain in layer ",S!.layer);
##  
##      # First create a few Schreier generators:
##      sgs := [];
##      i := 1;
##      j := 1;
##      o := S!.orb;
##      nrgens := Length(o!.gens);
##      MakeSchreierGens := function(n)
##          local sg;
##          Info(InfoGenSS,3,"Creating ",n," Schreier generators...");
##          while Length(sgs) < n and
##                i <= Length(o) do
##              sg := GENSS_CreateSchreierGenerator(S,i,j);
##              j := j + 1;
##              if j > nrgens then
##                  j := 1;
##                  i := i + 1;
##              fi;
##              if sg <> fail then
##                  Add(sgs,sg);
##              fi;
##          od;
##      end;
##  
##      nrschr := S!.opt.NumberSchreierGens;
##      MakeSchreierGens(nrschr);
##      f := FreeGroup(Length(strong));
##      gens := GeneratorsOfGroup(f);
##      gensi := List(gens,x->x^-1);
##      subgens := gens{List(strongbelow,i->Position(strong,i))};
##      Hrels := ResultOfStraightLineProgram(pres,subgens);
##      if S!.opt.Projective then
##          ords := List([1..nrgens],i->ProjectiveOrder(o!.gens[i]));
##      else
##          ords := List([1..nrgens],i->Order(o!.gens[i]));
##      fi;
##      Prels := List([1..nrgens],i->gens[i+sb]^ords[i]);
##      Grels := [];
##      cosetnrlimitfactor := 4;
##      while true do   # will be left by return eventually
##          for k in [Length(Grels)+1..Length(sgs)] do
##              Grels[k] := GENSS_Prod(gens,sgs[k][1]+sb) * gens[sgs[k][2]+sb] * 
##                          GENSS_Prod(gensi,sgs[k][3]+sb);
##              x := GENSS_Prod(o!.gens,sgs[k][1]) * o!.gens[sgs[k][2]] * 
##                   GENSS_Prod(o!.gensi,sgs[k][3]);
##              if S!.stab <> false then
##                  slp := SiftGroupElementSLP(S!.stab,x);
##                  if not(slp.isone) then
##                      return [fail,S!.layer];
##                  fi;
##                  Grels[k] := Grels[k] / ResultOfStraightLineProgram(slp.slp,
##                                                                     subgens);
##                  sgs[k][4] := slp.slp;
##              else
##                  if not(S!.IsOne(x)) then
##                      return [fail,S!.layer];
##                  fi;
##                  sgs[k][4] := false;
##              fi;
##          od;
##          Info(InfoGenSS,2,"Doing coset enumeration with limit ",
##               cosetnrlimitfactor*Length(o));
##          ct := CosetTableFromGensAndRels(gens,Concatenation(Prels,Hrels,Grels),
##                     subgens:max := cosetnrlimitfactor*Length(o),silent);
##          if ct = fail then   # did not close!
##              cosetnrlimitfactor := QuoInt(cosetnrlimitfactor*3,2);
##              Info(InfoGenSS,2,"Coset enumeration did not finish!");
##          #Error(1);
##              if nrschr > Length(sgs) # or
##                 # nrschr > S!.opt.MaxNumberSchreierGens 
##                 then   # we are done!
##                  # Something is wrong!
##                  return [fail, S!.layer];
##              fi;
##          else
##              Info(InfoGenSS,2,"Coset enumeration found ",Length(ct[1]),
##                   " cosets.");
##          #Error(2);
##              if Length(ct[1]) = Length(o) then
##                  # Verification is OK, now build a presentation:
##                  l := GeneratorsWithMemory(
##                         ListWithIdenticalEntries(S!.nrstrong,()));
##                  li := List(l,x->x^-1);
##                  newpres := ResultOfStraightLineProgram(pres,
##                                     l{[1..S!.strongbelow]});
##                  for k in [1..nrgens] do
##                      Add(newpres,l[k+sb]^ords[k]);
##                  od;
##                  for k in [1..Length(sgs)] do
##                      if sgs[k][4] <> false then
##                          Add(newpres,
##                              GENSS_Prod(l,sgs[k][1]+sb)*l[sgs[k][2]+sb]*
##                              GENSS_Prod(li,sgs[k][3]+sb)*
##                              ResultOfStraightLineProgram(sgs[k][4],l)^-1);
##                      else
##                          Add(newpres,GENSS_Prod(l,sgs[k][1]+sb)*l[sgs[k][2]+sb]*
##                                      GENSS_Prod(li,sgs[k][3]+sb));
##                      fi;
##                  od;
##                  Info(InfoGenSS,2,"Found presentation for layer ",S!.layer,
##                       " using ",Length(newpres)," relators.");
##                  return SLPOfElms(newpres);
##              fi;
##          fi;
##          nrschr := QuoInt(nrschr*4,2);
##          MakeSchreierGens(nrschr);
##      od;
##    end);


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


##  VerifyStabilizerChainTC5 := 
##    function( S )
##      local FindTwoWords,Grels,Hrels,Prels,allgens,cosetlimit,done,el,f,
##            gens,gensi,hom,k,l,li,newpres,newrel,nrcosets,nrgens,o,opgens,
##            ords,pres,r,rels,sb,stronggens,strongi,subgens,tc,words;
##  
##      if S!.stab <> false then
##          pres := VerifyStabilizerChainTC5(S!.stab);
##          if IsList(pres) then return pres; fi;
##      else
##          pres := StraightLineProgram([[]],0);
##      fi;
##      Info(InfoGenSS,1,"Verifying stabilizer chain in layer ",S!.layer);
##  
##      o := S!.orb;
##      nrgens := Length(o!.gens);
##      sb := S!.strongbelow;
##              
##      f := FreeGroup(sb+nrgens);
##      gens := GeneratorsOfGroup(f);
##      gensi := List(gens,x->x^-1);
##      allgens := 0*[1..2*Length(gens)];
##      allgens{[1,3..2*Length(gens)-1]} := gens;
##      allgens{[2,4..2*Length(gens)]} := gensi;
##      subgens := gens{[1..S!.strongbelow]};
##      Hrels := ResultOfStraightLineProgram(pres,subgens);
##      if S!.opt.Projective then
##          ords := List([1..nrgens],i->ProjectiveOrder(o!.gens[i]));
##      else
##          ords := List([1..nrgens],i->Order(o!.gens[i]));
##      fi;
##      Prels := List([1..nrgens],i->gens[i+sb]^ords[i]);
##      Grels := [];
##      stronggens := StrongGenerators(S);
##      strongi := List(stronggens,x->x^-1);
##      opgens := 0*[1..2*Length(stronggens)];
##      opgens{[1,3..2*Length(stronggens)-1]} := stronggens;
##      opgens{[2,4..2*Length(stronggens)]} := strongi;
##  
##      # Now start up a coset enumeration:
##      cosetlimit := Maximum(QuoInt(7 * Length(o),6),Length(o)+5);
##      Info(InfoGenSS,2,"Starting coset enumeration with limit ",
##           cosetlimit," and ",Length(Hrels),
##           "+",Length(Prels),"+",Length(Grels)," relations...");
##      rels := Concatenation(Hrels,Prels);
##      tc := GENSS_CosetTableFromGensAndRelsInit(gens,rels,subgens,cosetlimit);
##      done := GENSS_DoToddCoxeter(tc);
##          
##      FindTwoWords := function(o,opgens,table)
##          local TraceWord,cosets,cosetsrevtab,i,j,new,nrcosets,pt,pts,
##                ptsrev,schgen,schpt,w1,w2,x,y;
##          nrcosets := Length(table[1]);
##          cosets := [1];
##          cosetsrevtab := 0*[1..nrcosets];
##          cosetsrevtab[1] := 1;
##          pts := [o[1]];
##          ptsrev := 0*[1..Length(o)];
##          ptsrev[1] := 1;
##          schpt := [fail];    # the Schreier tree
##          schgen := [fail];
##          i := 1;
##          TraceWord := function(pos)
##              local w;
##              w := [];
##              while pos > 1 do
##                  Add(w,schgen[pos]);
##                  pos := schpt[pos];
##              od;
##              return Reversed(w);
##          end;
##          while i <= Length(cosets) do
##              for j in [1..Length(opgens)] do
##                  x := table[j][cosets[i]];
##                  if x <> 0 then   # image is defined:
##                      if cosetsrevtab[x] = 0 then   # not visited
##                          Add(cosets,x);
##                          new := Length(cosets);
##                          cosetsrevtab[x] := new;
##                          schpt[new] := i;
##                          schgen[new] := j;
##                          pt := o!.op(pts[i],opgens[j]);
##                          y := Position(o,pt);
##                          if ptsrev[y] = 0 then
##                              Add(pts,pt);
##                              ptsrev[y] := new;
##                          else
##                              # We have reached a new coset by a word that
##                              # maps the starting point of the orbit to the 
##                              # same point as the one of another coset!
##                              w1 := TraceWord(ptsrev[y]);
##                              w2 := TraceWord(new);
##                              return [w1,w2];
##                          fi;
##                      fi;
##                  fi;
##              od;
##              i := i + 1;
##          od;
##          Error("Bad, this should never have been reached!");
##          return fail;
##      end;
##  
##      while true do   # will be left by return eventually
##          if done = true then
##              nrcosets := Length(tc.table[1]);
##              Info(InfoGenSS,2,"Coset enumeration found ",nrcosets," cosets.");
##              if nrcosets = Length(o) then
##                  # Verification is OK, now build a presentation:
##                  l := GeneratorsWithMemory(
##                         ListWithIdenticalEntries(S!.nrstrong,()));
##                  newpres := ResultOfStraightLineProgram(pres,
##                                     l{[1..S!.strongbelow]});
##                  for k in [1..nrgens] do
##                      Add(newpres,l[k+sb]^ords[k]);
##                  od;
##                  hom := GroupHomomorphismByImagesNC(f,Group(l),gens,l);
##                  for k in Grels do
##                      Add(newpres,ImageElm(hom,k));
##                  od;
##                  Info(InfoGenSS,2,"Found presentation for layer ",S!.layer,
##                       " using ",Length(newpres)," relators.");
##                  return SLPOfElms(newpres);
##              elif nrcosets < Length(o) then
##                  Error("This cannot possibly have happened!");
##                  return [fail,S!.layer];
##              else   # nrcosets > Length(o)
##                  Info(InfoGenSS,2,"Too many cosets, we must have forgotten ",
##                       "another relation!");
##                  Error("This cannot possibly have happened2!");
##                  return [fail,S!.layer];
##              fi;
##          fi;
##          # Now we have to find another relation, we do a breadth-first
##          # search through the already defined cosets to find two cosets
##          # that are still different but ought to be equal because the
##          # corresponding orbit points are equal:
##          words := FindTwoWords(o,opgens,tc.table);
##          el := EvaluateWord(opgens,words[1])/EvaluateWord(opgens,words[2]);
##          r := SiftGroupElementSLP(S,el);
##          if not(r.isone) then
##              # Error, we found a new stabilizer element!
##              return [fail,S!.layer,el];
##          fi;
##          newrel := [(EvaluateWord(allgens,words[1])
##                      /EvaluateWord(allgens,words[2]))
##                     / ResultOfStraightLineProgram(r.slp,gens)];
##          GENSS_TCAddRelators(tc,newrel);
##          Add(Grels,newrel[1]);
##          if Length(words[1]) > 0 then
##              tc.app[10] := words[1][1];
##          else
##              # More difficult:
##              words := ExtRepOfObj(newrel[1]);
##              if words[2] > 0 then
##                  tc.app[10] := 2*words[1]-1;
##              else
##                  tc.app[10] := 2*words[1];
##              fi;
##          fi;
##          #Print("<\c");
##          #for i in [1..tc.limit] do
##          #    for j in [1..Length(allgens)] do
##          #        if tc.table[j][i] <> 0 then
##          #            tc.app[11] := i;
##          #            tc.app[10] := j;
##          #            tc.nrdel := tc.nrdel + MakeConsequences( tc.app );
##          #        fi;
##          #    od;
##          #od;
##          tc.app[11] := 1;
##          tc.nrdel := tc.nrdel + MakeConsequences( tc.app );
##          #Print("-\c");
##          done := GENSS_DoToddCoxeter(tc);
##          #Print(">\c");
##          if Length(Grels) mod 100 = 0 then
##              #Print("\n");
##              Info(InfoGenSS,2,"Currently using ",Length(Hrels),"+",
##                   Length(Prels),"+",Length(Grels)," relations.");
##          fi;
##      od;
##    end;
##  
##  VerifyStabilizerChainMax := 
##    function( S )
##      local bl,count,d,gens,i,invtab,j,k,o,p,r,res,s,w1,w2,x,xi,y,isone;
##  
##      isone := S!.IsOne;
##      if S!.stab <> false then
##          res := VerifyStabilizerChainMax(S!.stab);
##          if res <> true then return res; fi;
##      fi;
##      Info(InfoGenSS,1,"Verifying stabilizer chain in layer ",S!.layer,
##           " (orbit of length ",Length(S!.orb),")");
##  
##      count := 0;
##      gens := ShallowCopy(Set(S!.orb!.gens));
##      invtab := [];
##      k := Length(gens);
##      for i in [1..k] do
##          x := gens[i];
##          xi := x^-1;
##          p := Position(gens,xi);
##          if p = fail then
##              Add(gens,xi);
##              invtab[i] := Length(gens);
##              invtab[Length(gens)] := i;
##          else
##              invtab[i] := p;
##          fi;
##      od;
##      o := Orb(gens,S!.orb[1],S!.orb!.op,rec(schreier := true, log := true));
##      Enumerate(o);
##  
##      if Length(o) <> Length(S!.orb) then
##          Error("something is fishy, orbits do not have the same length");
##          return fail;
##      fi;
##  
##      # Now we have a nice breadth-first orbit such that all inverses of 
##      # generators are again generators.
##      # First check whether the generators fix the first point:
##      for x in gens do
##          if o!.op(o[1],x) = o[1] then
##              count := count + 1;
##              if S!.stab = false then
##                  r := rec( isone := isone(x) );
##              else
##                  r := SiftGroupElement(S!.stab,x);
##              fi;
##              if not(r.isone) then
##                  return [fail,S!.layer,x,"generator"];
##              fi;
##          fi;
##      od;
##  
##      bl := BlistList([1..Length(o)],[]);
##      for d in [1..o!.depth] do
##          Info(InfoGenSS,2,"Testing in depth ",d," (of ",o!.depth,") checking ",
##               o!.depthmarks[d+1]-o!.depthmarks[d]," points (of ",Length(o),")");
##          for i in [o!.depthmarks[d]..o!.depthmarks[d+1]-1] do
##              # These indices contain points in depth d
##              if bl[o!.schreierpos[i]] then
##                  # this subtree does not need to be done!
##                  bl[i] := true;   # mark subtree as done
##              else
##                  x := o[i];
##                  for j in [1..Length(gens)] do
##                      if j <> invtab[o!.schreiergen[i]] then
##                          y := o!.op(x,gens[j]);
##                          p := Position(o,y);
##                          # if p > i then this is a new point, do nothing
##                          if p <= i then
##                              count := count + 1;
##                              w1 := TraceSchreierTreeForward(o,i);
##                              w2 := TraceSchreierTreeForward(o,p);
##                              s := (EvaluateWord(gens,w1)*gens[j])
##                                   /EvaluateWord(gens,w2);
##                              #Print(Length(gens),w1,j,w2,"\n");
##                              if S!.stab = false then
##                                  r := rec( isone := isone(s) );
##                              else
##                                  r := SiftGroupElement(S!.stab,s);
##                              fi;
##                              if not(r.isone) then
##                                  return [fail,S!.layer,s,"schreier gen"];
##                              fi;
##                          fi;
##                          if p < o!.depthmarks[d] then
##                              # this goes up in the tree, this means, if this
##                              # Schreier generator is OK, we do not have to
##                              # look below!
##                              bl[i] := true;
##                              break;
##                          fi;
##                      fi;
##                  od;
##              fi;
##          od;
##      od;
##      Info(InfoGenSS,2,"Have sifted ",count," elements.");
##      return true;
##    end;

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

InstallMethod( BaseStabilizerChain,
  "generic method for stabilizer chains",
  [IsStabilizerChain],
  function( S )
    local SS,pts,ops;
    pts := [];
    ops := [];
    SS := S;
    while SS <> false do
        Add(pts,SS!.orb[1]);
        Add(ops,SS!.orb!.op);
        SS := SS!.stab;
    od;
    return rec( points := pts, ops := ops );
  end );

InstallMethod( SiftBaseImage,
  "generic method for genss stabilizer chains",
  [IsStabilizerChain, IsList],
  function(S,bi)
    local l, SS, i, o, pos;
    l := Length(bi);
    SS := S;
    i := 1;
    while i <= Length(bi) do
        o := SS!.orb;
        pos := Position(o,bi[i]);
        if pos = fail then
            return false;
        fi;
        while pos > 1 do
            if o!.memorygens then
                bi{[i..l]} := GENSS_MapBaseImage(bi{[i..l]},
                                o!.gensi[o!.schreiergen[pos]]!.el,SS);
            else
                bi{[i..l]} := GENSS_MapBaseImage(bi{[i..l]},
                                o!.gensi[o!.schreiergen[pos]],SS);
            fi;
            pos := o!.schreierpos[pos];
        od;
        if o[1] <> bi[i] then
            Error("this should not have happened, tell the authors");
        fi;
        i := i + 1;
        SS := SS!.stab;
    od;
    return true;
  end );

InstallMethod( IsProved, "for a stabilizer chain",
  [ IsStabilizerChain and IsStabilizerChainByOrb ],
  function( S )
    return S!.proof;
  end );

InstallMethod( StabChainOp, "for a permutation group and a stabilizer chain",
  [ IsPermGroup, IsStabilizerChain and IsStabilizerChainByOrb ],
  function( g, S )
    local base;
    if HasSize(g) and Size(g) <> Size(S) then
        Error("known size of group does not match stabiliser chain");
        return fail;
    else
        SetSize(g,Size(S));
    fi;
    base := BaseStabilizerChain(S);
    if not(ForAll(base.points,IsInt) and ForAll(base.ops,x->x=OnPoints)) then
        return StabChainOp(g);
    fi;
    return StabChainOp(g,rec( base := base.points, knownBase := base.points, 
                              size := Size(S) ));
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
        if S!.IsOne(el) then
            Info(InfoGenSS,3,"Found trivial stabilizer element.");
        else
            Add(stabgens,el);
            stab := GroupWithGenerators(stabgens);
            SS := StabilizerChain(stab,rec( Base := S, IsOne := S!.IsOne ));
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
            SS := StabilizerChain(stab,rec(Base := S, IsOne := S!.IsOne));
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
        if S!.IsOne(el) then
            Info(InfoGenSS,2,"Found trivial stabilizer element.");
        else
            Add(stabgens,el);
            stab := GroupWithGenerators(stabgens);
            SS := StabilizerChain(stab,rec( Base := S, IsOne := S!.IsOne ));
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
            SS := StabilizerChain(stab,rec(Base := S, IsOne := S!.IsOne));
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
        if S!.IsOne(el) then
            Info(InfoGenSS,3,"Found trivial stabilizer element.");
        else
            Add(stabgens,el);
            stab := GroupWithGenerators(stabgens);
            SS := StabilizerChain(stab,rec( Base := S, IsOne := S!.IsOne ));
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
       SS := StabilizerChain(stab,rec(Base := S, IsOne := S!.IsOne));
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

#############################################################################
# We can store a stabilizer chain in the attribute
# StoredStabilizerChain, then the following methods for groups apply:
#############################################################################

InstallMethod( SetStabilizerChain, "for a group and a stabilizer chain",
  [IsGroup, IsStabilizerChain],
  function(g,S)
    if IsIdenticalObj(S!.IsOne,IsOne) and
       HasSize(g) and Size(g) <> Size(S) then
      Error("you try to set a stabilizer chain for the wrong group");
      return;
    fi;
    SetStoredStabilizerChain(g,S);
  end );

InstallMethod( Size, "for a group with a stored stabilizer chain",
  [IsGroup and HasStoredStabilizerChain],
  function(g)
    return Size(StoredStabilizerChain(g));
  end );

InstallMethod( Size, "for a perm. group with a stored stabilizer chain",
  [IsPermGroup and HasStoredStabilizerChain],
  function(g)
    return Size(StoredStabilizerChain(g));
  end );

InstallMethod( Size, "for a group with a stored stabilizer chain",
  [IsGroup and HasStoredStabilizerChain and IsHandledByNiceMonomorphism],
  function(g)
    return Size(StoredStabilizerChain(g));
  end );

InstallMethod( \in, "for a group elm and a group with stored stabilizer chain",
  [IsPerm, IsPermGroup and HasStoredStabilizerChain], 1,
  function(x, g)
    local S,r;
    S := StoredStabilizerChain(g);
    r := SiftGroupElement(S,x);
    return r.isone;
  end );

InstallMethod( \in, "for a group elm and a group with stored stabilizer chain",
  [IsObject, IsMatrixGroup and HasStoredStabilizerChain and 
             IsHandledByNiceMonomorphism],
  function(x, g)
    local S,r;
    S := StoredStabilizerChain(g);
    r := SiftGroupElement(S,x);
    return r.isone;
  end );

InstallMethod( \in, "for a group elm and a group with stored stabilizer chain",
  [IsObject, IsGroup and HasStoredStabilizerChain],
  function(x, g)
    local S,r;
    S := StoredStabilizerChain(g);
    r := SiftGroupElement(S,x);
    return r.isone;
  end );

InstallOtherMethodWithRandomSource( Random, "for a random source and a stabilizer chain",
  [ IsRandomSource, IsStabilizerChainByOrb ],
  function(rs, S)
    local x,pos,w;
    x := One(S!.orb!.gens[1]);
    while S <> false do
        pos := Random(rs, 1,Length(S!.orb));
        w := TraceSchreierTreeForward(S!.orb,pos);
        x := GENSS_Prod(S!.orb!.gens,w) * x;
        S := S!.stab;
    od;
    return x;
  end );

InstallMethodWithRandomSource( Random, "for a random source and a group with a stored stabilizer chain",
  [ IsRandomSource, IsGroup and HasStoredStabilizerChain ],
  function(rs, g) return Random(rs, StoredStabilizerChain(g)); end );

InstallMethodWithRandomSource( Random, "for a random source and a permgroup with a stored stabilizer chain",
  [ IsRandomSource, IsPermGroup and HasStoredStabilizerChain ], 10,
  function(rs, g) return Random(rs, StoredStabilizerChain(g)); end );

InstallMethodWithRandomSource( Random, "for a random source and a matrixgroup with a stored stabilizer chain",
  [ IsRandomSource, IsMatrixGroup and IsHandledByNiceMonomorphism and 
    HasStoredStabilizerChain ],
  function(rs, g) return Random(rs, StoredStabilizerChain(g)); end );

InstallMethod( SizeMC, "for a group and an error bound",
  [IsGroup, IsRat],
  function( G, err )
    local S;
    S := StabilizerChain(G,rec( ErrorBound := err ));
    return Size(S);
  end );

InstallMethod( SizeMC, "for a permutation group and an error bound, genss",
  [IsGroup and IsPermGroup, IsRat], 1,
  function( G, err )
    local S;
    S := StabilizerChain(G,rec( ErrorBound := err ));
    return Size(S);
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
  function( g, h, images, opt1, opt2 )
    local Sg,Si,data,gm,im,slpstrongg,slpstrongi,strongg,stronggims,
          strongi,strongipre;
    gm := GroupWithMemory(g);
    Sg := StabilizerChain(gm,opt1);
    strongg := StrongGenerators(Sg);
    slpstrongg := SLPOfElms(strongg);
    ForgetMemory(Sg);
    stronggims := ResultOfStraightLineProgram(slpstrongg,images);
    im := GroupWithMemory(images);
    Si := StabilizerChain(im,opt2);
    strongi := StrongGenerators(Si);
    slpstrongi := SLPOfElms(strongi);
    ForgetMemory(Si);
    strongipre := ResultOfStraightLineProgram(slpstrongi,GeneratorsOfGroup(g));
    data := rec( Sg := Sg, stronggims := stronggims,
                 Si := Si, strongipre := strongipre,
                 slpstrongg := slpstrongg, slpstrongi := slpstrongi );
    return GroupHomByFuncWithData( g, h, GENSS_ImageElm, false,
                                   GENSS_PreImagesRepresentative, data );
  end );

InstallMethod( ORB_StabilizerChainKnownSize,
  "GENSS method for arbitrary groups",
  [IsGroup,IsPosInt],
  function(g,size)
    if HasStoredStabilizerChain(g) and
       size = Size(StoredStabilizerChain(g)) then
        return StoredStabilizerChain(g);
    fi;
    return StabilizerChain(g,rec( Size := size ));
  end );

InstallMethod( ORB_BaseStabilizerChain,
  "GENSS method for arbitary groups",
  [IsStabilizerChain],
  BaseStabilizerChain );

InstallMethod( ORB_StabilizerChainKnownBase,
  "GENSS method for arbitrary groups",
  [IsGroup,IsObject],
  function(g,base)
    if HasStoredStabilizerChain(g) then
        return StoredStabilizerChain(g);
    fi;
    if base = fail then
        return StabilizerChain( g );
    else
        return StabilizerChain( g, rec( Cand := ShallowCopy(base) ) );
    fi;
  end );

InstallMethod( ORB_SizeStabilizerChain,
  "GENSS method for arbitrary groups",
  [IsStabilizerChain], Size );

InstallMethod( ORB_IsWordInStabilizerChain,
  "GENSS method for arbitrary groups",
  [IsList,IsList,IsList,IsStabilizerChain],
  function(word, gens, gensi, S)
    local x, b, ops, i;
    if Size(S) = 1 then
        x := ORB_ApplyWord(gens[1]^0,word,gens,gensi,OnRight);
        return S!.IsOne(x);
    fi;
    b := BaseStabilizerChain(S);
    ops := b.ops;
    b := b.points;
    for i in [1..Length(word)] do
        if word[i] > 0 then
            b := List([1..Length(b)],j->ops[j](b[j],gens[word[i]]));
        else
            b := List([1..Length(b)],j->ops[j](b[j],gensi[-word[i]]));
        fi;
    od;
    return SiftBaseImage(S,b);
  end );

InstallMethod( ORB_IsElementInStabilizerChain,
  "GENSS method for arbitrary groups",
  [IsObject, IsStabilizerChain],
  function(el, S)
    return el in S;
  end );


#############################################################################
# The following methods are about methods to compute stabilisers:
#############################################################################

InstallMethod( Stab, "add empty options record",
  [IsGroup, IsObject, IsFunction],
  function( g, x, op )
    return Stab(g,x,op,rec());
  end );

InstallMethod( Stab, "add empty options record",
  [IsList, IsObject, IsFunction],
  function( l, x, op )
    return Stab(l,x,op,rec());
  end );

InstallMethod( Stab, "create group from list of generators",
  [IsList, IsObject, IsFunction, IsRecord],
  function( l, x, op, opt )
    return Stab(GroupWithGenerators(l),x,op,opt);
  end );

# Now the real one: 

InstallMethod( Stab, "by Orb orbit enumeration",
  [IsGroup, IsObject, IsFunction, IsRecord],
  function( g, x, op, opt )
    local S,count,el,errorprob,found,gens,i,j,limit,memperpt,nrrand,o,
          pat,pos,pr,res,stab,stabchain,stabel,stabgens,stabsizeest,w1,w2,y;
    GENSS_CopyDefaultOptions(GENSS,opt);
    if HasStoredStabilizerChain(g) then
      S := StoredStabilizerChain(g);
      if IsIdenticalObj(S!.orb!.op,op) and x in S!.orb then
        if S!.stab = false then
            y := rec( stab := TrivialSubgroup(g), size := 1, proof := true );
            y.stabilizerchain := StabilizerChain(y.stab);
            return y;
        fi;
        pos := Position(S!.orb,x);
        w1 := TraceSchreierTreeForward(S!.orb,pos);
        y := GENSS_Prod(S!.orb!.gens,w1);
        stab := Group(List(S!.stab!.orb!.gens,a->y^-1*a*y));
        return rec( stab := stab, size := Size(S!.stab),
                    stabilizerchain := StabilizerChain(stab,rec( Base := S)),
                    proof := true );
      fi;
    fi;
    # Now we have to do it ourselves:
    if not(IsBound(opt.ErrorBound)) then   # we are working deterministically
        if not(HasSize(g)) then
            # We are basically stuffed, unless we want to check all
            # Schreier generators of an orbit!
            Error("I do not want to check all Schreier generators, ",
                  "need size of group or ErrorBound option!");
            return fail;
        fi;
    else
        if opt.ErrorBound < opt.StabAssumeCompleteLimit then
            opt.StabAssumeCompleteLimit := opt.ErrorBound;
        fi;
    fi;
    # Now we either know the group size or work Monte Carlo:
    errorprob := 1;
    limit := opt.StabInitialLimit;
    pat := opt.StabInitialPatience;
    o := Orb(g,x,op,rec( report := opt.Report, 
                         treehashsize := opt.InitialHashSize,
                         schreier := true ) );
    stabgens := [];
    stabsizeest := 1;
    pr := ProductReplacer(GeneratorsOfGroup(g),
                rec( scramble := opt.StabScramble, 
                     scramblefactor := opt.StabScrambleFactor,
                     addslots := opt.StabAddSlots,
                     maxdepth := opt.StabMaxDepth ));
    gens := GeneratorsOfGroup(g);
    # Now go through a cycle of orbit enumeration and stabilizer generation:
    repeat
        if not(IsClosed(o)) then
            if HasSize(g) and limit > QuoInt(Size(g),stabsizeest*2)+1 then
                limit := QuoInt(Size(g),stabsizeest*2)+2;
            fi;
            Info(InfoGenSS,2,"Enumerating orbit with limit ",limit);
            Enumerate(o,limit);
            if IsClosed(o) then
                Info(InfoGenSS,2,"Orbit closed, size is ",Length(o));
            fi;
        fi;
        if errorprob < opt.StabAssumeCompleteLimit then 
            if HasSize(g) and 
               2 * Length(o) * stabsizeest > Size(g) then
                # Done!
                #stab := Subgroup(g,stabgens);
                if Length(stabgens) > 0 then
                    stab := Group(stabgens);
                    SetParent(stab,g);
                    SetSize(stab,stabsizeest);
                    SetStabilizerChain(stab,stabchain);
                else
                    stab := TrivialSubgroup(g);
                    stabchain := StabilizerChain(stab);
                fi;
                return rec( stab := stab, size := stabsizeest, 
                            stabilizerchain := stabchain,
                            proof := true );
            fi;
            limit := 2*limit;
            if not IsClosed(o) then continue; fi;
        fi;
        count := 0;
        found := false;
        repeat
            el := Next(pr);
            i := 1;
            while i <= Length(o) and not(found) do
                y := op(o[i],el);
                j := Position(o,y);
                if j <> fail then 
                    found := true; 
                else
                    i := i + 1;
                fi;
            od;
            count := count + 1;
        until found or count >= pat;
        if not(IsClosed(o)) then
            if not(found) then
                limit := QuoInt(3*limit,2);
            else
                if count < 10 then
                    limit := limit + 100;
                else
                    limit := QuoInt(3*limit,2);
                fi;
            fi;
        fi;
        pat := QuoInt(pat*5,4);
        if found then   # Have a potential stabiliser element:
            Info(InfoGenSS,3,"Found a potential stabilising element...");
            w1 := TraceSchreierTreeForward(o,i);
            w2 := TraceSchreierTreeForward(o,j);
            stabel := EvaluateWord(gens,w1)*el/EvaluateWord(gens,w2);
            if IsOne(stabel) then
                    errorprob := errorprob / 2;
                    Info(InfoGenSS,3,"... which was the identity.");
            else
                Add(stabgens,stabel);
                if Length(stabgens) < 2 then
                    Info(InfoGenSS,3,"Waiting for second stabiliser element.");
                else      # Length(stabgens) >= 2 then
                    if Length(stabgens) = 2 then
                        Info(InfoGenSS,2,"Estimating stabilizer...");
                        stabchain := StabilizerChain(Group(stabgens));
                        errorprob := 1;
                    else
                        res := AddGeneratorToStabilizerChain(stabchain,stabel);
                        if not(res) then
                            Remove(stabgens,Length(stabgens));
                            errorprob := errorprob / 2;
                            Info(InfoGenSS,2,"Error probablity now < ",
                                 errorprob);
                        else
                            Info(InfoGenSS,2,
                                 "Added generator to stabilizer...");
                            VerifyStabilizerChainMC(stabchain,10);
                            errorprob := 1;
                        fi;
                    fi;
                    if Size(stabchain) > stabsizeest then
                        stabsizeest := Size(stabchain);
                        Info(InfoGenSS,2,"New stabilizer estimate: ",
                             stabsizeest);
                    else
                        Info(InfoGenSS,2,"Stabiliser estimate unchanged.");
                    fi;
                    if HasSize(g) and
                       2 * Length(o) * stabsizeest > Size(g) then
                        # Done!
                        if Length(stabgens) > 0 then
                            stab := Group(stabgens);
                            SetParent(stab,g);
                            SetSize(stab,stabsizeest);
                            SetStabilizerChain(stab,stabchain);
                        else
                            stab := TrivialSubgroup(g);
                            stabchain := StabilizerChain(stab);
                        fi;
                        return rec( stab := stab, size := stabsizeest, 
                                    stabilizerchain := stabchain,
                                    proof := true );
                    fi;
                    if IsBound(opt.ErrorBound) and
                       errorprob < opt.ErrorBound then
                        #stab := Subgroup(g,stabgens);
                        if Length(stabgens) > 0 then
                            stab := Group(stabgens);
                            SetParent(stab,g);
                        else
                            stab := TrivialSubgroup(g);
                            stabchain := StabilizerChain(stab);
                        fi;
                        return rec( stab := stab, size := stabsizeest,
                                    stabilizerchain := stabchain,
                                    proof := false );
                    fi;
                fi;
            fi;
        else   # no element found
            Info(InfoGenSS,3,"No stabiliser element found!");
        fi;
    until Length(o) > opt.StabOrbitLimit;
    return fail;
  end );

InstallMethod( StabMC, "add empty options record",
  [IsGroup, IsObject, IsFunction],
  function( g, x, op )
    return StabMC(g,x,op,rec());
  end );

InstallMethod( StabMC, "add empty options record",
  [IsList, IsObject, IsFunction],
  function( l, x, op )
    return StabMC(l,x,op,rec());
  end );

InstallMethod( StabMC, "create group from list of generators",
  [IsList, IsObject, IsFunction, IsRecord],
  function( l, x, op, opt )
    return StabMC(GroupWithGenerators(l),x,op,opt);
  end );

InstallMethod( StabMC, "call Stab with errorbound",
  [IsGroup, IsObject, IsFunction, IsRecord],
  function( g, x, op, opt )
    if not(IsBound(opt.ErrorBound)) then
        opt.ErrorBound := 1/1025;
    fi;
    if not(IsBound(opt.DoEstimate)) then
        opt.DoEstimate := 1000;   # try for 1 second
    fi;
    return Stab(g,x,op,opt);
  end );


InstallGlobalFunction( BacktrackSearchStabilizerChainElement,
  function( S,     # a stabilizer chain
            P,     # a property function G -> {true,false}
            g,     # a group element
            pruner )
    # Let G be a group and U being the subgroup defined by S.
    # Does a backtrack search in U using S for an element x of U
    # for which P(xg) is true. g not necessarily in U!
    # pruner is either false or a function taking arguments as seen below.
    
    local cache,gen,i,ii,res,t,w,x,dotree;

    cache := WeakPointerObj([[S!.orb!.gens[1]^0,[]]]);
    for i in [1..Length(S!.orb)] do
        ii := S!.orb!.orbind[i];
        if ii > 1 then
            t := ElmWPObj(cache,S!.orb!.schreierpos[ii]);
            if t <> fail then
                gen := S!.orb!.schreiergen[ii];
                w := EmptyPlist(Length(t[2])+1);
                Append(w,t[2]);
                Add(w,gen);
                t := t[1] * S!.orb!.gens[gen];
            else
                w := TraceSchreierTreeForward(S!.orb,ii);
                t := Product(S!.orb!.gens{w});
            fi;
            SetElmWPObj(cache,ii,[t,w]);
            x := t*g;
        else
            w := [];
            t := S!.orb!.gens[1]^0;
            x := g;
        fi;
        if S!.stab = false then
            if P(x) then 
                return rec( elm := x, vec := [i], word := S!.layergens{w} );
            fi;
        else
            if pruner = false or pruner(S,ii,x,t,w) then
                res := BacktrackSearchStabilizerChainElement(S!.stab,P,x,
                                                             pruner);
                if res <> fail then
                    Add(res.vec,ii,1);
                    res.word := Concatenation(res.word,S!.layergens{w});
                    return res;
                fi;
            fi;
        fi;
    od;
    return fail;
  end );

InstallGlobalFunction( ComputeSuborbitsForStabilizerChain,
  function( S )
    local SS,gens;
    SS := S;
    while SS!.stab <> false do
        gens := StrongGenerators(S){SS!.stab!.layergens};
        SS!.suborbs := FindSuborbits(S!.orb,gens);
        SS := SS!.stab;
    od;
  end );

InstallGlobalFunction( BacktrackSearchStabilizerChainSubgroup,
  function(S,P,pruner) 
    # Let G be a group and U being the subgroup defined by S.
    # Does a backtrack search in U using S for the subgroup H of U
    # with H := {h \in U | P(h)=true}. P must be such that H is 
    # a subgroup of U.
    # pruner is either false or a function taking arguments as seen below.

    local SS,SSS,cache,dotree,gen,i,ii,newgens,o,res,t,w,T;

    cache := WeakPointerObj([[S!.orb!.gens[1]^0,[]]]);
    if S!.stab = false then    # lowest level
        newgens := [];
        o := false;
        for i in [2..Length(S!.orb)] do
            ii := S!.orb!.orbind[i];
            if o = false or not(S!.orb[ii] in o) then
                if ii > 1 then
                    t := ElmWPObj(cache,S!.orb!.schreierpos[ii]);
                    if t <> fail then
                        gen := S!.orb!.schreiergen[ii];
                        w := EmptyPlist(Length(t[2])+1);
                        Append(w,t[2]);
                        Add(w,gen);
                        t := t[1] * S!.orb!.gens[gen];
                    else
                        w := TraceSchreierTreeForward(S!.orb,ii);
                        t := Product(S!.orb!.gens{w});
                    fi;
                    SetElmWPObj(cache,ii,[t,w]);
                else
                    w := [];
                    t := S!.orb!.gens[1]^0;
                fi;
                if P(t) then
                    Add(newgens,t);
                    if Length(newgens) = 1 then
                        o := Orb(ShallowCopy(newgens),S!.orb[1],S!.orb!.op,
                                 rec( treehashsize := 200, schreier := true,
                                      log := true ));
                    else
                        AddGeneratorsToOrbit(o,[t]);
                    fi;
                    Enumerate(o);
                fi;
            fi;
        od;
        if o = false then   # No solution found
            newgens := [S!.orb!.gens[1]^0];
            o := Orb(newgens,S!.orb[1],S!.orb!.op,
                     rec(schreier := true, log := true));
            Enumerate(o);
        fi;
        SSS := rec( stab := false, size := Length(o), base := [o[1]],
                    opt := ShallowCopy(S!.opt), layer := S!.layer,
                    orb := o, stronggens := newgens, 
                    layergens := [1..Length(newgens)], randpool := [], 
                    IsOne := S!.opt.IsOne, proof := true,
                    parentS := false);
        Objectify( StabChainByOrbType, SSS );
        SSS!.opt.pr := ProductReplacer(SSS!.stronggens);
        SSS!.topS := SSS;   # this will be changed further up!
        o!.stabilizerchain := SSS;
        o!.gensi := List(o!.gens,x->x^-1);
        Info(InfoGenSS,2,"Layer ",SSS!.layer," completed, size ",Size(SSS));
        return SSS;
    else   # S!.stab <> false
        # First go down in tree, essentially trying coset rep 1 first:
        SS := BacktrackSearchStabilizerChainSubgroup(S!.stab,P,pruner);
        newgens := [];
        o := false;
        for i in [2..Length(S!.orb)] do
            ii := S!.orb!.orbind[i];
            if o = false or not(S!.orb[ii] in o) then
                if ii > 1 then
                    t := ElmWPObj(cache,S!.orb!.schreierpos[ii]);
                    if t <> fail then
                        gen := S!.orb!.schreiergen[ii];
                        w := EmptyPlist(Length(t[2])+1);
                        Append(w,t[2]);
                        Add(w,gen);
                        t := t[1] * S!.orb!.gens[gen];
                    else
                        w := TraceSchreierTreeForward(S!.orb,ii);
                        t := Product(S!.orb!.gens{w});
                    fi;
                    SetElmWPObj(cache,ii,[t,w]);
                else
                    w := [];
                    t := S!.orb!.gens[1]^0;
                fi;
                if pruner = false or pruner(S,ii,t,t,w) then
                    res := BacktrackSearchStabilizerChainElement(S!.stab,P,t,
                                                                 pruner);
                    if res <> fail then
                      Add(newgens,res.elm);
                      if Length(newgens) = 1 then
                        o := Orb(Concatenation(StrongGenerators(SS),newgens),
                                 S!.orb[1],S!.orb!.op,
                        rec( treehashsize := 
                                    Maximum(100,QuoInt(Length(S!.orb),3)), 
                             schreier := true, log := true ));
                      else
                        AddGeneratorsToOrbit(o,[res.elm]);
                      fi;
                      Enumerate(o);
                    fi;
                fi;
            fi;
        od;
        if o = false then   # No solution found
            o := Orb([S!.orb!.gens[1]^0],S!.orb[1],S!.orb!.op,
                     rec(schreier := true, log := true));
            Enumerate(o);
        fi;
        SSS := rec( stab := SS, size := Length(o)*Size(SS), 
                    base := SS!.base, opt := ShallowCopy(S!.opt), 
                    layer := S!.layer,
                    orb := o, stronggens := SS!.stronggens, 
                    layergens := [1..Length(SS!.stronggens)+Length(newgens)],
                    randpool := [], IsOne := SS!.opt.IsOne, proof := true,
                    parentS := false );
        Add(SSS.base,o[1],1);
        Append(SSS.stronggens,newgens);
        SSS.opt.pr := ProductReplacer(ShallowCopy(SSS.stronggens));
        Objectify( StabChainByOrbType, SSS );
        SS!.parentS := SSS;
        T := SSS;
        while T <> false do
            T!.topS := SSS;
            T := T!.stab;
        od;
        o!.stabilizerchain := SSS;
        o!.gensi := List(o!.gens,x->x^-1);
        Info(InfoGenSS,2,"Layer ",SSS!.layer," completed, size ",Size(SSS));
        return SSS;
    fi;
  end );

InstallMethod( FindShortGeneratorsOfSubgroup, 
  "without option rec or func, permutation group method",
  [ IsPermGroup, IsPermGroup ],
  function(G,U)
    local S; 
    S := StabilizerChain(U);
    return FindShortGeneratorsOfSubgroup(G,U,
             rec( membershiptest :=
                    function(a,b) 
                      return SiftGroupElement(S,a).isone; 
                    end,
                  sizetester :=
                    function(a) 
                      return Size(StabilizerChain(a)); 
                    end ));
  end );

InstallMethod( FindShortGeneratorsOfSubgroup, 
  "without option rec or func, matrix group method",
  [ IsMatrixGroup, IsMatrixGroup ],
  function(G,U)
    local S; 
    S := StabilizerChain(U);
    return FindShortGeneratorsOfSubgroup(G,U,
             rec( membershiptest :=
                    function(a,b) 
                      return SiftGroupElement(S,a).isone; 
                    end,
                  sizetester :=
                    function(a) 
                      return Size(StabilizerChain(a)); 
                    end ));
  end );

