GENSS.SchreierPatience := 200;
GENSS.SchreierPatienceBaseImage := infinity;
GENSS.VerificationElements := 10;  # this makes for error prob. 1/1024
GENSS.NumberSchreierGensWait := 5;

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

InstallGlobalFunction( GENSS_CompleteStabilizerChain,
  function( S )
    while S!.stab <> false do S := S!.stab; od;
    if Length(S!.stabgens) > 0 then
        GENSS_StartNewOrbitEnumeration(S);
    fi;
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
