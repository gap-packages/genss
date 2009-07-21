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

## InstallGlobalFunction( VerifyStabilizerChainTC2,
##   function( S )
##     local Grels,Hrels,MakeSchreierGens,Prels,ace,acenrGrels,
##           cosetnrlimitfactor,f,gens,gensi,i,j,k,l,li,max,newpres,
##           nrcosets,nrgens,nrschr,o,ords,pres,sb,sgs,slp,subgens,x;
##     if S!.stab <> false then
##         pres := VerifyStabilizerChainTC2(S!.stab);
##         if IsList(pres) then return pres; fi;
##     else
##         pres := StraightLineProgram([[]],0);
##     fi;
##     Info(InfoGenSS,1,"Verifying stabilizer chain in layer ",S!.layer);
##     # First create a few Schreier generators:
##     sgs := [];
##     i := 1;
##     j := 1;
##     o := S!.orb;
##     nrgens := Length(o!.gens);
##     sb := S!.strongbelow;
##     MakeSchreierGens := function(n)
##         local sg;
##         Info(InfoGenSS,3,"Creating ",n," Schreier generators...");
##         while Length(sgs) < n and
##               i <= Length(o) do
##             sg := GENSS_CreateSchreierGenerator(S,i,j);
##             j := j + 1;
##             if j > nrgens then
##                 j := 1;
##                 i := i + 1;
##             fi;
##             if sg <> fail then
##                 Add(sgs,sg);
##             fi;
##         od;
##     end;
## 
##     nrschr := S!.opt.NumberSchreierGens;
##     MakeSchreierGens(nrschr);
##     f := FreeGroup(S!.nrstrong);
##     gens := GeneratorsOfGroup(f);
##     gensi := List(gens,x->x^-1);
##     subgens := gens{[1..S!.strongbelow]};
##     Hrels := ResultOfStraightLineProgram(pres,subgens);
##     if S!.opt.Projective then
##         ords := List([1..nrgens],i->ProjectiveOrder(o!.gens[i]));
##     else
##         ords := List([1..nrgens],i->Order(o!.gens[i]));
##     fi;
##     Prels := List([1..nrgens],i->gens[i+sb]^ords[i]);
##     Grels := [];
##     cosetnrlimitfactor := 100;
##     ace := false;
##     while true do   # will be left by return eventually
##         for k in [Length(Grels)+1..Length(sgs)] do
##             Grels[k] := GENSS_Prod(gens,sgs[k][1]+sb) * gens[sgs[k][2]+sb] * 
##                         GENSS_Prod(gensi,sgs[k][3]+sb);
##             x := GENSS_Prod(o!.gens,sgs[k][1]) * o!.gens[sgs[k][2]] * 
##                  GENSS_Prod(o!.gensi,sgs[k][3]);
##             if S!.stab <> false then
##                 slp := SiftGroupElementSLP(S!.stab,x);
##                 if not(slp.isone) then
##                     if ace <> false then
##                         ACEQuit(ace);
##                     fi;
##                     return [fail,S!.layer];
##                 fi;
##                 Grels[k] := Grels[k] / ResultOfStraightLineProgram(slp.slp,
##                                                                    subgens);
##                 sgs[k][4] := slp.slp;
##             else
##                 if not(IsOne(x)) then
##                     return [fail,S!.layer];
##                 fi;
##                 sgs[k][4] := false;
##             fi;
##         od;
##         Info(InfoGenSS,2,"Doing ACE coset enumeration with limit ",
##              cosetnrlimitfactor*Length(o)," and ",Length(Hrels),
##              "+",Length(Prels),"+",Length(Grels)," relations...");
##         if ace = false then
##             ace := ACEStart(gens,Concatenation(Hrels,Prels,Grels),subgens:
##                             max := cosetnrlimitfactor * Length(o) );
##         else
##             ACEAddRelators(ace,Grels{[acenrGrels+1..Length(Grels)]});
##         fi;
##         acenrGrels := Length(Grels);
##         if not(IsCompleteACECosetTable(ace))  then   # did not close!
##             #cosetnrlimitfactor := QuoInt(cosetnrlimitfactor*3,2);
##             Info(InfoGenSS,2,"Coset enumeration did not finish!");
##             if nrschr > Length(sgs) # or
##                # nrschr > S!.opt.MaxNumberSchreierGens 
##                then   # we are done!
##                 # Something is wrong!
##                 Error();
##                 ACEQuit(ace);
##                 return [fail, S!.layer];
##             fi;
##         else
##             nrcosets := ACEStats(ace).index;
##             Info(InfoGenSS,2,"Coset enumeration found ",nrcosets," cosets.");
##             if nrcosets = Length(o) then
##                 ACEQuit(ace);
##                 # Verification is OK, now build a presentation:
##                 l := GeneratorsWithMemory(
##                        ListWithIdenticalEntries(S!.nrstrong,()));
##                 li := List(l,x->x^-1);
##                 newpres := ResultOfStraightLineProgram(pres,
##                                    l{[1..S!.strongbelow]});
##                 for k in [1..nrgens] do
##                     Add(newpres,l[k+sb]^ords[k]);
##                 od;
##                 for k in [1..Length(sgs)] do
##                     if sgs[k][4] <> false then
##                         Add(newpres,
##                             GENSS_Prod(l,sgs[k][1]+sb)*l[sgs[k][2]+sb]*
##                             GENSS_Prod(li,sgs[k][3]+sb)*
##                             ResultOfStraightLineProgram(sgs[k][4],l)^-1);
##                     else
##                         Add(newpres,GENSS_Prod(l,sgs[k][1]+sb)*l[sgs[k][2]+sb]*
##                                     GENSS_Prod(li,sgs[k][3]+sb));
##                     fi;
##                 od;
##                 Info(InfoGenSS,2,"Found presentation for layer ",S!.layer,
##                      " using ",Length(newpres)," relators.");
##                 return SLPOfElms(newpres);
##             fi;
##         fi;
##         # nrschr := nrschr + S!.opt.NumberSchreierGens;
##         nrschr := QuoInt(nrschr*4,3);
##         MakeSchreierGens(nrschr);
##     od;
##   end);

## VerifyStabilizerChainTC3 := 
##   function( S )
##     local Grels,Hrels,MakeSchreierGens,Prels,TestGenerationSGens,ace,
##           acenrGrels,cosetnrlimitfactor,dr,f,gens,gensi,i,j,k,l,li,max,
##           newpres,nrcosets,nrgens,nrschr,o,ok,ords,pres,sb,sgs,slp,subgens,x;
##     if S!.stab <> false then
##         pres := VerifyStabilizerChainTC3(S!.stab);
##         if IsList(pres) then return pres; fi;
##     else
##         pres := StraightLineProgram([[]],0);
##     fi;
##     Info(InfoGenSS,1,"Verifying stabilizer chain in layer ",S!.layer);
##     # First create a few Schreier generators:
##     sgs := [];
##     i := 1;
##     j := 1;
##     o := S!.orb;
##     nrgens := Length(o!.gens);
##     sb := S!.strongbelow;
##     MakeSchreierGens := function(n)
##         local sg;
##         Info(InfoGenSS,3,"Creating ",n," Schreier generators...");
##         while Length(sgs) < n and
##               i <= Length(o) do
##             sg := GENSS_CreateSchreierGenerator(S,i,j);
##             j := j + 1;
##             if j > nrgens then
##                 j := 1;
##                 i := i + 1;
##             fi;
##             if sg <> fail then
##                 Add(sgs,sg);
##             fi;
##         od;
##     end;
##     TestGenerationSGens := function()
##         local SS,g,k,l;
##         l := [];
##         for k in [1..Length(sgs)] do
##             Add(l,GENSS_Prod(o!.gens,sgs[k][1]) * o!.gens[sgs[k][2]] * 
##                   GENSS_Prod(o!.gensi,sgs[k][3]));
##         od;
##         g := Group(l);
##         SS := StabilizerChain(g,rec( Base := S!.stab, 
##                                      ErrorBound := 1/3,
##                                      Projective := S!.opt.Projective ));
##         Info( InfoGenSS,3,"Schreier gen test gave ",Size(SS)," of ",
##               Size(S!.stab) );
##         return Size(SS) = Size(S!.stab);
##     end;        
##             
##     nrschr := S!.opt.NumberSchreierGens;
##     while true do;
##         MakeSchreierGens(nrschr);
##         ok := S!.stab = false or TestGenerationSGens();
##         if ok then break; fi;
##         if i > Length(o) then
##             Error( "this cannot have happenend" );
##         fi;
##         nrschr := QuoInt(nrschr*4,3);
##     od;
##     # We are now pretty sure that the Schreier generators generate the
##     # stabilizer, just make some more and then start coset enumeration:
##     nrschr := nrschr * 2;
##     MakeSchreierGens(nrschr);
## 
##     f := FreeGroup(S!.nrstrong);
##     gens := GeneratorsOfGroup(f);
##     gensi := List(gens,x->x^-1);
##     subgens := gens{[1..S!.strongbelow]};
##     Hrels := ResultOfStraightLineProgram(pres,subgens);
##     if S!.opt.Projective then
##         ords := List([1..nrgens],i->ProjectiveOrder(o!.gens[i]));
##     else
##         ords := List([1..nrgens],i->Order(o!.gens[i]));
##     fi;
##     Prels := List([1..nrgens],i->gens[i+sb]^ords[i]);
##     Grels := [];
##     cosetnrlimitfactor := 4;
##     ace := false;
##     # Xrels := [];
##     while true do   # will be left by return eventually
##         for k in [Length(Grels)+1..Length(sgs)] do
##             Grels[k] := GENSS_Prod(gens,sgs[k][1]+sb) * gens[sgs[k][2]+sb] * 
##                         GENSS_Prod(gensi,sgs[k][3]+sb);
##             x := GENSS_Prod(o!.gens,sgs[k][1]) * o!.gens[sgs[k][2]] * 
##                  GENSS_Prod(o!.gensi,sgs[k][3]);
##             if S!.stab <> false then
##                 slp := SiftGroupElementSLP(S!.stab,x);
##                 if not(slp.isone) then
##                     if ace <> false then
##                         ACEQuit(ace);
##                     fi;
##                     return [fail,S!.layer];
##                 fi;
##                 Grels[k] := Grels[k] / ResultOfStraightLineProgram(slp.slp,
##                                                                    subgens);
##                 sgs[k][4] := slp.slp;
##             else
##                 if not(IsOne(x)) then
##                     return [fail,S!.layer];
##                 fi;
##                 sgs[k][4] := false;
##             fi;
##         od;
##         Info(InfoGenSS,2,"Doing ACE coset enumeration with limit ",
##              cosetnrlimitfactor*Length(o)," and ",Length(Hrels),
##              "+",Length(Prels),"+",Length(Grels)," relations...");
##         if ace = false then
##             ace := ACEStart(gens,Concatenation(Hrels,Prels,Grels),subgens:
##                             max := cosetnrlimitfactor * Length(o), hlt );
##             dr := ACEDataRecord(ace);
##             nrcosets := dr.stats.index;
##         else
##             # k := Random([Length(o)+1..nrcosts]).
##             if Length(Grels) > acenrGrels then
##                 ACEAddRelators(ace,Grels{[acenrGrels+1..Length(Grels)]}:
##                                max := cosetnrlimitfactor * Length(o), hlt );
##             else
##                 ACEContinue(ace : max := cosetnrlimitfactor * Length(o));
##             fi;
##             dr := ACEDataRecord(ace);
##             nrcosets := dr.stats.index;
##         fi;
##         acenrGrels := Length(Grels);
##         if not(IsCompleteACECosetTable(ace))  then   # did not close!
##             cosetnrlimitfactor := QuoInt(cosetnrlimitfactor*3,2);
##             Info(InfoGenSS,2,"Coset enumeration did not finish!");
##             #if nrschr > Length(sgs) # or
##             #   # nrschr > S!.opt.MaxNumberSchreierGens 
##             #   then   # we are done!
##             #    # Something is wrong!
##             #    Error("wrong1");
##             #    ACEQuit(ace);
##             #    return [fail, S!.layer];
##             #fi;
##         else
##             Info(InfoGenSS,2,"Coset enumeration found ",nrcosets," cosets.");
##             if nrcosets = Length(o) then
##                 ACEQuit(ace);
##                 # Verification is OK, now build a presentation:
##                 l := GeneratorsWithMemory(
##                        ListWithIdenticalEntries(S!.nrstrong,()));
##                 li := List(l,x->x^-1);
##                 newpres := ResultOfStraightLineProgram(pres,
##                                    l{[1..S!.strongbelow]});
##                 for k in [1..nrgens] do
##                     Add(newpres,l[k+sb]^ords[k]);
##                 od;
##                 for k in [1..Length(sgs)] do
##                     if sgs[k][4] <> false then
##                         Add(newpres,
##                             GENSS_Prod(l,sgs[k][1]+sb)*l[sgs[k][2]+sb]*
##                             GENSS_Prod(li,sgs[k][3]+sb)*
##                             ResultOfStraightLineProgram(sgs[k][4],l)^-1);
##                     else
##                         Add(newpres,GENSS_Prod(l,sgs[k][1]+sb)*l[sgs[k][2]+sb]*
##                                     GENSS_Prod(li,sgs[k][3]+sb));
##                     fi;
##                 od;
##                 Info(InfoGenSS,2,"Found presentation for layer ",S!.layer,
##                      " using ",Length(newpres)," relators.");
##                 return SLPOfElms(newpres);
##             fi;
##         fi;
##         # nrschr := nrschr + S!.opt.NumberSchreierGens;
##         nrschr := QuoInt(nrschr*4,3);
##         MakeSchreierGens(nrschr);
##     od;
##   end;


## VerifyStabilizerChainTC4 := 
##   function( S )
##     local Grels,Hrels,MakeSchreierGen,Prels,ace,cosetlimit,done,dr,f,
##           gens,gensi,guck1,guck2,hlt,i,j,k,l,li,max,newpres,nrcosets,
##           nrgens,o,ords,pres,sb,sg,sgs,slp,st,subgens,x,y,y1,y2;
## 
##     if S!.stab <> false then
##         pres := VerifyStabilizerChainTC4(S!.stab);
##         if IsList(pres) then return pres; fi;
##     else
##         pres := StraightLineProgram([[]],0);
##     fi;
##     Info(InfoGenSS,1,"Verifying stabilizer chain in layer ",S!.layer);
## 
##     # The following are global to "MakeSchreierGen":
##     i := 1;
##     j := 1;
##     MakeSchreierGen := function()
##         local sg;
##         Info(InfoGenSS,4,"Creating Schreier generator... i=",i," j=",j);
##         while i <= Length(o) do
##             sg := GENSS_CreateSchreierGenerator(S,i,j);
##             j := j + 1;
##             if j > nrgens then
##                 j := 1;
##                 i := i + 1;
##             fi;
##             if sg <> fail then return sg; fi;
##         od;
##         return fail;
##     end;
## 
##     o := S!.orb;
##     nrgens := Length(o!.gens);
##     sb := S!.strongbelow;
##             
##     if S!.nrstrong > 26 then
##         f := FreeGroup(List([1..S!.nrstrong],String));
##     else
##         gens := [];
##         st := "a";
##         for k in [1..S!.nrstrong] do
##             Add(gens,ShallowCopy(st));
##             st[1] := CHAR_INT(INT_CHAR(st[1])+1);
##         od;
##         f := FreeGroup(gens);
##     fi;
##     gens := GeneratorsOfGroup(f);
##     gensi := List(gens,x->x^-1);
##     subgens := gens{[1..S!.strongbelow]};
##     Hrels := ResultOfStraightLineProgram(pres,subgens);
##     if S!.opt.Projective then
##         ords := List([1..nrgens],i->ProjectiveOrder(o!.gens[i]));
##     else
##         ords := List([1..nrgens],i->Order(o!.gens[i]));
##     fi;
##     Prels := List([1..nrgens],i->gens[i+sb]^ords[i]);
##     sgs := [];
##     Grels := [];
## 
##     # Now start up a coset enumeration:
##     cosetlimit := QuoInt(5 * Length(o),4);
##     Info(InfoGenSS,2,"Starting ACE coset enumeration with limit ",
##          cosetlimit," and ",Length(Hrels),
##          "+",Length(Prels),"+",Length(Grels)," relations...");
##     ace := ACEStart(gens,Concatenation(Hrels,Prels,Grels),subgens:
##                     max := cosetlimit, hlt := true );
##     done := IsCompleteACECosetTable(ace);
##     if done then
##         dr := ACEDataRecord(ace);
##         nrcosets := dr.stats.index;
##     fi;
##         
##     while true do   # will be left by return eventually
##         if done then
##             Info(InfoGenSS,2,"Coset enumeration found ",nrcosets," cosets.");
##             if nrcosets = Length(o) then
##                 ACEQuit(ace);
##                 # Verification is OK, now build a presentation:
##                 l := GeneratorsWithMemory(
##                        ListWithIdenticalEntries(S!.nrstrong,()));
##                 li := List(l,x->x^-1);
##                 newpres := ResultOfStraightLineProgram(pres,
##                                    l{[1..S!.strongbelow]});
##                 for k in [1..nrgens] do
##                     Add(newpres,l[k+sb]^ords[k]);
##                 od;
##                 for k in [1..Length(sgs)] do
##                     if sgs[k][4] <> false then
##                         Add(newpres,
##                             GENSS_Prod(l,sgs[k][1]+sb)*l[sgs[k][2]+sb]*
##                             GENSS_Prod(li,sgs[k][3]+sb)*
##                             ResultOfStraightLineProgram(sgs[k][4],l)^-1);
##                     else
##                         Add(newpres,GENSS_Prod(l,sgs[k][1]+sb)*l[sgs[k][2]+sb]*
##                                     GENSS_Prod(li,sgs[k][3]+sb));
##                     fi;
##                 od;
##                 Info(InfoGenSS,2,"Found presentation for layer ",S!.layer,
##                      " using ",Length(newpres)," relators.");
##                 return SLPOfElms(newpres);
##             elif nrcosets < Length(o) then
##                 Error("This cannot possible have happened!");
##             else   # nrcosets > Length(o)
##                 Info(InfoGenSS,2,"Too many cosets, we must have forgotten ",
##                      "another relation!");
##                 done := false;
##             fi;
##         fi;
##         if not(done) then
##             while true do    # will be left by break
##                 sg := MakeSchreierGen();
##                 if sg = fail then
##                     Error("Something wrong, have processed all Schreier gens");
##                 fi;
##                 # Sift residue:
##                 x := GENSS_Prod(o!.gens,sg[1]) * o!.gens[sg[2]] * 
##                      GENSS_Prod(o!.gensi,sg[3]);
##                 if not(IsOne(x)) then
##                     if S!.stab <> false then
##                         slp := SiftGroupElementSLP(S!.stab,x);
##                         if not(slp.isone) then
##                             ACEQuit(ace);
##                             return [fail,S!.layer,x];
##                         fi;
##                     fi;
##                     sg[4] := slp.slp;
##                 else
##                     sg[4] := false;
##                 fi;
##                 y1 := GENSS_Prod(gens,sg[1]+sb) * gens[sg[2]+sb];
##                 y2 := GENSS_Prod(gens,Reversed(sg[3]+sb));
##                 # Now check with ACE:
##                 guck1 := ACETraceWord(ace,1,y1);
##                 guck2 := ACETraceWord(ace,1,y2);
##                 if guck1 = fail or guck2 = fail or guck1 <> guck2 then
##                     y := y1/y2;
##                     if sg[4] <> false then
##                         y := y / ResultOfStraightLineProgram(sg[4],subgens);
##                     fi;
##                     Add(sgs,sg);
##                     Add(Grels,y);
##                     break;
##                 fi;
##                 # Otherwise go to next Schreier generator.
##             od;
##             Info(InfoGenSS,2,"Redoing ACE coset enumeration with limit ",
##                  cosetlimit," and ",Length(Hrels),
##                  "+",Length(Prels),"+",Length(Grels)," relations...");
##             ACEAddRelators(ace,[y]);
##             done := IsCompleteACECosetTable(ace);
##             if done then
##                 dr := ACEDataRecord(ace);
##                 nrcosets := dr.stats.index;
##             fi;
##         fi;
##     od;
##   end;

