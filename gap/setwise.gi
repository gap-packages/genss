#############################################################################
##
##  setwise.gi              genss package           
##                                                           Max Neunhoeffer
##                                                              Felix Noeske
##
##  Copyright 2008 by the authors
##
##  Implementation stuff for setwise stabilizer code.
##
#############################################################################

InstallGlobalFunction( GENSS_FindElmMappingBaseSubsetIntoSet,
  function( G, op, S, M, N, NN, map, depth )
    # G is a group, op an action, S a stabilizer chain for G
    # M and N are sets of points in the domain of the action.
    # The base points of S must have been selected by successively
    # using elements of M. If some stabilizer has orbit length 1 on
    # the next point of M, then it is omitted. After M is exhausted
    # arbitrary new base points can be used in S. NN is a subset of N.
    # map is a list and depth gives the recursion depth.
    # Mapping point M[depth] to x is documented by map[depth] := x.
    # This function recursively searches for some element of G
    # that maps all points of M into N and the first basis point into
    # NN. It directly returns as soon as it has found one or returns
    # fail if none exists.
    # It returns the list map, filled with a (possibly partial) images
    # of the points in M.
    local Nti,i,j,orbpos,res,word;

    if depth > Length(M) then
        # Done!
        map[depth] := true;   # marks the end
        return map;
    fi;

    if S = false then   # stabchain exhausted, what do we do?
        # the only element left is the identity
        if M[depth] in NN and ForAll(M{[depth+1..Length(M)]},x -> x in N) then
            map[depth] := false;   # marks the end
            return map;
        else
            return fail;
        fi;
    fi;

    if S!.orb[1] = M[depth] then   # next point of M is in fact in next orbit
        for i in [1..Length(NN)] do
            orbpos := Position(S!.orb,NN[i]);
            if orbpos <> fail then
                map[depth] := NN[i];
                # Find transversal element from Schreier tree
                # trace N\{NN[i]} backwards through transversal element
                word := TraceSchreierTreeBack(S!.orb,orbpos);
                Nti := Filtered(N,x->x <> NN[i]);
                for j in [1..Length(word)] do
                    Nti := List(Nti,x->S!.orb!.op(x,S!.orb!.gensi[word[j]]));
                od;
                Sort(Nti);
                res := GENSS_FindElmMappingBaseSubsetIntoSet(
                   G, op, S!.stab, M, Nti, Nti, map, depth+1 );
                if res <> fail then return res; fi;
            fi;
        od;
    else
        # In this case we know that M[depth] is fixed by the pointwise
        # stabilizer of M[1..depth-1], thus we only have to check if
        # M[depth] is in NN or not.
        if M[depth] in NN then
            map[depth] := M[depth];
            res := GENSS_FindElmMappingBaseSubsetIntoSet(
              G, op, S, M, N, N, map, depth+1 );
            if res <> fail then return res; fi;
        fi;
    fi;
    return fail;
  end);

InstallGlobalFunction( GENSS_ComputeFoundElm,
function( G, op, S, M, map, depth )
  # See above, to be called upon completion to actually make the group
  # element found.
  local el,orbpos,word;
  if map = fail then
      return fail;
  elif depth > Length(M) then
      return One(G);
  elif S = false then
      return One(G);
  fi;
  if S!.orb[1] = M[depth] then
      orbpos := Position(S!.orb,map[depth]);
      word := TraceSchreierTreeForward(S!.orb,orbpos);
      if Length(word) = 0 then
          el := One(S!.orb!.gens[1]);
      else
          el := Product(S!.orb!.gens{word});
      fi;
      return GENSS_ComputeFoundElm( G, op, S!.stab, M, map, depth+1 ) * el;
  else
      return GENSS_ComputeFoundElm( G, op, S, M, map, depth+1 );
  fi;
end);

        
InstallMethod( SetwiseStabilizer, 
  "generic method",
  [IsGroup,IsFunction,IsList],
function( G, op, M )
  # G is a group, op an action, and M a set of points G is acting upon
  # with the action op. Returns a record with the following
  # components:
  #   setstab   : the setwise stabilizer
  #   S         : a stabilizer chain for G using (parts of) M as base
  local FindGens,S,SS,gens,i,m;
  
  Info( InfoGenSS, 1, "Computing adjusted stabilizer chain..." );
  S := StabilizerChain( G, 
          rec( Cand := rec( points := ShallowCopy(M),
                            ops := ListWithIdenticalEntries( Length(M), op ), 
                            used := 0 ),
               StrictlyUseCandidates := true ));
  gens := [];
  SS := S;
  for i in [1..Length(M)] do
      if SS = false then break; fi;
      if M[i] = SS!.orb[1] then SS := SS!.stab; fi;
  od;
  while SS <> false do
      Append(gens, SS!.orb!.gens);
      SS := SS!.stab;
  od;
  # These are now the generators of the pointwise stabilizer!
  # We add to them as we go.

  m := Length(M);

  FindGens := function(S)
      local MM,g,map,newgens,newgensp,o,tolook,x,depth;
      depth := 0;
      while true do  # this loop never terminates but there are "return" statements
          depth := depth + 1;
          # In depth <depth> this finds gens of the setwise stabilizer
          # that fix all points M{[1..depth-1]} and move M[depth] and determine
          # which points in M{[depth+1..Length(M)]} can be reached in this
          # way.
          if depth = m or S = false then return; fi;
          if M[depth] <> S!.orb[1] then
              continue;
          fi;
          o := [M[depth]];
          tolook := M{[depth+1..m]};
          newgens := [];
          g := Group(S!.orb!.gens);
          MM := M{[depth..m]};
          while true do
              map := GENSS_FindElmMappingBaseSubsetIntoSet(
                   g, S!.orb!.op, S, MM, MM, tolook, [], 1 );
              if map = fail then 
                  Info( InfoGenSS, 2, "No more to be found in depth ",depth );
                  break;  # no more to be found
              fi;
              # in SetwiseStabilizer
              x := GENSS_ComputeFoundElm( g, S!.orb!.op, S, MM, map, 1 );
              Info( InfoGenSS, 2, "Found new generator in depth ", depth );
              Add(newgens,x);
              Add(gens,x);   # make gen known outside
              o := Orb(newgens,M[depth],S!.orb!.op,rec( storenumbers := true ));
              Enumerate(o);
              tolook := Filtered(MM,x->not(x in o));
              if Length(tolook) = 0 then   # we have everything!
                  newgensp := ActionOnOrbit(o,newgens);
                  if IsNaturalSymmetricGroup(Group(newgensp)) then
                      # We now know that everything has been found!
                      Info( InfoGenSS, 2, "Found Sym in depth ",depth );
                      return;
                  fi;
                  break;
              fi;
          od;
          S := S!.stab;
      od;
  end;
  FindGens(S);
  if Length(gens) = 0 then Add(gens,One(G)); fi;
  return rec( setstab := GroupWithGenerators(gens), S := S );
end);

##  SetwiseStabilizer2 := function( G, op, M )
##    # G is a group, op an action, and M a set of points G is acting upon
##    # with the action op. Returns a record with the following
##    # components:
##    #   setstab   : the setwise stabilizer
##    #   S         : a stabilizer chain for G using (parts of) M as base
##    local FindGensInDepth,S,SS,gens,i,m;
##    
##    Info( InfoGenSS, 1, "Computing adjusted stabilizer chain..." );
##    S := StabilizerChain( G, 
##            rec( Cand := rec( points := ShallowCopy(M),
##                              ops := ListWithIdenticalEntries( Length(M), op ), 
##                              used := 0 ),
##                 StrictlyUseCandidates := true ));
##    gens := [];
##    SS := S;
##    for i in [1..Length(M)] do
##        if SS = false then break; fi;
##        if M[i] = SS!.orb[1] then SS := SS!.stab; fi;
##    od;
##    if SS = false then
##        gens := [];
##    else
##        gens := ShallowCopy(SS!.orb!.gens);
##    fi;
##    # These are now the generators of the pointwise stabilizer!
##    # We add to them as we go.
##  
##    m := Length(M);
##  
##    # We now proceed recursively (actually, this is tail recursion) again:
##    FindGensInDepth := function(S,depth)
##        local MM,done,g,map,newgens,newgensp,o,tolook,x;
##        # In depth <depth> this finds gens of the setwise stabilizer
##        # that fix all points M{[1..depth-1]} and move M[depth] and determine
##        # which points in M{[depth+1..Length(M)]} can be reached in this
##        # way.
##        if depth = m or S = false then return; fi;
##        if M[depth] <> S!.orb[1] then
##            FindGensInDepth(S,depth+1);
##            return;
##        fi;
##        FindGensInDepth(S!.stab,depth+1);
##        o := [M[depth]];
##        tolook := M{[depth+1..m]};
##        newgens := [];
##        g := Group(S!.orb!.gens);
##        MM := M{[depth..m]};
##        done := false;
##        repeat  # this loop never terminates but there are "return" statements
##            map := GENSS_FindElmMappingBaseSubsetIntoSet(
##                 g, S!.orb!.op, S, MM, MM, tolook, [], 1 );
##            if map <> fail then 
##                x := GENSS_ComputeFoundElm( g, S!.orb!.op, S, MM, map, 1 );
##                Info( InfoGenSS, 2, "Found new generator in depth ", depth );
##                Add(newgens,x);
##                Add(gens,x);   # make gen known outside
##                o := Orb(newgens,M[depth],S!.orb!.op,rec( storenumbers := true ));
##                Enumerate(o);
##                tolook := Filtered(MM,x->not(x in o));
##                #if Length(tolook) = 0 then   # we have everything!
##                #    done := true;
##                #    newgensp := ActionOnOrbit(o,newgens);
##                #    if IsNaturalSymmetricGroup(Group(newgensp)) then
##                #        # We now know that everything has been found!
##                #        Info( InfoGenSS, 2, "Found Sym in depth ",depth );
##                #        return;
##                #    fi;
##                #fi;
##            else
##                Info( InfoGenSS, 2, "No more to be found in depth ",depth );
##                done := true;  # no more to be found
##            fi;
##        until done;
##    end;
##    FindGensInDepth(S,1);
##    if Length(gens) = 0 then Add(gens,One(G)); fi;
##    return rec( setstab := GroupWithGenerators(gens), S := S );
##  end;

#st := "                                                          ";
InstallGlobalFunction(
  GENSS_FindElmMappingBaseSubsetIntoSetPartitionBacktrack,
  function( G, op, S, M, N, NN, map, depth )
    # G is a group, op an action, S a stabilizer chain for G
    # M and N are sets of points in the domain of the action.
    # The base points of S must have been selected by successively
    # using elements of M. If some stabilizer has orbit length 1 on
    # the next point of M, then it is omitted. After M is exhausted
    # arbitrary new base points can be used in S. NN is a subset of N.
    # map is a list and depth gives the recursion depth.
    # Mapping point M[depth] to x is documented by map[depth] := x.
    # This function recursively searches for some element of G
    # that maps all points of M into N and the first basis point into
    # NN. It directly returns as soon as it has found one or returns
    # fail if none exists.
    # It returns the list map, filled with a (possibly partial) images
    # of the points in M.
    local cellsM,cellsN,Nti,i,j,orbpos,pa,po,res,word,x;

    if depth > Length(M) then
        # Done!
        map[depth] := true;   # marks the end
        #Print("SUCCESS!\n");
        return map;
    fi;

    if S = false then   # stabchain exhausted, what do we do?
        # the only element left is the identity
        if M[depth] in NN and ForAll(M{[depth+1..Length(M)]},x -> x in N) then
            map[depth] := false;   # marks the end
            return map;
        else
            return fail;
        fi;
    fi;

    if IsBound( S!.mainorb ) then
        # We think about pruning the tree here, we check, whether the
        # cell distribution in M{[depth..]} is the same as in NN:
        # If not, this subtree can not be successful any more:
        pa := S!.parts;
        cellsM := S!.cellsM{[depth..Length(M)]};
        #cellsM := [];
        cellsN := [];
        #for i in [depth..Length(M)] do
        #    Add(cellsM,pa[Position(S!.mainorb,M[i])]);
        #od;
        Sort(cellsM);
        for i in [1..Length(N)] do
            x := pa[Position(S!.mainorb,N[i])];
            if x < cellsM[1] or x > cellsM[Length(cellsM)] then break; fi;
            Add(cellsN,x);
        od;
        Sort(cellsN);
        #if cellsM <> cellsN then Print(st{[1..depth]},"cells\n"); fi;
        #Error();
        if cellsM <> cellsN then return fail; fi;
    fi;

    if S!.orb[1] = M[depth] then   # next point of M is in fact in next orbit
        for i in [1..Length(NN)] do
            orbpos := Position(S!.orb,NN[i]);
            if orbpos <> fail then
                map[depth] := NN[i];
                #Print(st{[1..depth]},
                #      "Trying ",M[depth]," :-> ",NN[i],"\n");
                # Find transversal element from Schreier tree
                # trace N\{NN[i]} backwards through transversal element
                word := TraceSchreierTreeBack(S!.orb,orbpos);
                Nti := Filtered(N,x->x <> NN[i]);
                for j in [1..Length(word)] do
                    Nti := List(Nti,x->S!.orb!.op(x,S!.orb!.gensi[word[j]]));
                od;
                Sort(Nti);
                res := GENSS_FindElmMappingBaseSubsetIntoSetPartitionBacktrack(
                   G, op, S!.stab, M, Nti, Nti, map, depth+1 );
                if res <> fail then return res; fi;
                #Print(st{[1..depth]},"Back\n");
            fi;
        od;
    else
        # In this case we know that M[depth] is fixed by the pointwise
        # stabilizer of M[1..depth-1], thus we only have to check if
        # M[depth] is in NN or not.
        if M[depth] in NN then
            map[depth] := M[depth];
            Nti := Filtered(N,x->x <> M[depth]);
            res := GENSS_FindElmMappingBaseSubsetIntoSetPartitionBacktrack(
              G, op, S, M, Nti, Nti, map, depth+1 );
            if res <> fail then return res; fi;
        fi;
    fi;
    return fail;
  end);

GENSS_SETSIZELIMITPARTITIONS := 5;

##  InstallGlobalFunction( GENSS_FindSuborbitsLocal,
##  function( o, gens )
##    local i,j,k,l,oo,parts;
##    l := Length(o);
##    parts := 0*[1..l];
##    i := 1;
##    j := 1;
##    while i <= Length(parts) do
##        if parts[i] = 0 then 
##            oo := Orb(gens,o[i],o!.op);
##            Enumerate(oo);
##            for k in [1..Length(oo)] do
##                parts[Position(o,oo[k])] := j;
##            od;
##            j := j + 1;
##        fi;
##        i := i + 1;
##    od;
##    return rec( suborbnr := parts );
##  end);

InstallMethod( SetwiseStabilizerPartitionBacktrack,
  "generic method",
  [IsGroup,IsFunction,IsList],
function( G, op, M )
  # G is a group, op an action, and M a set of points G is acting upon
  # with the action op. Returns a record with the following
  # components:
  #   setstab   : the setwise stabilizer
  #   S         : a stabilizer chain for G using (parts of) M as base
  # We assume that all of M lies in the orbit of the first point in M.
  local FindGens,S,SS,gens,i,m,r,rt,rt2;
  
  rt := Runtime();
  Info( InfoGenSS, 1, "Computing adjusted stabilizer chain..." );
  S := StabilizerChain( G, 
          rec( Cand := rec( points := ShallowCopy(M),
                            ops := ListWithIdenticalEntries( Length(M), op ), 
                            used := 0 ),
               StrictlyUseCandidates := true ));
  Info( InfoGenSS, 1, "Time so far: ",Runtime()-rt);
  rt2 := Runtime();

  # If M is big enough, find the suborbit partitions for the stabilizers:
  if Length(M) >= GENSS_SETSIZELIMITPARTITIONS then
      SS := S!.stab;
      for i in [2..Length(M)] do
          if SS = false then break; fi;
          if M[i] = SS!.orb[1] then
              r := FindSuborbits(S!.orb,SS!.orb!.gens);
              SS!.mainorb := S!.orb;
              SS!.parts := r.suborbnr;
              Unbind(r);
              SS!.cellsM := List(M,x->SS!.parts[Position(S!.orb,x)]);
              SS := SS!.stab;
          fi;
      od;
  fi;

  Info( InfoGenSS, 1, "Time so far: ",Runtime()-rt," Lap: ",Runtime()-rt2);
  rt2 := Runtime();

  gens := [];
  SS := S;
  for i in [1..Length(M)] do
      if SS = false then break; fi;
      if M[i] = SS!.orb[1] then SS := SS!.stab; fi;
  od;
  while SS <> false do
      Append(gens, SS!.orb!.gens);
      SS := SS!.stab;
  od;
  # These are now the generators of the pointwise stabilizer!
  # We add to them as we go.

  m := Length(M);

  FindGens := function(S)
      local MM,g,map,newgens,newgensp,o,tolook,x,depth;
      depth := 0;
      while true do  # this loop never terminates but there are "return" statements
          depth := depth + 1;
          # In depth <depth> this finds gens of the setwise stabilizer
          # that fix all points M{[1..depth-1]} and move M[depth] and determine
          # which points in M{[depth+1..Length(M)]} can be reached in this
          # way.
          if depth = m or S = false then return; fi;
          if M[depth] <> S!.orb[1] then
              continue;
          fi;
          o := [M[depth]];
          tolook := M{[depth+1..m]};
          newgens := [];
          g := Group(S!.orb!.gens);
          MM := M{[depth..m]};
          while true do
              map := GENSS_FindElmMappingBaseSubsetIntoSetPartitionBacktrack(
                   g, S!.orb!.op, S, M, MM, tolook, M{[1..depth-1]}, depth );
              if map = fail then 
                  Info( InfoGenSS, 2, "No more to be found in depth ",depth );
                  break;  # no more to be found
              fi;
              # in SetwiseStabilizerPartitionBacktrack
              x := GENSS_ComputeFoundElm( g, S!.orb!.op, S, M, map, depth );
              Info( InfoGenSS, 2, "Found new generator in depth ", depth );
              Add(newgens,x);
              Add(gens,x);   # make gen known outside
              o := Orb(newgens,M[depth],S!.orb!.op,rec( storenumbers := true ));
              Enumerate(o);
              tolook := Filtered(MM,x->not(x in o));
              if Length(tolook) = 0 then   # we have everything!
                  newgensp := ActionOnOrbit(o,newgens);
                  if IsNaturalSymmetricGroup(Group(newgensp)) then
                      # We now know that everything has been found!
                      Info( InfoGenSS, 2, "Found Sym in depth ",depth );
                      return;
                  fi;
                  break;
              fi;
          od;
          S := S!.stab;
      od;
  end;
  FindGens(S);
  if Length(gens) = 0 then Add(gens,One(G)); fi;

  Info( InfoGenSS, 1, "Time so far: ",Runtime()-rt," Lap: ",Runtime()-rt2);
  return rec( setstab := GroupWithGenerators(gens), S := S );
end);


