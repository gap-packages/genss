# Test file for the genss package:

SetInfoLevel(InfoOrb,0);
SetInfoLevel(InfoGenSS,0);

errors := [];

TestGroup := function(g,size,opt)
  local S,i,pr,slp,x,y,success;
  success := "+"; 
  if opt = false then opt := rec(); fi;
  S := StabilizerChain(g,opt);
  if Size(S) <> size then
      Add(errors,[g,S,size,opt]);
      return "-";
  fi;
  pr := ProductReplacer(g);
  for i in [1..10] do
      x := Next(pr);
      slp := SiftGroupElementSLP(S,x);
      if slp.isone = false then
          Add(errors,[g,S,size,opt,x,slp]);
          success := "-";
      else
          y := ResultOfStraightLineProgram(slp.slp,StrongGenerators(S));
          if not(opt.IsOne(x*y^-1)) then
              Add(errors,[g,S,size,opt,x,y,slp]);
              success := "-";
          fi;
      fi;
  od;
  return success;
end;

GetAtlasGroup := function(name,rep,max)
  local ct,g,gens,m,s,size;
  gens := AtlasGenerators(name,rep);
  size := gens.size;
  gens := gens.generators;
  if max > 0 then
      s := AtlasStraightLineProgram(name,max).program;
      SlotUsagePattern(s);
      gens := ResultOfStraightLineProgram(s,gens);
      ct := CharacterTable(name);
      if ct <> fail then
          m := Maxes(ct);
          ct := CharacterTable(m[max]);
          if ct = fail then
              size := fail;
          else
              size := Size(ct);
          fi;
      else
          size := fail;
      fi;
  fi;
  g := Group(gens);
  if size <> fail then g!.mysize := size; fi;
  return g;
end;

name := "J1";
Print(name,": \c"); 
for r in [1,2,3,4,5,6,7,8,14] do
  for m in [0..7] do
    g := GetAtlasGroup(name,r,m); Print(TestGroup(g,g!.mysize,false),"\c");
  od;
od;
Print("\n");

name := "2.A10";
Print(name,": \c"); 
for r in [1..4] do
    g := GetAtlasGroup(name,r,0); Print(TestGroup(g,g!.mysize,false),"\c");
od;
Print("\n");

name := "M11";
Print(name,": \c"); 
for r in [1..25] do
  for m in [0..5] do
    g := GetAtlasGroup(name,r,m); Print(TestGroup(g,g!.mysize,false),"\c");
  od;
od;
Print("\n");

name := "M12";
Print(name,": \c");
for r in [1..26] do
  for m in [0..11] do
    g := GetAtlasGroup(name,r,m); Print(TestGroup(g,g!.mysize,false),"\c");
  od;
od;
Print("\n");

Print("GL(6,3): \c");
g := Group(GeneratorsOfGroup(GL(6,3)));
Print(TestGroup(g,Size(GL(6,3)),false),"\n");


