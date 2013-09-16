# Test file for the genss package:

LoadPackage("atlasrep");

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


Print("Diagonal: \c");
for p in [3,5,7] do
  g := Group(DiagonalMat([Z(p)^1, Z(p)^0, Z(p)^0]),
             DiagonalMat([Z(p)^0, Z(p)^1, Z(p)^0]),
             DiagonalMat([Z(p)^0, Z(p)^0, Z(p)^1]));
  Print(TestGroup(g,(p-1)^3,false),"\c");
  Print(TestGroup(g,(p-1)^2,rec( Projective := true)),"\c");
od;
Print("\n");


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
Print(TestGroup(g,Size(GL(6,3)),false),"\c");
Print(TestGroup(g,Size(GL(6,3))/2,rec( Projective := true)),"\n");

Print("SL(4,5): \c");
g := Group(GeneratorsOfGroup(SL(4,5)));
Print(TestGroup(g,Size(SL(4,5)),false),"\c");
Print(TestGroup(g,Size(SL(4,5))/4,rec( Projective := true)),"\n");

Print("SL(6,5): \c");
g := Group(GeneratorsOfGroup(SL(6,5)));
Print(TestGroup(g,Size(SL(6,5)),false),"\c");
Print(TestGroup(g,Size(SL(6,5))/2,rec( Projective := true)),"\n");

Print("SL(7,5): \c");
g := Group(GeneratorsOfGroup(SL(7,5)));
Print(TestGroup(g,Size(SL(7,5)),false),"\c");
Print(TestGroup(g,Size(SL(7,5)),rec( Projective := true)),"\n");

name := "J2";
Print(name,": \c"); 
for r in [1,2,3,4,5,6,7,8,11,16,17,19,20] do
  for m in [0..9] do
    g := GetAtlasGroup(name,r,m); Print(TestGroup(g,g!.mysize,false),"\c");
  od;
od;
Print("\n");

name := "J3";
Print(name,": \c"); 
for r in [1,2,4,9] do
  for m in [0..9] do
    g := GetAtlasGroup(name,r,m); Print(TestGroup(g,g!.mysize,false),"\c");
  od;
od;
Print("\n");

name := "Fi22";
Print(name,": \c"); 
for r in [1,2,4,6] do
  if r = 6 then
      g := GetAtlasGroup(name,r,0);
      h := GetAtlasGroup(name,r,1);
      m := GModuleByMats(GeneratorsOfGroup(h),FieldOfMatrixGroup(g));
      b := MutableCopyMat(MTX.BasesCompositionSeries(m)[2]);
      TriangulizeMat(b);
      opt := rec( Cand := rec(points := [b], 
                              ops := [OnSubspacesByCanonicalBasis] ) );
  else
      opt := false;
  fi;
  for m in [0..14] do
    g := GetAtlasGroup(name,r,m); Print(TestGroup(g,g!.mysize,opt),"\c");
  od;
od;
Print("\n");

name := "Fi23";
Print(name,": \c"); 
for r in [1,2,3] do    # 4 taken out since it used too much memory
  if r = 4 then
      g := GetAtlasGroup(name,r,0);
      h := GetAtlasGroup(name,r,1);
      m := GModuleByMats(GeneratorsOfGroup(h),FieldOfMatrixGroup(g));
      v := MutableCopyMat(MTX.BasesCompositionSeries(m)[2])[1];
      opt := rec( Cand := rec(points := [v], 
                              ops := [OnRight] ) );
  else
      opt := false;
  fi;
  for m in [0,1,2,3,4,5,6,9,10,13,14] do
    g := GetAtlasGroup(name,r,m); Print(TestGroup(g,g!.mysize,opt),"\c");
  od;
od;
Print("\n");


