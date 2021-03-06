LoadPackage("genss");
Print("GL(5,5) natural:\n\n");
g := GL(5,5);
gg := Group(GeneratorsOfGroup(g));
ti := Runtime();
S := StabilizerChain(gg);
ti2 := Runtime();
if Size(S) <> Size(g) then Error("wrong size"); fi;
ti3 := Runtime();
Print("\n");
S := StabilizerChain(g);
ti4 := Runtime();
Print("\nStabChain:\n");
ViewObj(S);
Print("\n\nTime with random verification: ",ti2-ti,"\n");
Print("Time with known size: ",ti4-ti3,"\n");
