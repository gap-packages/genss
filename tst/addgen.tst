gap> START_TEST("addgen.tst");

# AddGeneratorToStabilizerChain on a chain for the trivial group
gap> S := StabilizerChain(Group(()));;
gap> Size(S);
1
gap> AddGeneratorToStabilizerChain(S, (1,2,3));
true
gap> Size(S);
3
gap> AddGeneratorToStabilizerChain(S, (1,2,3)^2);
false
gap> AddGeneratorToStabilizerChain(S, (1,2));
true
gap> Size(S);
6

#
gap> S := StabilizerChain(Group(One(GL(3,5))), rec( Projective := true ));;
gap> Size(S);
1
gap> AddGeneratorToStabilizerChain(S, Z(5) * [[1,1,0],[0,1,0],[0,0,1]]);
true
gap> Size(S);
5

#
gap> STOP_TEST("addgen.tst");
