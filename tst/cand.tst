gap> START_TEST("cand.tst");

# Base point candidates given by the user, with actions other than the
# default one
gap> g := Group(GeneratorsOfGroup(GL(4,3)));;

#
gap> b := ImmutableMatrix(GF(3), [[1,0,0,0],[0,1,0,0]] * Z(3)^0);;
gap> S := StabilizerChain(g, rec( Cand := rec( points := [b],
>              ops := [OnSubspacesByCanonicalBasis] ) ));;
gap> Size(S) = Size(GL(4,3));
true
gap> BaseStabilizerChain(S).points[1] = b;
true
gap> Length(S!.orb);
130

#
gap> v := [1,0,0,0] * Z(3)^0;;
gap> S := StabilizerChain(g, rec( Cand := rec( points := [v],
>              ops := [OnRight] ) ));;
gap> Size(S) = Size(GL(4,3));
true
gap> BaseStabilizerChain(S).points[1] = v;
true
gap> Length(S!.orb);
80

#
gap> STOP_TEST("cand.tst");
