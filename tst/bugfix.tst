gap> START_TEST("bugfix.tst");

#
# See https://github.com/gap-packages/genss/issues/5
#
gap> for i in [1..100] do
>    g := PrimitiveGroup(9, 5);
>    actual := SetwiseStabilizer(g, OnPoints, [3]).setstab;
>    expected := Stabilizer(g, 3);
>    if actual <> expected then
>        Error("Not equal in iteration ", i);
>    fi;
> od;
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...
#I  Computing adjusted stabilizer chain...

#
gap> STOP_TEST("bugfix.tst", 1);
