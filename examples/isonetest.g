gens := [ (1,2,3,4,5,6,7,8,9,10), (1,2), (11,12) ];
g := Group(gens);
S := StabilizerChain(g);
if Size(S) <> Factorial(10)*2 then Error("wrong size"); fi;
f := function ( x )
    return OnTuples( [ 1 .. 10 ], x ) = [ 1 .. 10 ];
end;
S := StabilizerChain(g,rec(IsOne := f));
if Size(S) <> Factorial(10) then Error("wrong size"); fi;
