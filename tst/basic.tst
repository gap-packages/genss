#
gap> TestGroup := function(g,size,opt)
>   local S,i,pr,slp,x,y;
>   if opt = false then opt := rec(); fi;
>   S := StabilizerChain(g,opt);
>   Assert(0, Size(S) = size, "size mismatch");
>   pr := ProductReplacer(g);
>   for i in [1..10] do
>       x := Next(pr);
>       slp := SiftGroupElementSLP(S,x);
>       Assert(0, slp.isone = true, "slp.isone false");
>       y := ResultOfStraightLineProgram(slp.slp,StrongGenerators(S));
>       Assert(0, opt.IsOne(x*y^-1), "x*y^-1 is nontrivial");
>   od;
> end;;

#
# diagonal matrix group
#
gap> diag := p -> Group(DiagonalMat([Z(p)^1, Z(p)^0, Z(p)^0]),
>              DiagonalMat([Z(p)^0, Z(p)^1, Z(p)^0]),
>              DiagonalMat([Z(p)^0, Z(p)^0, Z(p)^1]));;

#
gap> p:=3;; g:=diag(p);;
gap> TestGroup(g,(p-1)^3,false);
gap> TestGroup(g,(p-1)^2,rec( Projective := true));

#
gap> p:=5;; g:=diag(p);;
gap> TestGroup(g,(p-1)^3,false);
gap> TestGroup(g,(p-1)^2,rec( Projective := true));

#
gap> p:=7;; g:=diag(p);;
gap> TestGroup(g,(p-1)^3,false);
gap> TestGroup(g,(p-1)^2,rec( Projective := true));

#
# various general and special linear groups
#

#
gap> g := Group(GeneratorsOfGroup(GL(6,3)));;
gap> TestGroup(g,Size(GL(6,3)),false);
gap> TestGroup(g,Size(GL(6,3))/2,rec( Projective := true));

#
gap> g := Group(GeneratorsOfGroup(SL(4,5)));;
gap> TestGroup(g,Size(SL(4,5)),false);
gap> TestGroup(g,Size(SL(4,5))/4,rec( Projective := true));

#
gap> g := Group(GeneratorsOfGroup(SL(6,5)));;
gap> TestGroup(g,Size(SL(6,5)),false);
gap> TestGroup(g,Size(SL(6,5))/2,rec( Projective := true));

#
gap> g := Group(GeneratorsOfGroup(SL(7,5)));;
gap> TestGroup(g,Size(SL(7,5)),false);
gap> TestGroup(g,Size(SL(7,5)),rec( Projective := true));
