gap> START_TEST("isone.tst");

#
gap> gens := [ (1,2,3,4,5,6,7,8,9,10), (1,2), (11,12) ];;
gap> g := Group(gens);
Group([ (1,2,3,4,5,6,7,8,9,10), (1,2), (11,12) ])

#
gap> S := StabilizerChain(g);
<stabchain size=7257600 orblen=10 layer=1 SchreierDepth=6>
 <stabchain size=725760 orblen=9 layer=2 SchreierDepth=2>
  <stabchain size=80640 orblen=8 layer=3 SchreierDepth=2>
   <stabchain size=10080 orblen=7 layer=4 SchreierDepth=3>
    <stabchain size=1440 orblen=6 layer=5 SchreierDepth=3>
     <stabchain size=240 orblen=5 layer=6 SchreierDepth=2>
      <stabchain size=48 orblen=4 layer=7 SchreierDepth=2>
       <stabchain size=12 orblen=3 layer=8 SchreierDepth=1>
        <stabchain size=4 orblen=2 layer=9 SchreierDepth=1>
         <stabchain size=2 orblen=2 layer=10 SchreierDepth=1>
gap> Size(S);
7257600

#
gap> f := x -> OnTuples( [ 1 .. 10 ], x ) = [ 1 .. 10 ];;
gap> S := StabilizerChain(g,rec(IsOne := f));
<stabchain size=3628800 orblen=10 layer=1 SchreierDepth=6>
 <stabchain size=362880 orblen=9 layer=2 SchreierDepth=3>
  <stabchain size=40320 orblen=8 layer=3 SchreierDepth=3>
   <stabchain size=5040 orblen=7 layer=4 SchreierDepth=3>
    <stabchain size=720 orblen=6 layer=5 SchreierDepth=2>
     <stabchain size=120 orblen=5 layer=6 SchreierDepth=3>
      <stabchain size=24 orblen=4 layer=7 SchreierDepth=1>
       <stabchain size=6 orblen=3 layer=8 SchreierDepth=2>
        <stabchain size=2 orblen=2 layer=9 SchreierDepth=1>
gap> Size(S);
3628800

#
gap> STOP_TEST("isone.tst", 1);
