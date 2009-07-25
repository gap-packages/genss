m := Group(AtlasGenerators("M24",1).generators);
S := StabilizerChain(m);
ComputeSuborbitsForStabilizerChain(S);
f := function(x)
    return OnTuples( [ 2, 3 ], x ) = [ 1, 5 ];
end;
# The following would take extremely long:
# BacktrackSearchStabilizerChainElement(S,f,(),false);
pos2 := Position(S!.orb,2);
pos3 := Position(S!.orb,3);
pruner := function(S,x,t,w)
    return S!.suborbs.suborbnr[Position(S!.topS!.orb,1/x)] = 
           S!.suborbs.suborbnr[pos2] and
           S!.suborbs.suborbnr[Position(S!.topS!.orb,5/x)] = 
           S!.suborbs.suborbnr[pos3];
end;
BacktrackSearchStabilizerChainElement(S,f,(),pruner);

f := function(x)
    return OnSets( [ 2, 3 ], x ) = [ 1, 5 ];
end;
# The following would take extremely long:
# BacktrackSearchStabilizerChainElement(S,f,(),false);
pos2 := Position(S!.orb,2);
pos3 := Position(S!.orb,3);
pruner := function(S,x,t,w)
    return Set(S!.suborbs.suborbnr{[Position(S!.topS!.orb,1/x),
                                    Position(S!.topS!.orb,5/x)]}) = 
           Set(S!.suborbs.suborbnr{[pos2,pos3]});
end;
BacktrackSearchStabilizerChainElement(S,f,(),pruner);

# Another example:
gens := AtlasGenerators("Fi23",1).generators;
g := Group(gens);
S := StabilizerChain(g);
ComputeSuborbitsForStabilizerChain(S);
f := function(x)
    return 15^x = 16988 and 17^x = 21754;
end;
# Takes too long:
#BacktrackSearchStabilizerChainElement(S,f,(),false);
pos15 := Position(S!.orb,15);
pos17 := Position(S!.orb,17);
pruner := function(S,x,t,w)
    return S!.suborbs.suborbnr[Position(S!.topS!.orb,16988/x)] = 
           S!.suborbs.suborbnr[pos15] and
           S!.suborbs.suborbnr[Position(S!.topS!.orb,21754/x)] = 
           S!.suborbs.suborbnr[pos17];
end;
BacktrackSearchStabilizerChainElement(S,f,(),pruner);

