
SizeScreen([1000,1000]);

for sz in [4..56] do
    Gs := AllSmallGroups(sz);;
    for G in Gs do
        desc := StructureDescription(G);
        Print(sz, " / ", desc, " / ");
H:=Image(IsomorphismFpGroup(G));
H:=SimplifiedFpGroup(Image(IsomorphismFpGroup(G)));
Print(H, " / ");
Print(RelatorsOfFpGroup(H), "\n");
    od;
od;

quit;

