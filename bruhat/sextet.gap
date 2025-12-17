#!/usr/bin/gap

SizeScreen([1000,1000]);
G := MathieuGroup(24);
Hs := MaximalSubgroupClassReps(G);;
H := Hs[5];
Js := ConjugateSubgroups(G, H);;
Print(Length(Js), "\n");

#for J in Js do Print(J, "\n"); od;


V := Subspace(GF(2)^24, Z(2)^0*[
    [0,0,0,0,1,0,0,0,0,0,1,1,0,0,1,0,1,0,0,1,0,0,1,1],
    [1,0,0,0,0,1,1,1,0,1,0,1,0,0,1,1,1,1,0,1,1,0,0,0],
    [0,0,1,0,0,0,1,0,0,0,0,1,1,0,1,1,1,1,0,1,1,0,1,1],
    [0,0,0,0,0,0,1,1,1,1,0,1,1,0,1,0,1,1,1,0,1,0,1,0],
    [0,0,0,0,0,1,1,0,0,1,1,0,1,1,0,0,0,0,0,0,1,0,0,1],
    [0,0,0,0,0,1,1,0,0,1,1,1,1,0,0,1,1,1,0,0,1,1,1,0],
    [0,0,0,1,0,1,0,1,1,0,0,1,0,0,0,1,0,0,0,1,1,0,0,0],
    [0,1,0,0,0,1,0,0,1,1,0,1,1,0,1,1,0,1,0,1,1,0,1,0]]);

#a := PermutationMat((1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23), 24);;
#b := PermutationMat((3,17,10,7,9)(4,13,14,19,5)(8,18,11,12,23)(15,20,22,21,16), 24);;
#c := PermutationMat((1,24)(2,23)(3,12)(4,16)(5,18)(6,10)(7,20)(8,14)(9,21)(11,17)(13,22)(15,19),24);;
#G := Group(a,b,c);

#G := SymmetricGroup(24);; # too big
G:= Group(List(GeneratorsOfGroup(G), g->PermutationMat(g,24)));;

H := Stabilizer(G,CanonicalBasis(V),OnSubspacesByCanonicalBasis);
#Print(H, "\n");
Print("|H| =", Order(H), "\n");
Print("GeneratorsOfGroup:", Length(GeneratorsOfGroup(H)));

Print(StructureDescription(H), "\n");
#Hs := MaximalSubgroupClassReps(G);;
for J in Js do 
    Print("J:", Length(J), "\n");
    Print("IsomorphicSubgroups:", Length(IsomorphicSubgroups(J,H)), "\n");
od;
Print("\n");





