#!/usr/bin/gap


w := E(3);

i := [[1,0,0],[0,1,0],[0,0,1]];
s := [[1,0,0],[0,1,0],[0,0,w]];
z := [[1,0,0],[0,w,0],[0,0,w^2]];
x := [[0,1,0],[0,0,1],[1,0,0]];
h := (1/(w^2-w))*[[1,1,1],[1,w,w^2],[1,w^2,w]];


Pauli1 := Group(x, z);
Cliff1 := Group(s, h);
Print("Pauli1: ", Order(Pauli1), "\n");
Print("Cliff1: ", Order(Cliff1), "\n");

Assert(0, w*i in Pauli1);
Assert(0, x in Cliff1);


xi := KroneckerProduct(x, i);;
ix := KroneckerProduct(i, x);;
zi := KroneckerProduct(z, i);;
iz := KroneckerProduct(i, z);;
si := KroneckerProduct(s, i);;
is := KroneckerProduct(i, s);;
hi := KroneckerProduct(h, i);;
ih := KroneckerProduct(i, h);;

cx := [
    [1, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 0]];

Pauli2 := Group(xi, ix, zi, iz);
Cliff2 := Group(si, is, hi, ih, cx);
Print("Pauli2: ", Order(Pauli2), "\n");
#Print("Cliff2: ", Order(Cliff2), "\n");

Print(Order(Center(Cliff2)), "\n");
#Print(IsomorphismGroups(Center(Cliff3), Group([[E(8)]]
center := Center(Cliff2);
#aff_sy := FactorGroup(Cliff2, Center(Cliff2));
aff_sy := FactorGroup(Cliff2, center);
Print(StructureDescription(aff_sy), "\n");

ASp := SemidirectProduct(Sp(4,3), GF(3)^4);
Print(StructureDescription(ASp), "\n");

Print(IsomorphismGroups(aff_sy, ASp), "\n");

#Print(Order(Sp(4,3)), "\n");

Print("Done.\n\n\n");

QUIT;

