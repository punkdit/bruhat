#!/usr/bin/gap


r2 := Sqrt(2);;
ir2 := 1/r2;;

i := [[1, 0], [0, 1]];;
w := [[E(4), 0], [0, E(4)]];;
w8 := [[E(8), 0], [0, E(8)]];;
x := [[0, 1], [1, 0]];;
z := [[1, 0], [0, -1]];;
s := [[1, 0], [0, E(4)]];;
h := [[ir2, ir2], [ir2, -ir2]];;




Pauli1 := Group(x, z);
Cliff1 := Group(s, h);
Print("Pauli1: ", Order(Pauli1), "\n");
Print("Cliff1: ", Order(Cliff1), "\n");


wi := KroneckerProduct(w, i);;
w8i := KroneckerProduct(w8, i);;
xi := KroneckerProduct(x, i);;
ix := KroneckerProduct(i, x);;
zi := KroneckerProduct(z, i);;
iz := KroneckerProduct(i, z);;
si := KroneckerProduct(s, i);;
is := KroneckerProduct(i, s);;
hi := KroneckerProduct(h, i);;
ih := KroneckerProduct(i, h);;
ii := KroneckerProduct(i, i);;


cz := [
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, -1]];

Print(si, "\n");
Pauli2 := Group(xi, ix, zi, iz, wi);
Cliff2 := Group(si, is, hi, ih, cz);
LCliff2 := Group(si, is, hi, ih); # already has w8i
Print("Pauli2: ", Order(Pauli2), "\n");
Print("Cliff2: ", Order(Cliff2), "\n");
Print("LCliff2: ", Order(LCliff2), "\n");
Print("IsNormal: ", IsNormal(Cliff2, Pauli2), "\n");
Print("IsNormal: ", IsNormal(Cliff2, LCliff2), "\n"); # nope

Print("Done.\n\n\n");

QUIT;

