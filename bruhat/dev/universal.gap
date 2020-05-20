#!/usr/bin/gap

# https://cp4space.wordpress.com/2020/05/10/minimalistic-quantum-computation/

Print("running universal.gap\n");;


r2 := Sqrt(2);;
ir2 := 1/r2;;

i := [[1, 0], [0, 1]];;
w := [[E(4), 0], [0, E(4)]];;
x := [[0, 1], [1, 0]];;
z := [[1, 0], [0, -1]];;
s := [[1, 0], [0, E(4)]];;
h := [[ir2, ir2], [ir2, -ir2]];;

Cliff1 := Group(w, s, h);; # Order 192
Pauli1 := Group(w, x, z);; # Order 32

xi := KroneckerProduct(x, i);;
ix := KroneckerProduct(i, x);;
zi := KroneckerProduct(z, i);;
iz := KroneckerProduct(i, z);;
si := KroneckerProduct(s, i);;
is := KroneckerProduct(i, s);;
hi := KroneckerProduct(h, i);;
ih := KroneckerProduct(i, h);;
wi := KroneckerProduct(w, i);;


cz := [
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, -1]];;

Cliff2 := Group(si, is, hi, ih, wi, cz);; # Order 92160
#for g in Cliff2 do Print(g, "\n"); od;
Pauli2 := Group(wi, xi, ix, zi, iz);;


xii := KroneckerProduct(xi, i);;
ixi := KroneckerProduct(i, xi);;
iix := KroneckerProduct(i, ix);;

zii := KroneckerProduct(zi, i);;
izi := KroneckerProduct(i, zi);;
iiz := KroneckerProduct(i, iz);;

sii := KroneckerProduct(si, i);;
isi := KroneckerProduct(i, si);;
iis := KroneckerProduct(i, is);;

hii := KroneckerProduct(hi, i);;
ihi := KroneckerProduct(i, hi);;
iih := KroneckerProduct(i, ih);;

wii := KroneckerProduct(wi, i);;

icz := KroneckerProduct(i, cz);;
czi := KroneckerProduct(cz, i);;


Tofolli := [
    [1, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 0, 0, 1, 0]];;


zhh := zii*ihi*iih;
hzh := izi*hii*iih;
hhz := iiz*hii*ihi;

xhh := xii*ihi*iih;
hxh := ixi*hii*iih;
hhx := iix*hii*ihi;

Gx := Group(Tofolli, xhh, hxh, hhx);
Gz := Group(Tofolli, zhh, hzh, hhz);

Print("Order(Gz) = ", Order(Gz), "\n");
Print("Order(Gx) = ", Order(Gx), "\n");

L7:= SimpleLieAlgebra("E", 7, Rationals);;
R7:= RootSystem(L7);;
W7:= WeylGroup(R7);
Print("Order(W7) = ", Order(W7), "\n");

f := IsomorphismGroups(Gz, W7);
Print("IsomorphismGroups(Gz, E7) = ", f, "\n");

L8:= SimpleLieAlgebra("E", 8, Rationals);;
R8:= RootSystem(L8);;
W8:= WeylGroup(R8);
Print("Order(W8) = ", Order(W8), "\n");

f := IsomorphismGroups(Gx, W8);
Print("IsomorphismGroups(Gx, E8) = ", f, "\n");





