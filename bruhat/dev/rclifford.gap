#!/usr/bin/gap

# time gap -o 8g clifford.gap

# Build clifford groups, up to 3 qubits, and
# search for T operators defined by field extensions......


Print("running clifford.gap\n");;

# See:
# https://www.mathstat.dal.ca/~selinger/papers/clifford.pdf

r2 := Sqrt(2);;
ir2 := 1/r2;;

i := [[1, 0], [0, 1]];;
w := [[E(4), 0], [0, E(4)]];;
w8 := [[E(8), 0], [0, E(8)]];;
x := [[0, 1], [1, 0]];;
z := [[1, 0], [0, -1]];;
s := [[1, 0], [0, E(4)]];;
h := [[ir2, ir2], [ir2, -ir2]];;

Phase  := Group(w8);;
Pauli1 := Group(w, x, z);; # Order 16

Cliff1 := Group(s, h);; # Order 192
RCliff1 := Group(z, h);; # Order 192

Print("Cliff1: ", Order(Cliff1), " ", StructureDescription(Cliff1),  "\n");
Print("RCliff1: ", Order(RCliff1),  " ", StructureDescription(RCliff1), "\n");

xi := KroneckerProduct(x, i);;
ix := KroneckerProduct(i, x);;
zi := KroneckerProduct(z, i);;
iz := KroneckerProduct(i, z);;
si := KroneckerProduct(s, i);;
is := KroneckerProduct(i, s);;
hi := KroneckerProduct(h, i);;
ih := KroneckerProduct(i, h);;
wi := KroneckerProduct(w, i);;

w8 := [[E(8), 0, 0, 0], [0, E(8), 0, 0], [0, 0, E(8), 0], [0, 0, 0, E(8)]];;

cz := [
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, -1]];;

Cliff2 := Group(si, is, hi, ih, cz);; # Order 92160
RCliff2 := Group(zi, iz, hi, ih, cz);; 

Print("Cliff2: ", Order(Cliff2), " ", StructureDescription(Cliff2), "\n");
Print("RCliff2: ", Order(RCliff2),  " ", StructureDescription(RCliff2),"\n");

#ASp := SemidirectProduct(Sp(4,2), GF(2)^4);
#Print(Order(ASp), "\n");
#Print(StructureDescription(ASp), "\n");


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


Cliff3 := Group(sii, isi, iis, hii, ihi, iih, icz, czi);; # Order 743178240
RCliff3 := Group(zii, izi, iiz, hii, ihi, iih, icz, czi);; # Order ?

#Print("Cliff3:\n");
#Print(Order(Cliff3), "\n");
#Print(StructureDescription(Cliff3), "\n");

Print("RCliff3:\n");
Print(Order(RCliff3), "\n");
Print(StructureDescription(RCliff3), "\n");

