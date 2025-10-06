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
Cliff1 := Group(w, s, h);; # Order 192
Pauli1 := Group(w, x, z);; # Order 16

#A := FactorGroup(Cliff1, Phase);
#A0 := SemidirectProduct(Sp(2,2), GF(2)^2);
#Print(Order(A), " ", Order(A0), "\n");
#Print(StructureDescription(A), " <-> ", StructureDescription(A0), "\n");
#quit;

for U in Cliff1 do
    found := false;;
    #Udag := Inverse(U);
    #Print("found? ");
    for g in Pauli1 do
        if (U*g*Inverse(U)*Inverse(g) in Pauli1)  then found:=true; break; fi;
    od;
    Assert(0, found);
od;

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

Cliff2 := Group(si, is, hi, ih, wi, cz);; # Order 92160
RCliff2 := Group(zi, iz, hi, ih, wi, cz);; 

Print(Order(Cliff2), "\n");
Print(Order(RCliff2), "\n");


ASp := SemidirectProduct(Sp(4,2), GF(2)^4);
Print(Order(ASp), "\n");
Print(StructureDescription(ASp), "\n");


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


Cliff3 := Group(sii, isi, iis, hii, ihi, iih, wii, icz, czi);; # Order 743178240
RCliff3 := Group(zii, izi, iiz, hii, ihi, iih, wii, icz, czi);; # Order ?

Print("Cliff3:\n");
Print(Order(Cliff3), "\n");
Print(StructureDescription(Cliff3), "\n");

Print("RCliff3:\n");
Print(Order(RCliff3), "\n");
Print(StructureDescription(RCliff3), "\n");

