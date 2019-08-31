
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
x := [[0, 1], [1, 0]];;
z := [[1, 0], [0, -1]];;
s := [[1, 0], [0, E(4)]];;
h := [[ir2, ir2], [ir2, -ir2]];;

Cliff1 := Group(w, s, h);; # Order 192
Pauli1 := Group(w, x, z);; # Order 32

for U in Cliff1 do
    found := false;;
    #Udag := Inverse(U);
    #Print("found? ");
    for g in Pauli1 do
        if (U*g*Inverse(U)*Inverse(g) in Pauli1)  then found:=true; break; fi;
    od;
    if not found then Print("Not found\n"); fi;
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


cz := [
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, -1]];;

Cliff2 := Group(si, is, hi, ih, wi, cz);; # Order 92160
Pauli2 := Group(wi, xi, ix, zi, iz);;

# Works:
for U in Cliff2 do
    found := false;;
    #Udag := Inverse(U);
    #Print("found? ");
    for g in Pauli2 do
        if (U*g*Inverse(U)*Inverse(g) in Pauli2)  then found:=true; break; fi;
    od;
    if not found then Print("Not found\n"); fi;
od;

a := [
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1],
    [-1, 0, 0, 0]];; # a in G2 = true

# Print(Order(G2));


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

ca := [
    [1, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, -1, 0, 0, 0]];;


Cliff3 := Group(sii, isi, iis, hii, ihi, iih, wii, icz, czi);; # Order 743178240
Pauli3 := Group(wii, xii, ixi, iix, zii, izi, iiz);; # Order 256

# Print(Order(Pauli3), "\n");;

# ca in Cliff3 = false

Cliff3_12 := Group(sii, isi, hii, ihi, wii, czi);;
Cliff3_13 := Group(sii, iis, hii, iih, wii);; # cz on 1&3 ??
Cliff3_23 := Group(isi, iis, ihi, iih, wii, icz);;

#SmCliff3 := Group(sii, isi, iis, hii, ihi, iih, icz, czi);; # Order 743178240

                      # WARNING:
Print(Order(Cliff3)); # need this line otherwise membership test eats all memory...!
Print("\n");

Print("Go:\n");
U := ca;
for g in Pauli3 do
    #if (U*g*Inverse(U) in Cliff3)  then Print("Found!\n"); break; fi;
    #Print(U*g*Inverse(U) in Cliff3); Print("\n");
    A := U*g*Inverse(U)*Inverse(g);
    if A in Cliff3_12 then continue; fi;
    if A in Cliff3_13 then continue; fi;
    if A in Cliff3_23 then continue; fi;
    if A in Cliff3 then continue; fi;
    #Print("...\n");
    #Print(A in Cliff3);
    #Print(A in SmCliff3);
    Print(A);
    Print("\n");
    Print("Not found!");
    Print("\n");
od;




Print("OK\n");


