#!/usr/bin/gap

# time gap -o 8g hierarchy.gap

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
Order(Cliff2);
#for g in Cliff2 do Print(g, "\n"); od;
Phase2 := Group(wi);;
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
    [-1, 0, 0, 0]];;

# Print(a in Cliff2, "\n"); # true

#hom:=EpimorphismFromFreeGroup(Cliff2:names:=["si", "is", "hi", "ih", "wi", "cz"]);
#Print(PreImagesRepresentative(hom, a), "\n"); # si^3*is^3*hi*si^2*ih*is^2*ih*cz*hi*si
#Print(PreImagesRepresentative(hom, a*a), "\n"); # (si^2*hi)^2
#Print(PreImagesRepresentative(hom, a*a*a), "\n"); # si^2*is^3*hi*(si^2*hi*si)^2*si*ih*is^2*ih*cz*hi*si

cls := ConjugacyClass(Cliff2, a);;
#Print(a,               a in cls, "\n");               # true
#Print(a*a,             a*a in cls, "\n");             # false
#Print(a*a*a,           a*a*a in cls, "\n");           # true
#Print(a*a*a*a,         a*a*a*a in cls, "\n");         # false
#Print(a*a*a*a*a,       a*a*a*a*a in cls, "\n");       # true
#Print(a*a*a*a*a*a,     a*a*a*a*a*a in cls, "\n");     # false
#Print(a*a*a*a*a*a*a,   a*a*a*a*a*a*a in cls, "\n");   # true
#Print(a*a*a*a*a*a*a*a, a*a*a*a*a*a*a*a in cls, "\n"); # false
#Print(cz, cz in cls, "\n"); # false


b := [
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 0, 1],
    [0, 0, -1, 0]];;
#Print(cz, cz in ConjugacyClass(Cliff2, b), "\n"); # false

#for b in ConjugacyClass(Cliff2, cz) do Print(b, "\n"); od;

quit;

#Print(a*a, "\n");
#Print(a*a*a, "\n");
#Print(a*a*a*a, "\n");

# Print(Order(G2));


U2 := [ 
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1]];



in_third_level_2 := function(U2)
    # Is U2 in the third level of the clifford hierarchy ?
    local A;
    for g in Pauli2 do
        A := U2*g*Inverse(U2)*Inverse(g);
        if A in Cliff2 then continue; fi;
        return false; # no
    od;
    return true; # yes
end;;

in_fourth_level_2 := function(U2)
    # Is U2 in the fourth level of the clifford hierarchy ?
    local A;
    for g in Pauli2 do
        A := U2*g*Inverse(U2)*Inverse(g);
        if in_third_level_2(U2) then continue; fi;
        return false; # no
    od;
    return true; # yes
end;;

in_fifth_level_2 := function(U2)
    # Is U2 in the fifth level of the clifford hierarchy ?
    local A;
    for g in Pauli2 do
        A := U2*g*Inverse(U2)*Inverse(g);
        if in_fourth_level_2(U2) then continue; fi;
        return false; # no
    od;
    return true; # yes
end;;

get_level_2 := function(U2)
    if U2 in Phase2 then return "0\c"; fi;
    if U2 in Pauli2 then return "1\c"; fi;
    if U2 in Cliff2 then return "2\c"; fi;
    if in_third_level_2(U2) then return "3\c"; fi;
    if in_fourth_level_2(U2) then return "4\c"; fi;
    if in_fifth_level_2(U2) then return "5\c"; fi;
    return ".\c";
end;;

Print("# Where do control 1-qubit Pauli gates live ?\n");
for g in Pauli1 do 
    #Print(g, "\n");
    U2{[3,4]}{[3,4]} := g;
    Print(get_level_2(U2));
od;
Print("\n");

Print("# Where do control 1-qubit clifford gates live ?\n");
for g in Cliff1 do 
    #Print(g, "\n");
    U2{[3,4]}{[3,4]} := g;
    Print(get_level_2(U2));
od;
Print("\n");

quit;

# ---------------------------------------------------
#
# 3 qubit gates
#


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

cb := [
    [1, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, -1, 0, 0, 0],
    [0, 0, 0, 0, 0, -1, 0, 0]];;

cc := [
    [1, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, -1, 0, 0, 0],
    [0, 0, 0, 0, 0, -1, 0, 0],
    [0, 0, 0, 0, 0, 0, -1, 0]];;


Tofolli := [
    [1, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 0, 0, 1, 0]];;


Cliff3 := Group(sii, isi, iis, hii, ihi, iih, wii, icz, czi);; # Order 743178240
Pauli3 := Group(wii, xii, ixi, iix, zii, izi, iiz);; # Order 256

# Print(Order(Pauli3), "\n");;

# ca in Cliff3 = false
# ca in third level of clifford hierarchy = true


Cliff3_12 := Group(sii, isi, hii, ihi, wii, czi);;
Cliff3_13 := Group(sii, iis, hii, iih, wii);; # cz on 1&3 ??
Cliff3_23 := Group(isi, iis, ihi, iih, wii, icz);;

#SmCliff3 := Group(sii, isi, iis, hii, ihi, iih, icz, czi);; # Order 743178240

Print("warming up...\n");
Order(Cliff3);; # need this line otherwise membership test eats all memory...!
Print("Ok\n");


in_first_level := function(U3)
    return U3 in Pauli3;
end;;

in_second_level := function(U3)
    if U3 in Cliff3_12 then return true; fi;
    if U3 in Cliff3_13 then return true; fi;
    if U3 in Cliff3_23 then return true; fi;
    if U3 in Cliff3 then return true; fi;
    return false;
end;;
    
#in_second_level := function(U3)
#    # Is U3 in the second level of the clifford hierarchy ?
#    local A;
#    for g in Pauli3 do
#        A := U3*g*Inverse(U3)*Inverse(g);
#        if in_first_level(A) then continue; fi;
#        return false; # no
#    od;
#    return true; # yes
#end;;


in_third_level := function(U3)
    # Is U3 in the third level of the clifford hierarchy ?
    local A;
    for g in Pauli3 do
        A := U3*g*Inverse(U3)*Inverse(g);
        if in_second_level(A) then continue; fi;
        return false; # no
    od;
    return true; # yes
end;;

in_fourth_level := function(U3)
    # Is U3 in the fourth level of the clifford hierarchy ?
    local A;
    for g in Pauli3 do
        A := U3*g*Inverse(U3)*Inverse(g);
        if in_third_level(A) then continue; fi;
        return false; # no
    od;
    return true; # yes
end;;


get_level := function(U3)
    if U3 in Pauli3 then return "1\c"; fi;
    if in_second_level(U3) then return "2\c"; fi;
    if in_third_level(U3) then return "3\c"; fi;
    if in_fourth_level(U3) then return "4\c"; fi;
    return ".\c";
end;;
    

Print("in_third_level(ca):", in_third_level(ca), "\n"); # true
Print("in_third_level(cb):", in_third_level(cb), "\n"); # true
Print("in_third_level(cc):", in_third_level(cc), "\n"); # true
Print("in_third_level(Tofolli):", in_third_level(Tofolli), "\n"); # true
Print("in_third_level(Tofolli*ca):", in_third_level(Tofolli*ca), "\n"); # false

U3 := [ 
    [1, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 1]];

Print("# Where do control control 1-qubit Pauli gates live ?\n");
for g in Pauli1 do 
    #Print(g, "\n");
    U3{[7,8]}{[7,8]} := g;
    Print(get_level(U3));
od;
Print("\n");

Print("# Where do control control 1-qubit clifford gates live ?\n");
for g in Cliff1 do 
    #Print(g, "\n");
    U3{[7,8]}{[7,8]} := g;
    Print(get_level(U3));
od;
Print("\n");

Print("# Where do control 2-qubit Pauli gates live ?\n");
for g in Pauli2 do 
    #Print(g, "\n");
    U3{[5,6,7,8]}{[5,6,7,8]} := g;
    Print(get_level(U3));
od;
Print("\n");

Print("# Where do control 2-qubit clifford gates live ?\n");
for g in Cliff2 do 
    #Print(g, "\n");
    U3{[5,6,7,8]}{[5,6,7,8]} := g;
    Print(get_level(U3));
od;

Print("\n");

Print("Done.\n");


