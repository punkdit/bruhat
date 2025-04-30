#!/usr/bin/gap

# time gap -o 8g clifford.gap

# Build clifford groups, up to 3 qubits, and

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

#for U in Cliff1 do
#    found := false;;
#    #Udag := Inverse(U);
#    #Print("found? ");
#    for g in Pauli1 do
#        if (U*g*Inverse(U)*Inverse(g) in Pauli1)  then found:=true; break; fi;
#    od;
#    Assert(0, found);
#od;

ii := KroneckerProduct(i, i);;
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

Assert(0, hi*si*si*hi*cz = cz*hi*is*si*is*si*hi);
Assert(0, w8 in Cliff2);
Assert(0, w8*cz*hi*cz = si*hi*cz*si*is*hi*si );

#for g in Center(Cliff2) do Print(g, "\n"); od;

#Print(Order(Center(Cliff2)), "\n");
#Print(IsomorphismGroups(Center(Cliff2), Group([[E(8)]])), "\n"); # yes

center := Group(w8);

G := FactorGroup(Cliff2, Group(E(8)*ii));

index := Order(G)/720;
Print(index,"\n");
for H in LowIndexSubgroups( G, index ) do
    #Print(H, "\n");
    #H := Representative(H);
    Print(Order(H), " ");
    if Order(H)=1152  then Print(StructureDescription(H),"\n"); fi;
od;
Print("\n");

#quit;
#
#Pauli2 := Group(wi, xi, ix, zi, iz);; # Order 64
#_Pauli2 := Group(E(8)*ii, xi, ix, zi, iz);; # Order 128
#
#Aut := AutomorphismGroup(Pauli2);
#Print(Order(Aut), "\n");
#
#count := 0;
#g := E(4)*ii;
#gens := [];
#for a in Aut do
#    #Print(xi^a, "\n");
#    if g^a=g then 
#        count:=count+1; 
#        Append(gens,[a]); 
#        A := Group(gens);
#        Print("\t", Order(A), "\n");
#        if Order(A) = 11520 then break; fi;
#    fi;
#od;
#
##A := AsGroup(gens);
##A := Group(gens);
#Print(Order(A), "\n");
#
#A := Image(ActionHomomorphism(A, Pauli2));
#
#G := FactorGroup(Cliff2, Group(E(8)*ii));
#
#Print(IsomorphismGroups(G,A), "\n"); # yes !
#
#quit;
#
##G := FactorGroup(Cliff2, Group(E(4)*ii));
##Print(Order(G), "\n");
##Print(IsomorphismGroups(G, Aut), "\n"); # not isomorphic
#k := Order(Aut);
#Print(StructureDescription(Aut), "\n");
#for H in LowIndexSubgroups( Cliff2, 4 ) do
#    if Order(H)=k  then 
#        Print(StructureDescription(H),"\n"); 
#        Print(IsomorphismGroups(H, Aut), "\n"); 
#    fi;
#od;
#
#quit;
#
#ASp := SemidirectProduct(Sp(4,2), GF(2)^4);
#Print(StructureDescription(ASp), "\n");

# is this the pseudo affine symplectic group... ? https://arxiv.org/abs/math/0007196
#APs := FactorGroup(Cliff2, Center(Cliff2));
#APs := FactorGroup(Cliff2, center);
#Print(StructureDescription(APs), "\n");
#Print(IsomorphismGroups(APs, ASp), "\n"); # fail

# APs is isomorphic to "ASp" group constructed in heisenberg.py

#quit;

# ------------------------------------------------------------------

#G := FactorGroup(Cliff2, _Pauli2);
#Print(Order(G), "\n");
#Print(IsomorphismGroups(G, Sp(4,2)), "\n"); # yes
#
#quit;


iii := KroneckerProduct(ii, i);;

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
Pauli3 := Group(wii, xii, ixi, iix, zii, izi, iiz);; # Order 256
_Pauli3 := Group(E(8)*iii, xii, ixi, iix, zii, izi, iiz);; # Order 512

G := FactorGroup(Cliff3, Group(E(4)*iii)); # no subgroup found
G := FactorGroup(Cliff3, Group(E(8)*iii)); # no subgroup found

G := SemidirectProduct(Sp(6,2), GF(2)^6); # has subgroup Sp(6,2)=GO(7,2)

Print(Order(G),"\n");
Print("LowIndexSubgroups:\n");
index := Order(G)/1451520; # looking for Order(Sp(6,2))
Print(index,"\n");
for H in LowIndexSubgroups( G, index ) do
    Print(Order(H), " ");
    if Order(H)=1451520  then Print(StructureDescription(H),"\n"); fi;
od;
Print("\n");

quit;

Print("Pauli3: ", Order(Pauli3), "\n");;
Print("_Pauli3: ", Order(_Pauli3), "\n");;
#Print("Cliff3: ", Order(Cliff3), "\n");; # slow ...

G := FactorGroup(Cliff3, _Pauli3);
#Print(Order(G), "\n");
#Print(Sp(6,2), "\n");

Print("IsomorphismGroups ?", "\n");
Print(IsomorphismGroups(G, Sp(6,2)), "\n");

quit;



# ca in Cliff3 = false
# ca in third level of clifford hierarchy = true


Cliff3_12 := Group(sii, isi, hii, ihi, wii, czi);;
Cliff3_13 := Group(sii, iis, hii, iih, wii);; # cz on 1&3 ??
Cliff3_23 := Group(isi, iis, ihi, iih, wii, icz);;

#SmCliff3 := Group(sii, isi, iis, hii, ihi, iih, icz, czi);; # Order 743178240

Print("warming up...\n");
Order(Cliff3);; # need this line otherwise membership test eats all memory...!
Print("Ok\n");

in_third_level := function(U)
    # Is U in the third level of the clifford hierarchy ?
    local A;
    for g in Pauli3 do
        A := U*g*Inverse(U)*Inverse(g);
        if A in Cliff3_12 then continue; fi;
        if A in Cliff3_13 then continue; fi;
        if A in Cliff3_23 then continue; fi;
        if A in Cliff3 then continue; fi;
        return false; # no
    od;
    return true; # yes
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
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0]];

Print("# Where do control control 1-qubit Pauli gates live ?\n");
for g in Pauli1 do 
    #Print(g, "\n");
    U3{[7,8]}{[7,8]} := g;
    if U3 in Pauli3 then Print("1\c"); continue; fi;
    if U3 in Cliff3 then Print("2\c"); continue; fi;
    if in_third_level(U3) then Print("3\c"); else Print(".\c"); fi;
od;
Print("\n");

Print("# Where do control control 1-qubit clifford gates live ?\n");
for g in Cliff1 do 
    #Print(g, "\n");
    U3{[7,8]}{[7,8]} := g;
    if U3 in Pauli3 then Print("1\c"); continue; fi;
    if U3 in Cliff3 then Print("2\c"); continue; fi;
    if in_third_level(U3) then Print("3\c"); else Print(".\c"); fi;
od;
Print("\n");

Print("# Where do control 2-qubit Pauli gates live ?\n");
for g in Pauli2 do 
    #Print(g, "\n");
    U3{[5,6,7,8]}{[5,6,7,8]} := g;
    if U3 in Pauli3 then Print("1\c"); continue; fi;
    if U3 in Cliff3 then Print("2\c"); continue; fi;
    if in_third_level(U3) then Print("3\c"); else Print(".\c"); fi;
od;
Print("\n");

Print("# Where do control 2-qubit clifford gates live ?\n");
for g in Cliff2 do 
    #Print(g, "\n");
    U3{[5,6,7,8]}{[5,6,7,8]} := g;
    if U3 in Pauli3 then Print("1\c"); continue; fi;
    if U3 in Cliff3 then Print("2\c"); continue; fi;
    if in_third_level(U3) then Print("3\c"); else Print(".\c"); fi;
od;

Print("\n");

Print("Done.\n");


