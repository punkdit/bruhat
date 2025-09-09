#!/usr/bin/gap

# time gap -o 8g hierarchy.gap

# Build clifford groups, up to 3 qubits, and
# search for T operators defined by field extensions......


Print("running hierarchy.gap\n");;

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

A := [[1,0], [0, E(3)]];
Print(Order(A), "\n");


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


in_third_level := function(U)
    # Is U in the third level of the clifford hierarchy ?
    local A;
    for g in Pauli1 do
        A := U*g*Inverse(U)*Inverse(g);
        if A in Cliff1 then continue; fi;
        return false; # no
    od;
    return true; # yes
end;;

in_fourth_level := function(U)
    # Is U in the fourth level of the clifford hierarchy ?
    local A;
    for g in Pauli1 do
        A := U*g*Inverse(U)*Inverse(g);
        if in_third_level(U) then continue; fi;
        return false; # no
    od;
    return true; # yes
end;;

in_fifth_level := function(U)
    # Is U in the fifth level of the clifford hierarchy ?
    local A;
    for g in Pauli1 do
        A := U*g*Inverse(U)*Inverse(g);
        if in_fourth_level(U) then continue; fi;
        return false; # no
    od;
    return true; # yes
end;;

in_sixth_level := function(U)
    # Is U in the sixth level of the clifford hierarchy ?
    local A;
    for g in Pauli1 do
        A := U*g*Inverse(U)*Inverse(g);
        if in_fourth_level(U) then continue; fi;
        return false; # no
    od;
    return true; # yes
end;;


Print("C3 ? ", in_third_level(A), "\n");
Print("C4 ? ", in_fourth_level(A), "\n");
Print("C5 ? ", in_fifth_level(A), "\n");
Print("C6 ? ", in_sixth_level(A), "\n");


