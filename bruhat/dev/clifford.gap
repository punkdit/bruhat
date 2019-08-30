
r2 := Sqrt(2);;
ir2 := 1/r2;;

i := [[1, 0], [0, 1]];;
w := [[E(4), 0], [0, E(4)]];;
x := [[0, 1], [1, 0]];;
z := [[1, 0], [0, -1]];;
s := [[1, 0], [0, E(4)]];;
h := [[ir2, ir2], [ir2, -ir2]];;

G1 := Group(w, x, z, s, h);; # Order 192


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

G2 := Group(xi, ix, zi, iz, si, is, hi, ih, wi, cz);; # Order 92160

a := [
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1],
    [-1, 0, 0, 0]];; # a in G2 = true

# Print(Order(G2));







