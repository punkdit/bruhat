
G := GL(2,3);;
Print(Order(G), "\n");;

Hs := AllSubgroups(G);;
C := CyclicGroup(8);;

for H in Hs do 
    if Order(H)=8 and IsAbelian(H) then
        break;
    fi;
od;

Print("IsomorphismGroups ", IsomorphismGroups(C,H), "\n");
Print(H, "\n");

#fn := InducedClassFunction( chi, H );
Print("Irr(G)\n");
irr := Irr(G);
for chi in irr do
    Print(chi, "\n");
od;

Print("Irr(H)\n");
for chi in Irr(H) do
    Print(chi, "\n");
od;

Print("InducedClassFunctions\n");

ind := InducedClassFunctions( Irr( H ), G );

#Print(ind, "\n");


for i in ind do
    Print(i, "\n");
    for j in irr do
        u := ScalarProduct(i,j);
        Print(u, " ");
    od;
    Print("\n");
od;



