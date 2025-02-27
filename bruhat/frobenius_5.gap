

G := GL(2,5);;
Print(Order(G), "\n");;

#Hs := AllSubgroups(G);;
#C := CyclicGroup(24);;
#
#for H in Hs do 
#    if Order(H)=24 and IsAbelian(H) then
#        Print("IsomorphismGroups ", IsomorphismGroups(C,H), "\n");
#        J := H;;
#    fi;
#od;
#
#Print(J, "\n");

#quit;

x := Indeterminate(GF(5), "x");
m := CompanionMatrix(x^2+2);
m in G;
H := Centralizer(G, m);
Order(H);


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



