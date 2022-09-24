#!/usr/bin/env gap

# See: arxiv 2101.09349 


LoadPackage("recog");
LoadPackage("LINS");
F := FreeGroup("l", "r", "t"); ;
AssignGeneratorVariables(F ); ;
G := F/[l^2,r^2,t^2,(l*r)^5,(l*t)^2,(r*t)^4];;
gr:=LowIndexNormalSubgroupsSearchForAll(G,200);
L := List(gr);
S := Grp(L[9]);
IsNormal(G,S);

#GeneratorsOfGroup(S);
#[ l*(r*l*t)^2*r*(t^-1*l^-1*r^-1)^2*t^-1, (l*t*r)^2*l*t*(r^-1*t^-1*l^-1)^2*r^-1,
#  (r*t*r*l)^2*(r^-1*t^-1*r^-1*l^-1)^2, (t*r*l*r)^2*(t^-1*r^-1*l^-1*r^-1)^2,
#  r*l*(r*l*t)^2*r*(t^-1*l^-1*r^-1)^2*t^-1*r^-1,
#  r*t*(r*l)^2*t*r*l*t^-1*r^-1*t^-1*(l^-1*r^-1)^2*t^-1 ]



