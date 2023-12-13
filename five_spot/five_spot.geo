//--------------------------------------------------------------------
//---------------- UNIVERSIDADE FEDERAL DE PERNAMBUCO ----------------
//---------------- CENTRO DE TECNOLOGIA E GEOCIENCIAS ----------------
//---------- PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL -----------
//--------------------------------------------------------------------

//Work developed by: Marcio Souza and Luiz E. Queiroz
//Adviser Professors: Paulo Lyra & Darlan Carvalho
//Create date: 2023/12/13;	hour: 16:59h

//--------------------------------------------------------------------
//This file has CAD parameters. It is related to building of domain

//"cl1" corresponds to element size attributed in "Start.dat";
cl1 = 0.031300;

Point(1) = {0.000000, 0.000000, 0.000000, cl1};
Point(2) = {1.000000, 0.000000, 0.000000, cl1};
Point(3) = {1.000000, 1.000000, 0.000000, cl1};
Point(4) = {0.000000, 1.000000, 0.000000, cl1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(6) = {2, 3, 4, 1};
Plane Surface(6) = {6};
Physical Point(101) = {1, 2, 3, 4};
Physical Line(101) = {1, 2, 3, 4};
Physical Surface(1) = {6};

