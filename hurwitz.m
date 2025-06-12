
(* 

Wolfram Mathematica version 13.3.1.0

This code reproduced figure 1, 2 and 3 in the summation formula paper. 
These figures shows the asymptotic behavior of S[x,rho], the weighted sum, given in theorem 6.4.

Figure 1 shows the behaviour of the ratio of S[x,rho] over the mainterm 
for rho = 0.50001,1,1.5,2
Figure 2 shows the behavior of the ratio of the (S[x,rho] - mainterm) over first error term 
for rho = 2,5,10

Figure 3 plots the asymptotic behavior of S[x,rho] for rho = 0,-0.1,-0.5 . 
These rho values lie beyond the region of convergence given in theorem 6.4 

This code was written by Olivia Beckwith and Kalani Thalagoda.
*)

(* basic number theory functions *)
IsDiscriminant[n_] := If[(Mod[n, 4] == 1) || Mod[n, 4] == 0, 1, 0]
SquarefreePart[n_] := 
 Times @@ Power @@@ ({#[[1]], Mod[#[[2]], 2]} & /@ FactorInteger[n])

FundamentalDiscriminant[n_] := 
 If[IsDiscriminant[n] == 1, 
  If[Mod[SquarefreePart[n], 4] == 1, SquarefreePart[n], 
   If[IntegerQ[Sqrt[n/4*SquarefreePart[n]]], 4*SquarefreePart[n], 
    8*SquarefreePart[n]]]]


wunits[D_] := 
 If[IsDiscriminant[D] == 1 && D < 0, 
  If[FundamentalDiscriminant[D] == -4, Return[2], 
   If[FundamentalDiscriminant[D] == -3, Return[3], Return[1]]]]

squarepart[N_] := 
 If[IsDiscriminant[N] == 1, Sqrt[N/FundamentalDiscriminant[N]]]
squareterm[f_, d_, g_] := 
 DivisorSigma[1, f/g]*JacobiSymbol[d, g]*MoebiusMu[g]

(* hurwitz class number*)
Hurwitz[N_] := 
 If[N == 0, -1/12, 
  If[IsDiscriminant[-N] == 
    1, (NumberFieldClassNumber[Sqrt[-N]]/wunits[-N])*
    DivisorSum[squarepart[-N], 
     squareterm[squarepart[-N], FundamentalDiscriminant[-N], #] &], 0]]


(* summation formula *)
level := 4;
k := 3/2

(* coefficients in theorem 6.4 *)
apzero := -1/12;
amzero := 1/(8*\[Pi]);
bpzero := -(1 + I)/(12*Sqrt[8]);
bmzero := (3*(1 + I))/(16*\[Pi]*Sqrt[2]);

(* weighted sum *)
S[x_, rho_] := 
  Sum[Hurwitz[n]*(x - n)^rho, {n, 1, Floor[x]}]/Gamma[rho + 1];

(* main term in theorem 6.4*)
mainterm[x_, rho_] := 
  x^(rho + k)*(bpzero*I^(3/2)*(2*\[Pi])^(3/2))/(level^(3/4)* Gamma[rho + 5/2]);

(* first error term in theorem 6.4*)
errorterm1[x_, rho_] := 
  x^(rho + 1)*(2*\[Pi]*level^(3/4 - 1)*bmzero*I^(3/2))/ Gamma[rho + 2];


(* asymptotics *)
mainasym[x_, rho_] := S[x, rho]/ mainterm[x, rho];
secondasym[x_, rho_] := (S[x, rho] - mainterm[x, rho])/
  errorterm1[x, rho]

(* figure 1 *)
table1 := Table[{x, mainasym[1.001*x, 0.50001]}, {x, 1, 100}]; 
table2 := Table[{x, mainasym[1.001*x, 1]}, {x, 1, 100}];
table3 := Table[{x, mainasym[1.001*x, 1.5]}, {x, 1, 100}];
table4 := Table[{x, mainasym[1.001*x, 2]}, {x, 1, 100}];
ListPlot[{table1, table2, table3, table4}, PlotLegends -> {"rho=5.00001", "rho=1", "rho=1.5", "rho=2"}]

(* figure 2 *)

table5 := Table[{x, secondasym[1.001*x, 2]}, {x, 1, 100}];
table6 := Table[{x, secondasym[1.001*x, 5]}, {x, 1, 100}];
table7 := Table[{x, secondasym[1.001*x, 10]}, {x, 1, 100}];
ListPlot[{table5, table6, table7}, PlotLegends -> {"rho=2", "rho=5", "rho=10"}]


(* figure 3 *)

table8 := Table[{x, secondasym[1.001*x, 0]}, {x, 1, 100}];
table9 := Table[{x, secondasym[1.001*x, -0.1]}, {x, 1, 100}];
table10 := Table[{x, secondasym[1.001*x, -0.5]}, {x, 1, 100}];
ListPlot[table8, PlotLegends -> {"rho=0"}]
ListPlot[table9, PlotLegends -> {"rho=-0.1"}]
ListPlot[table10, PlotLegends -> {"rho=-0.5"}]
