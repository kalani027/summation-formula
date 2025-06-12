(*

Wolfram Mathematica version 13.3.1.0

This code reproduced figure 4 in the summation formula paper.

We test the behavior of the ratio of S[x,rho] + G[x,rho], the weighted sum and second sum,
over the mainterm for values of rho = 1,3,5,10. 

From theorem 7.1, we prove this for rho > 3.5. We observed that the behviour presist for smaller values of rho. 
However, the our summation formula does not extend to these due to convergece issues. 

This code was written by Olivia Beckwith and Kalani Thalagoda.

*)

(* basic number theory functions *)

IsDiscriminant[n_] := If[(Mod[n, 4] == 1) || Mod[n, 4] == 0, 1, 0]

SquarefreePart[n_] := 
 Times @@ Power @@@ ({#[[1]], Mod[#[[2]], 2]} & /@ FactorInteger[n])

FundamentalDiscriminant[n_] :=
 If[IsDiscriminant[n] == 1,
  If[Mod[SquarefreePart[n], 4] == 1,
   SquarefreePart[n],
   If[IntegerQ[Sqrt[n/4*SquarefreePart[n]]], 4*SquarefreePart[n],
    8*SquarefreePart[n]
    ]
   ]
  ]

(* summation forumula *)
sfunction[v_, j_, n_] := 
  Sum[KroneckerSymbol[(-1)^v*2^j, l]*Exp[2*\Pi*I*n*l/2^j], {l, 1, 2^
    j}]; 

chi[n_, a_, k_] :=
  If[GCD[a, 4*SquarefreePart[n]] == 1,
   KroneckerSymbol[ (-1)^((k + 1)/2)*4*SquarefreePart[n], a],
   0
   ];

T[n_, k_] :=
  Sum[Sum[
    If[OddQ[a] && OddQ[b] && Divisible[Sqrt[n/SquarefreePart[n]], a*b],
             MoebiusMu[a]*chi[n, a, k]*a^((1 - k)/2)*(b)^(2 - k), 0],
    {b, 1, Sqrt[n/SquarefreePart[n]]}],
   {a, 1, Sqrt[n/SquarefreePart[n]]}];

(* Dirichlet L-series for the character chi[n,a,k] above *)
Lseries[n_, s_, prec_, k_] := Sum[chi[n, i, k]/i^s, {i, 1, prec}];

A[n_, k_] := (1 + I^-k)/2^k + 
   1/2*Sum[(1 - I^-k)*2^-(k*j/2)*
       sfunction[1, j, n] + (1 + I^-k)*2^-(k*j/2)*
       sfunction[2, j, n], {j, 2, Ceiling[Log[2, Abs[n]] + 3]}];


(* fourier coeffient of weight -1/2 eisenstein series *)
aplus[n_, k_] :=  (-2*I)^(
   2 - k/2)*\[Pi] *((Lseries[n, (k - 1)/2, 100, 
        k]*(1 - chi[n, 2, k]*2^((1 - k)/2)))/(Zeta[
        k - 1]*(1 - 2^(1 - k))))*T[n, k]*A[n, k];

aminus[n_, k_] :=  (-2*I)^(
   2 - k/2)*\[Pi] *((Lseries[-1*n, (k - 1)/2, 100, 
        k]*(1 - chi[-1*n, 2, k]*2^((1 - k)/2)))/(Zeta[
        k - 1]*(1 - 2^(1 - k))*Gamma[k/2 - 1]))*T[-1*n, k]*A[-1*n, k];

(* k = -1/2 to agree with the notation from main theorem!*)

(* equation 2.5 in theorem 2.8 *)
g[rho_, n_, x_, k_] := 
  If[2*n/(x + n) < 10^-10, (2*\[Pi] * I )/
    Gamma[rho + k]*(1 + n/x)^rho * 
    Integrate[v^-k *(1 - v)^(k + rho - 1), {v, 0, 1}], (
    2*\[Pi] * I )/Gamma[rho + k]*(1 + n/x)^rho * 
    Integrate[v^-k *(1 - v)^(k + rho - 1), {v, 2*n/(x + n), 1}]];

(* asymptotic in equation 6.18 for k = -1/2 or kappa = 5 *)
(* a_plus summation *)
sum1[x_, rho_] := 
  Sum[aplus[n, 5]*(x - n)^rho, {n, 1, Floor[x]}]/
  Gamma[rho + 1]; 

(* a_minus sum *)
sum2[x_, rho_] := 
  x^rho/(2*\[Pi] * I )*
   Sum[aminus[n, 5]*g[rho, n, x, 2 - 5/2], {n, 1, Floor[x]}];

bminuszero[x_, rho_, k_] := 
 I^(k/2) * 2^(1 - (k/2) )* Cos[\[Pi] * k/4]

mainterm[x_, rho_, k_] := (2 * \[Pi] * (4)^(k/2 - 1)*I^k)/
   Gamma[rho + 2] * bminuszero[x, rho, 4 - 2 * k]  * 
   x^(rho + 1);

asym1plus[x_, rho_] := Re[sum1[x, rho] + sum2[x, rho]]/
  Re[mainterm[x, rho, -1/2]];
