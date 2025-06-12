(* ::Package:: *)

(* ::Title:: *)
(*O2SSCalcSpinBoson*)


Get["TCL4Generator.mx"]
Get["TCL2Generator.mx"]
Get["TCL0Generator.mx"]
Get["TCL4SpinBosonFunctions.wl"]


(*===========================================================*)
(* 1.   removeIntegrals  \[Dash]   strip every Integrate wrapper   *)
(*===========================================================*)
ClearAll[removeIntegrals]
removeIntegrals[expr_] :=
  expr //. Integrate[integrand_, spec__] :> integrand   (*recursive //.  *)

(*===========================================================*)
(* 2.   splitByIntegrate  \[Dash]   build the 2-entry list         *)
(*===========================================================*)
ClearAll[splitByIntegrate]

splitByIntegrate[expr_] := Module[
  {terms, noInt, withInt},

  terms   = List @@ Expand[expr];            (* flatten \[CapitalSigma] into a term list *)

  noInt   =   Select[terms,  FreeQ[#, Integrate] &];                 (* \[CapitalSigma]\:2081 *)
  withInt =   Select[terms, !FreeQ[#, Integrate] &]       //          (* \[CapitalSigma]\:2082 \[Ellipsis] *)
              Map[removeIntegrals];                                  (* \[Ellipsis]strip *)

  {
    Total[noInt]   /. Plus[] :> 0,             (* 0 if list empty *)
    Total[withInt] /. Plus[] :> 0
  }
]


D[Tanh[x],x] /. x->0


(*
SpectralDensitySymmetryProperty = {J[0] -> 0, Derivative[1][f][0]-> 0, J[x_]-> f[x] Tanh[\[Beta] x/2], J'[x_]-> D[f[x] Tanh[\[Beta] x/2], x]}
*)

SpecZeroRule = {J[0] -> 0, Derivative[1][f][0]-> 0}

TCL4IntegralFormExpr[arr_]:= Module[{TCL4Generator2},
If[arr ==={0,0,0}, Return[0]];
TCL4Generator2 = arr //. SpecZeroRule;
With[{ZeroT = TCL4Generator2[[1]], integrandExpressionSingle = TCL4Generator2[[2]], 
integrandExpressionDouble = TCL4Generator2[[3]]},
  ZeroT + If[integrandExpressionSingle === 0, 0 ,HoldForm[Integrate[integrandExpressionSingle, {\[Omega], -Infinity, Infinity}]]] + 
  If[integrandExpressionDouble === 0, 0 ,HoldForm[Integrate[integrandExpressionDouble, {\[Omega]2, -Infinity, Infinity}, {\[Omega]1, -Infinity, Infinity}]]]
]
]


TCL0Generator;


TCL2Generator


L2 = Map[TCL4IntegralFormExpr, TCL2Generator, {2}]


L41 = Map[TCL4IntegralFormExpr, TCL4Generator, {2}];

L4 = modifyMatrix[L41, p3/p1];


L4[[2,1]];


L = TCL0Generator + \[Lambda]^2  L2 + \[Lambda]^4 L4


\[Rho]0={1,0,0,x0z}
\[Rho]2=\[Lambda]^2 {0,x2x,x2y,x2z}
\[Rho]4=\[Lambda]^4 {0,x4x,x4y,x4z}


\[Rho] = \[Rho]0 + \[Rho]2 + \[Rho]4
\[Rho] //MatrixForm


Final = L . \[Rho];


(* ::Chapter:: *)
(*Second order calculations*)


(* ::Section:: *)
(*Gibbs Population*)


(* ::Text:: *)
(*Gibbs population result comes from the second order eigenvalue equation. We subsequently use it to find the second order coherences.*)


G1 = Coefficient[Final[[4]], \[Lambda], 2]


G2 = G1 //. {J[-\[CapitalOmega]] -> -J[\[CapitalOmega]], f[x_] -> J[x] Coth[\[Beta] x/2]}


GSol = Solve[G2 == 0, x0z]


(* ::Section:: *)
(*Sigma x coherence*)


x2Dot = Final[[3]] //. GSol[[1]];


Solx2Dot = Coefficient[x2Dot, \[Lambda], 2]


ExpToTrig[Solx2Dot]


SigmaXSol = Solve[Solx2Dot == 0, x2x]


(* ::Text:: *)
(*There is an overall \[Omega] integral here.*)


xO2Val = x2x //. SigmaXSol[[1]]


xO2Val2 = xO2Val //. f[x_] -> J[x] Coth[\[Beta] x/2]


xO2Val3 = Collect[xO2Val2, J[\[Omega]], Simplify]


xO2Val4 = ReleaseHold[xO2Val3]


FullSimplify[xO2Val4]


xO2Val5 = splitByIntegrate[xO2Val4][[2]]


(* ::Subsection:: *)
(*Verifications*)


(* ::Text:: *)
(*Following is c4 from cresser anders. In the following, we do divide by a factor of 2 because the Cresser anders' intergral goes from 0 to Infinity whereas for us, it goes from -Infinity to Infinity.*)


CrAnx2 = SigmaXSpinBosonMFGS/2


CrAnx22 = CrAnx2 //. Sin[2 \[Theta]] -> 2 p1 p3


FullSimplify[ReleaseHold[xO2Val5 - CrAnx22]]


(* ::Text:: *)
(*Zero!*)


(* ::Section:: *)
(*Sigma  y  coherence*)


x1Dot = Final[[2]] //. Join[GSol[[1]], SigmaXSol[[1]]];


Solx1Dot = Coefficient[x1Dot, \[Lambda], 2]


SigmaYSol1 = Solve[Solx1Dot == 0, x2y]


(* ::Text:: *)
(*Just a constant?*)


SigmaYSol = FullSimplify[SigmaYSol1] //. {f[x_] -> J[x] Coth[\[Beta] x/2], J[-\[CapitalOmega]] -> - J[\[CapitalOmega]]}


(* ::Text:: *)
(*Consistent  with  Cresser  Anders!!!*)


(* ::Section:: *)
(*Diagonal population*)


(* ::Text:: *)
(*We use the second order coherences and the Gibbs population to find the second order population by going to the 4th order eigenvalue equation.*)


Clear[e4]


x3Dot = Final[[4]] //. Join[GSol[[1]], SigmaYSol[[1]], SigmaXSol[[1]]];


O4TermX3 = Coefficient[x3Dot, \[Lambda], 4]


X2zSol = Solve[O4TermX3==0, x2z]


x2zVal = x2z //. X2zSol[[1]]


(* ::Text:: *)
(*The following function separates the constant terms from the integrals by putting them in separate elements of an array. For the integral parts, we only save the integrand for ease of algebraic manipulation, since all the integral terms have the same variable and range of integration.*)


x2z0 = splitByIntegrate[x2zVal];


x2z0[[2]]


x2z1 = ReleaseHold[splitByIntegrate[x2zVal]]


(* ::Text:: *)
(*Since the contour of the integral is defined to go just above the real axis, we bring it back to the real axis by factoring in the residue contribution of the real axis poles. Why go back to the real axis? Because only then can we use the property like spectral density is an odd function to simplify the integral.*)


x2z2 = ResLoopSeg[x2z1[[2]] , \[Omega], z,\[Omega]]


(* ::Text:: *)
(*x2z2[[3]] below is the divergence test, which evaluates to zero, which is good. Divergence test is there just to make sure there are no fatal second order poles on the real line.*)


x2z2[[3]] //. {J[x_] -> f[x] Tanh[\[Beta] x/2], f[-\[CapitalOmega]] -> f[\[CapitalOmega]]}


(* ::Subsection:: *)
(*Constant Contribution*)


(* ::Text:: *)
(*Just subtracting the pole contribution above from the constant term in the 4th order eigenvalue equation for the populations.*)


ConstantCOntri = x2z1[[1]] - x2z2[[1]]


(* ::Text:: *)
(*Following assumption can now be made.*)


SpecSymRule = {J[x_] -> f[x] Tanh[\[Beta] x/2], J'[x_] -> D[f[x] Tanh[\[Beta] x/2], x], f[-\[Omega]+\[CapitalOmega]] -> f[\[Omega]-\[CapitalOmega]], f[-\[Omega]-\[CapitalOmega]] -> f[\[Omega]+\[CapitalOmega]], f[-\[CapitalOmega]] -> f[\[CapitalOmega]], f[-\[Omega]]-> f[\[Omega]]}


Simplify[ConstantCOntri /. Join[SpecSymRule, SpecZeroRule]]


(* ::Text:: *)
(*The constant term is zero! Now lets get back to the integral terms.*)


(* ::Subsection:: *)
(*Integral term*)


(* ::Text:: *)
(*Again, since the contour is now on the real axis, we can now assume that J[x] = - J[-x]. Note the factor of 1/2 I did below to compensate for the double counting during addition below.*)


IntContri = Simplify[Simplify[(x2z2[[2]] + (x2z1[[2]] /. \[Omega] -> -\[Omega]))/2] //. {J[x_] -> f[x] Tanh[\[Beta] x/2], 
							J'[x_] -> D[f[x] Tanh[\[Beta] x/2], x], f[-\[Omega]+\[CapitalOmega]] -> f[\[Omega]-\[CapitalOmega]], f[-\[Omega]-\[CapitalOmega]] -> f[\[Omega]+\[CapitalOmega]], 
							f[-\[CapitalOmega]] -> f[\[CapitalOmega]], f[-\[Omega]]-> f[\[Omega]]}]


FreeQ[IntContri, J]
FreeQ[IntContri, p3]


IntContri2 = IntContri //. {f'[-x_]-> -f'[x]}


(* ::Text:: *)
(*From the MFGS result, we know that the final result cannot depend upon p3. So lets verify whether the p3 dependent term disappears or not.*)


(* ::Subsubsection:: *)
(*p3 dependent part*)


P3Part = CoeffTerm[IntContri2, p3^2]


(* ::Text:: *)
(*There are 2 types of terms here, with factor f[\[Omega]] f[\[Omega]-\[CapitalOmega]] and f[\[Omega]] f[\[Omega]+\[CapitalOmega]]. But since the integral is from \[Omega] = -Infinity -> Infinity, a change in variable \[Omega] -> \[Omega] + \[CapitalOmega] for the first type of term will make both the terms identical and hopefully cancel. let's check this out.*)


P3Part2 = Collect[P3Part, {f[\[Omega]] f[\[Omega]-\[CapitalOmega]], f[\[Omega]] f[\[Omega]+\[CapitalOmega]]}]


OmShift[T_]:= T /. \[Omega] -> \[Omega] + \[CapitalOmega]

P3Part2Part1 = ConditionalOp[Expand[P3Part2], {f[\[Omega]] f[\[Omega]-\[CapitalOmega]]}, OmShift]


Simplify[P3Part2Part1]


(* ::Text:: *)
(*p3 dependent part cancels, as expected. Now only p1 dependent terms remain.*)


(* ::Subsubsection:: *)
(*P1 dependent part*)


P1Part = CoeffTerm[IntContri2, p1^2]


Simplify[P3Part + P1Part - IntContri2]


p1Pop2 = Simplify[P1Part /. f[x_]-> J[x] Coth[\[Beta] x/2]]


(* ::Subsection:: *)
(*Verification with MFGS result*)


CresserAndersPop = SigmaZSpinBosonMFGS/ (2 \[Lambda]^2) //. {Sin[\[Theta]]^2 -> p1^2, \[Omega]q -> \[CapitalOmega]}
diff = (CresserAndersPop - p1Pop2);


Simplify[diff]


(* ::Text:: *)
(*The difference with the MFGS result is zero, as expected!*)


SigmaZSol = {With[{expr = p1Pop2}, {x2z -> HoldForm[Integrate[expr, {\[Omega], -Infinity, Infinity}]]}]}


(* ::Chapter:: *)
(*4th order coherences*)


(* ::Text:: *)
(*Finally, let us calculate the 4th order sigmaX coherence.*)


x4Dot = Final[[3]] //. Join[ GSol[[1]], SigmaXSol[[1]], SigmaYSol[[1]], sigmaZSol[[1]]];


Solx4Dot = Coefficient[x4Dot, \[Lambda], 4];


sigmax4Sol = Solve[Solx4Dot == 0, x4x];


xO4Val = x4x //. sigmax4Sol[[1]];


xO4Val2 = ReleaseHold[xO4Val]


xO4Val3 = Simplify[xO4Val2]


(* ::Text:: *)
(*It is complicated. As for 4th order SigmaX coherence, as per argument given in our [earlier paper](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.111.115423), it has to be zero. Verifying this using the above expressions is left as an exercise for the readers.*)


SigmaX4Sol = {{x4x -> xO4Val3}};


SteadyStateSpinBosonResult = {SigmaXSol, SigmaYSol, SigmaZSol, SigmaX4Sol}


DumpVar[SteadyStateSpinBosonResult]



