(* ::Package:: *)

(* ::Title:: *)
(*TCL4PrecisionVerification*)


Get["MyFunctions.wl"]
Get["DrudeCutoffFunctions.wl"]
Get["TCL4SpinBosonFunctions.wl"]


LoadVar[TCL4Generator]


(*
GammaNumVal       = RandomReal[{0.1,1}];
LambdaNumVal      = RandomReal[{0.1,1}];
\[CapitalOmega]Val              = RandomReal[{0.1,1}];
BetaNumVal        = RandomReal[{0.1,1}];
\[Theta]NumVal           = RandomReal[{0.1,1}];
*)


{GammaNumVal, LambdaNumVal, \[CapitalOmega]Val, BetaNumVal, \[Theta]NumVal} = {0.431531423940297`,0.17796109372220548`,0.5270521973258038`,0.5997509586470544`,0.8272727677434009`}

p1NumVal          = Sin[\[Theta]NumVal];
p3NumVal          = Cos[\[Theta]NumVal];

PrecisionOfNumericalIntegration = 17;


DrudeRandomeNumReplace[expr_] := expr //. {
  p1 -> p1NumVal,
  p3 -> p3NumVal,
  \[CapitalLambda] -> LambdaNumVal,
  \[Beta] -> BetaNumVal,
  \[Gamma] -> GammaNumVal,
  \[Theta] -> \[Theta]NumVal,
  \[CapitalOmega] -> \[CapitalOmega]Val,
  f[0] -> Limit[fDrude[\[Omega]], \[Omega] -> 0], f'[0] -> 0,
  J -> Function[\[Omega], JDrude[\[Omega]]],
  J[0] -> Limit[JDrude[\[Omega]], {\[Omega] -> 0}],
  f -> Function[\[Omega], fDrude[\[Omega]]]
};


(* ::Chapter:: *)
(*Relative error evaluation*)


TCL4GeneratorNum = TCLNumericalEval[DrudeRandomeNumReplace[TCL4Generator], {LambdaNumVal, 2 \[Pi] /BetaNumVal}, DrudeRandomeNumReplace[p3/p1], PrecisionOfNumericalIntegration];


TCL4GeneratorNum // MatrixForm


LoadVar[TCL4DrudeGeneratorSimplified]


TCL4SimpleHDrude = TCL4DrudeGeneratorSimplified /. m-> m + 10^-15 I;


TCL4SimpleHDrude2 = DrudeRandomeNumReplace[TCL4SimpleHDrude //. SpecialSymbolVal];


TCL4SimpleHDrude2[[3,1]] // MatrixForm


AnaDrudeSum[ar1_, i_, j_]:=Module[{ar, SegTripleInt3, FinalResult},
ar = ar1[[i,j]];
SegTripleInt3 = Simplify[ar];
FinalResult = SegTripleInt3[[2,1]] + Sum[SegTripleInt3[[2,2]], {n,1, Infinity}];
FinalResult
]

SimpD31[i_,j_]:= AnaDrudeSum[TCL4DrudeGeneratorSimplified, i,j] //. SpecialSymbolVal

RelErrorOneSum[i_,j_]:= Abs[(DrudeRandomeNumReplace[SimpD31[i,j]] - TCL4GeneratorNum[[i, j]])/TCL4GeneratorNum[[i, j]]]


RelErrorOneSum[4,1]


RelErrorOneSum[3,1]


(* ::Chapter:: *)
(*Bugs encountered*)


(* ::Text:: *)
(*High precision summation is required for consistency. See:*)


Wiki[MathematicaSummationBug]
