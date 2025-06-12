(* ::Package:: *)

(* ::Title:: *)
(*TCL4DrudeGeneratorCalcV2*)


$RecursionLimit = 1024 * 1024;


Get["TCL4SpinBosonFunctions.wl"]
Get["DrudeCutoffFunctions.wl"]
Get["MyFunctions.wl"]


Get["TCLIntegrandCalcs.wl"]
TCL4Integrand2 = Simplify[K[4]];


For[i =1, i< 5, i++, TCL4Integrand2[[2,i]]=0];
TCL4Integrand2;


Dimensions[TCL4Integrand2]


(* ::Chapter:: *)
(*Main*)


(*
I am using ZeroTermShift here MAJORLY because for n,m=0, Exp[-nt] does not decay off.

Also because I have plans to do these double summations numerically.
So it will pay to define the summations from 1 to Ininifty rather than -Infinity to Infinity.
*)

ZeroTermShift[ar_] := Module[{newar},
  newar = ar;
  newar = ReplacePart[newar, {2, 1} -> (newar[[2,1]] + (newar[[2,2]] //. n -> 0) + (newar[[1,2]] //. {m -> 0, n -> 0}))];
  newar = ReplacePart[newar, {2, 2} -> (newar[[2,2]] + ((newar[[1,2]] //. n -> 0)//. m -> n) + (newar[[1,2]] //. m -> 0))];
  newar = ReplacePart[newar, {2, 2} -> 2 newar[[2,2]]];
  newar = ReplacePart[newar, {1, 2} -> 4 newar[[1,2]]];
  newar
]

tripleIntWrapper[expr_]:= TCL4TripleInt[expr, DrudeAssumptions]

SingleDrudeTermEval[expr_]:= Module[{IntegrandSeg,IntegrandSeg2, IntegrandSeg3, SegTripleInt, SegTripleInt2,NiceTerm},
IntegrandSeg = decomposeExpr[expr, \[Eta], \[Nu]];
Print["After decomposeExpr (IntegrandSeg): ", IntegrandSeg];
IntegrandSeg2 = (IntegrandSeg //. DrudeNuEtaReplace) //. SpecialSymbolValInverse;
Print["After DrudeNuEtaReplace (IntegrandSeg2): ", IntegrandSeg2];
IntegrandSeg3 = ZeroTermShift[IntegrandSeg2];
Print["After ZeroTermShift (IntegrandSeg3): ", IntegrandSeg3];
NiceTerm = ExprNiceForm[IntegrandSeg3];
Print["After ExprNiceForm (NiceTerm): ", NiceTerm];
SegTripleInt = Map[tripleIntWrapper, NiceTerm,{2}];
Print["After Map[tripleIntWrapper] (SegTripleInt): ", SegTripleInt];
SegTripleInt
]


(*
TCL4DrudeGeneratorSmallT = ParallelMap[SingleDrudeTermEval, TCL4Integrand2, {2}]
DumpVar[TCL4DrudeGeneratorSmallT]
*)


LoadVar[TCL4DrudeGeneratorSmallT]


FreeQ[TCL4DrudeGeneratorSmallT,\[Beta]]


Length[Dimensions[TCL4DrudeGeneratorSmallT]]


TCL4DrudeGeneratorSmallT2 = ParallelMapFinest[StabilizeDrudeTripleInt, TCL4DrudeGeneratorSmallT]


(* ::Text:: *)
(*Use either of the following 2:*)


TCL4DrudeGeneratorLargeTSimpler = ParallelMapFinest[LimitWrapper, TCL4DrudeGeneratorSmallT2]
(*
DumpVar[TCL4DrudeGeneratorLargeTSimple]
*)


(*
ArrayEvalCheckpoint[Simplify, TCL4DrudeGeneratorLargeTSimpler //. SpecialSymbolVal, TCL4DrudeGeneratorSimplified, 4]
*)


(*
ArrayEvalCheckpoint[LimitWrapper, TCL4DrudeGeneratorSmallT2, TCL4DrudeGeneratorLargeTSimpler, 4]
*)


LoadVar[TCL4DrudeGeneratorSimplified]
LoadVar[TCL4DrudeGeneratorNiceForm]


(*
randomValue[(TCL4DrudeGeneratorNiceForm - TCL4DrudeGeneratorLargeTSimpler) //. SpecialSymbolValInverse][[1]]// MatrixForm
*)


Simplify[(TCL4DrudeGeneratorNiceForm - TCL4DrudeGeneratorSimplified) //. SpecialSymbolVal]


(* ::Text:: *)
(*Output = {{{{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}}, {{{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}}, {{{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}}, {{{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}}}*)
