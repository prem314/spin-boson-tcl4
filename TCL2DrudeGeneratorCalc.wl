(* ::Package:: *)

(* ::Title:: *)
(*TCL2DrudeGeneratorCalc*)


Get["DrudeCutoffFunctions.wl"]
Get["MyFunctions.wl"]


Get["TCLIntegrandCalcs.wl"]
TCL2Integrand2 = Simplify[K[2]];


TCL2Integrand2[[4,2]]


ZeroTermShift[ar_] := Module[{newar},
  newar = ar;
  newar = ReplacePart[newar, {2, 1} -> (newar[[2,1]] + (newar[[1,2]] //. n -> 0))];
  newar = ReplacePart[newar, {1, 2} -> 2 newar[[1,2]]];
  newar
]


SingleDrudeTermEval[expr_]:= Module[{IntegrandSeg, SegTripleInt, IntegrandSeg2, IntegrandSeg1},
IntegrandSeg = decomposeExpr[expr, \[Eta], \[Nu]];
IntegrandSeg1 = (IntegrandSeg //. SingleDrudeNuEtaReplace) //. SpecialSymbolValInverse;
IntegrandSeg2 = ZeroTermShift[IntegrandSeg1];
Print["After decomposeExpr (IntegrandSeg): ", IntegrandSeg2];
SegTripleInt = Map[TCL2DrudeIntegral, IntegrandSeg2,{2}];
Print["After Map[tripleIntWrapper] (SegTripleInt): ", SegTripleInt];
SegTripleInt
]


decomposeExpr[ -4 p1 p3 \[Nu][t-t1], \[Eta], \[Nu]]


SingleDrudeTermEval[TCL2Integrand2[[4,2]]]


ar = {{a,n + 1},{b,c}}
ZeroTermShift[ar]


TCL2DrudeGeneratorSmallT = ParallelMap[SingleDrudeTermEval, TCL2Integrand2, {2}]


TCL2DrudeGeneratorSmallT[[4,2]]


TCL2DrudeGeneratorSmallT2 = ParallelMapFinest[StabilizeDrudeTripleInt, TCL2DrudeGeneratorSmallT]


TCL2DrudeGeneratorSmallT2[[4,2]]


TCL2DrudeGeneratorLargeTSimpler = ParallelMapFinest[LimitWrapper, TCL2DrudeGeneratorSmallT2]


TCL2Summation[ar_] := (ar[[2,1]] + Sum[ar[[1,2]],{n,1,Infinity},Assumptions->DrudeAssumptions]) //. SpecialSymbolVal

ArrayEvalCheckpoint[TCL2Summation, TCL2DrudeGeneratorLargeTSimpler, TCL2DrudeGenerator, 2]


LoadVar[TCL2DrudeGenerator]


TCL2DrudeGenerator2 = TCL2DrudeGenerator //. P1P3ToTheta;


Simplify[Tcl2drude - TCL2DrudeGenerator2]
