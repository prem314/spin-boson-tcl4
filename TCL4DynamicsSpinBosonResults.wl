(* ::Package:: *)

(* ::Title:: *)
(*TCL4DynamicsSpinBosonResults*)


Get["TCL4Generator.mx"]
Get["TCL2Generator.mx"]
Get["TCL4DrudeGeneratorSimplified.mx"]
Get["TCL2DrudeGenerator.mx"]


Get["DrudeCutoffFunctions.wl"]


(* ::Chapter:: *)
(*General result*)


LatexRule = {J[0] -> 0, Derivative[1][f][0]-> 0, p1 ->  Subscript[a, 1], p3 ->  Subscript[a, 3], \[Omega]1 ->  Subscript[\[Omega], 1], \[Omega]3 ->  Subscript[\[Omega], 3]};


TCL4IntegralForm[TCL4Generator_, i_,j_]:= Module[{TCL4Generator2},
TCL4Generator2 = TCL4Generator //. LatexRule;
With[{ZeroT = TCL4Generator2[[i,j]][[1]], integrandExpressionSingle = TCL4Generator2[[i,j]][[2]], 
integrandExpressionDouble = TCL4Generator2[[i,j]][[3]]},
  HoldForm[ZeroT + HoldForm[Integrate[integrandExpressionSingle, {\[Omega], -Infinity, Infinity}]]] + 
  If[integrandExpressionDouble === 0, 0 ,HoldForm[Integrate[integrandExpressionDouble, {\[Omega]2, -Infinity, Infinity}, {\[Omega]1, -Infinity, Infinity}]]]
]
]

TCL4IntegralForm[TCL4Generator, 3,1]


(* ::Section:: *)
(*TCL4*)


Do[Print[Subscript[Superscript[F,{4}], {i,j}], " = ",TCL4IntegralForm[TCL4Generator, i,j], "\n"], {i,3,4},{j,1,4}]


(* ::Section:: *)
(*TCL2*)


Do[Print[Subscript[Superscript[F,{2}], {i,j}], " = ",TCL4IntegralForm[TCL2Generator, i,j], "\n"], {i,3,4},{j,1,4}]


(* ::Chapter:: *)
(*Drude result*)


TCL4DrudeSummationForm[TCL4Generator_, i_,j_]:= Module[{TCL4Generator2},
TCL4Generator2 = TCL4Generator //. Join[LatexRule, SpecialSymbolVal];
With[{ZeroT = TCL4Generator2[[i,j]][[2,1]], integrandExpressionSingle = TCL4Generator2[[i,j]][[2,2]], 
integrandExpressionDouble = TCL4Generator2[[i,j]][[1,2]]},
  HoldForm[ZeroT + HoldForm[Sum[integrandExpressionSingle, {n, 1, Infinity}]]] + 
  If[integrandExpressionDouble === 0, 0 ,HoldForm[Sum[integrandExpressionDouble, {m, 1, Infinity}, {n, 1, Infinity}]]]
]
]

TCL4DrudeSummationForm[TCL4DrudeGeneratorSimplified, 3,1]


(* ::Section:: *)
(*TCL2*)


Do[Print[Subscript[Superscript[F,{2}], {i,j}], " = ",TCL2DrudeGenerator[[i,j]], "\n \n \n"], {i,3,4},{j,1,4}]


(* ::Section:: *)
(*TCL4*)


Do[Print[Subscript[Superscript[F,{4}], {i,j}], " = ",TCL4DrudeSummationForm[TCL4DrudeGeneratorSimplified, i,j], "\n \n \n"], {i,3,4},{j,1,4}]



(* ::Section:: *)
(*Coherences*)


Get["SteadyStateSpinBosonResult.mx"]


(* ::Text:: *)
(*The following expression takes a little time to load.*)


SteadyStateSpinBosonResult;


SteadyStateSpinBosonResult[[1]]


SteadyStateSpinBosonResult[[2]]


SteadyStateSpinBosonResult[[3]]


SteadyStateSpinBosonResult[[4]]
