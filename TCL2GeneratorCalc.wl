(* ::Package:: *)

(* ::Title:: *)
(*TCL2GeneratorCalc*)


Get["TCLIntegrandCalcs.wl"]
Get["TCL4SpinBosonFunctions.wl"]


(*
TCL0Generator = K[0]
*)


DumpVar[TCL0Generator]


IntegrandTCL2 = Simplify[K[2]];


IntegrandTCL2


Simplify[TCL2Integrand - IntegrandTCL2]


TCL2TimeIntegrated = ParallelMap[TCL2TimeIntegral, IntegrandTCL2, {2}]


ResLoopSegWrapper[expr_]:= ResLoopSeg[expr, \[Omega],z, \[Omega]]
TCL2OmegaIntegral = Map[ResLoopSegWrapper,TCL2TimeIntegrated,{2}]/2


TCL2GeneratorV2 = FullSimplify[TCL2OmegaIntegral]


TCL2GeneratorV2 //MatrixForm


LoadVar[TCL2Generator]


Simplify[TCL2Generator-TCL2GeneratorV2]


(*
DumpVar[TCL2Generator]
*)
