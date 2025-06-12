(* ::Package:: *)

(* ::Title:: *)
(*ProjectTCL4DynamicsForSpinBoson*)


(* ::Text:: *)
(*Following file has the better calculations of TCL4 dynamics for simple Hamiltonian. It also numerically evaluates it for Drude.*)


(* ::Chapter:: *)
(*General*)


Wiki[TCLIntegrandCalcs]


(* ::Section:: *)
(*TCL2*)


Wiki[TCL2GeneratorCalc]


(* ::Section:: *)
(*TCL4*)


(* ::Text:: *)
(*The idea is to do all the calcs again, neatly. From scratch.*)


Wiki[TCL4TripleTimeIntegral]
Wiki[TCL4OmegaIntegral]


Wiki[DQDTCL4GeneratorEvaluation]


Wiki[TCL4GeneratorEvaluationForDrude]


(* ::Chapter:: *)
(*Drude Calcs*)


(* ::Text:: *)
(*verify that the following result matches with our earlier results.*)


Wiki[TCL4DrudeGeneratorCalcV2]
Wiki[TCL2DrudeGeneratorCalc]


(* ::Chapter:: *)
(*Verification*)


Wiki[TCL4Verification]
Wiki[TCL4PrecisionVerification]
Wiki[TCL4VerificationPlotAgainstDrude]


(* ::Text:: *)
(*HEOM plot code is in this file.*)


Wiki[TCLVsHEOMFidelity]


(* ::Chapter:: *)
(*Important files*)


(* ::Section:: *)
(*Results*)


Get["TCL4Generator.mx"]
Get["TCL2Generator.mx"]
Get["TCL0Generator.mx"]

Get["TCL4DrudeGeneratorSimplified.mx"]
Get["TCL2DrudeGenerator.mx"]

Get["SteadyStateSpinBosonResult.mx"]


Wiki[TCL4DynamicsSpinBosonResults]
Wiki[O2SSCalcSpinBoson]

Wiki[NonMarkovianitySpinBoson]
Wiki[SBMProjectNMPlot]


(* ::Section:: *)
(*Functions*)


Wiki[MyFunctions]
Wiki[MathFuncMathematica]
Wiki[TCL4SpinBosonFunctions]
Wiki[DrudeCutoffFunctions]
Wiki[LatexExtract]


(* ::Section:: *)
(*Other useful files*)


Wiki[ArrayTCL4Verification]
Wiki[ArrayPlotTCLVsHEOMFidelity]
Wiki[TCLMarkovianityResonanceCondition]
