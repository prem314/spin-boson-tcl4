(* ::Package:: *)

(* ::Title:: *)
(*Outline*)


(* ::Text:: *)
(*Run the following code to enable hyperlinking between notebooks in this folder for easier navigation.*)


(* Set the Wiki function to hold its first argument unevaluated *)
SetAttributes[Wiki, HoldFirst]

(* Define the Wiki function to convert the symbol's name to a string and create the hyperlink *)
Wiki[symbol_Symbol] := Module[{name},
  name = SymbolName[symbol];
  Hyperlink[name, name <> ".wl"]
]


(* ::Chapter:: *)
(*General*)


Wiki[TCL4OperatorFormDerivation]


Wiki[TCLIntegrandCalcs]


(* ::Section:: *)
(*TCL2*)


Wiki[TCL2GeneratorCalc]


(* ::Section:: *)
(*TCL4*)


Wiki[TCL4TripleTimeIntegral]
Wiki[TCL4OmegaIntegral]


Wiki[DQDTCL4GeneratorEvaluation]


Wiki[TCL4GeneratorEvaluationForDrude]


(* ::Chapter:: *)
(*Drude Calcs*)


Wiki[TCL4DrudeGeneratorCalcV2]
Wiki[TCL2DrudeGeneratorCalc]


(* ::Chapter:: *)
(*Verification*)


Wiki[TCL4Verification]
Wiki[TCL4PrecisionVerification]
Wiki[TCL4VerificationPlotAgainstDrude]


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

Get["datasimpletcl2_gamma01_thetapi2_modifiedtcl2_400_antipodal.mx"]
Get["datasimpletcl4_gamma01_thetapi2_modifiedtcl4_400_antipodal.mx"]

Get["datasimpletcl2_gamma01_thetapi2_modifiedtcl2_4096.mx"]
Get["datasimpletcl4_gamma01_thetapi2_modifiedtcl4_4096.mx"]


Wiki[TCL4DynamicsSpinBosonResults]
Wiki[O2SSCalcSpinBoson]


Wiki[AntipodalTCL2NM]
Wiki[AntipodalTCL4NM]
Wiki[AntipodalSBMProjectNMPlot]


Wiki[TraceDistanceDynamicsPlot]


Wiki[GeneralTCL2NM]
Wiki[GeneralTCL4NM]
Wiki[SBMProjectNMPlot]


(* ::Section:: *)
(*Functions*)


Wiki[MyFunctions]
Wiki[MathFuncMathematica]
Wiki[TCL4SpinBosonFunctions]
Wiki[DrudeCutoffFunctions]


(* ::Section:: *)
(*Other useful files*)


Wiki[TCLMarkovianityResonanceCondition]
