(* ::Package:: *)

(* ::Title:: *)
(*TCL4VerificationPlotAgainstDrude*)


Get["MyFunctions.wl"]
Get["DrudeCutoffFunctions.wl"]

LoadVar[DrudeVsGeneralPlotData]


TCL4GenPlot = VerificationPlotCombined[DrudeVsGeneralPlotData, {{3,1},{3,2},{3,3},{3,4},{4,1},{4,2},{4,3},{4,4}}]


TCL4DynamicsFolder = "/home/premkr/Dropbox/work/projects/tcl4_dynamics"
fullPath = FileNameJoin[{TCL4DynamicsFolder,"TCL4Verification.pdf"}]
Export[fullPath, TCL4GenPlot, "PDF"];




NumVals = DrudeVsGeneralPlotData[[2]]


Get["LatexExtract.wl"]

DrudeVerificationParameters = {
  {\[Gamma], NumVals[[1]]},
  {\[CapitalLambda], NumVals[[2]]},                 (* Assuming tc represents t_c *)
  {\[CapitalOmega], NumVals[[3]]},
  {\[Beta], NumVals[[4]]},
  {\[Theta], NumVals[[5]]}
};

ParameterToLatex[DrudeVerificationParameters]

