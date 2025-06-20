(* ::Package:: *)

(* ::Title:: *)
(*DrudeCutoffFunctions*)


Get["MathFuncMathematica.wl"]


(* Define a function to save an array in Python list format (without assignment) *)
SaveArrayAsPythonList[array_, filename_String] := Module[
  {pythonList, filepath},
  
  (* Convert the Mathematica list to a string and replace curly braces with square brackets *)
  pythonList = StringReplace[ToString[array, InputForm], {"{" -> "[", "}" -> "]"}];
  
  (* Create the full path: file saved in the same folder as the current notebook *)
  filepath = FileNameJoin[{NotebookDirectory[], filename}];
  
  (* Export the Python-formatted list to a text file in that folder *)
  Export[filepath, pythonList, "Text"]
];


VerificationPlotCombined[SimpleHDrudeVsGeneralPlotData1_, idx_] :=
  Module[{dataToPlot, legends, SimpleHDrudeVsGeneralPlotData, PlotParameters},
  
  SimpleHDrudeVsGeneralPlotData = SimpleHDrudeVsGeneralPlotData1[[1]];
  
  
  
	PlotParameters = SimpleHDrudeVsGeneralPlotData1[[2]];

	Print["Generating plot for \[Gamma] = ", PlotParameters[[1]], ", \[CapitalLambda] = ", PlotParameters[[2]], ", \[CapitalOmega] = ", PlotParameters[[3]], ", \[Beta] = ", PlotParameters[[4]], ", \[Theta] = ", PlotParameters[[5]],"." ];

    dataToPlot = Table[
      Abs[SimpleHDrudeVsGeneralPlotData[[i]]],
      {i, Length@SimpleHDrudeVsGeneralPlotData}
    ];

    (* Prepare the list of legends directly as expressions *)
    legends = Table[
      Superscript[Subscript["F", idx[[i]] - {1,1} ], "(4)"],  (* REMOVED MakeBoxes[] wrapper *)
      {i, Length@SimpleHDrudeVsGeneralPlotData}
    ];

    plot = ListLinePlot[
      dataToPlot,
      PlotMarkers      -> Automatic,
      PlotLegends      -> legends,
      LabelStyle -> {Black, 24},               (* Provide the list of expressions *)
      AspectRatio->0.7,
      Frame            -> True,
      FrameLabel       -> {Style["No. of Matsubara Terms", 24], Style["Relative Difference", 24]},
      FrameStyle -> Directive[Black, AbsoluteThickness[1]],
      FrameTicks-> {{All,None},{{10,100,1000},None}},
      FrameTicksStyle -> Directive[Black, 20],
      PlotRange        -> All,
      ImageSize        -> Large,
      BaseStyle -> {FontSize -> 20}, 
      ScalingFunctions -> {"Log", "Log"}
    ];
    plot
  ]


(*\[Mu]0 = 2 \[Pi] kT;  *)

(*I removed an extra factor of 1/\[CapitalOmega] from here.*)

\[Mu][n_] := \[Mu]0 n
NuDrude[n_,t_] := \[Kappa] (\[CapitalLambda] Exp[-\[CapitalLambda] t] - \[Mu][n] Exp[- \[Mu][n] t])/(\[CapitalLambda]^2 - \[Mu][n]^2);(*I think there is a mod over n and m here inside the exponent. This wil resolve the divergence issue.*)
EtaDrude[t_]:= - \[Beta] \[Kappa]/2  Exp[- \[CapitalLambda] t]
(* - sign error in the paper??? No. The - factor should be there!*)



DrudeNuEtaReplace = {\[Nu][t_] \[Nu][t1_] -> NuDrude[n,t] NuDrude[m,t1], \[Nu][t_] \[Eta][t1_] -> NuDrude[n,t] EtaDrude[t1], \[Eta][t_] \[Eta][t1_] -> EtaDrude[t] EtaDrude[t1]};

SingleDrudeNuEtaReplace = {\[Nu][t_] -> NuDrude[n,t] , \[Eta][t1_] -> EtaDrude[t1]}


DrudeAssumptions = {t>0, p1 \[Element] Reals, p3 \[Element] Reals, \[Mu]0 > 0, m > 0, n > 0, \[CapitalLambda] > 0, \[CapitalOmega] > 0, \[Kappa] >0, \[Beta] > 0, \[Gamma] > 0} 

(* The assumption I made on \[CapitalOmega] here is necessary. Mathematica was fucking something up. See below*)


StabilizeDrudeTripleInt[expr_]:= ReleaseHold[Expand[TrigToExp[expr]] /. Exp[x_] -> Exp[HoldForm[Simplify[x, Assumptions->DrudeAssumptions]]]]


GenToDrudeRule = {f[0] -> Limit[fDrude[\[Omega]], \[Omega]->0], J[x_] -> JDrude[x], f[x_]-> fDrude[x], p1 -> Sin[\[Theta]], p3 -> Cos[\[Theta]], J'[x_] -> D[JDrude[x],x], f'[x_] -> D[fDrude[x],x]};


CothFunRes[Integrand_, \[Omega]_, n_]:= Residue[Integrand, {\[Omega], I 2 \[Pi] n /\[Beta]},Assumptions->{Element[n, Integers], \[Beta] >0 }]


DrudeInt[expr_]:= Module[{LContri, OmContri, CothContri, CothPole, Result},
	LContri = 2 \[Pi] I Rsidu[expr, \[Omega], {I \[CapitalLambda]}];
	OmContri = \[Pi] I Rsidu[expr, \[Omega], {\[CapitalOmega], -\[CapitalOmega]}];
	CothContri = 2 \[Pi] I Sum[CothFunRes[expr, \[Omega], n], {n, 1, Infinity}];
	
	Print["Coth contri = ", CothContri];
	
	Result = LContri + OmContri + CothContri;
	Result/2
]


(* ::Text:: *)
(*See the following for verification of the above expressions.*)


Wiki[CalcsDrudeCutoff]


SpecialSymbolVal = {\[Mu]0 -> 2 \[Pi]/\[Beta], \[Kappa] -> \[Pi] \[Gamma] \[CapitalLambda]^2/\[Beta]};
SpecialSymbolValInverse = {\[Beta] -> 2 \[Pi]/\[Mu]0};
(*Note the extra \[Pi] factor here as compared with Kappler. It is there because we are following Cresser's convention for the spectral density.*)

JDrude[\[Omega]_] := \[Gamma] \[Omega] \[CapitalLambda]^2/(\[CapitalLambda]^2 + \[Omega]^2)
fDrude[\[Omega]_] := JDrude[\[Omega]] Coth[\[Beta] \[Omega]/2]


LimitWrapper[expr_]:= Limit[expr, {t-> Infinity}, Assumptions->DrudeAssumptions]


TCL2DrudeIntegral[expr_]:= Integrate[expr , {t1, 0,t}, Assumptions -> DrudeAssumptions]


P1P3ToTheta = {p1 -> Sin[\[Theta]], p3 -> Cos[\[Theta]]};


(* ::Chapter:: *)
(*Simple Hamiltonian*)


(* ::Text:: *)
(*All the following calculations assumed a1 = 1/2 in the interaction Hamiltonian.*)


bz2Drude0 =-(\[Beta] \[Kappa] \[CapitalOmega])/(2 (\[CapitalLambda]^2+\[CapitalOmega]^2));

azz2Drude0 =-((\[Pi] \[Kappa] \[CapitalOmega] Coth[(\[Pi] \[CapitalOmega])/\[Mu]0])/(\[CapitalLambda]^2 \[Mu]0+\[Mu]0 \[CapitalOmega]^2));

azz4Drude0 =1/(2 \[CapitalLambda] \[Mu]0^3 (\[CapitalLambda]^2+\[CapitalOmega]^2)^3) \[Kappa]^2 \[CapitalOmega] (4 \[Mu]0^3 \[CapitalOmega]-I \[CapitalLambda] \[Mu]0 (\[CapitalLambda]^2-3 \[CapitalOmega]^2) PolyGamma[0,-((I \[CapitalOmega])/\[Mu]0)]^2-I \[Mu]0^2 (\[CapitalLambda]-I \[CapitalOmega]) (\[CapitalLambda]-3 I \[CapitalOmega]) PolyGamma[0,(I \[CapitalOmega])/\[Mu]0]+I \[CapitalLambda] \[Mu]0 (\[CapitalLambda]^2-3 \[CapitalOmega]^2) PolyGamma[0,(I \[CapitalOmega])/\[Mu]0]^2-I \[Mu]0 (\[CapitalLambda]-I \[CapitalOmega]) (\[CapitalLambda]+I \[CapitalOmega])^2 PolyGamma[1,-((I \[CapitalOmega])/\[Mu]0)]+PolyGamma[0,-((I \[CapitalOmega])/\[Mu]0)] (I \[Mu]0^2 (\[CapitalLambda]+I \[CapitalOmega]) (\[CapitalLambda]+3 I \[CapitalOmega])-2 \[CapitalLambda] \[CapitalOmega] (\[CapitalLambda]^2+\[CapitalOmega]^2) PolyGamma[1,-((I \[CapitalOmega])/\[Mu]0)])+\[Mu]0 (I \[CapitalLambda]+\[CapitalOmega]) (\[CapitalLambda]^2+\[CapitalOmega]^2) PolyGamma[1,(I \[CapitalOmega])/\[Mu]0]-2 \[CapitalLambda] \[CapitalOmega] (\[CapitalLambda]^2+\[CapitalOmega]^2) PolyGamma[0,(I \[CapitalOmega])/\[Mu]0] PolyGamma[1,(I \[CapitalOmega])/\[Mu]0]+\[Pi] \[CapitalLambda] Csch[(\[Pi] \[CapitalOmega])/\[Mu]0]^2 PolyGamma[0,\[CapitalLambda]/\[Mu]0] (-2 \[Pi] \[CapitalOmega] (\[CapitalLambda]^2+\[CapitalOmega]^2)+\[Mu]0 (\[CapitalLambda]^2-3 \[CapitalOmega]^2) Sinh[(2 \[Pi] \[CapitalOmega])/\[Mu]0]));

bz4Drude0 =1/(2 \[CapitalLambda] \[Mu]0^2 (\[CapitalLambda]^2+\[CapitalOmega]^2)^3) \[Beta] \[Kappa]^2 \[CapitalOmega] (2 \[Mu]0^2 \[CapitalOmega]^2+2 \[CapitalLambda] \[Mu]0 (\[CapitalLambda]^2-2 \[CapitalOmega]^2) PolyGamma[0,(\[CapitalLambda]+\[Mu]0)/\[Mu]0]+\[CapitalLambda] (\[Mu]0 (-\[CapitalLambda]^2+I \[CapitalLambda] \[CapitalOmega]+2 \[CapitalOmega]^2) PolyGamma[0,1-(I \[CapitalOmega])/\[Mu]0]-\[Mu]0 (\[CapitalLambda]^2+I \[CapitalLambda] \[CapitalOmega]-2 \[CapitalOmega]^2) PolyGamma[0,1+(I \[CapitalOmega])/\[Mu]0]+I \[CapitalOmega] (\[CapitalLambda]^2+\[CapitalOmega]^2) (PolyGamma[1,1-(I \[CapitalOmega])/\[Mu]0]-PolyGamma[1,1+(I \[CapitalOmega])/\[Mu]0])));


bz2Drude = bz2Drude0 //. SpecialSymbolVal;

azz2Drude = azz2Drude0 //. SpecialSymbolVal;

azz4Drude = azz4Drude0 //. SpecialSymbolVal;

bz4Drude = bz4Drude0 //. SpecialSymbolVal;


(*
IntegrandTerm1 = (\[Omega]n Coth[\[Beta] \[Omega]/2])/(\[Omega]^2 - \[Omega]n^2);
Integrandterm2 = \[Omega]/(\[Omega]^2 - \[Omega]n^2);
Integrand = JDrude[\[Omega]](IntegrandTerm1+Integrandterm2)
*)

AbetaExpr = (A \[CapitalOmega]^2)/(\[CapitalOmega]^2+\[Omega]n^2) ( (\[Pi] \[CapitalOmega])/2 + \[Omega]n (PolyGamma[0,\[Beta] \[CapitalOmega]/(2\[Pi])]-(1/2)(PolyGamma[0,I \[Beta] \[Omega]n/(2\[Pi])] + PolyGamma[0,-I \[Beta] \[Omega]n/(2\[Pi])])+\[Pi]/(\[Beta]*\[CapitalOmega]))) //.{\[CapitalOmega]->\[CapitalLambda],A->\[Gamma]};

(*Part 1 is term with coth. part 2 is the one without.*)

AbetaExprPart2 = (\[Pi] \[Gamma] \[CapitalLambda]^3)/(2 (\[CapitalLambda]^2+\[Omega]n^2))

AbetaExprPart1 = AbetaExpr - AbetaExprPart2

Abeta[\[Omega]n_]:= Evaluate[AbetaExpr];



AbetaPart1[\[Omega]n_] := Evaluate[AbetaExprPart1];

AbetaPart2[\[Omega]n_] := Evaluate[AbetaExprPart2];

(*
Hyperlink["Abeta Integral for Drude cutoff","CresserAndersABetaIntegralForDrudeCutoff.wl"]
*)


(* ::Section:: *)
(*Integrand*)


(*Do not depend on the following code. Instead use TCL4SpinBosonFunctions.
This code is causingconflict with that.*)


(*
Azz4IntTerm1 = Sin[\[CapitalOmega] (t1 - t2)] Sin[\[CapitalOmega] (t - t3)] \[Nu][t - t2] \[Nu][t1 - t3];
       
Azz4IntTerm2 = Sin[\[CapitalOmega] (t1 - t3)] Sin[\[CapitalOmega] (t - t2)] \[Nu][t - t3] \[Nu][t1 - t2];

Azz4TotalIntegrand = Azz4IntTerm1 + Azz4IntTerm2
*)


(* ::Chapter:: *)
(*General Hamiltonian*)


TwoPtCorrDrude = {\[Nu][t1_] \[Nu][t2_] -> NuDrude[n, t1] NuDrude[m,t2], \[Eta][t1_] \[Eta][t2_] -> EtaDrude[t1] EtaDrude[t2]}


Bz4TripleIntegralResult = 1/(\[CapitalLambda]^2-n^2 \[Mu]0^2) 2 I E^(-t (\[CapitalLambda]+n \[Mu]0+I \[CapitalOmega])) p1^2 \[Beta] \[Kappa]^2 (-((E^(t (-\[CapitalLambda]+n \[Mu]0+3 I \[CapitalOmega])) p1^2 \[CapitalLambda])/(\[CapitalLambda]-I \[CapitalOmega])^3)+(E^(-t (\[CapitalLambda]-n \[Mu]0+I \[CapitalOmega])) p1^2 \[CapitalLambda])/(\[CapitalLambda]+I \[CapitalOmega])^3+(2 E^(t (-\[CapitalLambda]+n \[Mu]0+I \[CapitalOmega])) p3^2)/(\[CapitalLambda] (\[CapitalLambda]+I \[CapitalOmega]))-(2 E^(t (-\[CapitalLambda]+n \[Mu]0+2 I \[CapitalOmega])) p3^2 (\[CapitalLambda]+I \[CapitalOmega]))/(\[CapitalLambda] (\[CapitalLambda]-I \[CapitalOmega]) (2 \[CapitalLambda]-I \[CapitalOmega]))+(2 E^(-t \[CapitalLambda]+n t \[Mu]0) p3^2 (2 \[CapitalLambda]-I \[CapitalOmega]))/(\[CapitalLambda] (\[CapitalLambda]+I \[CapitalOmega]) (2 \[CapitalLambda]+I \[CapitalOmega]))+(E^(t (-\[CapitalLambda]+n \[Mu]0+I \[CapitalOmega])) p1^2 (\[CapitalLambda]-2 I \[CapitalOmega]))/((\[CapitalLambda]-I \[CapitalOmega]) (\[CapitalLambda]+I \[CapitalOmega])^2)+(E^(-I t \[CapitalOmega]) n p1^2 \[Mu]0)/((\[CapitalLambda]+I \[CapitalOmega]) (n \[Mu]0+I \[CapitalOmega]) (\[CapitalLambda]+n \[Mu]0+2 I \[CapitalOmega]))+(E^(-I t \[CapitalOmega]) n p1^2 \[Mu]0 (2 \[CapitalLambda]+n \[Mu]0+3 I \[CapitalOmega]))/((n \[Mu]0+I \[CapitalOmega]) (\[CapitalLambda]+n \[Mu]0+2 I \[CapitalOmega]) (-I \[CapitalLambda]+\[CapitalOmega])^2)-(E^(3 I t \[CapitalOmega]) n p1^2 \[Mu]0)/((n \[Mu]0-I \[CapitalOmega]) (I \[CapitalLambda]+\[CapitalOmega])^2)+(2 E^(t (-\[CapitalLambda]+n \[Mu]0+2 I \[CapitalOmega])) p3^2)/((I \[CapitalLambda]+\[CapitalOmega]) (2 I \[CapitalLambda]+\[CapitalOmega]))-(2 E^(n t \[Mu]0+2 I t \[CapitalOmega]) p3^2)/((I \[CapitalLambda]+\[CapitalOmega]) (2 I \[CapitalLambda]+\[CapitalOmega]))+(2 E^(n t \[Mu]0+2 I t \[CapitalOmega]) n p1^2 \[Mu]0 \[CapitalOmega])/((\[CapitalLambda]-n \[Mu]0) (\[CapitalLambda]-I \[CapitalOmega]) (n \[Mu]0+I \[CapitalOmega]) (I n \[Mu]0+\[CapitalOmega]))+(2 E^(n t \[Mu]0) n p1^2 \[Mu]0 \[CapitalOmega])/((\[CapitalLambda]-n \[Mu]0) (\[CapitalLambda]+I \[CapitalOmega]) (n \[Mu]0+I \[CapitalOmega]) (I n \[Mu]0+\[CapitalOmega]))+(E^(t (\[CapitalLambda]+2 I \[CapitalOmega])) n \[Mu]0 (p1^2 \[CapitalLambda] (\[CapitalLambda]+n \[Mu]0-I \[CapitalOmega])+2 p3^2 (\[CapitalLambda]-I \[CapitalOmega]) (\[CapitalLambda]+n \[Mu]0-2 I \[CapitalOmega])))/(\[CapitalLambda] (\[CapitalLambda]-I \[CapitalOmega])^2 (I n \[Mu]0+\[CapitalOmega]) (I \[CapitalLambda]+I n \[Mu]0+\[CapitalOmega]))-(2 E^(t \[CapitalLambda]+I t \[CapitalOmega]) p3^2)/(\[CapitalLambda]^2-I \[CapitalLambda] \[CapitalOmega])+(2 E^(t \[CapitalLambda]+I t \[CapitalOmega]) p3^2)/(\[CapitalLambda]^2+I \[CapitalLambda] \[CapitalOmega])-(2 E^(n t \[Mu]0+2 I t \[CapitalOmega]) p3^2)/(\[CapitalLambda]^2+I \[CapitalLambda] \[CapitalOmega])-(2 I E^(t (\[CapitalLambda]+n \[Mu]0+I \[CapitalOmega])) \[CapitalOmega] (n^2 \[Mu]0^2 (p1^2 \[CapitalLambda]+p3^2 (\[CapitalLambda]+I \[CapitalOmega]))+n \[Mu]0 (p1^2 \[CapitalLambda]+p3^2 (\[CapitalLambda]+I \[CapitalOmega])) (\[CapitalLambda]+I \[CapitalOmega])+I p3^2 \[CapitalLambda] (\[CapitalLambda]+I \[CapitalOmega]) \[CapitalOmega]))/(\[CapitalLambda] (\[CapitalLambda]+n \[Mu]0) (\[CapitalLambda]+I \[CapitalOmega]) (n \[Mu]0+I \[CapitalOmega])^2 (\[CapitalLambda]+n \[Mu]0+I \[CapitalOmega]))-(2 p3^2 (n^3 \[Mu]0^3+n^2 \[Mu]0^2 (\[CapitalLambda]+2 I \[CapitalOmega])+\[CapitalLambda] \[CapitalOmega] (-I \[CapitalLambda]+\[CapitalOmega])))/(n \[CapitalLambda] \[Mu]0 (\[CapitalLambda]+I \[CapitalOmega]) (n \[Mu]0+I \[CapitalOmega]) (\[CapitalLambda]+n \[Mu]0+I \[CapitalOmega]))+(2 E^(2 I t \[CapitalOmega]) p3^2 (n^3 \[Mu]0^3+n^2 \[Mu]0^2 (\[CapitalLambda]-2 I \[CapitalOmega])+\[CapitalLambda] \[CapitalOmega] (I \[CapitalLambda]+\[CapitalOmega])))/(n \[CapitalLambda] \[Mu]0 (\[CapitalLambda]-I \[CapitalOmega]) (n \[Mu]0-I \[CapitalOmega]) (\[CapitalLambda]+n \[Mu]0-I \[CapitalOmega]))-(E^(t (-\[CapitalLambda]+n \[Mu]0+I \[CapitalOmega])) (p1^2 \[CapitalLambda] (\[CapitalLambda]+2 I \[CapitalOmega])+2 p3^2 (\[CapitalLambda]^2+\[CapitalOmega]^2)))/(\[CapitalLambda] (\[CapitalLambda]-I \[CapitalOmega])^2 (\[CapitalLambda]+I \[CapitalOmega]))-(2 E^(t (n \[Mu]0+I \[CapitalOmega])) p3^2 (3 n^3 \[Mu]0^3-4 n^2 \[Mu]0^2 (\[CapitalLambda]-I \[CapitalOmega])+I \[CapitalOmega] (\[CapitalLambda]^2+\[CapitalOmega]^2)+n \[Mu]0 (\[CapitalLambda]^2-2 I \[CapitalLambda] \[CapitalOmega]+2 \[CapitalOmega]^2)))/(\[CapitalLambda] (n \[Mu]0+I \[CapitalOmega]) (-\[CapitalLambda]+n \[Mu]0+I \[CapitalOmega]) (I \[CapitalLambda]+\[CapitalOmega]) (I n \[Mu]0+\[CapitalOmega]))+(2 E^(t (n \[Mu]0+I \[CapitalOmega])) p3^2 (3 n^3 \[Mu]0^3-4 n^2 \[Mu]0^2 (\[CapitalLambda]+I \[CapitalOmega])-I \[CapitalOmega] (\[CapitalLambda]^2+\[CapitalOmega]^2)+n \[Mu]0 (\[CapitalLambda]^2+2 I \[CapitalLambda] \[CapitalOmega]+2 \[CapitalOmega]^2)))/(\[CapitalLambda] (\[CapitalLambda]+I \[CapitalOmega]) (\[CapitalLambda]-n \[Mu]0+I \[CapitalOmega]) (-I n \[Mu]0+\[CapitalOmega]) (I n \[Mu]0+\[CapitalOmega]))-(2 E^(t (\[CapitalLambda]+n \[Mu]0+I \[CapitalOmega])) \[CapitalOmega] (p3^2 \[CapitalLambda] (\[CapitalLambda]-I \[CapitalOmega]) \[CapitalOmega]+n \[Mu]0 (p1^2 \[CapitalLambda]+p3^2 (\[CapitalLambda]-I \[CapitalOmega])) (I \[CapitalLambda]+\[CapitalOmega])+n^2 \[Mu]0^2 (I p1^2 \[CapitalLambda]+p3^2 (I \[CapitalLambda]+\[CapitalOmega]))))/(\[CapitalLambda] (\[CapitalLambda]+n \[Mu]0) (\[CapitalLambda]-I \[CapitalOmega]) (n \[Mu]0-I \[CapitalOmega])^2 (\[CapitalLambda]+n \[Mu]0-I \[CapitalOmega]))-(2 I E^(n t \[Mu]0+2 I t \[CapitalOmega]) (3 n^4 p3^2 \[Mu]0^4-n^2 p3^2 \[Mu]0^2 (\[CapitalLambda]-I \[CapitalOmega])^2+2 n^3 p3^2 \[Mu]0^3 (\[CapitalLambda]+2 I \[CapitalOmega])+p3^2 \[CapitalLambda] (\[CapitalLambda]+2 I \[CapitalOmega]) \[CapitalOmega]^2+n \[Mu]0 \[CapitalOmega] (p1^2 \[CapitalLambda] (I \[CapitalLambda]+\[CapitalOmega])+2 p3^2 (-I \[CapitalLambda]^2+\[CapitalLambda] \[CapitalOmega]+I \[CapitalOmega]^2))))/(\[CapitalLambda] (\[CapitalLambda]+n \[Mu]0) (\[CapitalLambda]-I \[CapitalOmega]) (n \[Mu]0+I \[CapitalOmega])^2 (I n \[Mu]0+\[CapitalOmega]))+(2 E^(t (\[CapitalLambda]+n \[Mu]0+I \[CapitalOmega])) \[CapitalOmega] (n \[Mu]0 (p1^2 \[CapitalLambda]^3+p3^2 (\[CapitalLambda]-I \[CapitalOmega])^2 (\[CapitalLambda]+I \[CapitalOmega])) (2 I \[CapitalLambda]+\[CapitalOmega])+n^2 \[Mu]0^2 (2 I p3^2 (\[CapitalLambda]-I \[CapitalOmega])^2 (\[CapitalLambda]+I \[CapitalOmega])+p1^2 \[CapitalLambda]^2 (2 I \[CapitalLambda]+\[CapitalOmega]))-\[CapitalLambda] (\[CapitalLambda]-I \[CapitalOmega]) \[CapitalOmega] (p1^2 \[CapitalLambda] (2 \[CapitalLambda]-I \[CapitalOmega])-p3^2 (\[CapitalLambda]^2+\[CapitalOmega]^2))))/(\[CapitalLambda] (\[CapitalLambda]+n \[Mu]0) (\[CapitalLambda]-I \[CapitalOmega])^3 (2 \[CapitalLambda]-I \[CapitalOmega]) (\[CapitalLambda]+n \[Mu]0-I \[CapitalOmega]) (\[CapitalLambda]+I \[CapitalOmega]))+(2 I E^(t (\[CapitalLambda]+n \[Mu]0+I \[CapitalOmega])) \[CapitalOmega] (n^2 \[Mu]0^2 (2 p3^2 (\[CapitalLambda]-I \[CapitalOmega]) (\[CapitalLambda]+I \[CapitalOmega])^2+p1^2 \[CapitalLambda]^2 (2 \[CapitalLambda]+I \[CapitalOmega]))+n \[Mu]0 (p1^2 \[CapitalLambda]^3+p3^2 (\[CapitalLambda]-I \[CapitalOmega]) (\[CapitalLambda]+I \[CapitalOmega])^2) (2 \[CapitalLambda]+I \[CapitalOmega])+\[CapitalLambda] \[CapitalOmega] (-I \[CapitalLambda]+\[CapitalOmega]) (p1^2 \[CapitalLambda] (2 \[CapitalLambda]+I \[CapitalOmega])-p3^2 (\[CapitalLambda]^2+\[CapitalOmega]^2))))/(\[CapitalLambda] (\[CapitalLambda]+n \[Mu]0) (\[CapitalLambda]-I \[CapitalOmega]) (\[CapitalLambda]+I \[CapitalOmega])^3 (2 \[CapitalLambda]+I \[CapitalOmega]) (\[CapitalLambda]+n \[Mu]0+I \[CapitalOmega]))+(2 I E^(I t \[CapitalOmega]) \[CapitalOmega] (n^3 p1^2 \[CapitalLambda] \[Mu]0^3 (-\[CapitalLambda]^2+\[CapitalOmega]^2)+2 p3^2 \[CapitalOmega]^4 (\[CapitalLambda]^2+\[CapitalOmega]^2)+n p1^2 \[CapitalLambda] \[Mu]0 \[CapitalOmega]^2 (\[CapitalLambda]^2+3 \[CapitalOmega]^2)+2 n^4 \[Mu]0^4 (p1^2 \[CapitalLambda]^2+p3^2 (\[CapitalLambda]^2+\[CapitalOmega]^2))+2 n^2 \[Mu]0^2 (2 p3^2 \[CapitalOmega]^2 (\[CapitalLambda]^2+\[CapitalOmega]^2)+p1^2 (\[CapitalLambda]^4+2 \[CapitalLambda]^2 \[CapitalOmega]^2))))/(\[CapitalLambda] (\[CapitalLambda]^2+\[CapitalOmega]^2)^2 (n^2 \[Mu]0^2+\[CapitalOmega]^2)^2)+(E^(t (\[CapitalLambda]+2 I \[CapitalOmega])) (-n^4 \[CapitalLambda] \[Mu]0^4 (p1^2 (\[CapitalLambda]-I \[CapitalOmega])+4 p3^2 (\[CapitalLambda]+I \[CapitalOmega]))+4 I n^3 p3^2 \[CapitalLambda] \[Mu]0^3 (\[CapitalLambda]+I \[CapitalOmega]) \[CapitalOmega]+2 p3^2 \[CapitalLambda] (\[CapitalLambda]+I \[CapitalOmega])^3 \[CapitalOmega]^2+n \[Mu]0 \[CapitalOmega] (-I \[CapitalLambda]+\[CapitalOmega]) (p1^2 \[CapitalLambda] (\[CapitalLambda]^2+I \[CapitalLambda] \[CapitalOmega]+2 \[CapitalOmega]^2)-2 p3^2 (\[CapitalLambda]^3+I \[CapitalLambda]^2 \[CapitalOmega]+3 \[CapitalLambda] \[CapitalOmega]^2+I \[CapitalOmega]^3))+n^2 \[Mu]0^2 (p1^2 (\[CapitalLambda]^4+2 I \[CapitalLambda]^3 \[CapitalOmega]+3 I \[CapitalLambda] \[CapitalOmega]^3)+2 p3^2 (2 \[CapitalLambda]^4+5 I \[CapitalLambda]^3 \[CapitalOmega]-5 \[CapitalLambda]^2 \[CapitalOmega]^2-I \[CapitalLambda] \[CapitalOmega]^3-\[CapitalOmega]^4))))/(\[CapitalLambda] (\[CapitalLambda]+n \[Mu]0) (\[CapitalLambda]-I \[CapitalOmega]) (\[CapitalLambda]+I \[CapitalOmega])^2 (\[CapitalLambda]-n \[Mu]0+I \[CapitalOmega]) (-I n \[Mu]0+\[CapitalOmega]) (I n \[Mu]0+\[CapitalOmega]))+(2 I E^(t \[CapitalLambda]) n \[Mu]0 (n^2 \[Mu]0^2 (\[CapitalLambda]-I \[CapitalOmega]) (-I p3^2 (\[CapitalLambda]+I \[CapitalOmega])^2 (2 \[CapitalLambda]-3 I \[CapitalOmega])-2 p1^2 \[CapitalLambda] \[CapitalOmega] (-2 \[CapitalLambda]+t \[CapitalLambda]^2+t \[CapitalOmega]^2))+n^3 \[Mu]0^3 (p3^2 (\[CapitalLambda]+I \[CapitalOmega])^2 (I \[CapitalLambda]+\[CapitalOmega])+p1^2 \[CapitalLambda] \[CapitalOmega] (-2 \[CapitalLambda]+t \[CapitalLambda]^2+t \[CapitalOmega]^2))+\[CapitalLambda] \[CapitalOmega] (I \[CapitalLambda]+\[CapitalOmega]) (I p3^2 (\[CapitalLambda]^3+3 \[CapitalLambda] \[CapitalOmega]^2+2 I \[CapitalOmega]^3)+p1^2 (-\[CapitalLambda]^2 \[CapitalOmega]+\[CapitalOmega]^3+\[CapitalLambda]^3 (-I+t \[CapitalOmega])+\[CapitalLambda] \[CapitalOmega]^2 (-I+t \[CapitalOmega])))+n \[Mu]0 (p3^2 (\[CapitalLambda]+I \[CapitalOmega])^2 (I \[CapitalLambda]^3+5 \[CapitalLambda]^2 \[CapitalOmega]-6 I \[CapitalLambda] \[CapitalOmega]^2-2 \[CapitalOmega]^3)+p1^2 \[CapitalLambda] \[CapitalOmega] (t \[CapitalLambda]^4+5 I \[CapitalLambda]^2 \[CapitalOmega]+\[CapitalLambda]^3 (-3-3 I t \[CapitalOmega])+\[CapitalLambda] \[CapitalOmega]^2 (1-3 I t \[CapitalOmega])-\[CapitalOmega]^3 (I+t \[CapitalOmega])))))/(\[CapitalLambda] (\[CapitalLambda]-n \[Mu]0) (\[CapitalLambda]-n \[Mu]0-I \[CapitalOmega]) (\[CapitalLambda]+I \[CapitalOmega])^2 (I \[CapitalLambda]+\[CapitalOmega])^2 (-I n \[Mu]0+\[CapitalOmega])^2)-(2 E^(t (\[CapitalLambda]+2 I \[CapitalOmega])) \[CapitalOmega] (-n^6 p1^2 t \[CapitalLambda] \[Mu]0^6-n^5 p1^2 t \[CapitalLambda] \[Mu]0^5 (\[CapitalLambda]-I \[CapitalOmega])+n p1^2 t \[CapitalLambda]^3 \[Mu]0 (\[CapitalLambda]-I \[CapitalOmega]) \[CapitalOmega]^2+p3^2 \[CapitalLambda]^2 \[CapitalOmega]^2 (\[CapitalLambda]^2+\[CapitalOmega]^2)+n^3 p1^2 \[CapitalLambda] \[Mu]0^3 (\[CapitalLambda]-I \[CapitalOmega]) (2 \[CapitalLambda]+t \[CapitalLambda]^2+\[CapitalOmega] (2 I-t \[CapitalOmega]))-n^4 \[Mu]0^4 (p1^2 t \[CapitalLambda] (-\[CapitalLambda]^2+\[CapitalOmega]^2)+p3^2 (\[CapitalLambda]^2+\[CapitalOmega]^2))+n^2 \[Mu]0^2 (p3^2 (\[CapitalLambda]^4-\[CapitalOmega]^4)+p1^2 \[CapitalLambda] (2 \[CapitalLambda]^3+2 \[CapitalLambda] \[CapitalOmega]^2-2 I \[CapitalOmega]^3+\[CapitalLambda]^2 \[CapitalOmega] (-2 I+t \[CapitalOmega])))))/(\[CapitalLambda] (\[CapitalLambda]-n \[Mu]0) (\[CapitalLambda]+n \[Mu]0) (\[CapitalLambda]-I \[CapitalOmega]) (\[CapitalLambda]+n \[Mu]0-I \[CapitalOmega]) (\[CapitalLambda]+I \[CapitalOmega]) (-I n \[Mu]0+\[CapitalOmega]) (I n \[Mu]0+\[CapitalOmega])^2)+(E^(n t \[Mu]0+2 I t \[CapitalOmega]) (2 p3^2 \[CapitalLambda]^2 (\[CapitalLambda]+I \[CapitalOmega])^2 \[CapitalOmega]^2 (2 I \[CapitalLambda]^2+3 \[CapitalLambda] \[CapitalOmega]-I \[CapitalOmega]^2)+n^3 \[Mu]0^3 (2 p3^2 (\[CapitalLambda]+I \[CapitalOmega])^2 (I \[CapitalLambda]^3+2 I \[CapitalLambda] \[CapitalOmega]^2+\[CapitalOmega]^3)+p1^2 \[CapitalLambda] (2 \[CapitalLambda]-I \[CapitalOmega]) \[CapitalOmega] (4 t \[CapitalLambda]^3-I t \[CapitalOmega]^3+\[CapitalLambda]^2 (-6+3 I t \[CapitalOmega])))+n^2 \[CapitalLambda] \[Mu]0^2 (-2 I p3^2 (\[CapitalLambda]-I \[CapitalOmega]) (\[CapitalLambda]+I \[CapitalOmega])^2 \[CapitalOmega]^2-p1^2 (2 \[CapitalLambda]-I \[CapitalOmega]) (-2 I \[CapitalOmega]^4+4 \[CapitalLambda] \[CapitalOmega]^3 (1-I t \[CapitalOmega])+\[CapitalLambda]^4 (-I+t \[CapitalOmega])-\[CapitalLambda]^2 \[CapitalOmega]^2 (I+3 t \[CapitalOmega])))+n^4 \[Mu]0^4 (2 I p3^2 (\[CapitalLambda]^2+\[CapitalOmega]^2)^2+p1^2 \[CapitalLambda] (2 \[CapitalLambda]-I \[CapitalOmega]) (4 I \[CapitalLambda] \[CapitalOmega] (I+t \[CapitalOmega])+\[CapitalOmega]^2 (I+t \[CapitalOmega])+\[CapitalLambda]^2 (-I+5 t \[CapitalOmega])))+n \[CapitalLambda] \[Mu]0 \[CapitalOmega] (2 p3^2 (\[CapitalLambda]+I \[CapitalOmega])^2 (2 \[CapitalLambda]^3-2 I \[CapitalLambda]^2 \[CapitalOmega]+\[CapitalLambda] \[CapitalOmega]^2-I \[CapitalOmega]^3)+p1^2 \[CapitalLambda]^2 \[CapitalOmega] (2 I \[CapitalLambda]+\[CapitalOmega]) (6 I \[CapitalOmega]+t (\[CapitalLambda]^2-4 I \[CapitalLambda] \[CapitalOmega]+5 \[CapitalOmega]^2)))))/(n \[CapitalLambda] \[Mu]0 (\[CapitalLambda]+n \[Mu]0) (n \[Mu]0-I \[CapitalOmega]) (n \[Mu]0+I \[CapitalOmega]) (2 I \[CapitalLambda]+\[CapitalOmega]) (\[CapitalLambda]^2+\[CapitalOmega]^2)^2)+1/(n \[CapitalLambda] \[Mu]0 (\[CapitalLambda]-I \[CapitalOmega])^2 (-I \[CapitalLambda]+\[CapitalOmega])^2 (I n \[Mu]0+\[CapitalOmega])^2) E^(n t \[Mu]0) (2 p3^2 \[CapitalLambda] (\[CapitalLambda]-I \[CapitalOmega])^2 (\[CapitalLambda]+I \[CapitalOmega]) \[CapitalOmega]^2+n^2 \[Mu]0^2 (-2 p3^2 (\[CapitalLambda]^4+5 I \[CapitalLambda]^3 \[CapitalOmega]+3 \[CapitalLambda]^2 \[CapitalOmega]^2+5 I \[CapitalLambda] \[CapitalOmega]^3+2 \[CapitalOmega]^4)+p1^2 \[CapitalLambda] (-9 t \[CapitalLambda]^2 \[CapitalOmega]^2+\[CapitalOmega]^3 (2 I-t \[CapitalOmega])+\[CapitalLambda]^3 (-1+I t \[CapitalOmega])+9 \[CapitalLambda] \[CapitalOmega]^2 (1+I t \[CapitalOmega])))+n^3 \[Mu]0^3 (2 p3^2 (3 \[CapitalLambda]^3-I \[CapitalLambda]^2 \[CapitalOmega]+3 \[CapitalLambda] \[CapitalOmega]^2-I \[CapitalOmega]^3)+p1^2 \[CapitalLambda] (\[CapitalOmega]^2 (-1-I t \[CapitalOmega])+\[CapitalLambda]^2 (1-5 I t \[CapitalOmega])-4 \[CapitalLambda] \[CapitalOmega] (-I+t \[CapitalOmega])))+n \[Mu]0 \[CapitalOmega] (2 p3^2 (\[CapitalLambda]+I \[CapitalOmega])^2 (2 I \[CapitalLambda]^2+3 \[CapitalLambda] \[CapitalOmega]-I \[CapitalOmega]^2)+p1^2 \[CapitalLambda]^2 \[CapitalOmega] (-6 I \[CapitalOmega]+t (\[CapitalLambda]^2+4 I \[CapitalLambda] \[CapitalOmega]+5 \[CapitalOmega]^2)))));


Bz4GenHDrude = 1/(\[Beta] (\[CapitalLambda]^2+\[CapitalOmega]^2)^3 (4 \[CapitalLambda]^2+\[CapitalOmega]^2)) 2 p1^2 \[Gamma]^2 \[CapitalLambda] \[CapitalOmega] (4 \[Pi]^2 \[CapitalOmega]^2 (2 p1^2 \[CapitalLambda]^2 (4 \[CapitalLambda]^2+\[CapitalOmega]^2)+p3^2 (-5 \[CapitalLambda]^4-4 \[CapitalLambda]^2 \[CapitalOmega]^2+\[CapitalOmega]^4))+4 \[Pi] \[Beta] \[CapitalLambda]^3 (6 p3^2 (\[CapitalLambda]^2+\[CapitalOmega]^2)^2+p1^2 (\[CapitalLambda]^2-2 \[CapitalOmega]^2) (4 \[CapitalLambda]^2+\[CapitalOmega]^2)) PolyGamma[0,1+(\[Beta] \[CapitalLambda])/(2 \[Pi])]-2 p3^2 \[Pi] \[Beta] \[CapitalLambda] (\[CapitalLambda]-I \[CapitalOmega])^2 (\[CapitalLambda]+I \[CapitalOmega])^3 (2 \[CapitalLambda]+I \[CapitalOmega]) PolyGamma[0,1+(\[Beta] (\[CapitalLambda]-I \[CapitalOmega]))/(2 \[Pi])]+\[Beta] \[CapitalLambda] (2 \[CapitalLambda]-I \[CapitalOmega]) (-2 p3^2 \[Pi] (\[CapitalLambda]-I \[CapitalOmega])^3 (\[CapitalLambda]+I \[CapitalOmega])^2 PolyGamma[0,1+(\[Beta] (\[CapitalLambda]+I \[CapitalOmega]))/(2 \[Pi])]-(2 \[CapitalLambda]+I \[CapitalOmega]) (2 \[Pi] (p1^2 \[CapitalLambda]^2 (\[CapitalLambda]^2-I \[CapitalLambda] \[CapitalOmega]-2 \[CapitalOmega]^2)+p3^2 (\[CapitalLambda]^2+\[CapitalOmega]^2)^2) PolyGamma[0,1-(I \[Beta] \[CapitalOmega])/(2 \[Pi])]+2 \[Pi] (p1^2 \[CapitalLambda]^2 (\[CapitalLambda]^2+I \[CapitalLambda] \[CapitalOmega]-2 \[CapitalOmega]^2)+p3^2 (\[CapitalLambda]^2+\[CapitalOmega]^2)^2) PolyGamma[0,1+(I \[Beta] \[CapitalOmega])/(2 \[Pi])]-I p1^2 \[Beta] \[CapitalLambda]^2 \[CapitalOmega] (\[CapitalLambda]^2+\[CapitalOmega]^2) (PolyGamma[1,1-(I \[Beta] \[CapitalOmega])/(2 \[Pi])]-PolyGamma[1,1+(I \[Beta] \[CapitalOmega])/(2 \[Pi])]))));


Tcl2drude= {{0,0,0,0},{-((\[Pi] \[Gamma] \[CapitalLambda]^2 \[CapitalOmega] Sin[2 \[Theta]])/(\[CapitalLambda]^2+\[CapitalOmega]^2)),-((4 \[Pi] \[Gamma] Cos[\[Theta]]^2)/\[Beta]),0,-((4 \[Pi]^2 \[Gamma] \[CapitalLambda]^2 \[CapitalOmega] Cos[\[Theta]] Coth[(\[Beta] \[CapitalOmega])/2] Sin[\[Theta]])/(\[Beta] ((2 \[Pi] \[CapitalLambda]^2)/\[Beta]+(2 \[Pi] \[CapitalOmega]^2)/\[Beta])))},{-((\[Pi] \[Gamma] \[CapitalLambda] \[CapitalOmega]^2 Sin[2 \[Theta]])/(\[CapitalLambda]^2+\[CapitalOmega]^2)),(2 \[Gamma] \[CapitalLambda] \[CapitalOmega] ((2 \[Pi])/\[Beta]-2 \[CapitalLambda] PolyGamma[0,(\[Beta] ((2 \[Pi])/\[Beta]+\[CapitalLambda]))/(2 \[Pi])]+\[CapitalLambda] PolyGamma[0,1-(I \[Beta] \[CapitalOmega])/(2 \[Pi])]+\[CapitalLambda] PolyGamma[0,1+(I \[Beta] \[CapitalOmega])/(2 \[Pi])]) Sin[\[Theta]]^2)/(\[CapitalLambda]^2+\[CapitalOmega]^2),-((4 \[Pi] \[Gamma] Cos[\[Theta]]^2)/\[Beta])-(2 \[Pi] \[Gamma] \[CapitalLambda]^2 \[CapitalOmega] Coth[(\[Beta] \[CapitalOmega])/2] Sin[\[Theta]]^2)/(\[CapitalLambda]^2+\[CapitalOmega]^2),-((\[Gamma] \[CapitalLambda] \[CapitalOmega] ((2 \[Pi])/\[Beta]-2 \[CapitalLambda] PolyGamma[0,(\[Beta] ((2 \[Pi])/\[Beta]+\[CapitalLambda]))/(2 \[Pi])]+\[CapitalLambda] PolyGamma[0,1-(I \[Beta] \[CapitalOmega])/(2 \[Pi])]+\[CapitalLambda] PolyGamma[0,1+(I \[Beta] \[CapitalOmega])/(2 \[Pi])]) Sin[2 \[Theta]])/(\[CapitalLambda]^2+\[CapitalOmega]^2))},{-((2 \[Pi] \[Gamma] \[CapitalLambda]^2 \[CapitalOmega] Sin[\[Theta]]^2)/(\[CapitalLambda]^2+\[CapitalOmega]^2)),-((4 \[Pi] \[Gamma] Cos[\[Theta]] Sin[\[Theta]])/\[Beta]),0,-((2 \[Pi] \[Gamma] \[CapitalLambda]^2 \[CapitalOmega] Coth[(\[Beta] \[CapitalOmega])/2] Sin[\[Theta]]^2)/(\[CapitalLambda]^2+\[CapitalOmega]^2))}};
