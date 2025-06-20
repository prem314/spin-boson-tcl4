(* ::Package:: *)

(* ::Title:: *)
(*AntiPodalSBMProjectNMPlot*)


(* ::Chapter:: *)
(*Data load*)


Get["TCL4SpinBosonFunctions.wl"]


(* ::Input:: *)
(*(*\[Gamma]=0.01;\[CapitalOmega]=1.6;\[Theta]=\[Pi]/2; StateDt=0.125*)*)


(* ::Input:: *)
(**)


lambdaList=Range[0.2,8,0.1]; 
temperatureList=Range[0.2,8,0.1]; 

Loadtcl2=Get[FileNameJoin[{NotebookDirectory[],"datasimpletcl2_gamma01_thetapi2_modifiedtcl2_400_antipodal.mx"}]];
Loadtcl4=Get[FileNameJoin[{NotebookDirectory[],"datasimpletcl4_gamma01_thetapi2_modifiedtcl4_400_antipodal.mx"}]];

tcl2Name = "AntiPtcl2lamdaTplot.pdf"
tcl4Name = "AntiPtcl4lamdaTplot.pdf"
diffPlotName = "AntiPredifftcl.pdf"



Length[lambdaList]


(*
OldlambdaList=Range[0.2,8,0.2]; 
OldtemperatureList=Range[0.2,8,0.2]; 

OldLoadtcl2=Get[FileNameJoin[{NotebookDirectory[],"datasimpletcl2_gamma01_thetapi2_modifiedtcl2_4096.mx"}]];
OldLoadtcl4=Get[FileNameJoin[{NotebookDirectory[],"datasimpletcl4_gamma01_thetapi2_modifiedtcl4_4096.mx"}]];

Loadtcl2 = OldLoadtcl2[[3 ;;, 9 ;;]];
Loadtcl4 = OldLoadtcl4[[3 ;;, 9 ;;]];

lambdaList=Range[2,8,0.2]
temperatureList=Range[0.6,8,0.2]

tcl2Name = "tcl2lamdaTplot.pdf"
tcl4Name = "tcl4lamdaTplot.pdf"
diffPlotName = "redifftcl.pdf"
*)


(*
lambdaList=Range[0.5,15,0.5]; 
temperatureList=Range[0.5,15,0.5]; 

Loadtcl2=Get[FileNameJoin[{NotebookDirectory[],"datasimpletcl2_gamma01_thetapi2_modifiedtcl2_100_antipodal_l15.mx"}]];
Loadtcl4=Get[FileNameJoin[{NotebookDirectory[],"datasimpletcl4_gamma01_thetapi2_modifiedtcl4_100_antipodal_l15.mx"}]];
*)


(*
lambdaList=Range[3,50,0.5]; 
temperatureList=Range[0.5,50,0.5]; 

OldLoadtcl2=Get[FileNameJoin[{NotebookDirectory[],"datasimpletcl2_gamma01_thetapi2_modifiedtcl2_100_antipodal_lambdaT50.mx"}]];

OldLoadtcl4=Get[FileNameJoin[{NotebookDirectory[],"datasimpletcl4_gamma01_thetapi2_modifiedtcl4_100_antipodal_lambdaT50.mx"}]];

Loadtcl2 = OldLoadtcl2[[1 ;;, 6 ;;]];
Loadtcl4 = OldLoadtcl4[[1 ;;, 6 ;;]];
*)


(* ::Input:: *)
(**)


Wiki[TCLMarkovianityResonanceCondition]


LambdaCurve[T_] := (\[CapitalOmega] Sqrt[2 Coth[\[CapitalOmega]/(2 T)]+(\[CapitalOmega] Csch[\[CapitalOmega]/(2 T)]^2)/T])/Sqrt[2 Coth[\[CapitalOmega]/(2 T)]-(\[CapitalOmega] Csch[\[CapitalOmega]/(2 T)]^2)/T]


(* ::Chapter:: *)
(*Plots*)


\[CapitalOmega] = 1.6


(* ::Section:: *)
(*TCL2*)


Font1 = 20
Font2 = 16


TCl2Normalized = Log10[Loadtcl2 + 10^-6];

(* Get the min and max of your data just once *)
dataRange = MinMax[TCl2Normalized];
{dmin, dmax} = dataRange;

(* This new function manually scales any input value from your data's range 
  to the [0, 1] range that GrayLevel expects.
*)
(*
universalColorFunction = (GrayLevel[1 - ((# - dmin)/(dmax - dmin))]) &;
*)

universalColorFunction = (Blend[{White, Darker[Darker[Purple]]}, (# - dmin)/(dmax - dmin)]) &;

(* 1. THE PLOT *)
figtcl2 = ListDensityPlot[TCl2Normalized,
   FrameLabel -> {Style["\[CapitalLambda]", Font1], Style["T", Font1]},
   FrameStyle -> Directive[Black, AbsoluteThickness[1]],
   FrameTicksStyle -> Directive[Black, Font2],
   ColorFunction -> universalColorFunction, (* Use the new function *)
   ColorFunctionScaling -> False, (* IMPORTANT: Turn off automatic scaling *)
   DataRange -> {MinMax[lambdaList], MinMax[temperatureList]},
   PlotRange -> All,
   ImageSize -> 280
   ];

(* 2. THE OVERLAID LINE *)
(* A black line is a good choice for a mostly white background *)
linePlot = ContourPlot[
   \[Lambda] == LambdaCurve[T],
   {\[Lambda], MinMax[lambdaList][[1]], MinMax[lambdaList][[2]]},
   {T, MinMax[temperatureList][[1]], MinMax[temperatureList][[2]]},
   ContourStyle -> {Black, Thick}
   ];


(* 3. COMBINE AND ADD THE LEGEND *)

(* Create a list of integer tick positions from your data's range *)
tickPositions = Range[Ceiling[dmin], Floor[dmax]];

(* Create the corresponding labels in the format 10^x using Superscript *)
tickLabels = Superscript[10, #] & /@ tickPositions;

(* Combine the positions and labels into a list for the Ticks option *)
customTicks = Transpose[{tickPositions, tickLabels}];


finalPlotTCL2 = Legended[
   Show[figtcl2, linePlot],
   (* The legend now uses the same function and range, ensuring a perfect match *)
   BarLegend[{universalColorFunction, dataRange},
    LegendLayout -> "Column", 
    LegendMarkerSize -> 220,
    Ticks -> customTicks,
    LabelStyle -> {Black, Font2}
    ]
   ];

(* Display the final, correct plot *)
finalPlotTCL2


Export[FileNameJoin[{TCL4DynamicsFolder,tcl2Name}],finalPlotTCL2,ImageResolution->300];


(* ::Section:: *)
(*TCL4*)


TCL4Normalized = Log10[Loadtcl4 + 10^-6];

(* 1. THE PLOT *)
figtcl4 = ListDensityPlot[TCL4Normalized,
   FrameLabel -> {Style["\[CapitalLambda]", Font1], Style["T", Font1]},
   FrameStyle -> Directive[Black, AbsoluteThickness[1]],
   FrameTicksStyle -> Directive[Black, Font2],
   ColorFunction -> universalColorFunction, (* Use the new function *)
   ColorFunctionScaling -> False, (* IMPORTANT: Turn off automatic scaling *)
   DataRange -> {MinMax[lambdaList], MinMax[temperatureList]},
   PlotRange -> All,
   ImageSize -> 280
   ];
(* 3. COMBINE AND ADD THE LEGEND *)
finalPlotTCL4 = Legended[
   Show[figtcl4, linePlot],
   (* The legend now uses the same function and range, ensuring a perfect match *)
   BarLegend[{universalColorFunction, dataRange},
    LegendLayout -> "Column", 
    LegendMarkerSize -> 220,
    Ticks -> customTicks,
    LabelStyle -> {Black, Font2}
    ]
   ]


Export[FileNameJoin[{TCL4DynamicsFolder,tcl4Name}],finalPlotTCL4,ImageResolution->300];


(* ::Section:: *)
(*TCL2 and TCL4 difference*)


(* ::Input:: *)
(*difftcl2tcl4=Loadtcl4-Loadtcl2;*)
(**)
(*DiffMax = Max[difftcl2tcl4]*)
(*DiffMin = Abs[Min[difftcl2tcl4]]*)
(**)
(**)
(*rediffplot=Module[{values=Flatten[difftcl2tcl4],posMax,negMax,linThresh,symlog,symlogPosMax,symlogNegMax,normalizeValue,(*The new function to map data to[-1,1]*)plotColorFunction,legendColorFunction,legendTicks},(*---1. Define Asymmetric Range and Logarithmic Scaling---*)posMax=DiffMax;*)
(*negMax=DiffMin;*)
(*(*Base the linear threshold on the smaller positive scale to capture its detail*)linThresh=0.1*posMax;*)
(*If[linThresh==0,linThresh=1];*)
(*symlog=Function[x,Sign[x]*Log[1+Abs[x]/linThresh]];*)
(*symlogPosMax=symlog[posMax];*)
(*symlogNegMax=Abs[symlog[-negMax]];*)
(*If[symlogPosMax==0,symlogPosMax=1];*)
(*If[symlogNegMax==0,symlogNegMax=1];*)
(*(*---2. Create the Normalization Function---*)(*This is the key step.It maps any data value to a point in the[-1,1] range.*)normalizeValue=Function[val,Which[val>0,symlog[val]/symlogPosMax,(*Maps (0,posMax] to (0,1]*)val<0,-(Abs[symlog[val]]/symlogNegMax),(*Maps[-negMax,0) to[-1,0)*)True,0]];*)
(*(*---3. Define Color Functions---*)(*The plot's color function now uses the normalized value*)plotColorFunction=Function[val,Module[{norm=normalizeValue[val]},Which[norm>0,Blend[{White,Red},norm],(*norm is already 0..1*)norm<0,Blend[{White,Blue},-norm],(*-norm makes it 0..1*)True,White]]];*)
(*(*A simpler version for the legend,which works directly on the[-1,1] scale*)legendColorFunction=Function[norm,Which[norm>0,Blend[{White,Red},norm],norm<0,Blend[{White,Blue},-norm],True,White]];*)
(**)
(*(*---4. Create Ticks for the Symmetric Legend---*)(*The real-world values you want to label on the bar*)*)
(*posTicksData={10^-4,posMax};*)
(*negTicksData={-negMax,-10^-1, -10^-2,-10^-3};*)
(*allTicksData=Sort@DeleteDuplicates@Join[negTicksData,{0},posTicksData];*)
(**)
(*(*Create {position,label} pairs.The position is the value's location on the normalized[-1,1] bar.*)legendTicks=Table[{normalizeValue[t],NumberForm[t,2,ExponentFunction->(If[-3<#<4,Null,#]&)]},{t,allTicksData}];*)
(**)
(*(*---5. Generate the Plot---*)*)
(*ListDensityPlot[difftcl2tcl4,ColorFunction->plotColorFunction,*)
(*ColorFunctionScaling->False,*)
(*DataRange->{MinMax[lambdaList],*)
(*MinMax[temperatureList]},*)
(*PlotLegends->BarLegend[*)
(*{legendColorFunction,{-1,1}},*)
(*Ticks->legendTicks,*)
(*LegendMarkerSize->250,*)
(*LabelStyle->{Black, Font2}],*)
(*FrameLabel->{Style["\[CapitalLambda]",Font1],*)
(*Style["T",Font1]},*)
(*PlotRange->Full,*)
(*ImageSize->280,*)
(*FrameStyle->Directive[Black,AbsoluteThickness[1]],*)
(*FrameTicksStyle->Directive[Black,Font2]*)
(*]];*)
(**)
(*finalPlotDiff = Show[rediffplot, linePlot]*)


Export[FileNameJoin[{TCL4DynamicsFolder,diffPlotName}],finalPlotDiff,ImageResolution->300];


temperatureList


lambdaList


Length[lambdaList]


difftcl2tcl4[[9,39]]
TCl2Normalized[[9,39]]
TCL4Normalized[[9,39]]


difftcl2tcl4 // MatrixForm


rules = ArrayRules[difftcl2tcl4 /. 0 -> 10^-16];
(* Sort the rules by their value (from smallest to largest) *)
sortedRules = SortBy[rules, Last];

redRules2 = Take[sortedRules, -1]

(* Extract just the positions from the rules *)
TopRedRegion2 = redRules2[[All, 1]]


BlueRules2 = Take[sortedRules, 1]

(* Extract just the positions from the rules *)
TopBlueRegion2 = BlueRules2[[All, 1]]










