(* ::Package:: *)

(* ::Title:: *)
(*SBMProjectNMPlot*)


(* ::Chapter:: *)
(*Data load*)


Get["TCL4SpinBosonFunctions.wl"]


(* ::Input:: *)
(*(*\[Gamma]=0.01;\[CapitalOmega]=1.6;\[Theta]=\[Pi]/2; StateDt=0.125*)*)


(* ::Input:: *)
(*lambdaList=Range[0.2,8,0.2]; *)
(*temperatureList=Range[0.2,8,0.2]; *)


(* ::Input:: *)
(*Loadtcl2=Get[FileNameJoin[{NotebookDirectory[],"datasimpletcl2_gamma01_thetapi2_modifiedtcl2_4096.mx"}]];*)


(* ::Input:: *)
(*Loadtcl4=Get[FileNameJoin[{NotebookDirectory[],"datasimpletcl4_gamma01_thetapi2_modifiedtcl4_4096.mx"}]];*)


(* ::Input:: *)
(*lambdaList*)


Get["TCLMarkovianityResonanceCondition.wl"]


(* ::Chapter:: *)
(*Plots*)


\[CapitalOmega] = 1.6


(* ::Section:: *)
(*TCL2*)


(* ::Input:: *)
(*ColorScheme = "ReverseGrayTones"*)


TCl2Normalized = Log10[Loadtcl2 + 10^-6];

(* Get the min and max of your data just once *)
dataRange = MinMax[TCl2Normalized];
{dmin, dmax} = dataRange;

(* This new function manually scales any input value from your data's range 
  to the [0, 1] range that GrayLevel expects.
*)
universalColorFunction = (GrayLevel[1 - ((# - dmin)/(dmax - dmin))]) &;

(* 1. THE PLOT *)
figtcl2 = ListDensityPlot[TCl2Normalized,
   FrameLabel -> {"\[CapitalLambda]", "T"},
   FrameStyle -> Directive[Black],
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
finalPlotTCL2 = Legended[
   Show[figtcl2, linePlot],
   (* The legend now uses the same function and range, ensuring a perfect match *)
   BarLegend[{universalColorFunction, dataRange},
    LegendLayout -> "Column", LegendMarkerSize -> 220]
   ];

(* Display the final, correct plot *)
finalPlotTCL2


Export[FileNameJoin[{TCL4DynamicsFolder,"tcl2lamdaTplot.pdf"}],finalPlotTCL2,ImageResolution->300];


(* ::Section:: *)
(*TCL4*)


TCL4Normalized = Log10[Loadtcl4 + 10^-6];

(* 1. THE PLOT *)
figtcl4 = ListDensityPlot[TCL4Normalized,
   FrameLabel -> {"\[CapitalLambda]", "T"},
   FrameStyle -> Directive[Black],
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
    LegendLayout -> "Column", LegendMarkerSize -> 220]
   ]


Export[FileNameJoin[{TCL4DynamicsFolder,"tcl4lamdaTplot.pdf"}],finalPlotTCL4,ImageResolution->300];


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
(*(*---4. Create Ticks for the Symmetric Legend---*)(*The real-world values you want to label on the bar*)posTicksData={0.0001,0.0003,posMax};*)
(*negTicksData={-negMax,-0.01, -0.0001};*)
(*allTicksData=Sort@DeleteDuplicates@Join[negTicksData,{0},posTicksData];*)
(*(*Create {position,label} pairs.The position is the value's location on the normalized[-1,1] bar.*)legendTicks=Table[{normalizeValue[t],NumberForm[t,2,ExponentFunction->(If[-3<#<4,Null,#]&)]},{t,allTicksData}];*)
(*(*---5. Generate the Plot---*)ListDensityPlot[difftcl2tcl4,ColorFunction->plotColorFunction,ColorFunctionScaling->False,DataRange->{MinMax[lambdaList],MinMax[temperatureList]},PlotLegends->BarLegend[(*Use the simple legend function and a SYMMETRIC[-1,1] range*){legendColorFunction,{-1,1}},Ticks->legendTicks,LegendMarkerSize->250],FrameLabel->{"\[CapitalLambda]","T"},PlotRange->Full,ImageSize->280,FrameStyle->Directive[Black]]];*)


finalPlotDiff = Show[rediffplot, linePlot]


Export[FileNameJoin[{TCL4DynamicsFolder,"redifftcl.pdf"}],finalPlotDiff,ImageResolution->300];


difftcl2tcl4 // MatrixForm


difftcl2tcl4[[5,-6]]
