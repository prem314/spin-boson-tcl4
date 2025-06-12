(* ::Package:: *)

(* ::Title:: *)
(*NonMarkovianitySpinBoson*)


(* ::Section:: *)
(*Initialization*)


Wiki[DecayTimeCalcSBMProject]


Get["MyFunctions.wl"]
Get["DrudeCutoffFunctions.wl"]
Get["TCL4SpinBosonFunctions.wl"]


LoadVar[TCL0Generator]

LoadVar[TCL4Generator]

LoadVar[TCL2Generator]


GammaNumVal       = 0.01;
(*LambdaNumVal      = 0.5;*)
\[CapitalOmega]Val = 1.6;
(*BetaNumVal        = 2;*)

\[Theta]NumVal = N[\[Pi]/2,16];

p1NumVal          = Sin[\[Theta]NumVal];
p3NumVal          = Cos[\[Theta]NumVal];


PrecVal = 16

GammaIntVal = IntegerPart[1000 GammaNumVal];


NumValList = {GammaNumVal, \[CapitalOmega]Val, \[Theta]NumVal, PrecVal}


(*
GammaNumVal    = RandomReal[{0.005, 0.05}];
LambdaNumVal     = RandomReal[{0.5, 1}];
\[CapitalOmega]Val = RandomReal[{0.5, 1}];
BetaNumVal       = RandomReal[{0.5, 1}];
p1NumVal         = RandomReal[{0.5, 1}];
p3NumVal         = 0;
*)


NumReplace[expr_] := expr //. {
  p1 -> p1NumVal,
  p3 -> p3NumVal,
  \[Gamma] -> GammaNumVal,
  \[CapitalOmega] -> \[CapitalOmega]Val,
  f[0] -> Limit[fDrude[\[Omega]], \[Omega] -> 0], f'[0] -> 0,
  J -> Function[\[Omega], JDrude[\[Omega]]],
  f -> Function[\[Omega], fDrude[\[Omega]]]
};


L0gen = NumReplace[TCL0Generator];


(* ::Section:: *)
(*Something else*)


BetaNumVal = 1
LambdaNumVal = 1


NumReplace[TCL2Generator]


TCL2SampleGen = TCLNumericalEval[NumReplace[TCL2Generator], {LambdaNumVal, 2 \[Pi] /BetaNumVal}]


TCL2SampleGen


(*
(* How many sample points? *)
n = 5000;                       (* increase for a denser cloud *)

(* Generate {x,y,z} triples *)
pts = Table[
   Module[{u = RandomReal[], v = RandomReal[]},
     {xcomp[u, v], ycomp[u, v], zcomp[v]}
   ],
   {n}
];

(* 3-D scatter plot on top of a faint sphere for context *)
Graphics3D[
 {
   Opacity[.15], EdgeForm[None], Sphere[],      (* unit sphere \[OpenCurlyDoubleQuote]ghost\[CloseCurlyDoubleQuote] *)
   Directive[Red, PointSize[Medium]], Point[pts]
 },
 BoxRatios -> 1, Axes -> False, Lighting -> "Neutral",
 ViewPoint -> {2, -2, 1}
]
*)


(*
(* --- single random sample and plot --- *)
Module[{u = RandomReal[], v = RandomReal[], p1, p2},
  p1 = {xcomp[u, v],  ycomp[u, v],  zcomp[v]};
  p2 = {xcompa[u, v], ycompa[u, v], zcompa[v]};
  
  Graphics3D[
    {
      {Opacity[.15], Sphere[]},           (* faint reference sphere   *)
      {Red,  PointSize[Large],  Point[p1]},  (* first point: red       *)
      {Blue, PointSize[Large],  Point[p2]}   (* second point: blue     *)
    },
    Axes      -> True,
    Boxed     -> False,
    PlotRange -> 1.1,                     (* show full sphere & pts   *)
    ImageSize -> 400
  ]
]
*)


(* ::Chapter::Closed:: *)
(*Long Calc*)


(* ::Section:: *)
(*TCL4 and TCL2 combined evaluation*)


lambdaList=Range[2,8,0.5]
temperatureList=Range[0.6,8,0.5]

numPoint = 10




(* ::Text:: *)
(*ParallelTable gave error below.*)


T1 = AbsoluteTime[];

TotalData = Table[NonMarkovianityEval[\[CapitalLambda],T,numPoint],{T,temperatureList},{\[CapitalLambda],lambdaList}];

T2 = AbsoluteTime[];
\[CapitalDelta]T = (T2 - T1)/60;



Print["Took ",\[CapitalDelta]T, " minutes."];


NMTraceDistanceFullData = {TotalData, lambdaList, temperatureList}

DumpVar[NMTraceDistanceFullData]


(* ::Chapter:: *)
(*Plot*)


Get[FileNameJoin[{TCL4DynamicsFolder, "NMTraceDistanceFullData.mx"}]]

{TotalData, lambdaList, temperatureList} = NMTraceDistanceFullData


lambdaList=Range[2,8,0.5]
temperatureList=Range[0.6,8,0.5]



Dimensions[TotalData]


datasimpletc4 = TotalData[[All, All, 1,1]]


TotalData[[3,2,1]];
datasimpletc4[[3,2]];


datasimpletc2 = TotalData[[All, All, 2 , 1]]


TotalData[[3,2,2]];
datasimpletc2[[3,2]];


Print["Took ", \[CapitalDelta]T, " minutes."]


datasimpletcl4r=datasimpletc4[[Reverse[Range[Dimensions[datasimpletc4][[1]]]],All]];


(*Save["datasimpletcl4_gamma02_thetapi2.mx",datasimpletc4]*)


EmitSound[Sound[SoundNote[]]]


(* ::Section:: *)
(*Data Load*)


(*
Get["NomMarkovTLambdaData2016.mx"]
{datasimpletc4, datasimpletc2, NumValList} = NomMarkovTLambdaData2016;
*)
flatArray = Flatten[datasimpletc2];
smallestNonZeroTCL2 = Min[DeleteCases[flatArray, 0 | 0.]]


(* ::Subsection:: *)
(*Relative Data*)


relativeDiffPercent = (datasimpletc4 - datasimpletc2);


(* ::Subsection:: *)
(*Plot*)


TCL4DynamicsFolder


SaveMyPlot[plotVariable_] :=
  Export[
    FileNameJoin[{
      TCL4DynamicsFolder,
      StringJoin[{SymbolName[Unevaluated[plotVariable]], ToString[GammaIntVal], ToString[PrecVal],".pdf"}]
    }],
    plotVariable,
    "PDF"
  ];
SetAttributes[SaveMyPlot, HoldFirst];


TCL4NonMarkovTLambdaPlot = ListDensityPlot[Log10[datasimpletc4 + smallestNonZeroTCL2/10],
  FrameLabel -> {"\[CapitalLambda]", "T"},
  ColorFunction -> "SunsetColors", (* or "TemperatureMap", "SunsetColors", etc. *)
  PlotLegends -> Automatic,
  DataRange -> {
    MinMax[lambdaList], (* X-axis range *)
    MinMax[temperatureList] (* Y-axis range *)
  },PlotRange->Full
]

SaveMyPlot[TCL4NonMarkovTLambdaPlot]


TCL2NonMarkovTLambdaPlot = ListDensityPlot[Log10[datasimpletc2 + smallestNonZeroTCL2/10],
  FrameLabel -> {"\[CapitalLambda]", "T"},
  ColorFunction -> "SunsetColors", (* or "TemperatureMap", "SunsetColors", etc. *)
  PlotLegends -> Automatic,
  DataRange -> {
    MinMax[lambdaList], (* X-axis range *)
    MinMax[temperatureList] (* Y-axis range *)
  },PlotRange->Full
]

SaveMyPlot[TCL2NonMarkovTLambdaPlot]


difftcl2tcl4 = datasimpletc4 - datasimpletc2;


DiffMax = Max[difftcl2tcl4]
DiffMin = Abs[Min[difftcl2tcl4]]


RelativeNonMarkovTLambdaPlot = Module[{values=Flatten[difftcl2tcl4],posMax,negMax,linThresh,symlog,symlogPosMax,symlogNegMax,normalizeValue,(*The new function to map data to[-1,1]*)plotColorFunction,legendColorFunction,legendTicks},(*---1. Define Asymmetric Range and Logarithmic Scaling---*)posMax=DiffMax;
negMax=DiffMin;
(*Base the linear threshold on the smaller positive scale to capture its detail*)linThresh=0.1*posMax;
If[linThresh==0,linThresh=1];
symlog=Function[x,Sign[x]*Log[1+Abs[x]/linThresh]];
symlogPosMax=symlog[posMax];
symlogNegMax=Abs[symlog[-negMax]];
If[symlogPosMax==0,symlogPosMax=1];
If[symlogNegMax==0,symlogNegMax=1];
(*---2. Create the Normalization Function---*)(*This is the key step.It maps any data value to a point in the[-1,1] range.*)normalizeValue=Function[val,Which[val>0,symlog[val]/symlogPosMax,(*Maps (0,posMax] to (0,1]*)val<0,-(Abs[symlog[val]]/symlogNegMax),(*Maps[-negMax,0) to[-1,0)*)True,0]];
(*---3. Define Color Functions---*)(*The plot's color function now uses the normalized value*)plotColorFunction=Function[val,Module[{norm=normalizeValue[val]},Which[norm>0,Blend[{White,Red},norm],(*norm is already 0..1*)norm<0,Blend[{White,Blue},-norm],(*-norm makes it 0..1*)True,White]]];
(*A simpler version for the legend,which works directly on the[-1,1] scale*)legendColorFunction=Function[norm,Which[norm>0,Blend[{White,Red},norm],norm<0,Blend[{White,Blue},-norm],True,White]];
(*---4. Create Ticks for the Symmetric Legend---*)(*The real-world values you want to label on the bar*)posTicksData={posMax/100,posMax/10,posMax};
negTicksData={-negMax,-negMax/10, -negMax/100, -negMax/1000};
allTicksData=Sort@DeleteDuplicates@Join[negTicksData,{0},posTicksData];
(*Create {position,label} pairs.The position is the value's location on the normalized[-1,1] bar.*)legendTicks=Table[{normalizeValue[t],NumberForm[t,2,ExponentFunction->(If[-3<#<4,Null,#]&)]},{t,allTicksData}];
(*---5. Generate the Plot---*)ListDensityPlot[difftcl2tcl4,ColorFunction->plotColorFunction,ColorFunctionScaling->False,DataRange->{MinMax[lambdaList],MinMax[temperatureList]},PlotLegends->BarLegend[(*Use the simple legend function and a SYMMETRIC[-1,1] range*){legendColorFunction,{-1,1}},Ticks->legendTicks,LegendMarkerSize->250],FrameLabel->{"\[CapitalLambda]","T"},PlotRange->Full,ImageSize->280,FrameStyle->Directive[Black]]]

SaveMyPlot[RelativeNonMarkovTLambdaPlot]


relativeDiffPercent // MatrixForm


(* ::Section:: *)
(*Analysis of different regions*)


TraceDistancePairPlot[pairs_]:=Module[{plotData,legendData},
plotData = Flatten[
    Table[
      {i, j} = k;
      {TotalData[[i, j, 1, 2]], TotalData[[i, j, 2, 2]]}
      , {k, pairs}
    ], 1];

legendData = Flatten[
    Table[
      {i, j} = pair;
      temp = temperatureList[[i]];
      lambda = lambdaList[[j]];
      {StringJoin["TCL4 (\[CapitalLambda]=", ToString[lambda], ", T=", ToString[temp], ")"],
       StringJoin["TCL2 (\[CapitalLambda]=", ToString[lambda], ", T=", ToString[temp], ")"]}
      , {pair, pairs}
    ], 1];
     
ListPlot[
    plotData,
    Joined -> True,
    ImageSize -> 1000,
    PlotLegends -> Placed[
      LineLegend[
        legendData,
        LegendMarkers -> Graphics[{Thick, Line[{{-1, 0}, {1, 0}}]}],
        LegendMarkerSize -> 30
        ],
      Below],
    PlotStyle -> Thickness[0.003],
    PlotRange -> Full,
    AxesLabel -> {"Time", "Trace Distance"},
    PlotLabel -> "Trace Distance vs. Time for various T and \[CapitalLambda]"
    ]
]


Length[Position[relativeDiffPercent, x_ /; x > 0]]


(* ::Text:: *)
(*7 Positive values.*)


(* ::Section:: *)
(*Red region*)


(* Get a list of rules: {i, j} -> value *)
rules = ArrayRules[relativeDiffPercent /. 0 -> 10^-16];

(* Sort the rules by their value (from smallest to largest) *)
sortedRules = SortBy[rules, Last];


redRules2 = Take[sortedRules, -1]

(* Extract just the positions from the rules *)
TopRedRegion2 = redRules2[[All, 1]]


TraceDistancePairPlot[TopRedRegion2]


(* ::Section:: *)
(*Blue region*)


top5Rules2 = Take[sortedRules, 1];

(* Extract just the positions from the rules *)
TopBlueRegion2 = top5Rules2[[All, 1]]


TraceDistancePairPlot[TopBlueRegion2]


(* ::Section:: *)
(*White region*)


explicitZeroIndices = Position[relativeDiffPercent, 0]


TraceDistancePairPlot[explicitZeroIndices]


(* ::Section:: *)
(*higher and lower values together*)


numCurves = 2

(* Take the last 5 rules (the ones with the largest values) *)
top5Rules = Take[sortedRules, numCurves];

(* Extract just the positions from the rules *)
TopBlueRegion = top5Rules[[All, 1]]


bottom5Rules = Take[sortedRules, -numCurves];

TopRedRegion = bottom5Rules[[All, 1]]


TraceDistancePairPlot[Join[TopBlueRegion, TopRedRegion]]


(* ::Section:: *)
(*Specific values*)


lambdaList
temperatureList


relativeDiffPercent[[9,6]]


TraceDistancePairPlot[{{9,6}}]


relativeDiffPercent[[1,1]]
TraceDistancePairPlot[{{1,1}}]


(* ::Chapter:: *)
(*Rough*)


NegValues = Position[relativeDiffPercent, x_ /; x < 0]

NegValuesSmall = RandomSample[NegValues, 20]

TraceDistancePairPlot[NegValuesSmall]


Join[TopBlueRegion, TopRedRegion]


PosValues = Position[relativeDiffPercent, x_ /; x > 0]

TraceDistancePairPlot[PosValues]


(*
TraceDistancePlot[T_, \[CapitalLambda]_]:=ListPlot[{TotalData[[T, \[CapitalLambda],1,2]], TotalData[[T, \[CapitalLambda],2,2]]},
  Joined -> True, (* This ensures you are plotting lines *)
  ImageSize -> 1000,
  PlotLegends -> {
  StringJoin["TCL4 (", ToString[lambdaList[[\[CapitalLambda]]]],",", ToString[temperatureList[[\[CapitalLambda]]]],")"],
  StringJoin["TCL2 (", ToString[lambdaList[[\[CapitalLambda]]]],",", ToString[temperatureList[[\[CapitalLambda]]]],")"]
},
  PlotStyle -> Thickness[0.002],
  PlotRange->Full
]



TraceDistancePlot[15,13]
*)


TotalData[[15,13]];


Dimensions[TotalData]


TraceDistanceMultiPlot[tempIndices_List, lambdaIndices_List] :=
 ListPlot[
  Flatten[
   Table[{TotalData[[i, j, 1, 2]], TotalData[[i, j, 2, 2]]}, {i, 
     tempIndices}, {j, lambdaIndices}],
   2
   ],
  Joined -> True,
  ImageSize -> 1000,
  PlotLegends -> Placed[
    LineLegend[
     Flatten[
      Table[{StringJoin["TCL4 (\[CapitalLambda]=", 
          ToString[lambdaList[[j]]], ", T=", 
          ToString[temperatureList[[i]]], ")"], 
         StringJoin["TCL2 (\[CapitalLambda]=", 
          ToString[lambdaList[[j]]], ", T=", 
          ToString[temperatureList[[i]]], ")"]}, {i, tempIndices}, {j, 
        lambdaIndices}],
      2
      ],
     LegendMarkers -> Graphics[{Thick, Line[{{-1, 0}, {1, 0}}]}],
     LegendMarkerSize -> 30
     ],
    Below],
  PlotStyle -> Thickness[0.003],
  PlotRange -> Full,
  AxesLabel -> {"Time", "Trace Distance"},
  PlotLabel -> 
   "Trace Distance vs. Time for various \[CapitalLambda] and T"
  ]
  
TraceDistanceMultiPlot[{1,4,8,15},{1}]


TraceDistanceMultiPlot[{10},{1,5,10}]


TraceDistanceMultiPlot[{2}, {1, 6,13}]


TraceDistanceMultiPlot[{1,2,6,15},{13}]


(* ::Chapter:: *)
(*Specific plot*)


numPoint = 20;

Get["TCL4SpinBosonFunctions.wl"]

T1 = AbsoluteTime[];

Vals = NonMarkovianityEval[7,1,numPoint, 0.1];

Print["TCL4 NM = ",Vals[[1,1]]]
Print["TCL2 NM = ",Vals[[2,1]]]

Print["TCL4 - TCL2 = ",Vals[[1,1]] - Vals[[2,1]]]

ListPlot[{Vals[[1, 2]], Vals[[2, 2]]},
  Joined -> True, (* This ensures you are plotting lines *)
  ImageSize -> 1000,
  PlotLegends -> {"TCL4", "TCL2"},
  PlotStyle -> Thickness[0.002],
  PlotRange->Full
]

T2 = AbsoluteTime[];

(T2 - T1)/60


(* ::Text:: *)
(*Checking the red region at \[CapitalLambda] = 2 and T = 0.6 at higher precision (dt = 0.01, trace distance tolerance = 0.001) gives TCL4 - TCL2 = 0.0000335524 for numPoint = 10.*)
(**)
(*TCL4 - TCL2 = 0.0000490002 for dt = 0.01, trace distance tolerance = 0.0001 and numPoint = 10.*)
(**)
(*TCL4 - TCL2 = 0000437882 for dt = 0.01, trace distance tolerance = 0.0001 and numPoint = 10. Too 7 minutes to calculate.*)


(* ::Chapter:: *)
(*Precision test*)


(* ::Text:: *)
(*TCL4 Value at precision = 16 is 0.0421467*)
(**)
(*TCL2 val = 0.225345*)
(**)
(*Setting precision to 17 did not change anything. very precise result.*)


(*
TCL4NMEval[1,7,18]
TCL2NMEval[1,7,18]
*)


(* ::Text:: *)
(*  - TCL2 is sufficient for the markovian regime.*)
(*  - red and blue region means TCL2 undershoots or overshoots against the exact result.*)
(*  - Antipodal vs general plot.*)
