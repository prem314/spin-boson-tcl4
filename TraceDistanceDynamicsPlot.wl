(* ::Package:: *)

(* ::Title:: *)
(* Trace distance  dynamics  using TCL2 and TCL4 master equation for spin boson model*)


(* ::Input:: *)
(*(* Initial satate *)*)
(*vec1={{1},{1},{0},{0}};*)
(*vec2={{1},{-1},{0},{0}};*)
(*(*Time interval over which we analyse the trace distance dynamics *)*)
(*tvals=Range[0,40,0.1];*)


(* ::Input:: *)
(*(*Import TCL2 and TCL4 generator*)*)


Get["MyFunctions.wl"]
Get["DrudeCutoffFunctions.wl"]
Get["TCL4SpinBosonFunctions.wl"]


LoadVar[TCL4Generator]
LoadVar[TCL2Generator]
LoadVar[TCL0Generator]


GammaNumVal       = 0.01;
(*LambdaNumVal      = 0.5;*)
\[CapitalOmega]Val = 1.6;
(*BetaNumVal        = 2;*)

\[Theta]NumVal = N[\[Pi]/2,16];

p1NumVal          = Sin[\[Theta]NumVal];
p3NumVal          = Cos[\[Theta]NumVal];


NumReplace[expr_] := expr //. {
  p1 -> p1NumVal,
  p3 -> p3NumVal,
  \[Gamma] -> GammaNumVal,
  \[CapitalOmega] -> \[CapitalOmega]Val,
  f[0] -> Limit[fDrude[\[Omega]], \[Omega] -> 0], f'[0] -> 0,
  J -> Function[\[Omega], JDrude[\[Omega]]],
  f -> Function[\[Omega], fDrude[\[Omega]]]
};


ClearAll[NumericalIntegralTCL];
NumericalIntegralTCL[gen : {_List ..}, lambdaVal_, TVal_] := Module[
  {t, nTerms, cons, single, dbl, total},
  (* 1) substitute numerical values *)
  t      = Map[NumReplace, gen, {2}] //. {\[CapitalLambda] -> lambdaVal, \[Beta] -> 1/TVal};
  
  (* 2) inspect how many slots each entry has *)
  nTerms = Length@t[[1, 1]];
  (* 3) constant and single always present *)
  cons   = t[[All, All, 1]];
  single = Map[ContourIntegralBetter[#, lambdaVal, 1/TVal] & , t[[All, All, 2]], {2}];
  (* 4) only do double if there are at least 4 elements *)
  If[nTerms >= 4,
    dbl   = Map[DoubleContourIntegralBetter[#, lambdaVal, 1/TVal] &, t[[All, All, 3]], {2}];
    total = cons + single + dbl,
    (* else skip double *)
    dbl   = ConstantArray[0, Dimensions[cons]];
    total = cons + single
  ];
  (* 5) return an Association *)
  <|
    "Const"  -> cons,
    "Single" -> single,
    "Double" -> dbl,
    "Total"  -> total
  |>
];




L0gen = NumReplace[TCL0Generator];


(* ::Input:: *)
(*(*Calculating TCL2 generator numerically*)*)


(* ::Input:: *)
(*computeLtotalfinalTCL2[lambdaVal_,TVal_]:=Module[{gen2,TCL2Final,LtotalModifiedTCL2},*)
(*gen2=NumericalIntegralTCL[TCL2Generator,lambdaVal,TVal];*)
(*TCL2Final=gen2["Total"];*)
(*LtotalModifiedTCL2=L0gen+TCL2Final;*)
(*N[LtotalModifiedTCL2] ];*)
(**)


(* ::Input:: *)
(*tracedataOptimizedTCL2[lambdaVal_,TVal_]:=Module[{LtotalfinalTCL2,evolvedTCL21,evolvedTCL22,blochTCL21,blochTCL22},*)
(*LtotalfinalTCL2=computeLtotalfinalTCL2[lambdaVal,TVal];*)
(*  Table[evolvedTCL21=MatrixExp[LtotalfinalTCL2 t] . vec1;*)
(*evolvedTCL22=MatrixExp[LtotalfinalTCL2 t] . vec2;*)
(*blochTCL21=Flatten[Rest[evolvedTCL21]];*)
(*blochTCL22=Flatten[Rest[evolvedTCL22]];*)
(*{t,trdist[blochTCL21,blochTCL22]},{t,tvals}]];*)


(* ::Input:: *)
(*(*Precompute Ltotalfinal for TCL4 for fixed lambdaVal and TVal*)*)
(*     computeLtotalfinal[lambdaVal_,TVal_]:=Module[{gen2,TCL2Final,gen4,raw4,TCL4GenFinal,LtotalModified},*)
(*gen4=NumericalIntegralTCL[TCL4Generator,lambdaVal,TVal];*)
(*raw4=gen4["Total"];*)
(*TCL4GenFinal=modifyMatrix[raw4,p3NumVal/p1NumVal];*)
(*gen2=NumericalIntegralTCL[TCL2Generator,lambdaVal,TVal];*)
(*TCL2Final=gen2["Total"];*)
(*LtotalModified=L0gen+TCL2Final+TCL4GenFinal;*)
(*N[LtotalModified] ];*)
(**)


(* ::Input:: *)
(*(*Generate tracedata for tcl4 using precomputed Ltotalfinal*)*)
(*     tracedataOptimizedTCL4[lambdaVal_,TVal_]:=Module[{Ltotalfinal,evolved1,evolved2,bloch1,bloch2},*)
(*Ltotalfinal=computeLtotalfinal[lambdaVal,TVal];*)
(*Table[evolved1=MatrixExp[Ltotalfinal t] . vec1;*)
(*evolved2=MatrixExp[Ltotalfinal t] . vec2;*)
(*bloch1=Flatten[Rest[evolved1]];*)
(*bloch2=Flatten[Rest[evolved2]];*)
(*{t,trdist[bloch1,bloch2]},{t,tvals}]];*)


PlotValues = {tracedataOptimizedTCL2[8,1.4],tracedataOptimizedTCL4[8,1.4],tracedataOptimizedTCL2[1,8],tracedataOptimizedTCL4[1,8]};


MoreDataPoints = {tracedataOptimizedTCL2[8,4],tracedataOptimizedTCL4[8,4], tracedataOptimizedTCL2[2,0.6],tracedataOptimizedTCL4[2,0.6]};


WhiteCurve = {tracedataOptimizedTCL2[4,1],tracedataOptimizedTCL4[4,1]};


Font1 = 20
Font2 = 18
Font3 = 12


(* ::Input:: *)
(*figtracedis=ListLinePlot[Join[PlotValues,MoreDataPoints, WhiteCurve],*)
(*PlotRange->All,*)
(*LabelStyle->Directive[Black,{Font2,GrayLevel[0.1]}],*)
(*Frame->True,*)
(*FrameTicks->True,*)
(*FrameTicksStyle->Black,*)
(*AspectRatio-> 1.1,*)
(*PlotStyle->{*)
(*Directive[Red,Thick,Line],*)
(*Directive[Darker[Red],Thick,Dashed],*)
(*Directive[Blue,Thick,Line],*)
(*Directive[Darker[Blue],Thick,Dashed],*)
(*Directive[Green,Thick,Line],*)
(*Directive[Darker[Green],Thick,Dashed],*)
(*Directive[Black,Thick,Line],*)
(*Directive[Darker[Black],Thick,Dashed],*)
(*Directive[Brown,Thick,Line],*)
(*Directive[Darker[Brown],Thick,Dashed]*)
(*},*)
(*FrameLabel->{Style["Time", Font1],Style["Trace distance",  Font1]},*)
(*RotateLabel->True,*)
(*PlotLegends->Placed[LineLegend[{*)
(*"TCL2: \[CapitalLambda] = 8, T = 1.4","TCL4: \[CapitalLambda] = 8, T = 1.4",*)
(*"TCL2: \[CapitalLambda] = 1, T = 8","TCL4: \[CapitalLambda] = 1, T = 8", *)
(*"TCL2: \[CapitalLambda] = 8, T = 4","TCL4: \[CapitalLambda] = 8, T = 4",*)
(*"TCL2: \[CapitalLambda] = 2, T = 0.6","TCL4: \[CapitalLambda] = 2, T = 0.6",*)
(*"TCL2: \[CapitalLambda] = 4, T = 1","TCL4: \[CapitalLambda] = 4, T = 1"},*)
(*LabelStyle->Directive[Black,Font2],*)
(*Background->Opacity[0.75,White]],*)
(*Scaled[{0.77,0.7}]],*)
(*ImageSize->600]*)
(**)


(* ::Input:: *)
(*Export[FileNameJoin[{TCL4DynamicsFolder,"tdistancelambda.pdf"}],figtracedis]*)
