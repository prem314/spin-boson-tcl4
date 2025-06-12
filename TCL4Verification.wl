(* ::Package:: *)

(* ::Title:: *)
(*TCL4Verification*)


Get["MyFunctions.wl"]
Get["DrudeCutoffFunctions.wl"]
Get["TCL4SpinBosonFunctions.wl"]

LoadVar[TCL4Generator];
LoadVar[TCL2Generator];

(*
GammaNumVal       = RandomReal[{0.001,0.01}];
LambdaNumVal      = RandomReal[{0.1,1}];
\[CapitalOmega]Val              = RandomReal[{0.1,1}];
BetaNumVal        = RandomReal[{0.1,1}];
\[Theta]NumVal           = RandomReal[{0.1,1}];
*)

(*
GammaNumVal       =  0.00358;
LambdaNumVal      = 0.207;
\[CapitalOmega]Val              = 0.191;
BetaNumVal        = 0.315;
\[Theta]NumVal           = 0.183;

p1NumVal          = Sin[\[Theta]NumVal];
p3NumVal          = Cos[\[Theta]NumVal];
*)

LambdaNumVal      = 7;
TNumVal = 1


GammaNumVal       = 0.01;

\[CapitalOmega]Val = 1.6;

BetaNumVal        = 1/TNumVal;

\[Theta]NumVal = N[\[Pi]/2]

p1NumVal          = Sin[\[Theta]NumVal]
p3NumVal          = Cos[\[Theta]NumVal]

NumReplace[expr_] := expr //. {
  p1 -> p1NumVal,
  p3 -> p3NumVal,
  \[CapitalLambda] -> LambdaNumVal,
  \[Beta] -> BetaNumVal,
  \[Gamma] -> GammaNumVal,
  \[Theta] -> \[Theta]NumVal,
  \[CapitalOmega] -> \[CapitalOmega]Val,
  f[0] -> Limit[fDrude[\[Omega]], \[Omega] -> 0], f'[0] -> 0,
  J -> Function[\[Omega], JDrude[\[Omega]]],
  J[0] -> Limit[JDrude[\[Omega]], {\[Omega] -> 0}],
  f -> Function[\[Omega], fDrude[\[Omega]]]
};

TCL2GeneratorNum = TCLNumericalEval[NumReplace[TCL2Generator], {LambdaNumVal, 2 \[Pi] /BetaNumVal}];

TCL4GeneratorNum = TCLNumericalEval[NumReplace[TCL4Generator], {LambdaNumVal, 2 \[Pi] /BetaNumVal} , p3NumVal/p1NumVal];

TCL0GeneratorNum = NumReplace[TCL0Gen];

Re[TCL0GeneratorNum] // MatrixForm
Re[TCL2GeneratorNum] // MatrixForm
Re[TCL4GeneratorNum] // MatrixForm


LoadVar[TCL4DrudeGeneratorSimplified]


(NumReplace[Tcl2drude] - TCL2GeneratorNum)// MatrixForm


(* ::Text:: *)
(*TCL2 result is consistent.*)


(* ::Chapter:: *)
(*TCL4*)


TCL4SimpleHDrude = TCL4DrudeGeneratorSimplified /. m-> m + 10^-15 I;


TCL4SimpleHDrude2 = NumReplace[TCL4SimpleHDrude //. SpecialSymbolVal];


NumDrudeSum[ar_, UpLim_]:=Module[{SegTripleInt3, FinalResult},
SegTripleInt3 = Simplify[ar];
FinalResult = (SegTripleInt3[[2,1]] + ParallelSum[SegTripleInt3[[2,2]], {n,1, UpLim}] + ParallelSum[SegTripleInt3[[1,2]], {n,1, UpLim},{m,1,UpLim}]);
FinalResult
]


(* ::Section:: *)
(*All terms*)


(* ::Text:: *)
(*Bag = fast, mutable accumulator held internally by the kernel.*)
(*StuffBag = push an item into that accumulator.*)
(**)
(*https://chatgpt.com/share/680cdbf1-b6c4-800d-8da6-6fad7329bae0*)


T1 = AbsoluteTime[];

minUp  = 2;
maxUp  = 2^9;
base   = 2;


upVals = Table[base^i, {i, Log[base, minUp], Log[base,maxUp]}];


upValQ = AssociationThread[upVals -> True];              


idxs   = {{3,1},{3, 2}, {3,3}, {3,4},{4,1},{4,2},{4,3},{4,4}};


(* ------------------------------------------------------------------ *)

ClearAll[buildSeriesOptimized];

buildSeriesOptimized[{i_, j_}] := Module[

  {cell   = TCL4SimpleHDrude2[[i, j]],

   exact  = N @ Re@TCL4GeneratorNum[[i, j]],

   c0, f1, f2,

   k,

   s1 = 0.0,   (* running \[Sum] f1  *)

   s2 = 0.0,   (* running \[Sum] f2  *)

   bag = Internal`Bag[],

   relVals

  },

  {c0, f1, f2} = {N@cell[[2, 1]], cell[[2, 2]], cell[[1, 2]]};


  Do[

    (* update 1-D part *)

    s1 += N[f1 /. n -> k];


    (* update 2-D part: row k and (k-1) entries of column k *)

    (*Two single summation to capture the effect of incrimented value of k.
    One runs till k and the other till k-1 to avoid duplicate value.
    *)

    s2 += Sum[N[f2 /. {n -> k, m -> c}], {c, 1, k}] +

          Sum[N[f2 /. {n -> r, m -> k}], {r, 1, k - 1}];


    (* NEW: store the approximation iff k is in upVals *)

    If[KeyExistsQ[upValQ, k],

      Internal`StuffBag[bag, c0 + s1 + s2]

    ],

    {k, 1, maxUp}

  ];


  (* convert Bag \[RightArrow] List and compute relative errors *)

  relVals = Map[

    If[exact == 0,

       If[# == 0, 0, \[Infinity]],

       (exact - #)/exact] &,

    Internal`BagPart[bag, All]

  ];
  Transpose[{upVals, relVals}]

];


(* ------------------------------------------------------------------ *)

DistributeDefinitions[buildSeriesOptimized, idxs, maxUp, upVals, upValQ];


DrudeVsGeneralPlotData2 = ParallelTable[

   buildSeriesOptimized[idx],

   {idx, idxs}
];

T2 = AbsoluteTime[];

DT = T2 - T1;


(* ::Text:: *)
(*6: 0.0024584728*)
(*7: 0.0098762744*)
(*8: 0.0400908325*)
(*9: 0.1695223033*)


Print["Took ", DT/3600, " hours."]


TimeReq[n_]:= 0.0098762744 * 2^(2(n-7))


TimeReq[11]*60


RelAccuracy[n_]:= 10^-4 * 2^-(n-9)


TimeVsReward[n_]:= {TimeReq[n], RelAccuracy[n]}


DrudeVsGeneralPlotData


Dimensions[DrudeVsGeneralPlotData]


PlotParameters = {GammaNumVal, LambdaNumVal, \[CapitalOmega]Val,BetaNumVal, \[Theta]NumVal}


DrudeVsGeneralPlotData = {DrudeVsGeneralPlotData2, PlotParameters};


Dimensions[DrudeVsGeneralPlotData]


TCL4GenPlot = VerificationPlotCombined[DrudeVsGeneralPlotData, {{3,1},{3,2},{3,3},{3,4},{4,1},{4,2},{4,3},{4,4}}]

plotFileName = FileNameJoin[{TCL4DynamicsFolder,"VerificationPlot_Lambda" <> ToString[LambdaNumVal] <> "_T" <> ToString[TNumVal] <> ".png"}];
Export[plotFileName, TCL4GenPlot];


(*
DumpVar[DrudeVsGeneralPlotData]
*)


Re[TCL0GeneratorNum] // MatrixForm
Re[TCL2GeneratorNum] // MatrixForm
Re[TCL4GeneratorNum] // MatrixForm



