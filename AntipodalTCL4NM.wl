(* ::Package:: *)

(* ::Title:: *)
(*AntipodalTCL4NM*)


start = AbsoluteTime[];



Get["MyFunctions.wl"]
Get["DrudeCutoffFunctions.wl"]
Get["TCL4SpinBosonFunctions.wl"]


LoadVar[TCL4Generator]

LoadVar[TCL2Generator]


GammaNumVal       = 0.01;
(*LambdaNumVal      = 0.5;*)
\[CapitalOmega]Val = 1.6;
(*BetaNumVal        = 2;*)

\[Theta]NumVal = N[\[Pi]/2,16];

p1NumVal          = Sin[\[Theta]NumVal];
p3NumVal          = Cos[\[Theta]NumVal];


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



modifyMatrix[A_?MatrixQ, c_] := Module[{B = A},
  B[[2, All]] = c\[ThinSpace]B[[4, All]];
  B
];


Get["TCL4SpinBosonFunctions.wl"]



(* \[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash] *)
(* 2.  NumericalIntegralTCL with auto\[Hyphen]double detection                  *)
(* \[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash] *)
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




L0gen={{0,0,0,0},{0,0,-\[CapitalOmega],0},{0,\[CapitalOmega],0,0},{0,0,0,0}}/.{\[CapitalOmega]->\[CapitalOmega]Val};


(* Density matrix definition*)
den[{x_,y_,z_}]:=0.5 (IdentityMatrix[2]+x PauliMatrix[1]+y PauliMatrix[2]+z PauliMatrix[3])


(*trace distance definition*)
trdist[{x1_,y1_,z1_},{x2_,y2_,z2_}]:= 0.5 *Re[Sqrt[(x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2]];


(* parametrising the arbitrary state on the bloch sphere*)
xcomp[u_,v_]:=Cos[2Pi u]Sin[ArcCos[2 v-1]]
ycomp[u_,v_]:=Sin[2Pi u]Sin[ArcCos[2v-1]]
zcomp[v_]:=Cos[ArcCos[2v-1]]


(*Antipodal state for given state *)
xcompa[u_,v_]:=Cos[Pi+2Pi u]Sin[Pi-ArcCos[-1+2 v]]
ycompa[u_,v_]:=Sin[Pi+2Pi u]Sin[Pi-ArcCos[-1+2 v]]
zcompa[v_]:=Cos[Pi-ArcCos[-1+2 v]]


StateDt = 0.1


TrDist[x_,y_]:= Re[0.5 Sqrt[Total[(x - y)^2]]]

StateEvolve[t_, vin_, L2_]:= Re[MatrixExp[L2 t] . vin]

FindThresholdTime[nul2_, vin_, L2_] := Module[{t = 0},
  While[TrDist[nul2, StateEvolve[t, vin, L2]] >= 0.001,
   t = t + 20
   ];
   
   (*print["The trace distance = ", TrDist[nul2, StateEvolve[t, vin, L2]]];*)
   
     Return[t]
]

TMaxCalc[L2_]:= Module[{nul1, nul2, antinul},
nul1 = NullSpace[L2][[1]];

nul2 = Re[nul1/nul1[[1]]];
antinul = {1,0,0, - Sign[nul2[[4]]]};

(*print["vectors = ", nul2, antinul];*)

FindThresholdTime[nul2, antinul, L2]
]


computeN[lambdaVal_,TVal_]:=Module[{gen2, TCL2Final,L2total,Ltotalff,map,vec,vec1,allTrace,backflowIntegrals,maxBackflow,tmin,tmax,
dt,tvals,evolved,evolved1,x,y,z,x1,y1,z1,GammaNumVal,LambdaNumVal,OmegaVal,BetaNumVal,ThetaNumVal,gen4,raw4,TCL4GenFinal,LtotalModified},
(*Parameter setup*)

(*Process generator matrix*)

gen4=NumericalIntegralTCL[TCL4Generator, lambdaVal, TVal];
raw4=gen4["Total"];

TCL4GenFinal=modifyMatrix[raw4, p3NumVal/p1NumVal];


gen2      = NumericalIntegralTCL[TCL2Generator, lambdaVal, TVal ];
TCL2Final = gen2["Total"];

L2total=L0gen + TCL2Final;

(*Combine matrices*)LtotalModified = L0gen + TCL2Final + TCL4GenFinal;
Ltotalff = N[LtotalModified];

(*Time evolution setup*)
tmin=0;
tmax=TMaxCalc[L2total];
Print["\[CapitalLambda],T and tmax = ",lambdaVal,",",TVal,",",tmax];

dt=0.1;
tvals=Range[tmin,tmax,dt];
map[t_]:=MatrixExp[Ltotalff*t];
(*State vectors*)
vec={{1},{xcomp[u,v]},{ycomp[u,v]},{zcomp[v]}};
vec1={{1},{xcompa[u,v]},{ycompa[u,v]},{zcompa[v]}};
(*Trace distance calculation*)allTrace=Table[evolved=map[t] . vec;

evolved1=map[t] . vec1;
{x,y,z}=Flatten[evolved[[2;;4]]];
{x1,y1,z1}=Flatten[evolved1[[2;;4]]];
trdist[{x,y,z},{x1,y1,z1}],{u,0,1,StateDt},{v,0,1,StateDt},{t,tvals}];
(*Backflow analysis*)
backflowIntegrals=Table[traceSeries=allTrace[[i,j]];
derivatives=Differences[traceSeries]/dt;
derivatives=Join[derivatives,{0}];
positives=Select[derivatives,#>0&];
Total[positives]*dt,{i,Length[allTrace]},{j,Length[allTrace[[1]]]}];
maxBackflow=Max[Flatten[backflowIntegrals]];
maxBackflow]
 



(* ::Chapter:: *)
(*Calc*)


(*Parameter ranges*)
lambdaList=Range[0.5,50,0.5]; 
temperatureList=Range[0.5,50,0.5];


datasimpletc4=ParallelTable[computeN[\[CapitalLambda],T],{T,temperatureList},{\[CapitalLambda],lambdaList}]


datasimpletcl4r=datasimpletc4[[Reverse[Range[Dimensions[datasimpletc4][[1]]]],All]];



(*Save["datasimpletcl4_gamma01_thetapi2_modifiedtcl4_100_antipodal_lambdaT50.mx",datasimpletc4]*)


figtcl4=ListDensityPlot[Log10[datasimpletc4 + 10^-8],
  FrameLabel -> {"\[CapitalLambda]", "T"},
  ColorFunction -> "Rainbow", (* or "TemperatureMap", "SunsetColors", etc. *)
  PlotLegends -> Automatic,
  DataRange -> {
    MinMax[lambdaList], (* X-axis range *)
    MinMax[temperatureList] (* Y-axis range *)
  },PlotRange->Full
]


end = AbsoluteTime[];



executionTime = end - start


 


Dimensions[datasimpletc4]

