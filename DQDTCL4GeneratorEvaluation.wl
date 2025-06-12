(* ::Package:: *)

(* ::Title:: *)
(*TCL4GeneratorEvaluationForDrude*)


Get["MyFunctions.wl"]
Get["DrudeCutoffFunctions.wl"]
Get["TCL4SpinBosonFunctions.wl"]
Get["LatexExtract.wl"]


LoadVar[TCL4Generator]
LoadVar[TCL2Generator]


(* ::Chapter:: *)
(*Spectral densities*)


(*
GammaNumVal       = 0.01;
LambdaNumVal      = 1;
\[CapitalOmega]Val = 1.6;

BetaNumVal        = 0.2;
\[Omega]maxVal=1; 
\[Theta]NumVal = N[\[Pi]/2,16];

p1NumVal          = Sin[\[Theta]NumVal];
p3NumVal          = Cos[\[Theta]NumVal];
*)



(*DQD parameters*)

(*Set \[Epsilon] initially to zero to reproduce the fig 2a first.*)

\[Epsilon]DQD = 1;
tc = 0.5;
GammaNumVal   = 0.4;
BetaNumVal    = 1;
\[Omega]maxVal       = 1;
\[CapitalOmega]Val = (\[Epsilon]DQD^2 + 4 tc^2)^0.5

p1NumVal          = 2 tc/\[CapitalOmega]Val
p3NumVal          = \[Epsilon]DQD/\[CapitalOmega]Val;

LambdaNumVal = 2 \[Pi]/BetaNumVal;

(* Constructing Data_list *)
DQDParameters = {
  {\[Epsilon], \[Epsilon]DQD},
  {Subscript[t, c], tc},                 (* Assuming tc represents t_c *)
  {\[Lambda]^2 \[Gamma], GammaNumVal},
  {\[Beta], BetaNumVal},
  {Subscript[\[Omega], "max"], \[Omega]maxVal}
};

ParameterToLatex[DQDParameters]


N[\[CapitalOmega]Val,4]


(*
GammaNumVal       = 0.05;
LambdaNumVal      = 0.05;
\[CapitalOmega]Val = 0.5;
BetaNumVal        = 20;
\[Omega]maxVal = 8
\[Theta]NumVal = 3.14/4;

p1NumVal          = Sin[\[Theta]NumVal];
p3NumVal          = Cos[\[Theta]NumVal];
*)


(* ::Text:: *)
(*I tried decreasing the working precision here, but that messes up the numerics.*)


(*
GammaNumVal       = RandomReal[{0.001,0.01}];
LambdaNumVal      = RandomReal[{0.1,1}];
\[CapitalOmega]Val              = RandomReal[{0.1,1}];
BetaNumVal        = RandomReal[{0.1,1}];
\[Theta]NumVal           = RandomReal[{0.1,1}];

p1NumVal          = Sin[\[Theta]NumVal];
p3NumVal          = Cos[\[Theta]NumVal];
*)


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


JDrude[\[Omega]]


JDQD[\[Omega]_]:=\[Gamma] \[Omega] (1-Sin[\[Omega]]/\[Omega]) Exp[-\[Omega]^2 /(2 \[Omega]max^2)];
fDQD[\[Omega]_]:=JDQD[\[Omega]] Coth[\[Beta] \[Omega]/2]

NumReplaceDQD[expr_] := expr //. {
  p1 -> p1NumVal,
  p3 -> p3NumVal,
  \[CapitalLambda] -> LambdaNumVal,
  \[Beta] -> BetaNumVal,
  \[Gamma] -> GammaNumVal,
  \[Omega]max -> \[Omega]maxVal,
  \[Theta] -> \[Theta]NumVal,
  \[CapitalOmega] -> \[CapitalOmega]Val,
  f[0] -> Limit[fDQD[\[Omega]], \[Omega] -> 0], 
  f'[0] -> 0, 
  J -> Function[\[Omega], JDQD[\[Omega]]],
  J[0] -> Limit[JDQD[\[Omega]], {\[Omega] -> 0}],
  J'[0]->0,
  f -> Function[\[Omega], fDQD[\[Omega]]]
};


Plot[NumReplaceDQD[JDQD[\[Omega]]],{\[Omega],0,40}, PlotRange->Full]


Plot[NumReplace[JDrude[\[Omega]]],{\[Omega],0,40}, PlotRange->Full]


NumReplaceFinal[expr_]:= NumReplaceDQD[expr]


(* ::Chapter:: *)
(*Integral calcs*)


TCL2GeneratorNum = TCLNumericalEval[NumReplaceFinal[TCL2Generator], {LambdaNumVal, 2 \[Pi] /BetaNumVal}];


TCL4GeneratorNum = TCLNumericalEval[NumReplaceFinal[TCL4Generator], {LambdaNumVal, 2 \[Pi] /BetaNumVal}, p3NumVal/p1NumVal];


TCL0GeneratorNum = NumReplaceFinal[TCL0Gen];


Re[TCL0GeneratorNum] // MatrixForm
Re[TCL2GeneratorNum] // MatrixForm
Re[TCL4GeneratorNum] // MatrixForm


(*
DQDTCL0GeneratorNum = {TCL0GeneratorNum, DQDParameters};
DQDTCL2GeneratorNum = {TCL2GeneratorNum, DQDParameters};
DQDTCL4GeneratorNum = {TCL4GeneratorNum, DQDParameters};

DumpVar[DQDTCL0GeneratorNum]
DumpVar[DQDTCL2GeneratorNum]
DumpVar[DQDTCL4GeneratorNum]
*)


(* ::Chapter:: *)
(*Plot*)


Get["MyFunctions.wl"]
Get["DrudeCutoffFunctions.wl"]
Get["TCL4SpinBosonFunctions.wl"]
Get["LatexExtract.wl"]


(*
LoadVar[DQDTCL0GeneratorNum]
LoadVar[DQDTCL2GeneratorNum]
LoadVar[DQDTCL4GeneratorNum]

TCL0GeneratorNum = DQDTCL0GeneratorNum[[1]]
TCL2GeneratorNum = DQDTCL2GeneratorNum[[1]]
TCL4GeneratorNum = DQDTCL4GeneratorNum[[1]]
*)


(*
GammaTarget = 0.4

\[Lambda]2 = GammaTarget/GammaNumVal
*)





\[Lambda]2 = 0.001



L2Num = TCL0GeneratorNum + \[Lambda]2  TCL2GeneratorNum;
L4Num = TCL0GeneratorNum + \[Lambda]2 TCL2GeneratorNum + \[Lambda]2^2 TCL4GeneratorNum;


(*
z0 = RandomReal[{0,1}]
x0 = RandomReal[{0, Sqrt[1-z0^2]}]
y0 = RandomReal[{0, Sqrt[1-(z0^2 + x0^2)]}]

x0^2 + y0^2 + z0^2

\[Rho]0 = {1,x0,y0,z0}
*)


\[Rho]0 = {1,0,0,-0.5}


InitialStateParameters = {
{AngleBracket[Subscript[\[Sigma],1][0]], \[Rho]0[[2]]},
{AngleBracket[Subscript[\[Sigma],2][0]], \[Rho]0[[3]]},
{AngleBracket[Subscript[\[Sigma],3][0]], \[Rho]0[[4]]}
};

ParameterToLatex[InitialStateParameters]




TCL2Rho[t_, \[Rho]0_] := MatrixExp[Re[L2Num] t] . \[Rho]0
TCL4Rho[t_, \[Rho]0_] := MatrixExp[Re[L4Num] t] . \[Rho]0


(* Ensure necessary variables are defined before calling the function *)
(* Example definitions (replace with your actual data/paths): *)
(* TCL2Rho[t_] := Table[Sin[k t / 10], {k, 1, 5}]; *)
(* TCL4Rho[t_] := Table[Cos[k t / 20], {k, 1, 5}]; *)
(* TCL4DynamicsFolder = FileNameJoin[{$HomeDirectory, "Desktop", "PlotOutput"}]; *)


PlotFunc[n_, range_:Full] := Module[{plot, fullPath, yLabel, legendLabels, plotStyles},

  (* Define plot styles for consistency *)
  plotStyles = {
    Directive[Thickness[0.006], Darker[Blue]],  (* Thicker line, distinct color 1 *)
    Directive[Thickness[0.006], Darker[Red]]   (* Thicker line, distinct color 2 *)
   };

  (* Construct a meaningful y-axis label using TraditionalForm for math symbols *)
  (* Using HoldForm to prevent evaluation of AngleBracket, Subscript, etc. *)
  yLabel = Style[TraditionalForm@HoldForm@AngleBracket[Subscript["\[Sigma]", n - 1]][t], FontSize -> 16, Bold, Black];

  (* Construct legend labels using TraditionalForm and Style *)
  legendLabels = {
    Style[Row[{TraditionalForm@AngleBracket[Subscript["\[Sigma]", n - 1]], "(t) (TCL2)"}], FontSize -> 14],
    Style[Row[{TraditionalForm@AngleBracket[Subscript["\[Sigma]", n - 1]], "(t) (TCL4)"}], FontSize -> 14]
   };

  plot = Plot[
    {TCL2Rho[t,\[Rho]0][[n]], TCL4Rho[t,\[Rho]0][[n]]}, (* Functions to plot *)
    {t, 0, 40},                        (* Plot range for t *)

    (* --- Style Enhancements --- *)
    PlotTheme -> "Scientific",         (* Use a theme designed for scientific plots - adds frame, adjusts ticks *)
    PlotStyle -> plotStyles,            (* Apply the defined thicker lines and colors *)
    ImageSize -> Large,                (* Increase the overall size of the plot graphic *)

    (* --- Labeling Enhancements --- *)
    AxesLabel -> {
      Style["t", FontSize -> 16, Bold, Black], (* X-axis label: Larger, Bold *)
      yLabel                                   (* Y-axis label: Using pre-defined styled label *)
     },
    LabelStyle -> Directive[Black, FontSize -> 14], (* Default style for labels (like ticks) *)
    (* TicksStyle -> Directive[FontSize -> 12], (* Uncomment to manually set tick font size if needed *) *)


(* --- Legend Enhancements --- *)
    PlotLegends -> Placed[              (* Use Placed for better legend position control *)
      LineLegend[
        plotStyles,                    (* Match legend colors/styles to PlotStyle *)
        legendLabels,                  (* Use the pre-defined styled labels *)
        LegendFunction -> (Framed[#, RoundingRadius -> 4, FrameStyle -> LightGray] &), (* Optional: Add a subtle frame *)
        LegendMargins -> {{5, 5}, {5, 5}} (* Optional: Adjust spacing around legend content *)
      ],
      (* Corrected Position: Use Scaled coordinates for inside placement *)
      (* Scaled[{x, y}] where x,y range from 0 to 1 relative to plot dimensions *)
      (* Example: Scaled[{0.8, 0.8}] places the anchor point 80% across and 80% up *)
      Scaled[{0.8, 0.8}]
      (* Other options: *)
      (* Scaled[{0.1, 0.9}] (* Near top-left corner *) *)
      (* Scaled[{0.9, 0.9}] (* Near top-right corner *) *)
      (* Scaled[{0.5, 0.1}] (* Near bottom-center *) *)
    ],


    (* --- Other Options --- *)
    PlotRange -> range,                  (* Ensure the entire range of the functions is shown *)
    AxesStyle -> Directive[Black, Thickness[0.002]] (* Make axes lines slightly thicker *)
    (* FrameTicksStyle -> Directive[FontSize -> 12] (* Alternative/addition to TicksStyle with Frame *) *)
  ];

  (* --- Exporting --- *)
  (* Ensure the output directory exists *)
  If[Not@DirectoryQ[TCL4DynamicsFolder],
    CreateDirectory[TCL4DynamicsFolder, CreateIntermediateDirectories -> True];
    Print["Created directory: ", TCL4DynamicsFolder];
  ];

  fullPath = FileNameJoin[{TCL4DynamicsFolder, StringJoin[{"Sigma", ToString[n - 1], "Plot.pdf"}]}];

  (* Export the plot as a high-quality PDF (vector format) *)
  Export[fullPath, plot, "PDF"];
  (* For raster formats like PNG, you might add ImageResolution: *)
  (* Export[FileNameJoin[{TCL4DynamicsFolder, StringJoin[{"Sigma", ToString[n-1], "Plot.png"}]}], plot, "PNG", ImageResolution -> 300]; *)

  plot (* Return the plot object so it displays in the notebook *)
];

(* Example Usage: *)
(* PlotFunc[3] *)
PlotFunc[2]

PlotFunc[3]

PlotFunc[4]



(* ::Chapter:: *)
(*Rough*)


CorrRhoBell = {{1, 0,  0, 0},  (* Coeffs for I\[Tensor]I, I\[Tensor]X, I\[Tensor]Y, I\[Tensor]Z *)
               {0, 1,  0, 0},  (* Coeffs for X\[Tensor]I, X\[Tensor]X, X\[Tensor]Y, X\[Tensor]Z *)
               {0, 0, -1, 0},  (* Coeffs for Y\[Tensor]I, Y\[Tensor]X, Y\[Tensor]Y, Y\[Tensor]Z *)
               {0, 0,  0, 1}}; (* Coeffs for Z\[Tensor]I, Z\[Tensor]X, Z\[Tensor]Y, Z\[Tensor]Z *)


FirstQubitEvolve[t_, CorrRho_, TCL2Rho_]:= Module[{ReducedRho, i, j, EvolvedReducedState, EvolvedRho},
EvolvedRho = ConstantArray[0, {4,4}];
For[j=1, j<5, j++, 
ReducedRho = CorrRho[[All, j]];
EvolvedReducedState = TCL2Rho[t, ReducedRho];
For[i=1,i<5,i++, EvolvedRho[[i,j]] = EvolvedReducedState[[i]]]];
EvolvedRho
]

FirstQubitEvolve[1, CorrRhoBell, TCL2Rho] //MatrixForm


(* Define the single-qubit Pauli matrices *)
(* SigmaMatrices[[1]] = Identity, SigmaMatrices[[2]] = PauliX, etc. *)
SigmaMatrices = {IdentityMatrix[2], PauliMatrix[1], PauliMatrix[2], PauliMatrix[3]};

(* Function to convert a 4x4 density matrix to its 4x4 Pauli coefficients matrix *)
DensityMatrixToPauliCoefficients[rho_?MatrixQ] :=
  Module[{dimCheck = Dimensions[rho]},
   If[dimCheck != {4, 4},
    Print["Error: Density matrix (rho) must be 4x4. Dimensions received: ", dimCheck];
    Return[$Failed];
   ];
   Table[
    Tr[rho . KroneckerProduct[SigmaMatrices[[k]], SigmaMatrices[[l]]]],
    {k, 1, 4}, {l, 1, 4}
   ]
  ];

(* Function to convert a 4x4 Pauli coefficients matrix back to a 4x4 density matrix *)
PauliCoefficientsToDensityMatrix[coeffs_?MatrixQ] :=
  Module[{dimCheck = Dimensions[coeffs]},
   If[dimCheck != {4, 4},
    Print["Error: Pauli coefficients matrix (coeffs) must be 4x4. Dimensions received: ", dimCheck];
    Return[$Failed];
   ];
   (1/4) * Sum[
     coeffs[[k, l]] * KroneckerProduct[SigmaMatrices[[k]], SigmaMatrices[[l]]],
     {k, 1, 4}, {l, 1, 4}
   ]
  ];


(* --- Example Usage --- *)

(* 1. Create an example two-qubit density matrix (e.g., a Bell state |Phi+>) *)
(* |Phi+> = (1/Sqrt[2]) * (|00> + |11>) *)
(* ket00 = {1,0,0,0}; ket11 = {0,0,0,1}; *)
(* psiPlus = (1/Sqrt[2]) (ket00 + ket11); *)
(* rhoPhiPlus = Outer[Times, psiPlus, Conjugate[psiPlus]]; *)
(* More directly: *)
rhoPhiPlus = {{1/2, 0, 0, 1/2}, {0, 0, 0, 0}, {0, 0, 0, 0}, {1/2, 0, 0, 1/2}};

Print["Original Density Matrix (rhoPhiPlus):"];
Print[MatrixForm[rhoPhiPlus]];

(* 2. Convert it to Pauli basis representation *)
pauliCoeffsPhiPlus = DensityMatrixToPauliCoefficients[rhoPhiPlus];

Print["\nPauli Coefficients Matrix (P):"];
If[pauliCoeffsPhiPlus =!= $Failed, Print[MatrixForm[N[pauliCoeffsPhiPlus]]]];
(* Expected for Phi+: p_II=1, p_XX=1, p_YY=-1, p_ZZ=1, others 0 *)
(* P = {{1,0,0,0},{0,1,0,0},{0,0,-1,0},{0,0,0,1}} *)

(* 3. Convert the Pauli coefficients back to the density matrix *)
rhoReconstructed = PauliCoefficientsToDensityMatrix[pauliCoeffsPhiPlus];

Print["\nReconstructed Density Matrix (from P):"];
If[rhoReconstructed =!= $Failed, Print[MatrixForm[Chop[N[rhoReconstructed]]]]];

(* 4. Verify the reconstruction *)
If[rhoReconstructed =!= $Failed,
 Print["\nVerification: Is reconstructed rho close to original rho? ",
  Chop[N[rhoReconstructed]] == Chop[N[rhoPhiPlus]]]
 ];

(* Example with a product state: |0><0| tensor |+><+| *)
(* rhoA = {{1,0},{0,0}}; *) (* |0><0| *)
(* rhoB = {{1/2, 1/2},{1/2, 1/2}}; *) (* |+><+| *)
(* rhoProduct = KroneckerProduct[rhoA, rhoB]; *)
rhoProduct = {{1/2, 1/2, 0, 0}, {1/2, 1/2, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};

Print["\n\n--- Example with Product State rhoA x rhoB ---"];
Print["Original Product Density Matrix (rhoProduct):"];
Print[MatrixForm[rhoProduct]];

pauliCoeffsProduct = DensityMatrixToPauliCoefficients[rhoProduct];
Print["\nPauli Coefficients Matrix for Product State (P_product):"];
If[pauliCoeffsProduct =!= $Failed, Print[MatrixForm[N[pauliCoeffsProduct]]]];
(* Expected for |0><0| (A) x |+><+| (B):
   Pauli vec for A (I,Z): {1,0,0,1}
   Pauli vec for B (I,X): {1,1,0,0}
   P_product should be OuterProduct[{1,0,0,1}, {1,1,0,0}]
   P_product = {{1,1,0,0},{0,0,0,0},{0,0,0,0},{1,1,0,0}}
*)

rhoReconstructedProduct = PauliCoefficientsToDensityMatrix[pauliCoeffsProduct];
Print["\nReconstructed Product Density Matrix:"];
If[rhoReconstructedProduct =!= $Failed, Print[MatrixForm[Chop[N[rhoReconstructedProduct]]]]];

If[rhoReconstructedProduct =!= $Failed,
 Print["\nVerification for product state: ",
  Chop[N[rhoReconstructedProduct]] == Chop[N[rhoProduct]]]
 ];


NegativityEval[t_, CorrRho_, TCL2Rho_]:= Total[Abs[Eigenvalues[PauliCoefficientsToDensityMatrix[FirstQubitEvolve[t,CorrRho, TCL2Rho]]]]] -1

NegativityEval[1, CorrRhoBell, TCL2Rho]


NegativityEval[1, CorrRhoBell, TCL4Rho]


Plot[
  {NegativityEval[t, CorrRhoBell, TCL2Rho], NegativityEval[t, CorrRhoBell, TCL4Rho]},
  {t, 0, 10},
  PlotLegends -> {"TCL2", "TCL4"}
]
