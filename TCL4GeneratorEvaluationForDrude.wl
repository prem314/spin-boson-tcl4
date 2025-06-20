(* ::Package:: *)

(* ::Title:: *)
(*TCL4GeneratorEvaluationForDrude*)


Get["MyFunctions.wl"]
Get["DrudeCutoffFunctions.wl"]
Get["TCL4SpinBosonFunctions.wl"]
Get["LatexExtract.wl"]


LoadVar[TCL4Generator]
LoadVar[TCL2Generator]
LoadVar[TCL0Generator]


(* ::Chapter:: *)
(*Spectral densities*)


LambdaNumVal      = 1;
TNumVal = 8;


GammaNumVal       = 0.01;


\[CapitalOmega]Val = 1.6;
BetaNumVal        = 1/TNumVal;

\[Theta]NumVal = N[\[Pi]/2,16];

p1NumVal          = Sin[\[Theta]NumVal];
p3NumVal          = Cos[\[Theta]NumVal];


(*
GammaNumVal       = 0.05;
LambdaNumVal      = 0.05;
\[CapitalOmega]Val = 0.5;
BetaNumVal        = 20;
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


Plot[NumReplace[JDrude[\[Omega]]],{\[Omega],0,40}, PlotRange->Full]


NumReplaceFinal[expr_]:= NumReplace[expr]


(* ::Chapter:: *)
(*Integral calcs*)


TCL2GeneratorNum = TCLNumericalEval[NumReplaceFinal[TCL2Generator], {LambdaNumVal, 2 \[Pi] /BetaNumVal}];


TCL4GeneratorNum = TCLNumericalEval[NumReplaceFinal[TCL4Generator], {LambdaNumVal, 2 \[Pi] /BetaNumVal}, p3NumVal/p1NumVal];


TCL0GeneratorNum = NumReplaceFinal[TCL0Generator];


Re[TCL0GeneratorNum] // MatrixForm
Re[TCL2GeneratorNum] // MatrixForm
Re[TCL4GeneratorNum] // MatrixForm


Re[TCL2GeneratorNum]
Re[TCL4GeneratorNum]


(* ::Text:: *)
(*T = 8, Lambda = 0.2*)
(*{{0, 0, 0, 0}, {0., 0., 0, 0.}, {0., 0.123725, -0.0155178, 0.}, {-0.00154663, 0., 0, -0.0155178}}*)
(*{{0, 0, 0, 0}, {0., 0., 0., 0.}, {0., -0.00478193, -0.0382674, 0.}, {0.00023552, 0., 0., 0.00235902}*)


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





\[Lambda]2 = 1

(\[Lambda]2 TCL2GeneratorNum/(\[Lambda]2^2 TCL4GeneratorNum)) // MatrixForm


L2Num = TCL0GeneratorNum + \[Lambda]2  TCL2GeneratorNum;
L4Num = TCL0GeneratorNum + \[Lambda]2 TCL2GeneratorNum + \[Lambda]2^2 TCL4GeneratorNum;


zeroEigenstates = Re[ NullSpace[L2Num][[1]]/NullSpace[L2Num][[1,1]]]


(* ::Section:: *)
(*Plot*)


z0 = RandomReal[{0,1}]
x0 = RandomReal[{0, Sqrt[1-z0^2]}]
y0 = RandomReal[{0, Sqrt[1-(z0^2 + x0^2)]}]

x0^2 + y0^2 + z0^2

\[Rho]0 = {1,x0,y0,z0}



(*
\[Rho]0 = {1,0,0,-0.5}
*)


InitialStateParameters = {
{AngleBracket[Subscript[\[Sigma],1][0]], \[Rho]0[[2]]},
{AngleBracket[Subscript[\[Sigma],2][0]], \[Rho]0[[3]]},
{AngleBracket[Subscript[\[Sigma],3][0]], \[Rho]0[[4]]}
};

ParameterToLatex[InitialStateParameters]




TCl2Map[t_]:= MatrixExp[Re[L2Num] t]
TCl4Map[t_]:= MatrixExp[Re[L4Num] t]


TCL2Rho[t_, \[Rho]0_] := TCl2Map[t] . \[Rho]0
TCL4Rho[t_, \[Rho]0_] := TCl4Map[t] . \[Rho]0


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
    {t, 0, 200},                        (* Plot range for t *)

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

