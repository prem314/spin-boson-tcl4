(* ::Package:: *)

(* Generates individual DQD sigma_i(t) plots using cached numerical generators.
   The caches are the DQDTCL*GeneratorNum.mx files written by
   DQDTCL4GeneratorEvaluation.wl via DumpVar. *)

ClearAll["Global`*"];

scriptDir = If[
   StringQ[$InputFileName] && StringLength[$InputFileName] > 0,
   DirectoryName[ExpandFileName[$InputFileName]],
   Directory[]
];
SetDirectory[scriptDir];

oldScriptPath = FileNameJoin[{scriptDir, "DQDTCL4GeneratorEvaluation.wl"}];
If[! FileExistsQ[oldScriptPath],
  Print["Missing old DQD script: ", oldScriptPath];
  Exit[1];
];
Print["Using cached generators produced by: ", oldScriptPath];

cacheSearchDirs = DeleteDuplicates @ Select[
   ExpandFileName /@ {
     scriptDir,
     "/home/premkr/git/papers/mathematica_wiki",
     "/home/premkr/Dropbox/work/configurations/git/mymodule/mathematica_wiki",
     "/home/premkr/Dropbox/master_eqn_project/prem_mathematica_files"
   },
   DirectoryQ
];

ClearAll[findCacheFile, loadCacheSymbol];
findCacheFile[symbolName_String] := Module[{fileName, matches},
  fileName = symbolName <> ".mx";
  matches = Select[FileNameJoin[{#, fileName}] & /@ cacheSearchDirs, FileExistsQ];
  If[matches === {},
    Print["Missing cached numerical generator: ", fileName];
    Print["Searched directories: ", cacheSearchDirs];
    Exit[1];
  ];
  First[matches]
];

SetAttributes[loadCacheSymbol, HoldFirst];
loadCacheSymbol[sym_Symbol] := Module[{symbolName, file},
  symbolName = SymbolName[Unevaluated[sym]];
  file = findCacheFile[symbolName];
  Get[file];
  Print["Loaded ", symbolName, " from ", file];
];

loadCacheSymbol[DQDTCL0GeneratorNum];
loadCacheSymbol[DQDTCL2GeneratorNum];
loadCacheSymbol[DQDTCL4GeneratorNum];

If[! And @@ (MatchQ[#, {_?MatrixQ, _List}] & /@
      {DQDTCL0GeneratorNum, DQDTCL2GeneratorNum, DQDTCL4GeneratorNum}),
  Print["One or more cached generator files did not contain {matrix, parameters}."];
  Exit[1];
];

TCL0GeneratorNum = DQDTCL0GeneratorNum[[1]];
TCL2GeneratorNum = DQDTCL2GeneratorNum[[1]];
TCL4GeneratorNum = DQDTCL4GeneratorNum[[1]];

If[! And @@ ((Dimensions[#] === {4, 4}) & /@
      {TCL0GeneratorNum, TCL2GeneratorNum, TCL4GeneratorNum}),
  Print["Expected all cached generators to be 4 x 4 matrices."];
  Print["Dimensions: ", Dimensions /@ {TCL0GeneratorNum, TCL2GeneratorNum, TCL4GeneratorNum}];
  Exit[1];
];

lambda2 = 1;
rho0 = {1, 0, 0, -0.5};
tmax = 40;
gibbsPopulation = -Tanh[
   DQDTCL0GeneratorNum[[2, 4, 2]] Abs[TCL0GeneratorNum[[3, 2]]] / 2
];

curves = {
   <|
     "Tag" -> "TCL0",
     "Label" -> "TCL0",
     "Matrix" -> Re[TCL0GeneratorNum],
     "Style" -> Directive[Thickness[0.006], Darker[Blue]]
   |>,
   <|
     "Tag" -> "TCL0_TCL2",
     "Label" -> "TCL2",
     "Matrix" -> Re[TCL0GeneratorNum + lambda2 TCL2GeneratorNum],
     "Style" -> Directive[Thickness[0.006], Darker[Red]]
   |>,
   <|
     "Tag" -> "TCL0_TCL2_TCL4",
     "Label" -> "TCL4",
     "Matrix" -> Re[TCL0GeneratorNum + lambda2 TCL2GeneratorNum + lambda2^2 TCL4GeneratorNum],
     "Style" -> Directive[Thickness[0.006], Darker[Green]]
   |>
};

plotGroups = {
   <|
     "Tag" -> "TCL0",
     "Curves" -> curves[[{1}]]
   |>,
   <|
     "Tag" -> "TCL2",
     "Curves" -> curves[[{1, 2}]]
   |>,
   <|
     "Tag" -> "TCL4",
     "Curves" -> curves[[{1, 2, 3}]]
   |>
};

outputDir = FileNameJoin[{scriptDir, "dqd_sigma_plots"}];
If[! DirectoryQ[outputDir],
  CreateDirectory[outputDir, CreateIntermediateDirectories -> True];
];

font1 = 18;
font2 = 16;

subfigureLabels = {"(a)", "(b)", "(c)"};

ClearAll[sigmaPanel];
sigmaPanel[group_Association, sigmaIndex_Integer] := Module[
  {componentIndex, panelTitle, groupCurves, curveExpressions, plotStyles, legendLabels},
  componentIndex = sigmaIndex + 1;
  panelTitle = Row[{
     subfigureLabels[[sigmaIndex]],
     "  ",
     TraditionalForm @ AngleBracket[Subscript["\[Sigma]", sigmaIndex]][t]
   }];
  groupCurves = group["Curves"];
  curveExpressions = ((MatrixExp[#["Matrix"] t] . rho0)[[componentIndex]] &) /@ groupCurves;
  plotStyles = (#["Style"] &) /@ groupCurves;
  legendLabels = (Style[#["Label"], FontSize -> font2] &) /@ groupCurves;
  If[MemberQ[{"TCL2", "TCL4"}, group["Tag"]] && sigmaIndex === 3,
    curveExpressions = Append[curveExpressions, gibbsPopulation];
    plotStyles = Append[plotStyles, Directive[Black, Dotted, Thickness[0.005]]];
    legendLabels = Append[
      legendLabels,
      Style["Gibbs", FontSize -> font2]
    ];
  ];
  Plot[
    Evaluate[curveExpressions],
    {t, 0, tmax},
    PlotTheme -> "Scientific",
    PlotStyle -> plotStyles,
    ImageSize -> 360,
    AspectRatio -> 1,
    Frame -> True,
    FrameLabel -> {
      Style["t", FontSize -> font1, Black],
      None
    },
    PlotLabel -> Style[panelTitle, FontSize -> font1, Black],
    LabelStyle -> Directive[Black, FontSize -> font2],
    PlotLegends -> Placed[
      LineLegend[
        plotStyles,
        legendLabels,
        LegendFunction -> (
          Framed[
            #,
            RoundingRadius -> 4,
            FrameStyle -> LightGray,
            Background -> Opacity[0.9, White]
          ] &
        ),
        LegendMargins -> {{5, 5}, {5, 5}}
      ],
      Scaled[{0.5, 0.25}]
    ],
    ImagePadding -> {{58, 16}, {54, 46}},
    ImageMargins -> 0,
    PlotRange -> All,
    PlotRangePadding -> Scaled[0.03],
    FrameStyle -> Directive[Black, Thickness[0.002]]
  ]
];

ClearAll[combinedSigmaPlot];
combinedSigmaPlot[group_Association] := Module[
  {panels, figure, file},
  panels = Table[sigmaPanel[group, sigmaIndex], {sigmaIndex, 1, 3}];
  figure = GraphicsRow[
    panels,
    Spacings -> 0,
    ImageSize -> 1080
  ];
  file = FileNameJoin[
    {outputDir, "DQDSigma_" <> group["Tag"] <> ".pdf"}
  ];
  Export[file, figure, "PDF"];
  Print["Wrote ", file];
  file
];

exportedFiles = combinedSigmaPlot /@ plotGroups;

Print["Finished. Exported ", Length[exportedFiles], " PDF files."];
