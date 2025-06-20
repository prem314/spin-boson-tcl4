(* ::Package:: *)

(* ::Title:: *)
(*TCLVsHEOMFidelity*)


Get["MyFunctions.wl"]
Get["DrudeCutoffFunctions.wl"]
Get["TCL4SpinBosonFunctions.wl"]
Get["LatexExtract.wl"]


LoadVar[TCL0Generator]
LoadVar[TCL4Generator]
LoadVar[TCL2Generator]

TCL4SimpleGen2 = TCL4Generator;
TCL2Integrand4 = TCL2Generator;


(* \[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash] *)
(* 1.  Load & parameters                                        *)
(* \[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash] *)

Get["TCL4SimpleHGenerator.mx"];
Get["TCL2Generator.mx"];


(*
GammaNumVal       = 0.02;
LambdaNumVal      = 1;
\[CapitalOmega]Val = 1.6;
BetaNumVal        = 1;

\[Theta]NumVal = N[\[Pi]/2];

p1NumVal          = Sin[\[Theta]NumVal];
p3NumVal          = Cos[\[Theta]NumVal];
*)


(*
NonMarkovianity red region parameter
*)


(*
GammaNumVal       = 0.01;
LambdaNumVal      = 2;
\[CapitalOmega]Val = 1.6;
TNumVal = 0.6
BetaNumVal        = 1/TNumVal;

\[Theta]NumVal = N[\[Pi]/2]

p1NumVal          = Sin[\[Theta]NumVal]
p3NumVal          = Cos[\[Theta]NumVal]
*)



HEOMFidelityParameter = {
  {\[Lambda]^2 \[Gamma], GammaNumVal},
  {\[CapitalLambda], LambdaNumVal},                 (* Assuming tc represents t_c *)
  {\[CapitalOmega], \[CapitalOmega]Val},
  {\[Beta], BetaNumVal},
  {\[Theta], \[Theta]NumVal}
};

ParameterToLatex[HEOMFidelityParameter]


(*
GammaNumVal       = 0.01;
LambdaNumVal      = 0.5;
\[CapitalOmega]Val = 0.5;
BetaNumVal        = 2;

\[Theta]NumVal = 3.14/4;

p1NumVal          = Sin[\[Theta]NumVal];
p3NumVal          = Cos[\[Theta]NumVal];
*)




GammaNumVal    = RandomReal[{0.001, 0.01}];
LambdaNumVal     = RandomReal[{0.5, 1.5}];
\[CapitalOmega]Val = RandomReal[{1, 2}];
BetaNumVal       = RandomReal[{0.1, 0.3}];
\[Theta]NumVal = RandomReal[{3, 4}]/4;

p1NumVal          = Sin[\[Theta]NumVal];
p3NumVal          = Cos[\[Theta]NumVal];



NumReplace[expr_] := expr //. {
  p1 -> p1NumVal,
  p3 -> p3NumVal,
  \[CapitalLambda] -> LambdaNumVal,
  \[Beta] -> BetaNumVal,
  \[Gamma] -> GammaNumVal,
  \[CapitalOmega] -> \[CapitalOmega]Val,
  f[0] -> Limit[fDrude[\[Omega]], \[Omega] -> 0], f'[0] -> 0,
  J -> Function[\[Omega], JDrude[\[Omega]]],
  f -> Function[\[Omega], fDrude[\[Omega]]]
};




TCL4GenFinal = TCLNumericalEval[NumReplace[TCL4Generator], {LambdaNumVal, 2 \[Pi] /BetaNumVal}, p3NumVal/p1NumVal];


TCL2Final = TCLNumericalEval[NumReplace[TCL2Generator], {LambdaNumVal, 2 \[Pi] /BetaNumVal}];


TCL2Final // MatrixForm


TCL4GenFinal // MatrixForm


(* ::Package:: *)
(**)


(* Define a function to save a list of Mathematica objects to a file,
   one object per line, in Python-readable format.
   Overwrites the file if it exists. *)

(* 2. Specify the output filename *)
(*
outputFilename = FileNameJoin[{NotebookDirectory[],"python_objects.txt"}];
*)
outputFilename = "/home/premkr/Dropbox/work/projects/tcl4_dynamics/python_objects.txt";

ExportToPythonReadableLines[objectsList_List, filename_String] := Module[
  {pythonStrings, fileContent, filepath},

  (* Convert each object in the input list to its Python-like string representation *)
  (* Uses InputForm for standard representation and replaces {} with [] for lists *)
  pythonStrings = Map[
    StringReplace[ToString[#, InputForm], {"{" -> "[", "}" -> "]"}] &,
    objectsList
  ];

  (* Join the individual object strings with newline characters *)
  fileContent = StringRiffle[pythonStrings, "\n"];

  (* Create the full path: file saved in the same folder as the current notebook *)
  filepath = outputFilename;

  (* Export the combined string (with newlines) to the text file, overwriting existing content *)
  (* "Text" format ensures plain text output and overwriting *)
  Export[filepath, fileContent, "Text"];

  (* Optional: Return the filepath *)
  Return[filepath];
];

(* Example Usage: *)

(* 1. Define some Mathematica objects (scalars, lists/arrays, strings) *)
myObjects = {GammaNumVal, LambdaNumVal, \[CapitalOmega]Val, BetaNumVal, \[Theta]NumVal,
Re[TCL2Final],
Re[TCL4GenFinal]                    (* Null - Note: Mathematica Null -> Python None *)
};


Re[TCL2Final]


GammaNumVal


(* 3. Call the function to export the objects *)
exportedFile = ExportToPythonReadableLines[myObjects, outputFilename];

(* 4. Optional: Print confirmation *)
Print["Mathematica objects exported to: ", exportedFile];



(*\[LongDash]\:2011 adjust both paths if necessary \[LongDash]\:2011*)
cmd  = "source ~/anaconda3/etc/profile.d/conda.sh ; \
		cd /home/premkr/Dropbox/work/projects/tcl4_dynamics
        conda activate tcl && \
        python tcl_vs_heom.py";

result = RunProcess[{"bash", "-lc", cmd}];   (* no 2nd argument *)

(* examine what came back *)
result["ExitCode"]
result["StandardOutput"]
result["StandardError"]


parameterLaTeXInline[]


FidelityTCLData = Join[myObjects, {NumReplace[TCL0Generator]}]
(*
DumpVar[FidelityTCLData]
*)

