(* ::Package:: *)

(* ::Title:: *)
(*MyFunctions*)


(*\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]*)
(*  Check\:2011pointed evaluator that works for ANY array rank and ANY Map level *)
(*\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]*)

SetAttributes[ArrayEvalCheckpoint, HoldAll];   (* keep the output symbol unevaluated *)

ArrayEvalCheckpoint[f_, arr_, sym_Symbol, level_Integer : 1] := Module[
  {file, todo, tag},
  
  (* Where to store the checkpoint (one file per output symbol) *)
  file = FileNameJoin @ {NotebookDirectory[], SymbolName[Unevaluated[sym]] <> ".mx"};
  
  (* 1 \:25ba  Load previous run, or create a fresh \[OpenCurlyDoubleQuote]Null skeleton\[CloseCurlyDoubleQuote] on first use *)
  If[FileExistsQ[file],
    Get[file];
    Print["\[RightArrow] loaded \[OpenCurlyDoubleQuote]", SymbolName[Unevaluated[sym]], "\[CloseCurlyDoubleQuote] from ", file],
    sym = Map[Null &, arr, {level}];
    Print["\[RightArrow] initialized \[OpenCurlyDoubleQuote]", SymbolName[Unevaluated[sym]], "\[CloseCurlyDoubleQuote] with Nulls at level ", level]
  ];
  
  (* 2 \:25ba  Positions at the requested level (Heads\[RightArrow]False = ignore List heads) *)
  todo = Position[arr, _, {level}, Heads -> False];
  
  (* 3 \:25ba  Walk over those positions; only compute if still Null *)
  Do[
    If[Extract[sym, pos] === Null,
      Print["  \[Bullet] computing at position ", pos];
      sym = ReplacePart[sym, pos -> f[Extract[arr, pos]]];
      DumpSave[file, sym];
      Print["    \[Dash] progress saved"]
    ],
    {pos, todo}
  ];
  
  Print["\:2714 finished \[Dash] \[OpenCurlyDoubleQuote]", SymbolName[Unevaluated[sym]], "\[CloseCurlyDoubleQuote] is now up\:2011to\:2011date."]
]



ParallelMapFinest[f_,expr_] := ParallelMap[f,expr, Length[Dimensions[expr]]]
MapFinest[f_,expr_] := Map[f,expr, Length[Dimensions[expr]]]


(* --- Start of DumpVar Definition --- *)

(* Clear any previous definitions and attributes for DumpVar to ensure a clean setup *)
ClearAll[DumpVar];

(* Define any messages DumpVar might use *)
DumpVar::novalue = "Symbol `1` has no value to dump.";
DumpVar::invalidarg = "Argument `1` to DumpVar is not a symbol or a string name.";

(* Set HoldFirst attribute for DumpVar. This is beneficial when symbols are passed directly,
   allowing SymbolName[Unevaluated[var]] to get the name without evaluating the symbol.
   It does not negatively affect the string case, as strings evaluate to themselves. *)
SetAttributes[DumpVar, HoldFirst];

(* Definition for when DumpVar is called with a STRING (the variable's name) *)
DumpVar[varName_String?StringQ] := Module[{path, symbolToTest},
  path = FileNameJoin[{NotebookDirectory[], varName <> ".mx"}];
  
  (* We need to get the actual symbol to check its ValueQ *)
  symbolToTest = Symbol[varName]; 
  
  If[!ValueQ[symbolToTest],
    Message[DumpVar::novalue, varName];
    Return[$Failed];
  ];
  
  (* DumpSave can accept the symbol's name as a string directly *)
  DumpSave[path, varName]; 
  Print["Wrote (via string name): ", varName, " \[RightArrow] ", path];
  Return[path]; (* Optionally return the path *)
];

(* Definition for when DumpVar is called with a SYMBOL directly *)
DumpVar[var_Symbol] := Module[{actualSymbolName, path},
  (* Because var is held (due to HoldFirst), Unevaluated[var] refers to the symbol itself *)
  actualSymbolName = SymbolName[Unevaluated[var]];
  path = FileNameJoin[{NotebookDirectory[], actualSymbolName <> ".mx"}];
  
  If[!ValueQ[Unevaluated[var]],
    Message[DumpVar::novalue, HoldForm[Unevaluated[var]]]; (* HoldForm for better message formatting *)
    Return[$Failed];
  ];
  
  (* DumpSave needs the symbol itself, which Unevaluated[var] provides here *)
  DumpSave[path, Unevaluated[var]];
  Print["Wrote (via symbol): ", actualSymbolName, " \[RightArrow] ", path];
  Return[path]; (* Optionally return the path *)
];

(* Optional: A catch-all definition for arguments that are neither strings nor symbols *)
DumpVar[other_] := (
  Message[DumpVar::invalidarg, HoldForm[other]];
  Return[$Failed];
);





(* --- End of DumpVar Definition --- *)
SetAttributes[LoadVar, HoldFirst];

LoadVar[var_Symbol] := Module[{path},
  
  path = FileNameJoin[{ NotebookDirectory[], SymbolName[Unevaluated[var]] <> ".mx" }];
  
  If[! FileExistsQ[path],
    Message[LoadVar::nofile, path];
    Return[$Failed];
  ];
  
  Get[path];
  Print["Loaded ", SymbolName[Unevaluated[var]], " \[LeftArrow] ", path];
];

LoadVar::nofile = "MX file `1` not found.";



HomeFolder = "/home/premkr/git/mymodule/mathematica_wiki"


AllSymbolsInExpr[expr_] := Module[{symbols},
  symbols = Union@Cases[expr, s_Symbol :> s, {0, Infinity}];
  Print[symbols];
]


CoeffTerm[expr_, coeff_, n_:1] := Coefficient[expr, coeff, n] coeff^n


(* Converts a expression into a list of its terms. Also works when there is just one term in the expression, hence the if condition*)
TermsAsList[expr_]:= If[Head[expr]===Plus, Apply[List,expr],{expr}]


SaveData[expr_]:= DumpSave[NotebookFileName[] <> ".mx", expr ]; 
SetAttributes[SaveData, HoldFirst];


ReadData[] := Module[{}, Print["Reading Data from ", NotebookFileName[] <> ".mx"];Get[NotebookFileName[] <> ".mx"]];


(* Define the commutator function *)
Commutator[A_, B_] := A . B - B . A

(* Define the anticommutator function *)
AntiCommutator[A_, B_] := A . B + B . A


ZeroCheck[expr_]:= If[Simplify[expr] == 0, Print["The result is consistent"], Print["WARNING: Something is wrong!!!"]]


Clear[randomValue]

randomValue[expr_, exemptVars_ : {}] := Module[
  {syms, varsToReplace, vals, numValue},
  
  (* Extract all unique global symbols from the expression *)
  syms = Union[
    Cases[Hold[expr], s_Symbol /; Context[s] === "Global`", {0, Infinity}]
  ];
  
  (* Determine which symbols to replace by excluding the exempt variables *)
  varsToReplace = Complement[syms, exemptVars];
  
  (* Assign random real values between 0 and 1 to the non-exempt symbols *)
  If[Length[varsToReplace] > 0,
    vals = Thread[varsToReplace -> RandomReal[{0, 1}, Length[varsToReplace]]],
    vals = {}
  ];
  
  (* Print the assigned variables and their values, if any *)
  If[varsToReplace =!= {},
    Print["Assigned Values:"];
    Do[
      Print[varsToReplace[[i]], " -> ", vals[[i, 2]]],
      {i, Length[varsToReplace]}
    ];
  ];
  
  (* Compute the numerical value with replacements *)
  numValue = N[expr /. vals];
  
  (* Return the numerical value and the list of assignments *)
  {numValue, vals}
]

(*
(* Define the expression *)
expr = a + b^2 + c*d;

(* Assign random values, exempting 'a' *)
result = randomValue[expr, {a, c}]
*)

(*
 Sample Output:
Assigned Values:
b -> 0.537
c -> 0.812
d -> 0.123

result =
{a + 0.537^2 + 0.812*0.123, {b -> 0.537, c -> 0.812, d -> 0.123}}
*)


randVal2[expr_, vals_] := Module[{syms, numValue},   
   (* Compute the numerical value *)
   numValue = N[expr /. vals];
   numValue
]

(*
expr = Sin[a] + b^2 + c*d;
randomValue[expr]
*)



(*Collects terms with same argument together*)

CollectTerms[cterm2_,f_]:= Module[{transformedExpr, collectedExpr, CollectedSinresult},
(*Define a replacement rule that replaces Sin[u_] with a new variable sin[u] *)
transformedExpr = cterm2 /. f[u_] :> f1[u];

(* Collect terms based on the new variable sin[u] and simplify each coefficient *)
collectedExpr = Collect[transformedExpr, _f1, FullSimplify];

(* Replace back the sin[u] with Sin[u] *)
CollectedSinresult = collectedExpr /. f1[u_] :> f[u];

CollectedSinresult
]


(*Collects t in the argument of f*)
collectCoefficientT[expr_, f_]:=expr/. f[arg_]:>f[FactorTerms[arg,t]];


(*swaps for each element in the ExprList*)
RepeatedSwap[cterm6_, ExprList_, \[Omega]1_, \[Omega]2_]:= Module[{e1, e2, e3, e4, x1, x2, expr},
expr = cterm6;
Do[ (Clear[e1,e2,e3,e4];
e1 = element Coefficient[expr, element];
e2 = e1 //. {\[Omega]1 -> x1, \[Omega]2 -> x2};
e3 = e2 //. {x1 -> \[Omega]2, x2 -> \[Omega]1};
e4 = expr - e1 + e3;
Print[e1, " swapped with ", e3];
Clear[expr];
expr = e4;
),{element, ExprList}
];
CollectTerms[Expand[expr],Head[First[ExprList]]]
]

VarSwap[expr_, x_, y_]:= Module[{x1, y1, e1, e2},
e1 = expr //. {x -> x1, y-> y1};
e2 = e1 //. {x1 -> y, y1 -> x};
Print[expr, " replaced with ", e2];
e2
]


(*
The following code applies the given operation to all terms that have an expression from  ExprList, one by one.
Careful that the operation is not applied twice unintentionally.
*)

ConditionalOp[cterm6_, ExprList_, OpFunction_, args___] := Module[
  {terms, newTerms},
  
  (* Split cterm6 into a list of terms if it's a sum, otherwise make it a single-element list *)
  terms = TermsAsList[cterm6];
  
  (* Process each term *)
  newTerms = Map[
    Function[term,
      (* Check if the term contains any expression from ExprList *)
      If[
        (Or @@ (Not[FreeQ[term, #]] & /@ ExprList)) || Length[ExprList] == 0,
        (
          OpFunction[term, args]
        ),
        (
          term
        )
      ]
    ],
    terms
  ];
  
  (* Recombine the terms into a single expression *)
  Total[newTerms]
]



PullFactor[expr_, factor_]:= Module[{func},
	func[expr2_]:= Simplify[expr2 /factor];
	factor ConditionalOp[expr, {}, func]
]


SelectFactorsByVariable[expr_, var_Symbol] := Module[{factors, ouput},
  (* Ensure the expression is treated as a product *)
  factors = If[Head[expr] === Times, List @@ expr, {expr}];
  
  (* Select factors that contain the specified variable *)
  ouput = Times @@ Select[factors, ! FreeQ[#, var] &];
  {ouput, Simplify[expr/ouput]}
]


(* ::Chapter:: *)
(*Latex*)


LatexRules = {t1 -> Subscript[t, 1], t2 -> Subscript[t, 2], t3 -> Subscript[t, 3], \[Omega]1 -> Subscript[\[Omega], 1],
\[Omega]2 -> Subscript[\[Omega], 2], J[x_] Coth[\[Beta] x/2]-> f[x], \[Omega][2] -> Subscript[\[Omega], 2], \[Omega][1] -> Subscript[\[Omega], 1],
p1 -> Subscript[a, 1], p3 -> Subscript[a, 3]}

Latex[expr_]:= Module[{latexExpr, TransRule},
TransRule = expr //. LatexRules;
latexExpr = TeXForm[TransRule];
CopyToClipboard[latexExpr];
TransRule // MatrixForm
]

Latex2[expr_]:= Module[{latexExpr, TransRule},
TransRule = expr //. LatexRules;
latexExpr = TeXForm[TransRule];
CopyToClipboard[latexExpr];
latexExpr
]


(*To use the following code in mathematica, select the output text, right click, copy as, plain text.*)

GenerateLaTeXMatrixElements[F_] := Module[{N, equations, elemLaTeX, result},
  
  N = Length[F]; (* Assuming F is a square matrix *)
  
  equations = Flatten[
    Table[
      elemLaTeX = ToString[Latex2[F[[i + 1, j + 1]]]];
      "\\begin{align}\n" <>
      "F_{" <> ToString[i] <> ToString[j] <> "}^{(2)} = " <> elemLaTeX <> ",\n" <>
      "\\end{align}",
      {i, 0, N - 1}, {j, 0, N - 1}
    ]
  ];
  
  result = StringRiffle[equations, "\n"];
  CopyToClipboard[result];
  result
];


GenerateLaTeXFunctionValues[func_, Factor_, N_Integer, fname_ : "f"] := 
 Module[{equations, elemLaTeX, result},
  equations = Table[
    elemLaTeX = ToString[Latex2[func[i]]];
    (* Trim leading spaces and check if the first character is not '-' *)
    If[!StringStartsQ[StringTrim[elemLaTeX], "-"],
     elemLaTeX = "+" <> elemLaTeX
    ];
    If[i == 1,
     fname <> " &= " <> ToString[Latex2[Factor]] <> "\\bigg\\{" <> elemLaTeX <> " \\notag\\\\",
     "&\\quad \\quad " <> elemLaTeX <> If[i == N, "\\bigg\\}", " \\\\"]],
    {i, 1, N}
  ];
  result = "\\begin{align}\n" <> StringRiffle[equations, "\n"] <> "\n\\end{align}";
  CopyToClipboard[result];
  result
];
