(* ::Package:: *)

(* ::Title:: *)
(*MathFuncMathematica*)


(*
(* Define the first function using Limit *)
RsiduO1Limit[expr_, z_, c_,n_] := 1/Factorial[n-1] Limit[D[Simplify[(z - c)^n expr],{z,n-1}], z -> c]

(* Define the second function using direct substitution.
This one comes handy when you have symbolic function and
mathematica is afraid to set the limit for the argument
of the function through direct substitution. *)
RsiduO1Direct[expr_, z_, c_,n_] := 1/Factorial[n-1] D[Simplify[(z - c)^n expr],{z,n-1}] /. z -> c

(* Define the combined function *)
RsiduO1[expr_, z_, c_, n_:1] := Module[{result},
  
  (* Attempt to evaluate using RsiduO1 *)
  result = Quiet[RsiduO1Limit[expr, z, c, n]];
  
  (* Check if the result still contains Limit, indicating it wasn't evaluated *)
  If[FreeQ[result, Limit],
    (* If Limit was successful, return the result *)
    result,
    
    (* If Limit wasn't successful, use RsiduO1Direct and print a warning *)
    (
      Print["Warning: Limit did not evaluate. Using direct evaluation."];
      RsiduO1Direct[expr, z, c , n]
    )
  ]
]

Rsidu[expr_, z_, PoleList_, n_: 1] := 
  Total[ParallelMap[RsiduO1[expr, z, #, n] &, PoleList]]
*)

Rsidu[expr_, z_, PoleList_, assumptions_: {}] := Module[{UniquePoleList},
    UniquePoleList = DeleteDuplicates[PoleList];
    Total[ParallelMap[
      Residue[expr, {z, #}, Assumptions-> assumptions] &,
      UniquePoleList
    ]]
]


  
O2Res[expr_, z_, PoleList_] := Total[ParallelMap[Residue[(z-#) expr, {z, #}] &, DeleteDuplicates[PoleList]]]


(* ::Text:: *)
(*https://en.wikipedia.org/wiki/Residue_(complex_analysis)#Limit_formula_for_higher-order_poles*)


Wiki[Contour]


ToScientificForm[number_] := ScientificForm[N[number]]


(* ::Chapter:: *)
(*TrignometricFunctions*)


(* Reduces product of trig functions into trig function*)
TrigFactorCollect = {Cos[a_]Cos[b_]:> 1/2 (Cos[a+b]+Cos[a-b]), Cos[a_]Sin[b_]:> 1/2 (Sin[a+b]-Sin[a-b]),Sin[a_]Sin[b_]:>1/2 (Cos[a-b]-Cos[a+b]), 
Sin[a_]^2 :>1/2 (1 - Cos[2 a]), Cos[a_]^2:> 1/2 (Cos[2 a] + 1)};
