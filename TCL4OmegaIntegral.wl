(* ::Package:: *)

(* ::Title:: *)
(*TCL4OmegaIntegral*)


(* ::Text:: *)
(*Do not delete the numerical and verification code from here. Wiki[TCL4And2DrudeEval] is just for getting the result fast, not for record, verifications and tests. This notebook is there for all of these.*)


Wiki[TCL4And2DrudeEval]
Wiki[GenHTCL4GenCalc]


Get["TCL4SpinBosonFunctions.wl"]
Get["DrudeCutoffFunctions.wl"]
Get["MyFunctions.wl"]


Get["TCL4TripleInt.mx"]
TripleIntR3 = (TCL4TripleInt/4);


TCL4TripleInt[[2,1]]


(* ::Text:: *)
(*At the time of doing these integrals, I did not calculate the [[2,i]] kind of terms.*)


(* ::Chapter:: *)
(*TCL4 Calc*)


(* ::Text:: *)
(*In the following, I am getting rid of the term without any zero exponential because in the cases it is nonzero, its integral evaluates to zero by symmetry. See derivation in the subsequent sections.*)


(* ::Text:: *)
(*The Divergent term was coming out to be non-zero because I ignored the zero zero term. Including it fixed the issue.*)


(*
ArrayEvalCheckpoint[AllResSeg, TripleIntR3, TCL4Generator ,2];
*)


LoadVar[TCL4Generator]


TCL4Generator[[4,1]]


(*
TCL4Generator2 = FullSimplify[TCL4Generator]

TCL4Generator = TCL4Generator2;

DumpVar[TCL4Generator]
*)
