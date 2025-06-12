(* ::Package:: *)

(* ::Title:: *)
(*TCL4TripleTimeIntegral*)


(* ::Chapter:: *)
(*Integrand*)


Get["TCLIntegrandCalcs.wl"]


SimpIntegrand = Simplify[K[4]];


(* ::Chapter:: *)
(*Triple Time integral*)


(* Initialize TripleIntR as a 4x4 table filled with Null *)
TripleIntR = Table[Null, {i, 4}, {j, 4}];

(* Populate TripleIntR for i = 3 to 4 and j = 1 to 4 using ParallelTable *)
TripleIntR[[3 ;; 4, 1 ;; 4]] = ParallelTable[
  CleanInt[SimpIntegrand[[i, j]]],
  {i, 3, 4},
  {j, 1, 4}
];


TripleIntR2 = TripleIntR //. Null -> 0;


Dimensions[TripleIntR2]


TripleIntR3 = ParallelMap[ExprNiceForm, TripleIntR2,{2}]


Get["TCL4TripleInt.mx"]


Simplify[TripleIntR3 - TCL4TripleInt]


TCL4TripleIntV2 = TripleIntR3;


DumpVar[TCL4TripleIntV2]
