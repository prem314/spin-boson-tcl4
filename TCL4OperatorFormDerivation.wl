(* ::Package:: *)

(* ::Title:: *)
(*TCL4OperatorFormDerivation*)


(* ::Text:: *)
(*This notebook derives eqn 9.61 and 9.62 of BR02 starting from eqn 9.45 and 9.51.*)


(* This allows signs and scalars to "leak out" of your Commutators automatically *)
Comm[a_ + b_, c_] := Comm[a, c] + Comm[b, c];
Comm[k_ a_, b_] := k Comm[a, b] /; NumberQ[k] || FreeQ[k, _Dot];

(* Linearity: Distribute over addition *)
AntiComm[a_ + b_, c_] := AntiComm[a, c] + AntiComm[b, c];
AntiComm[c_, a_ + b_] := AntiComm[c, a] + AntiComm[c, b];

(* Scalar Homogeneity: Pull out constants/scalars *)
AntiComm[k_ a_, b_] := k AntiComm[a, b] /; NumberQ[k] || FreeQ[k, _Dot | Comm | AntiComm];
AntiComm[a_, k_ b_] := k AntiComm[a, b] /; NumberQ[k] || FreeQ[k, _Dot | Comm | AntiComm];


(* ::Chapter:: *)
(*TCL2*)


(* Your Definitions *)
Commutator[A_, B_] := A . B - B . A
P[Rho_] := TensorProduct[TrB[Rho], RhoB]
Q[Rho_] := Rho - P[Rho]
HamI[t_] := TensorProduct[A[t], B[t]]
L[t_, Rho_] := -I Commutator[HamI[t], Rho]

(* Note: I added t1 as an argument or global symbol for this to work *)
K2[t_, Rho_] := TrB[P[L[t, L[t1, P[Rho]]]]]

(* --- THE FIX --- *)

(* Define a comprehensive set of simplification rules *)
TensorSimplification = {
   (* 1. Distribute Dot over Plus/Subtract (Left and Right) *)
   Dot[x_, y_ + z_] :> Dot[x, y] + Dot[x, z],
   Dot[x_ + y_, z_] :> Dot[x, z] + Dot[y, z],

   (* 2. Pull numeric scalars (like -I) out of Dot *)
   Dot[x_, s_?NumericQ * y_] :> s * Dot[x, y],
   Dot[s_?NumericQ * x_, y_] :> s * Dot[x, y],
   
   (* 3. Handle simple negation which Mathematica treats differently than -1* *)
   Dot[x_, -y_] :> -Dot[x, y],
   Dot[-x_, y_] :> -Dot[x, y],

   (* 4. Your Tensor Product Rule *)
   TensorProduct[A_, B_] . TensorProduct[C_, D_] :> TensorProduct[A . C, B . D],
   
   (* 5. (Optional) Linearity of TrB to clean up the final result *)
   TrB[x_ + y_] :> TrB[x] + TrB[y],
   TrB[s_?NumericQ * x_] :> s * TrB[x]
};


(* Apply the rules repeatedly *)
result = K2[t, Rho] //. TensorSimplification

PartialTraceRule = {TrB[TensorProduct[A_, B_]] -> TensorProduct[A, TrB[B]], TrB[Rho] -> RhoS};

r2 = result //. PartialTraceRule

r3 = r2 //. TensorProduct[X_, TrB[RhoB]] -> X


(* Define the Cyclic Rule anchored to RhoB *)
CyclicTraceRule = {
    (* 1. Ensure linearity is maintained (distribute TrB over Plus) *)
    TrB[x_ + y_] :> TrB[x] + TrB[y],
    TrB[s_?NumericQ * x_] :> s * TrB[x],

    (* 2. The Cyclic Shift Rule *)
    (* Matches: [Anything Before] . RhoB . [Anything After] *)
    (* Replaces with: [Anything After] . [Anything Before] . RhoB *)
    TrB[pre___ . RhoB . post__] :> TrB[post . pre . RhoB]
};

(* Apply the rule *)
r4 = r3 //. CyclicTraceRule

r5 = r4 //. TensorProduct[X_,TrB[B[t_] . B[t1_] . RhoB]] -> (\[Nu][t-t1] + I \[Eta][t-t1]) X

r6 = r5 //. {\[Eta][-t+t1_] -> -\[Eta][t-t1], \[Nu][-t+t1_] -> \[Nu][t-t1]}

r7 = Simplify[Collect[r6, {\[Nu][t-t1], \[Eta][t-t1]}]]


e8 = TensorExpand[r7]


e9 = e8 //. {
   c1_. A_ . B_ + c2_. B_ . A_ :> 
     Which[
       c1 == -c2, c1 Comm[A, B],
       c1 == c2,  c1 AntiComm[A, B] ,
       True, c1 A . B + c2 B . A 
     ]
}


e10 = e9 //. Comm[A[t], X_] :> - Comm[X,A[t]]


e11 = e10 //. {
   (* Case 1: Coefficients are equal -> AntiCommutator *)
   c1_. Comm[A_ . B_, C_] + c2_. Comm[B_ . A_, C_] /; (c1 == c2) :> 
    c1 Comm[AntiComm[A, B], C],
   
   (* Case 2: Coefficients are opposites -> Commutator *)
   c1_. Comm[A_ . B_, C_] + c2_. Comm[B_ . A_, C_] /; (c1 == -c2) :> 
    c1 Comm[Comm[A, B], C]
}


(* ::Chapter:: *)
(*TCL4*)


K4[t_, t1_, t2_, t3_, Rho_] := 
  TrB[P[L[t, L[t1, L[t2, L[t3, P[Rho]]]]]] - 
  P[L[t, L[t1, P[L[t2, L[t3, P[Rho]]]]]]] - 
  P[L[t, L[t2, P[L[t1, L[t3, P[Rho]]]]]]] - 
  P[L[t, L[t3, P[L[t1, L[t2, P[Rho]]]]]]]]


(* ::Text:: *)
(*Note that there is no clear time ordering in the expression above. This will also reflect in the final result.*)


(* Apply the rules repeatedly *)
result4 = K4[t[0], t[1], t[2], t[3], Rho] //. TensorSimplification;


PartialTraceRule


(* ::Text:: *)
(*The double tensor product below is causing me trouble.*)


r24 = result4 //. PartialTraceRule;

r24 = TensorExpand[r24];


r34 = r24 //. TensorProduct[X_, TrB[RhoB]] -> X;

r34 = TensorExpand[r34];


(* ::Text:: *)
(*Following moves RhoB to the right most part within the TrB*)


r44 = r34 //. CyclicTraceRule;

r44 = TensorExpand[r44];


(* ::Text:: *)
(*Wicks theorem needs to be applied here!*)


w1 = r44 //. TrB[B[t_] . B[t1_] . B[t2_] . B[t3_] . RhoB] :> 
TrB[B[t] . B[t1] . RhoB] TrB[B[t2] . B[t3] . RhoB] + 
TrB[B[t] . B[t2] . RhoB] TrB[B[t1] . B[t3] . RhoB] + 
TrB[B[t] . B[t3] . RhoB] TrB[B[t1] . B[t2] . RhoB];


w2 = TensorExpand[w1];


w3 = w2 //. TrB[B[t[n_]] . B[t[m_]] . RhoB] :> (\[Nu][t[n]-t[m]] + I \[Eta][t[n]-t[m]]);


r54 = w3 //. TensorProduct[X_, Y_] :> X Y;


FreeQ[r54, TrB]


VariableOrderRule = {\[Eta][t[n_]-t[m_]] :> -\[Eta][t[m]-t[n]] /; m < n , \[Nu][t[n_]-t[m_]] :> \[Nu][t[m]-t[n]] /; m < n};
expr = ((A[t[2]] . RhoS . A[t[1]] (I \[Eta][t[1]-t[2]]+\[Nu][t[1]-t[2]])) . A[t[3]] . A[t[0]]-(RhoS . A[t[2]] . A[t[1]] (I \[Eta][-t[1]+t[2]]+\[Nu][-t[1]+t[2]])) . A[t[3]] . A[t[0]]);
expr //. VariableOrderRule


r64 = r54 //. VariableOrderRule;


(* ::Text:: *)
(*I need to separate the nu and eta terms. Find a way to do it.*)


expr = A[t[0]] . (RhoS . A[t[3]] . A[t[2]] (-I \[Eta][t[2] - t[3]] + \[Nu][t[2] - t[3]])) . A[t[1]] (-I \[Eta][t[0] - t[1]] + \[Nu][t[0] - t[1]]) - A[t[0]] . (A[t[2]] . RhoS . A[t[3]] (-I \[Eta][t[2] - t[3]] + \[Nu][t[2] - t[3]])) . A[t[1]] (-I \[Eta][t[0] - t[1]] + \[Nu][t[0] - t[1]]);

expr //. Dot[a___, s_ * b_, c___] /; ! FreeQ[s, \[Eta] | \[Nu]] :> s * Dot[a, b, c]


r74 = r64 //. Dot[a___, s_ * b_, c___] /; ! FreeQ[s, \[Eta] | \[Nu]] :> s * Dot[a, b, c];


Wiki[PullingScalarOutOfOrderedOperators]


r84 = Collect[r74, {\[Nu][t[n_]-t[m_]], \[Eta][t[n_]-t[m_]]}];


(* ::Text:: *)
(*What is the best strategy to break these terms in terms of commutators and anticommutators?*)


r84


e94 = r84 //. {
   c1_. A_ . B_ + c2_. B_ . A_ :> 
     Which[
       c1 == -c2, c1 Comm[A, B],
       c1 == c2,  c1 Anti[A, B] ,
       True, c1 A . B + c2 B . A 
     ]
}


OperatorReorder[expr_,n_] := expr //. {Comm[A[t[n]], X__ ] :> - Comm[X, A[t[n]]], Anti[A[t[n]], X__] :> Anti[X, A[t[n]]]}

e104 = OperatorReorder[e94, 0]


FreeQ[e104, Comm[A[t[0]],A[t[1]] . RhoS . A[t[3]] . A[t[2]]]]
FreeQ[e94, Comm[A[t[0]],A[t[1]] . RhoS . A[t[3]] . A[t[2]]]]


CommAntiRule = {
   (* Case 1: Coefficients are equal -> Antiutator *)
   c1_. X_[A_ . B_, C_] + c2_. X_[B_ . A_, C_] /; (c1 == c2) :> 
    c1 X[Anti[A, B], C],
   
   (* Case 2: Coefficients are opposites -> Commutator *)
   c1_. X_[A_ . B_, C_] + c2_. X_[B_ . A_, C_] /; (c1 == -c2) :> 
    c1 X[Comm[A, B], C]
    
};

RearrangeCommAnti[expr_]:= expr //. CommAntiRule;

e114 = RearrangeCommAnti[e104]


CommAntiRule = {
   (* Case 1: Coefficients are equal -> Antiutator *)
   c1_. X_[A_ . B_, C_] + c2_. X_[B_ . A_, C_] /; (c1 == c2) :> 
    c1 X[Anti[A, B], C],
   
   (* Case 2: Coefficients are opposites -> Commutator *)
   c1_. X_[A_ . B_, C_] + c2_. X_[B_ . A_, C_] /; (c1 == -c2) :> 
    c1 X[Comm[A, B], C],
    (* Case 1: Coefficients are equal -> Antiutator *)
   c1_. X_[C_, A_ . B_] + c2_. X_[C_, B_ . A_] /; (c1 == c2) :> 
    c1 X[C, Anti[A, B]],
   
   (* Case 2: Coefficients are opposites -> Commutator *)
   c1_. X_[C_, A_ . B_] + c2_. X_[C_, B_ . A_] /; (c1 == -c2) :> 
    c1 X[C, Comm[A, B]]
};

RearrangeCommAnti[expr_]:= expr //. CommAntiRule;

e114 = RearrangeCommAnti[e104]


expandedTarget


(* ::Chapter:: *)
(*Verification*)


term1=\[Nu][t-t2]  \[Nu][t1 - t3]comm[f0,comm[ comm[f1, f2],comm[f3, \[Rho]tS]]];
term2=I \[Nu][t-t2] \[Eta][t1 - t3]comm[f0,comm[ comm[f1, f2],anticomm[f3, \[Rho]tS]]];
term3=I \[Eta][t-t2] \[Nu][t1 - t3]comm[f0,anticomm[ comm[f1, f2],comm[f3, \[Rho]tS]]];
term4 = - \[Eta][t-t2] \[Eta][t1 - t3]comm[f0,anticomm[ comm[f1, f2],anticomm[f3, \[Rho]tS]]];
term5=\[Nu][t - t3]\[Nu][t1-t2]comm[f0,comm[comm[f1, f3], comm[f2, \[Rho]tS]]];
term6=I \[Nu][t - t3] \[Eta][t1-t2] comm[f0,comm[comm[f1, f3],anticomm[f2, \[Rho]tS]]];
term7=I \[Eta][t - t3]\[Nu][t1-t2]comm[f0,anticomm[comm[f1, f3],comm[f2, \[Rho]tS]]];
term8=- \[Eta][t - t3]\[Eta][t1-t2]comm[f0,anticomm[comm[f1, f3],anticomm[f2, \[Rho]tS]]];
term9=\[Nu][t-t3]\[Nu][t1-t2] comm[f0,comm[f1,comm[comm[f2,f3],\[Rho]tS]]];
term10=- \[Eta][t-t3]\[Eta][t1-t2] comm[f0,comm[f1,comm[comm[f2,f3],\[Rho]tS]]];
term11=I \[Nu][t-t3]\[Eta][t1-t2] comm[f0,comm[f1,anticomm[comm[f2,f3],\[Rho]tS]]];
term12=I \[Eta][t-t3]\[Nu][t1-t2] comm[f0,comm[f1,anticomm[comm[f2,f3],\[Rho]tS]]];


dot\[Rho]TCL4 = term1+term2+term3+term4+term5+term6+term7+term8+term9+term10+term11+term12;


comm[a_, b_] := a . b - b . a
anticomm[a_, b_] := a . b + b . a


(* ::Text:: *)
(*Matches!*)


term1


OperatorMapping = {
  f0 -> A[t[0]], f1 -> A[t[1]], f2 -> A[t[2]], f3 -> A[t[3]],
  \[Rho]tS -> RhoS
};


TimeMapping = {t -> t[0], t1 -> t[1], t2 -> t[2], t3 -> t[3]}


TCLReplace[expr_] := (expr /. TimeMapping) /. OperatorMapping


storedRes = TCLReplace[dot\[Rho]TCL4]


Difference = r84 - storedRes;


Simplify[TensorExpand[r84 - storedRes]]


(* ::Text:: *)
(*Not putting TensorExpand above does not give the result.*)


dot\[Rho]TCL2 = -( \[Nu][t-t1] * ( comm[f0,comm[f1,\[Rho]tS]] ) + I*\[Eta][t-t1] * (comm[f0,anticomm[f1,\[Rho]tS]]) );


StoredResTCL2 = TCLReplace[dot\[Rho]TCL2]


TensorExpand[TCLReplace[e8] - StoredResTCL2]


(* ::Text:: *)
(*Hence, the results match as they are used in*)


Wiki[TCLIntegrandCalcs]
