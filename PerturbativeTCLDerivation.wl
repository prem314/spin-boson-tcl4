(* ::Package:: *)

(* ::Title:: *)
(*TCL2 and TCL4 generator from operator equations*)


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

PartialTraceRule = {TrB[TensorProduct[A_, B_]] -> TensorProduct[A, TrB[B]], TrB[Rho] -> RhoS}

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


(* 1. Define the full expression *)
expr = I (RhoS . A[t1] . A[t] - A[t] . RhoS . A[t1] - A[t] . A[t1] . RhoS + 
       A[t1] . RhoS . A[t]) \[Eta][t - t1] - (RhoS . A[t1] . A[t] - 
       A[t] . RhoS . A[t1] + A[t] . A[t1] . RhoS - 
       A[t1] . RhoS . A[t]) \[Nu][t - t1];

(* 2. Robust Grouping Function *)
(* This handles scalars (I, -1, etc.) correctly *)
RobustGroup[exp_, op_] := Module[{terms, leftGroup, rightGroup},
  terms = List @@ Expand[exp];
  
  (* Helper to check if a term starts or ends with op *)
  (* pattern: optional scalar (c_.) times Dot product *)
  
  leftGroup = Total @ Cases[terms, c_. * Dot[op, rest__] :> c * Dot[rest]];
  rightGroup = Total @ Cases[terms, c_. * Dot[rest__, op] :> c * Dot[rest]];
  
  op . leftGroup + rightGroup . op
]

(* 3. Separate Eta and Nu parts *)
coeffEta = Coefficient[expr, \[Eta][t - t1]];
coeffNu = Coefficient[expr, \[Nu][t - t1]];

(* 4. Apply the robust grouping *)
groupedEta = RobustGroup[coeffEta, A[t]];
groupedNu = RobustGroup[coeffNu, A[t]];

(* 5. Display Result *)
FinalResult = groupedEta * \[Eta][t - t1] + groupedNu * \[Nu][t - t1]


(* ::Text:: *)
(*I just realize that I will be better off doing this calcs In sympy.*)
