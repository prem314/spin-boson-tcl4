(* ::Package:: *)

(* ::Title:: *)
(*TCLIntegrandCalcs*)


Get["TCL4SpinBosonFunctions.wl"]
Wiki[TCL4DynamicsLatex]


(* ::Chapter:: *)
(*Initialization*)


For[i=0,i<4,i++,\[Rho]1[i] = Vec\[Rho][i]*PauliMatrix[i]];

(* General system density matrix in the Schrodinger picture *)
\[Rho]tS=ConstantArray[0,{2,2}];
For[i=0,i<4,i++,\[Rho]tS += 1/2 \[Rho]1[i]];  (*Added an extra 1/2 factor so that now we have \[Rho] = 1/2 \[Sum] ri \[Sigma]i. This mean \[Rho]1[0] = 1.*)

\[Rho]tS // MatrixForm;


(* ::Chapter:: *)
(*TCL0*)


dot\[Rho]TCL = -I*comm[ha1,\[Rho]tS]; (*Bare schrodinger equation*)
dot\[Rho]TCL // MatrixForm


K[0] = VecGen[dot\[Rho]TCL]
K[0] // MatrixForm


(* ::Chapter:: *)
(*TCL2*)


f0=DiffIntPic[t, \[Sigma]1];

f1=DiffIntPic[t1, \[Sigma]1];

(* decay2 is K_2(t)\[Rho](t) in schrodinger picture. Eqn 31 of Kappler *)
dot\[Rho]TCL2 = -( \[Nu][t-t1] * ( comm[f0,comm[f1,\[Rho]tS]] ) + I*\[Eta][t-t1] * (comm[f0,anticomm[f1,\[Rho]tS]]) );


K[2] = VecGen[dot\[Rho]TCL2];
ExpToTrig[K[2]] // MatrixForm


(* ::Chapter:: *)
(*TCL4*)


f2=DiffIntPic[t2, \[Sigma]1];

f3=DiffIntPic[t3, \[Sigma]1];
(* Fourth order generator K_4(t) is defined *)

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

(*
I verified that the above matches with A.1 of Kappler.
*)


(* decay4 is K_4(t)\[Rho](t) in the interaction picture *)
dot\[Rho]TCL4 = term1+term2+term3+term4+term5+term6+term7+term8+term9+term10+term11+term12;

(*VecDot\[Rho]tcl4 is the vectorized form of K_4(t)\[Rho](t*)
K[4] = VecGen[dot\[Rho]TCL4];


(* ::Chapter:: *)
(*Verification*)


(*
Simplify[K[0] - TCL0Gen]
Simplify[K[2] - TCL2Integrand]
Simplify[K[4] - TCL4Integrand]
*)


(* ::Text:: *)
(*Matches!*)
