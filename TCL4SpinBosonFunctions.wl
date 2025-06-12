(* ::Package:: *)

(* ::Title:: *)
(*TCL4SpinBosonFunctions*)


TCL4DynamicsFolder = "/home/premkr/Dropbox/work/projects/tcl4_dynamics";


Get["MyFunctions.wl"]
Get["MathFuncMathematica.wl"]


(* ::Text:: *)
(*Based  on*)


Wiki[TheTimeConvolutionlessProjectionOperatorTechniqueInTheQuantumTheoryOfDissipationAndDecoherenceHeinzpeterBreuerBerndKapplerFrancescoPetruccionePaper]
Wiki[TCL4ForSpinBosonDerivation]


(* ::Chapter:: *)
(*Initialization functions*)


ha1=\[CapitalOmega] PauliMatrix[3]/2; (* System free Hamiltonian, eqn 101 *)

U[t_] := MatrixExp[-I ha1 t] (*The unitary matrix*)

\[Sigma]1 = p3 PauliMatrix[3] - p1 PauliMatrix[1]; (*The system side interaction operator, eqn 103*)

DiffIntPic[t1_, op_]:= U[t -t1] . op . U[-(t -t1) ] (*Going to the interaction picture*)

(* Define commutator and anti-commulator *)
comm[a_,b_]:=a . b - b . a;
anticomm[a_,b_]:=a . b + b . a;

VecGen[expr_] := Module[{i, j, VecDotRhoLocal, a},
  
  (* Step 1: Compute VecDotRhoLocal[i] for i = 0 to 3 *)
  VecDotRhoLocal = Table[
    Tr[expr . PauliMatrix[i]],
    {i, 0, 3}
  ];
  
  (* Step 2 and 3: Define and return function a[i, j] *)
  a = Table[
    Coefficient[VecDotRhoLocal[[i + 1]], Vec\[Rho][j]], {i, 0, 3}, {j, 0, 3}
  ];
  
  Return[a];
]


(* ::Chapter:: *)
(*Spectral Density Symmetry*)


AntiSymmetricJRule = {J[-\[Omega]+\[CapitalOmega]] -> - J[\[Omega] - \[CapitalOmega]], J[-\[Omega] - \[CapitalOmega]] -> - J[\[Omega] + \[CapitalOmega]], f[-\[Omega]+\[CapitalOmega]] -> f[\[Omega] - \[CapitalOmega]], f[-\[Omega] - \[CapitalOmega]] -> f[\[Omega] + \[CapitalOmega]]}
P2SinRule = {p1 -> Sin[\[Theta]], p3 -> Cos[\[Theta]]}
Sin2PRule = {Sin[\[Theta]] -> p1, Cos[\[Theta]] -> p3}


(* ::Chapter:: *)
(*General hamiltonian results*)


EtaIntegrand[\[Omega]_,t_] := - J[\[Omega]] Sin[\[Omega] t]
NuIntegrand[\[Omega]_,t_] := f[\[Omega]] Cos[\[Omega] t]


(* ::Section:: *)
(*TCL0 Expressions*)


SpinBosonHamiltonian = \[CapitalOmega]*PauliMatrix[3]/2;
GibbsDiagPop = Tr[PauliMatrix[3] . MatrixExp[- \[Beta] SpinBosonHamiltonian ]/Tr[MatrixExp[- \[Beta] SpinBosonHamiltonian]]]


TCL0Gen = {{0,0,0,0},{0,0,-\[CapitalOmega],0},{0,\[CapitalOmega],0,0},{0,0,0,0}}


StateZero = {{1}, {0}, {0}, {HoldForm[Subscript[v, 3]^(0)]}}


(* ::Section:: *)
(*TCL2*)


Wiki[O2SSGenSpinBoson]


NicerNuEtaRule = {nu[t_] -> \[Nu][t], eta[t_]-> \[Eta][t]};


TCL2Integrand = {{0,0,0,0},{4 p1 p3 Sin[t \[CapitalOmega]-t1 \[CapitalOmega]] \[Eta][t-t1],-4 p3^2 \[Nu][t-t1],0,-4 p1 p3 Cos[t \[CapitalOmega]-t1 \[CapitalOmega]] \[Nu][t-t1]},{4 p1 p3 \[Eta][t-t1]-4 p1 p3 Cos[t \[CapitalOmega]-t1 \[CapitalOmega]] \[Eta][t-t1],4 p1^2 Sin[t \[CapitalOmega]-t1 \[CapitalOmega]] \[Nu][t-t1],-4 p3^2 \[Nu][t-t1]-4 p1^2 Cos[t \[CapitalOmega]-t1 \[CapitalOmega]] \[Nu][t-t1],-4 p1 p3 Sin[t \[CapitalOmega]-t1 \[CapitalOmega]] \[Nu][t-t1]},{4 p1^2 Sin[t \[CapitalOmega]-t1 \[CapitalOmega]] \[Eta][t-t1],-4 p1 p3 \[Nu][t-t1],0,-4 p1^2 Cos[t \[CapitalOmega]-t1 \[CapitalOmega]] \[Nu][t-t1]}}


TCL2Gen = {{0,0,0,0},{-2 p1 p3 \[Pi] \[Lambda]^2 J[\[CapitalOmega]],-2 p3^2 \[Pi] \[Lambda]^2 f[0],0,-2 p1 p3 \[Pi] \[Lambda]^2 f[\[CapitalOmega]]},{-((4 p1 p3 \[Lambda]^2 \[CapitalOmega]^2 J[\[Omega]])/(\[Omega] (-\[Omega]^2+\[CapitalOmega]^2))),(4 p1^2 \[Lambda]^2 \[CapitalOmega] f[\[Omega]])/(-\[Omega]^2+\[CapitalOmega]^2),\[Lambda]^2 (-2 p3^2 \[Pi] f[0]-2 p1^2 \[Pi] f[\[CapitalOmega]]),-((4 p1 p3 \[Lambda]^2 \[CapitalOmega] f[\[Omega]])/(-\[Omega]^2+\[CapitalOmega]^2))},{-2 p1^2 \[Pi] \[Lambda]^2 J[\[CapitalOmega]],-2 p1 p3 \[Pi] \[Lambda]^2 f[0],0,-2 p1^2 \[Pi] \[Lambda]^2 f[\[CapitalOmega]]}};
TCL2Gen // MatrixForm


SigmaZO2 = 1/((1+E^(\[Beta] \[CapitalOmega])) (\[Omega]^3-\[Omega] \[CapitalOmega]^2)^2 J[\[CapitalOmega]]) 2 E^((\[Beta] \[CapitalOmega])/2) Tanh[(\[Beta] \[CapitalOmega])/2] (-4 p1^2 \[Omega]^2 \[CapitalOmega]^2 Cosh[(\[Beta] \[CapitalOmega])/2] Coth[(\[Beta] \[Omega])/2] J[\[Omega]] J[\[CapitalOmega]]-4 p1^2 \[Omega]^3 \[CapitalOmega] Cosh[(\[Beta] \[CapitalOmega])/2] Coth[(\[Beta] \[CapitalOmega])/2] J[\[Omega]] J[\[CapitalOmega]]+Cosh[(\[Beta] \[CapitalOmega])/2] Coth[(\[Beta] \[Omega])/2] J[\[Omega]] (2 p1^2 \[Omega]^4 J[\[CapitalOmega]]+6 p1^2 \[Omega]^2 \[CapitalOmega]^2 J[\[CapitalOmega]]+2 p1^2 \[Omega]^2 (\[Omega]-\[CapitalOmega]) \[CapitalOmega] (\[Omega]+\[CapitalOmega]) Derivative[1][J][\[CapitalOmega]])-2 p1^2 \[Omega]^2 (\[Omega]-\[CapitalOmega]) \[CapitalOmega] (\[Omega]+\[CapitalOmega]) Coth[(\[Beta] \[Omega])/2] J[\[Omega]] Sinh[(\[Beta] \[CapitalOmega])/2] (-(1/2) \[Beta] Csch[(\[Beta] \[CapitalOmega])/2]^2 J[\[CapitalOmega]]+Coth[(\[Beta] \[CapitalOmega])/2] Derivative[1][J][\[CapitalOmega]]))

SigmaXO2 = -((4 p1 p3 ((1+E^(\[Beta] \[CapitalOmega])) \[CapitalOmega]-(-1+E^(\[Beta] \[CapitalOmega])) \[Omega] Coth[(\[Beta] \[Omega])/2]) J[\[Omega]])/((1+E^(\[Beta] \[CapitalOmega])) \[Omega] (\[Omega]^2-\[CapitalOmega]^2)))


(* ::Section:: *)
(*TCL4*)


TCL4Integrand = {{0,0,0,0},{16 p1 p3 ((p3^2 (-Sin[(t-t1) \[CapitalOmega]]+Sin[(t-t2) \[CapitalOmega]])+p1^2 Cos[(t-t3) \[CapitalOmega]] Sin[(t1-t2) \[CapitalOmega]]) \[Eta][t1-t3] \[Nu][t-t2]+(p3^2 (Sin[(t-t2) \[CapitalOmega]]-Sin[(t-t3) \[CapitalOmega]])-p1^2 Cos[(t-t1) \[CapitalOmega]] Sin[(t2-t3) \[CapitalOmega]]) \[Eta][t-t3] \[Nu][t1-t2]+1/2 (2 p3^2 (-Sin[(t-t1) \[CapitalOmega]]+Sin[(t-t2) \[CapitalOmega]])+p1^2 (Sin[(t+t1-t2-t3) \[CapitalOmega]]-Sin[(t-t1+t2-t3) \[CapitalOmega]])) \[Eta][t1-t2] \[Nu][t-t3]),16 p1^2 p3^2 (Cos[(2 t-t1-t2) \[CapitalOmega]] \[Nu][t1-t2] \[Nu][t-t3]+Cos[(t-t3) \[CapitalOmega]] ((Cos[(t-t1) \[CapitalOmega]]-Cos[(t-t2) \[CapitalOmega]]) (\[Eta][t1-t2] \[Eta][t-t3]+\[Eta][t-t2] \[Eta][t1-t3])-Cos[(t-t1) \[CapitalOmega]] \[Nu][t1-t2] \[Nu][t-t3])+Sin[(t-t3) \[CapitalOmega]] (Sin[(t-t2) \[CapitalOmega]] \[Nu][t1-t2] \[Nu][t-t3]+(-Sin[(t-t1) \[CapitalOmega]]+Sin[(t-t2) \[CapitalOmega]]) \[Nu][t-t2] \[Nu][t1-t3])),8 p1^2 p3^2 ((Sin[(2 t-t1-t3) \[CapitalOmega]]+Sin[(t1-t3) \[CapitalOmega]]-Sin[(2 t-t2-t3) \[CapitalOmega]]-Sin[(t2-t3) \[CapitalOmega]]) \[Eta][t1-t2] \[Eta][t-t3]+(Sin[(2 t-t1-t3) \[CapitalOmega]]+Sin[(t1-t3) \[CapitalOmega]]-Sin[(2 t-t2-t3) \[CapitalOmega]]-Sin[(t2-t3) \[CapitalOmega]]) \[Eta][t-t2] \[Eta][t1-t3]+(2 Sin[(2 t-t1-t2) \[CapitalOmega]]-Sin[(2 t-t1-t3) \[CapitalOmega]]+Sin[(t1-t3) \[CapitalOmega]]-Sin[(2 t-t2-t3) \[CapitalOmega]]+Sin[(t2-t3) \[CapitalOmega]]) \[Nu][t1-t2] \[Nu][t-t3]+(2 Sin[(t1-t2) \[CapitalOmega]]+Sin[(2 t-t1-t3) \[CapitalOmega]]-Sin[(t1-t3) \[CapitalOmega]]-Sin[(2 t-t2-t3) \[CapitalOmega]]+Sin[(t2-t3) \[CapitalOmega]]) \[Nu][t-t2] \[Nu][t1-t3]),16 p1 p3 (p3^2 (-Cos[(t-t1) \[CapitalOmega]]+Cos[(t-t2) \[CapitalOmega]]) \[Eta][t1-t2] \[Eta][t-t3]+p3^2 (-Cos[(t-t1) \[CapitalOmega]]+Cos[(t-t2) \[CapitalOmega]]) \[Eta][t-t2] \[Eta][t1-t3]+(p3^2 (-Cos[(t-t2) \[CapitalOmega]]+Cos[(t-t3) \[CapitalOmega]])+p1^2 Sin[(t-t2) \[CapitalOmega]] Sin[(t1-t3) \[CapitalOmega]]) \[Nu][t1-t2] \[Nu][t-t3]+p1^2 Sin[(t1-t2) \[CapitalOmega]] Sin[(t-t3) \[CapitalOmega]] \[Nu][t-t2] \[Nu][t1-t3])},{8 p1 p3 ((2 p3^2 (Cos[(t-t1) \[CapitalOmega]]-Cos[(t-t2) \[CapitalOmega]])+p1^2 (2 Cos[(t1-t3) \[CapitalOmega]]-Cos[(t+t1-t2-t3) \[CapitalOmega]]-2 Cos[(t2-t3) \[CapitalOmega]]+Cos[(t-t1+t2-t3) \[CapitalOmega]])) \[Eta][t1-t3] \[Nu][t-t2]+(2 p3^2 (-Cos[(t-t2) \[CapitalOmega]]+Cos[(t-t3) \[CapitalOmega]])+p1^2 (-2 Cos[(t1-t2) \[CapitalOmega]]+2 Cos[(t1-t3) \[CapitalOmega]]+Cos[(t-t1+t2-t3) \[CapitalOmega]]-Cos[(t-t1-t2+t3) \[CapitalOmega]])) \[Eta][t-t3] \[Nu][t1-t2]+(2 p3^2 (Cos[(t-t1) \[CapitalOmega]]-Cos[(t-t2) \[CapitalOmega]])+p1^2 (2 Cos[(t1-t3) \[CapitalOmega]]-Cos[(t+t1-t2-t3) \[CapitalOmega]]-2 Cos[(t2-t3) \[CapitalOmega]]+Cos[(t-t1+t2-t3) \[CapitalOmega]])) \[Eta][t1-t2] \[Nu][t-t3]),16 p1^2 ((p3^2 Cos[(t-t3) \[CapitalOmega]] Sin[(t-t1) \[CapitalOmega]]-Cos[(t-t2) \[CapitalOmega]] (p3^2 Sin[(t-t3) \[CapitalOmega]]+p1^2 Sin[(t1-t3) \[CapitalOmega]])+(p3^2+p1^2 Cos[(t-t1) \[CapitalOmega]]) Sin[(t2-t3) \[CapitalOmega]]) \[Eta][t1-t2] \[Eta][t-t3]+Cos[(t-t3) \[CapitalOmega]] (p3^2 (Sin[(t-t1) \[CapitalOmega]]-Sin[(t-t2) \[CapitalOmega]])-p1^2 Sin[(t1-t2) \[CapitalOmega]]) \[Eta][t-t2] \[Eta][t1-t3]-1/2 (p3^2 (2 Sin[(t-t1) \[CapitalOmega]]-2 Sin[(2 t-t1-t2) \[CapitalOmega]]-2 Sin[(t-t3) \[CapitalOmega]]+Sin[(2 t-t1-t3) \[CapitalOmega]]+Sin[(t1-t3) \[CapitalOmega]]+Sin[(2 t-t2-t3) \[CapitalOmega]]+Sin[(t2-t3) \[CapitalOmega]])+p1^2 (Sin[(t-t1+t2-t3) \[CapitalOmega]]-Sin[(t-t1-t2+t3) \[CapitalOmega]])) \[Nu][t1-t2] \[Nu][t-t3]-p3^2 (Sin[(t-t1) \[CapitalOmega]]-Sin[(t-t2) \[CapitalOmega]]+Sin[(t1-t2) \[CapitalOmega]]+(-Cos[(t-t1) \[CapitalOmega]]+Cos[(t-t2) \[CapitalOmega]]) Sin[(t-t3) \[CapitalOmega]]) \[Nu][t-t2] \[Nu][t1-t3]),16 p1^2 ((-Sin[(t-t2) \[CapitalOmega]] (p3^2 Sin[(t-t3) \[CapitalOmega]]+p1^2 Sin[(t1-t3) \[CapitalOmega]])+Sin[(t-t1) \[CapitalOmega]] (p3^2 Sin[(t-t3) \[CapitalOmega]]+p1^2 Sin[(t2-t3) \[CapitalOmega]])) \[Eta][t1-t2] \[Eta][t-t3]+(p3^2 (Sin[(t-t1) \[CapitalOmega]]-Sin[(t-t2) \[CapitalOmega]])-p1^2 Sin[(t1-t2) \[CapitalOmega]]) Sin[(t-t3) \[CapitalOmega]] \[Eta][t-t2] \[Eta][t1-t3]+1/2 (p3^2 (2 Cos[(t-t1) \[CapitalOmega]]-2 Cos[(2 t-t1-t2) \[CapitalOmega]]-2 Cos[(t-t3) \[CapitalOmega]]+Cos[(2 t-t1-t3) \[CapitalOmega]]+Cos[(2 t-t2-t3) \[CapitalOmega]]-Cos[(-t1+t3) \[CapitalOmega]]+Cos[(-t2+t3) \[CapitalOmega]])-2 p1^2 Sin[(t-t1) \[CapitalOmega]] Sin[(t2-t3) \[CapitalOmega]]) \[Nu][t1-t2] \[Nu][t-t3]-p3^2 (Cos[(t-t1) \[CapitalOmega]]-Cos[(t-t2) \[CapitalOmega]]) (-1+Cos[(t-t3) \[CapitalOmega]]) \[Nu][t-t2] \[Nu][t1-t3]),8 p1 p3 (2 (p3^2 (-Sin[(t-t1) \[CapitalOmega]]+Sin[(t-t2) \[CapitalOmega]])+p1^2 Sin[(t1-t2) \[CapitalOmega]]) \[Eta][t1-t2] \[Eta][t-t3]+2 (p3^2 (-Sin[(t-t1) \[CapitalOmega]]+Sin[(t-t2) \[CapitalOmega]])+p1^2 Sin[(t1-t2) \[CapitalOmega]]) \[Eta][t-t2] \[Eta][t1-t3]+(2 p3^2 (-Sin[(t-t2) \[CapitalOmega]]+Sin[(t-t3) \[CapitalOmega]])+p1^2 (2 Sin[(t1-t3) \[CapitalOmega]]-Sin[(t+t1-t2-t3) \[CapitalOmega]]+2 Sin[(t2-t3) \[CapitalOmega]]+Sin[(t-t1-t2+t3) \[CapitalOmega]])) \[Nu][t1-t2] \[Nu][t-t3]+p1^2 (2 Sin[(t1-t3) \[CapitalOmega]]-Sin[(t+t1-t2-t3) \[CapitalOmega]]-2 Sin[(t2-t3) \[CapitalOmega]]+Sin[(t-t1+t2-t3) \[CapitalOmega]]) \[Nu][t-t2] \[Nu][t1-t3])},{16 p1^2 ((p3^2 (-Sin[(t-t1) \[CapitalOmega]]+Sin[(t-t2) \[CapitalOmega]])+p1^2 Cos[(t-t3) \[CapitalOmega]] Sin[(t1-t2) \[CapitalOmega]]) \[Eta][t1-t3] \[Nu][t-t2]+(p3^2 (Sin[(t-t2) \[CapitalOmega]]-Sin[(t-t3) \[CapitalOmega]])-p1^2 Cos[(t-t1) \[CapitalOmega]] Sin[(t2-t3) \[CapitalOmega]]) \[Eta][t-t3] \[Nu][t1-t2]+1/2 (2 p3^2 (-Sin[(t-t1) \[CapitalOmega]]+Sin[(t-t2) \[CapitalOmega]])+p1^2 (Sin[(t+t1-t2-t3) \[CapitalOmega]]-Sin[(t-t1+t2-t3) \[CapitalOmega]])) \[Eta][t1-t2] \[Nu][t-t3]),16 p1^3 p3 (Cos[(2 t-t1-t2) \[CapitalOmega]] \[Nu][t1-t2] \[Nu][t-t3]+Cos[(t-t3) \[CapitalOmega]] ((Cos[(t-t1) \[CapitalOmega]]-Cos[(t-t2) \[CapitalOmega]]) (\[Eta][t1-t2] \[Eta][t-t3]+\[Eta][t-t2] \[Eta][t1-t3])-Cos[(t-t1) \[CapitalOmega]] \[Nu][t1-t2] \[Nu][t-t3])+Sin[(t-t3) \[CapitalOmega]] (Sin[(t-t2) \[CapitalOmega]] \[Nu][t1-t2] \[Nu][t-t3]+(-Sin[(t-t1) \[CapitalOmega]]+Sin[(t-t2) \[CapitalOmega]]) \[Nu][t-t2] \[Nu][t1-t3])),-16 p1^3 p3 ((Cos[(t-t3) \[CapitalOmega]] Sin[(t-t2) \[CapitalOmega]]-Cos[(t-t1) \[CapitalOmega]] Sin[(t-t3) \[CapitalOmega]]+Sin[(t2-t3) \[CapitalOmega]]) \[Eta][t1-t2] \[Eta][t-t3]+(-Cos[(t-t1) \[CapitalOmega]]+Cos[(t-t2) \[CapitalOmega]]) Sin[(t-t3) \[CapitalOmega]] \[Eta][t-t2] \[Eta][t1-t3]-1/2 (2 Sin[(2 t-t1-t2) \[CapitalOmega]]-Sin[(2 t-t1-t3) \[CapitalOmega]]+Sin[(t1-t3) \[CapitalOmega]]-Sin[(2 t-t2-t3) \[CapitalOmega]]+Sin[(t2-t3) \[CapitalOmega]]) \[Nu][t1-t2] \[Nu][t-t3]-(Cos[(t-t3) \[CapitalOmega]] (Sin[(t-t1) \[CapitalOmega]]-Sin[(t-t2) \[CapitalOmega]])+Sin[(t1-t2) \[CapitalOmega]]) \[Nu][t-t2] \[Nu][t1-t3]),16 p1^2 (p3^2 (-Cos[(t-t1) \[CapitalOmega]]+Cos[(t-t2) \[CapitalOmega]]) \[Eta][t1-t2] \[Eta][t-t3]+p3^2 (-Cos[(t-t1) \[CapitalOmega]]+Cos[(t-t2) \[CapitalOmega]]) \[Eta][t-t2] \[Eta][t1-t3]+(p3^2 (-Cos[(t-t2) \[CapitalOmega]]+Cos[(t-t3) \[CapitalOmega]])+p1^2 Sin[(t-t2) \[CapitalOmega]] Sin[(t1-t3) \[CapitalOmega]]) \[Nu][t1-t2] \[Nu][t-t3]+p1^2 Sin[(t1-t2) \[CapitalOmega]] Sin[(t-t3) \[CapitalOmega]] \[Nu][t-t2] \[Nu][t1-t3])}};


(* ::Text:: *)
(*The omega integrals for the full range version of these variables is assumed to run from -Inf to +Inf.*)


TCL4IntegrandFullRange = TCL4Integrand/4;


TCL4GenRule = {e4[4,1] -> 1/(\[Omega]^3-\[Omega] \[CapitalOmega]^2)^2 2 p1^2 \[Pi] (-4 \[Omega] \[CapitalOmega] (p3^2 (\[Omega]-\[CapitalOmega]) (\[Omega]+\[CapitalOmega]) f[0]+p1^2 \[Omega]^2 f[\[CapitalOmega]]) J[\[Omega]]+f[\[Omega]] (p3^2 (\[Omega]-\[CapitalOmega])^2 \[CapitalOmega]^2 J[-\[Omega]-\[CapitalOmega]]+p3^2 \[CapitalOmega]^2 (\[Omega]+\[CapitalOmega])^2 J[\[Omega]-\[CapitalOmega]]+2 p1^2 \[Omega]^4 J[\[CapitalOmega]]+4 p3^2 \[Omega]^4 J[\[CapitalOmega]]+6 p1^2 \[Omega]^2 \[CapitalOmega]^2 J[\[CapitalOmega]]-8 p3^2 \[Omega]^2 \[CapitalOmega]^2 J[\[CapitalOmega]]+4 p3^2 \[CapitalOmega]^4 J[\[CapitalOmega]]-p3^2 \[Omega]^2 \[CapitalOmega]^2 J[-\[Omega]+\[CapitalOmega]]-2 p3^2 \[Omega] \[CapitalOmega]^3 J[-\[Omega]+\[CapitalOmega]]-p3^2 \[CapitalOmega]^4 J[-\[Omega]+\[CapitalOmega]]-p3^2 \[Omega]^2 \[CapitalOmega]^2 J[\[Omega]+\[CapitalOmega]]+2 p3^2 \[Omega] \[CapitalOmega]^3 J[\[Omega]+\[CapitalOmega]]-p3^2 \[CapitalOmega]^4 J[\[Omega]+\[CapitalOmega]]+2 p1^2 \[Omega]^2 (\[Omega]-\[CapitalOmega]) \[CapitalOmega] (\[Omega]+\[CapitalOmega]) Derivative[1][J][\[CapitalOmega]])),
e4[4,4] -> 1/(\[Omega]^3-\[Omega] \[CapitalOmega]^2)^2 2 p1^2 \[Pi] (-p3^2 (\[Omega]-\[CapitalOmega]) \[CapitalOmega] (\[Omega]+\[CapitalOmega]) J[\[Omega]] ((\[Omega]-\[CapitalOmega]) J[-\[Omega]-\[CapitalOmega]]+(\[Omega]+\[CapitalOmega]) (J[\[Omega]-\[CapitalOmega]]-J[-\[Omega]+\[CapitalOmega]])+(-\[Omega]+\[CapitalOmega]) J[\[Omega]+\[CapitalOmega]])-f[\[Omega]] (4 p3^2 \[Omega]^2 (\[Omega]-\[CapitalOmega]) (\[Omega]+\[CapitalOmega]) f[0]-p3^2 \[Omega] (\[Omega]-\[CapitalOmega])^2 \[CapitalOmega] f[-\[Omega]-\[CapitalOmega]]+p3^2 \[Omega] \[CapitalOmega] (\[Omega]+\[CapitalOmega])^2 f[\[Omega]-\[CapitalOmega]]-4 p3^2 \[Omega]^4 f[\[CapitalOmega]]-4 p1^2 \[Omega]^2 \[CapitalOmega]^2 f[\[CapitalOmega]]+8 p3^2 \[Omega]^2 \[CapitalOmega]^2 f[\[CapitalOmega]]-4 p3^2 \[CapitalOmega]^4 f[\[CapitalOmega]]+p3^2 \[Omega]^3 \[CapitalOmega] f[-\[Omega]+\[CapitalOmega]]+2 p3^2 \[Omega]^2 \[CapitalOmega]^2 f[-\[Omega]+\[CapitalOmega]]+p3^2 \[Omega] \[CapitalOmega]^3 f[-\[Omega]+\[CapitalOmega]]-p3^2 \[Omega]^3 \[CapitalOmega] f[\[Omega]+\[CapitalOmega]]+2 p3^2 \[Omega]^2 \[CapitalOmega]^2 f[\[Omega]+\[CapitalOmega]]-p3^2 \[Omega] \[CapitalOmega]^3 f[\[Omega]+\[CapitalOmega]]-2 p1^2 \[Omega]^2 (\[Omega]-\[CapitalOmega]) \[CapitalOmega] (\[Omega]+\[CapitalOmega]) Derivative[1][f][\[CapitalOmega]]))};


F334Calc = e4[4,4] //. TCL4GenRule;


F334 = Simplify[F334Calc //. AntiSymmetricJRule]


(* ::Text:: *)
(*In the following, the Full range version of these variables is *)


F334FullRange = F334/2


F304Calc = e4[4,1] //. TCL4GenRule;


F304 = Simplify[F304Calc //. AntiSymmetricJRule]


F304FullRange = F304/2


(* ::Subsection:: *)
(*Triple Integral saved*)


F334Int = -((I E^(2 I t \[CapitalOmega]-I t \[Omega]1-I t \[Omega]2) p1^4 f[\[Omega]1] f[\[Omega]2])/((\[CapitalOmega]-\[Omega]1)^2 (\[CapitalOmega]-\[Omega]2)))+(I E^(-2 I t \[CapitalOmega]+I t \[Omega]1+I t \[Omega]2) p1^4 f[\[Omega]1] f[\[Omega]2])/((\[CapitalOmega]-\[Omega]1)^2 (\[CapitalOmega]-\[Omega]2))-(I E^(2 I t \[CapitalOmega]+I t \[Omega]1-I t \[Omega]2) p1^4 f[\[Omega]1] f[\[Omega]2])/((\[CapitalOmega]+\[Omega]1)^2 (\[CapitalOmega]-\[Omega]2))+(I E^(-2 I t \[CapitalOmega]-I t \[Omega]1+I t \[Omega]2) p1^4 f[\[Omega]1] f[\[Omega]2])/((\[CapitalOmega]+\[Omega]1)^2 (\[CapitalOmega]-\[Omega]2))+(I E^(-2 I t \[CapitalOmega]+I t \[Omega]1-I t \[Omega]2) p1^4 f[\[Omega]1] f[\[Omega]2])/((\[CapitalOmega]-\[Omega]1)^2 (\[CapitalOmega]+\[Omega]2))-(I E^(2 I t \[CapitalOmega]-I t \[Omega]1+I t \[Omega]2) p1^4 f[\[Omega]1] f[\[Omega]2])/((\[CapitalOmega]-\[Omega]1)^2 (\[CapitalOmega]+\[Omega]2))+(I E^(-2 I t \[CapitalOmega]-I t \[Omega]1-I t \[Omega]2) p1^4 f[\[Omega]1] f[\[Omega]2])/((\[CapitalOmega]+\[Omega]1)^2 (\[CapitalOmega]+\[Omega]2))-(I E^(2 I t \[CapitalOmega]+I t \[Omega]1+I t \[Omega]2) p1^4 f[\[Omega]1] f[\[Omega]2])/((\[CapitalOmega]+\[Omega]1)^2 (\[CapitalOmega]+\[Omega]2))-(2 I E^(-I t \[CapitalOmega]+I t \[Omega]1-I t \[Omega]2) p1^2 p3^2 (\[CapitalOmega] \[Omega]2 f[\[Omega]1] f[\[Omega]2]+\[Omega]1 (\[Omega]1-\[Omega]2) J[\[Omega]1] J[\[Omega]2]))/((\[CapitalOmega]-\[Omega]1) \[Omega]1^2 \[Omega]2 (\[CapitalOmega]-\[Omega]1+\[Omega]2))+(2 I E^(I t \[CapitalOmega]-I t \[Omega]1+I t \[Omega]2) p1^2 p3^2 (\[CapitalOmega] \[Omega]2 f[\[Omega]1] f[\[Omega]2]+\[Omega]1 (\[Omega]1-\[Omega]2) J[\[Omega]1] J[\[Omega]2]))/((\[CapitalOmega]-\[Omega]1) \[Omega]1^2 \[Omega]2 (\[CapitalOmega]-\[Omega]1+\[Omega]2))+(2 I E^(I t \[CapitalOmega]+I t \[Omega]1-I t \[Omega]2) p1^2 p3^2 (\[CapitalOmega] \[Omega]2 f[\[Omega]1] f[\[Omega]2]+\[Omega]1 (-\[Omega]1+\[Omega]2) J[\[Omega]1] J[\[Omega]2]))/(\[Omega]1^2 (\[CapitalOmega]+\[Omega]1) (\[CapitalOmega]+\[Omega]1-\[Omega]2) \[Omega]2)-(2 I E^(-I t \[CapitalOmega]-I t \[Omega]1+I t \[Omega]2) p1^2 p3^2 (\[CapitalOmega] \[Omega]2 f[\[Omega]1] f[\[Omega]2]+\[Omega]1 (-\[Omega]1+\[Omega]2) J[\[Omega]1] J[\[Omega]2]))/(\[Omega]1^2 (\[CapitalOmega]+\[Omega]1) (\[CapitalOmega]+\[Omega]1-\[Omega]2) \[Omega]2)+(2 I E^(-I t \[CapitalOmega]-I t \[Omega]1-I t \[Omega]2) p1^2 p3^2 (-\[CapitalOmega] \[Omega]2 f[\[Omega]1] f[\[Omega]2]+\[Omega]1 (\[Omega]1+\[Omega]2) J[\[Omega]1] J[\[Omega]2]))/(\[Omega]1^2 (\[CapitalOmega]+\[Omega]1) \[Omega]2 (\[CapitalOmega]+\[Omega]1+\[Omega]2))-(2 I E^(I t \[CapitalOmega]+I t \[Omega]1+I t \[Omega]2) p1^2 p3^2 (-\[CapitalOmega] \[Omega]2 f[\[Omega]1] f[\[Omega]2]+\[Omega]1 (\[Omega]1+\[Omega]2) J[\[Omega]1] J[\[Omega]2]))/(\[Omega]1^2 (\[CapitalOmega]+\[Omega]1) \[Omega]2 (\[CapitalOmega]+\[Omega]1+\[Omega]2))+(2 I E^(I t \[CapitalOmega]-I t \[Omega]1-I t \[Omega]2) p1^2 p3^2 (\[CapitalOmega] \[Omega]2 f[\[Omega]1] f[\[Omega]2]+\[Omega]1 (\[Omega]1+\[Omega]2) J[\[Omega]1] J[\[Omega]2]))/((\[CapitalOmega]-\[Omega]1) \[Omega]1^2 (\[CapitalOmega]-\[Omega]1-\[Omega]2) \[Omega]2)-(2 I E^(-I t \[CapitalOmega]+I t \[Omega]1+I t \[Omega]2) p1^2 p3^2 (\[CapitalOmega] \[Omega]2 f[\[Omega]1] f[\[Omega]2]+\[Omega]1 (\[Omega]1+\[Omega]2) J[\[Omega]1] J[\[Omega]2]))/((\[CapitalOmega]-\[Omega]1) \[Omega]1^2 (\[CapitalOmega]-\[Omega]1-\[Omega]2) \[Omega]2)-(2 I E^(I t \[CapitalOmega]-I t \[Omega]2) p1^2 ((2 p3^2 (\[CapitalOmega]-\[Omega]1)^2 (\[CapitalOmega]+\[Omega]1)^2 (\[CapitalOmega]-\[Omega]2) (\[Omega]1-\[Omega]2) (\[Omega]1+\[Omega]2)-I p1^2 \[CapitalOmega] \[Omega]1^2 (\[CapitalOmega]^2 (-I+t \[CapitalOmega]) \[Omega]1^2-(I+t \[CapitalOmega]) \[Omega]1^4+(-2 I \[CapitalOmega]^3+\[CapitalOmega] (4 I-t \[CapitalOmega]) \[Omega]1^2+t \[Omega]1^4) \[Omega]2+(\[CapitalOmega]^2 (3 I-t \[CapitalOmega])+(-I+t \[CapitalOmega]) \[Omega]1^2) \[Omega]2^2+(\[CapitalOmega] (-2 I+t \[CapitalOmega])-t \[Omega]1^2) \[Omega]2^3)) f[\[Omega]1] f[\[Omega]2]+2 p3^2 \[CapitalOmega] (\[CapitalOmega]-\[Omega]1) \[Omega]1 (\[CapitalOmega]+\[Omega]1) (\[CapitalOmega]-\[Omega]2) (\[Omega]1-\[Omega]2) (\[Omega]1+\[Omega]2) J[\[Omega]1] J[\[Omega]2]))/((\[CapitalOmega]-\[Omega]1)^2 \[Omega]1^2 (\[CapitalOmega]+\[Omega]1)^2 (\[CapitalOmega]-\[Omega]2)^2 (\[Omega]1-\[Omega]2) (\[Omega]1+\[Omega]2))+(2 I E^(-I t \[CapitalOmega]+I t \[Omega]2) p1^2 ((2 p3^2 (\[CapitalOmega]-\[Omega]1)^2 (\[CapitalOmega]+\[Omega]1)^2 (\[CapitalOmega]-\[Omega]2) (\[Omega]1-\[Omega]2) (\[Omega]1+\[Omega]2)+I p1^2 \[CapitalOmega] \[Omega]1^2 (\[CapitalOmega]^2 (I+t \[CapitalOmega]) \[Omega]1^2+(I-t \[CapitalOmega]) \[Omega]1^4+(2 I \[CapitalOmega]^3-\[CapitalOmega] (4 I+t \[CapitalOmega]) \[Omega]1^2+t \[Omega]1^4) \[Omega]2+(-\[CapitalOmega]^2 (3 I+t \[CapitalOmega])+(I+t \[CapitalOmega]) \[Omega]1^2) \[Omega]2^2+(\[CapitalOmega] (2 I+t \[CapitalOmega])-t \[Omega]1^2) \[Omega]2^3)) f[\[Omega]1] f[\[Omega]2]+2 p3^2 \[CapitalOmega] (\[CapitalOmega]-\[Omega]1) \[Omega]1 (\[CapitalOmega]+\[Omega]1) (\[CapitalOmega]-\[Omega]2) (\[Omega]1-\[Omega]2) (\[Omega]1+\[Omega]2) J[\[Omega]1] J[\[Omega]2]))/((\[CapitalOmega]-\[Omega]1)^2 \[Omega]1^2 (\[CapitalOmega]+\[Omega]1)^2 (\[CapitalOmega]-\[Omega]2)^2 (\[Omega]1-\[Omega]2) (\[Omega]1+\[Omega]2))-(2 I E^(I t \[CapitalOmega]+I t \[Omega]2) p1^2 ((2 p3^2 (\[CapitalOmega]-\[Omega]1)^2 (\[CapitalOmega]+\[Omega]1)^2 (\[Omega]1-\[Omega]2) (\[CapitalOmega]+\[Omega]2) (\[Omega]1+\[Omega]2)-I p1^2 \[CapitalOmega] \[Omega]1^2 (\[CapitalOmega]^2 (-I+t \[CapitalOmega]) \[Omega]1^2-(I+t \[CapitalOmega]) \[Omega]1^4+(2 I \[CapitalOmega]^3+\[CapitalOmega] (-4 I+t \[CapitalOmega]) \[Omega]1^2-t \[Omega]1^4) \[Omega]2+(\[CapitalOmega]^2 (3 I-t \[CapitalOmega])+(-I+t \[CapitalOmega]) \[Omega]1^2) \[Omega]2^2+(2 I \[CapitalOmega]-t \[CapitalOmega]^2+t \[Omega]1^2) \[Omega]2^3)) f[\[Omega]1] f[\[Omega]2]-2 p3^2 \[CapitalOmega] (\[CapitalOmega]-\[Omega]1) \[Omega]1 (\[CapitalOmega]+\[Omega]1) (\[Omega]1-\[Omega]2) (\[CapitalOmega]+\[Omega]2) (\[Omega]1+\[Omega]2) J[\[Omega]1] J[\[Omega]2]))/((\[CapitalOmega]-\[Omega]1)^2 \[Omega]1^2 (\[CapitalOmega]+\[Omega]1)^2 (\[Omega]1-\[Omega]2) (\[CapitalOmega]+\[Omega]2)^2 (\[Omega]1+\[Omega]2))+(2 I E^(-I t \[CapitalOmega]-I t \[Omega]2) p1^2 ((2 p3^2 (\[CapitalOmega]-\[Omega]1)^2 (\[CapitalOmega]+\[Omega]1)^2 (\[Omega]1-\[Omega]2) (\[CapitalOmega]+\[Omega]2) (\[Omega]1+\[Omega]2)+I p1^2 \[CapitalOmega] \[Omega]1^2 (\[CapitalOmega]^2 (I+t \[CapitalOmega]) \[Omega]1^2+(I-t \[CapitalOmega]) \[Omega]1^4+(-2 I \[CapitalOmega]^3+\[CapitalOmega] (4 I+t \[CapitalOmega]) \[Omega]1^2-t \[Omega]1^4) \[Omega]2+(-\[CapitalOmega]^2 (3 I+t \[CapitalOmega])+(I+t \[CapitalOmega]) \[Omega]1^2) \[Omega]2^2+(-\[CapitalOmega] (2 I+t \[CapitalOmega])+t \[Omega]1^2) \[Omega]2^3)) f[\[Omega]1] f[\[Omega]2]-2 p3^2 \[CapitalOmega] (\[CapitalOmega]-\[Omega]1) \[Omega]1 (\[CapitalOmega]+\[Omega]1) (\[Omega]1-\[Omega]2) (\[CapitalOmega]+\[Omega]2) (\[Omega]1+\[Omega]2) J[\[Omega]1] J[\[Omega]2]))/((\[CapitalOmega]-\[Omega]1)^2 \[Omega]1^2 (\[CapitalOmega]+\[Omega]1)^2 (\[Omega]1-\[Omega]2) (\[CapitalOmega]+\[Omega]2)^2 (\[Omega]1+\[Omega]2))-(2 I E^(-I t \[Omega]1-I t \[Omega]2) (p1^4 \[Omega]1 \[Omega]2 f[\[Omega]1] f[\[Omega]2]-2 p1^2 p3^2 (\[CapitalOmega]^2-\[Omega]1 \[Omega]2) J[\[Omega]1] J[\[Omega]2]))/((\[CapitalOmega]-\[Omega]1) \[Omega]1 (\[CapitalOmega]+\[Omega]1) (\[CapitalOmega]-\[Omega]2) (\[CapitalOmega]+\[Omega]2))+(2 I E^(I t \[Omega]1+I t \[Omega]2) (p1^4 \[Omega]1 \[Omega]2 f[\[Omega]1] f[\[Omega]2]-2 p1^2 p3^2 (\[CapitalOmega]^2-\[Omega]1 \[Omega]2) J[\[Omega]1] J[\[Omega]2]))/((\[CapitalOmega]-\[Omega]1) \[Omega]1 (\[CapitalOmega]+\[Omega]1) (\[CapitalOmega]-\[Omega]2) (\[CapitalOmega]+\[Omega]2))+(2 I E^(-I t \[Omega]1+I t \[Omega]2) (p1^4 \[Omega]1 \[Omega]2 f[\[Omega]1] f[\[Omega]2]-2 p1^2 p3^2 (\[CapitalOmega]^2+\[Omega]1 \[Omega]2) J[\[Omega]1] J[\[Omega]2]))/((\[CapitalOmega]-\[Omega]1) \[Omega]1 (\[CapitalOmega]+\[Omega]1) (\[CapitalOmega]-\[Omega]2) (\[CapitalOmega]+\[Omega]2))+(2 I E^(I t \[Omega]1-I t \[Omega]2) (-p1^4 \[Omega]1 \[Omega]2 f[\[Omega]1] f[\[Omega]2]+2 p1^2 p3^2 (\[CapitalOmega]^2+\[Omega]1 \[Omega]2) J[\[Omega]1] J[\[Omega]2]))/((\[CapitalOmega]-\[Omega]1) \[Omega]1 (\[CapitalOmega]+\[Omega]1) (\[CapitalOmega]-\[Omega]2) (\[CapitalOmega]+\[Omega]2))-(4 I E^(I t \[CapitalOmega]-I t \[Omega]1) p1^2 \[CapitalOmega] (p1^2 \[Omega]1 (-\[CapitalOmega]+\[Omega]1-\[Omega]2) \[Omega]2 (-\[CapitalOmega]+\[Omega]1+\[Omega]2) f[\[Omega]1] f[\[Omega]2]+p3^2 (\[Omega]1-\[Omega]2) (\[Omega]1+\[Omega]2) (\[CapitalOmega]^2-\[CapitalOmega] \[Omega]1+\[Omega]2^2) J[\[Omega]1] J[\[Omega]2]))/((\[CapitalOmega]-\[Omega]1) (\[CapitalOmega]-\[Omega]2) (\[CapitalOmega]-\[Omega]1-\[Omega]2) (\[Omega]1-\[Omega]2) \[Omega]2 (\[CapitalOmega]+\[Omega]2) (\[CapitalOmega]-\[Omega]1+\[Omega]2) (\[Omega]1+\[Omega]2))+(4 I E^(-I t \[CapitalOmega]+I t \[Omega]1) p1^2 \[CapitalOmega] (p1^2 \[Omega]1 (-\[CapitalOmega]+\[Omega]1-\[Omega]2) \[Omega]2 (-\[CapitalOmega]+\[Omega]1+\[Omega]2) f[\[Omega]1] f[\[Omega]2]+p3^2 (\[Omega]1-\[Omega]2) (\[Omega]1+\[Omega]2) (\[CapitalOmega]^2-\[CapitalOmega] \[Omega]1+\[Omega]2^2) J[\[Omega]1] J[\[Omega]2]))/((\[CapitalOmega]-\[Omega]1) (\[CapitalOmega]-\[Omega]2) (\[CapitalOmega]-\[Omega]1-\[Omega]2) (\[Omega]1-\[Omega]2) \[Omega]2 (\[CapitalOmega]+\[Omega]2) (\[CapitalOmega]-\[Omega]1+\[Omega]2) (\[Omega]1+\[Omega]2))+(8 I E^(-I t \[Omega]2) p1^2 p3^2 ((f[\[Omega]1] f[\[Omega]2])/\[Omega]2+(\[CapitalOmega]^2 (-\[CapitalOmega]^2-3 \[Omega]1^2+\[Omega]2^2) J[\[Omega]1] J[\[Omega]2])/(\[Omega]1 (-\[CapitalOmega]+\[Omega]1-\[Omega]2) (\[CapitalOmega]+\[Omega]1-\[Omega]2) (-\[CapitalOmega]+\[Omega]1+\[Omega]2) (\[CapitalOmega]+\[Omega]1+\[Omega]2))))/((\[CapitalOmega]-\[Omega]1) (\[CapitalOmega]+\[Omega]1))-(8 I E^(I t \[Omega]2) p1^2 p3^2 ((f[\[Omega]1] f[\[Omega]2])/\[Omega]2+(\[CapitalOmega]^2 (-\[CapitalOmega]^2-3 \[Omega]1^2+\[Omega]2^2) J[\[Omega]1] J[\[Omega]2])/(\[Omega]1 (-\[CapitalOmega]+\[Omega]1-\[Omega]2) (\[CapitalOmega]+\[Omega]1-\[Omega]2) (-\[CapitalOmega]+\[Omega]1+\[Omega]2) (\[CapitalOmega]+\[Omega]1+\[Omega]2))))/((\[CapitalOmega]-\[Omega]1) (\[CapitalOmega]+\[Omega]1))-(4 I E^(-I t \[CapitalOmega]-I t \[Omega]1) p1^2 \[CapitalOmega] (p1^2 \[Omega]1 (\[CapitalOmega]+\[Omega]1-\[Omega]2) \[Omega]2 (\[CapitalOmega]+\[Omega]1+\[Omega]2) f[\[Omega]1] f[\[Omega]2]+p3^2 (\[Omega]1-\[Omega]2) (\[Omega]1+\[Omega]2) (\[CapitalOmega] (\[CapitalOmega]+\[Omega]1)+\[Omega]2^2) J[\[Omega]1] J[\[Omega]2]))/((\[CapitalOmega]+\[Omega]1) (\[CapitalOmega]-\[Omega]2) (\[Omega]1-\[Omega]2) (\[CapitalOmega]+\[Omega]1-\[Omega]2) \[Omega]2 (\[CapitalOmega]+\[Omega]2) (\[Omega]1+\[Omega]2) (\[CapitalOmega]+\[Omega]1+\[Omega]2))+(4 I E^(I t \[CapitalOmega]+I t \[Omega]1) p1^2 \[CapitalOmega] (p1^2 \[Omega]1 (\[CapitalOmega]+\[Omega]1-\[Omega]2) \[Omega]2 (\[CapitalOmega]+\[Omega]1+\[Omega]2) f[\[Omega]1] f[\[Omega]2]+p3^2 (\[Omega]1-\[Omega]2) (\[Omega]1+\[Omega]2) (\[CapitalOmega] (\[CapitalOmega]+\[Omega]1)+\[Omega]2^2) J[\[Omega]1] J[\[Omega]2]))/((\[CapitalOmega]+\[Omega]1) (\[CapitalOmega]-\[Omega]2) (\[Omega]1-\[Omega]2) (\[CapitalOmega]+\[Omega]1-\[Omega]2) \[Omega]2 (\[CapitalOmega]+\[Omega]2) (\[Omega]1+\[Omega]2) (\[CapitalOmega]+\[Omega]1+\[Omega]2));


F304Int =-((I E^(2 I t \[CapitalOmega]-I t \[Omega]1-I t \[Omega]2) p1^4 f[\[Omega]2] J[\[Omega]1])/((\[CapitalOmega]-\[Omega]1)^2 (\[CapitalOmega]-\[Omega]2)))+(I E^(-2 I t \[CapitalOmega]+I t \[Omega]1+I t \[Omega]2) p1^4 f[\[Omega]2] J[\[Omega]1])/((\[CapitalOmega]-\[Omega]1)^2 (\[CapitalOmega]-\[Omega]2))+(I E^(2 I t \[CapitalOmega]+I t \[Omega]1-I t \[Omega]2) p1^4 f[\[Omega]2] J[\[Omega]1])/((\[CapitalOmega]+\[Omega]1)^2 (\[CapitalOmega]-\[Omega]2))-(I E^(-2 I t \[CapitalOmega]-I t \[Omega]1+I t \[Omega]2) p1^4 f[\[Omega]2] J[\[Omega]1])/((\[CapitalOmega]+\[Omega]1)^2 (\[CapitalOmega]-\[Omega]2))+(I E^(-2 I t \[CapitalOmega]+I t \[Omega]1-I t \[Omega]2) p1^4 f[\[Omega]2] J[\[Omega]1])/((\[CapitalOmega]-\[Omega]1)^2 (\[CapitalOmega]+\[Omega]2))-(I E^(2 I t \[CapitalOmega]-I t \[Omega]1+I t \[Omega]2) p1^4 f[\[Omega]2] J[\[Omega]1])/((\[CapitalOmega]-\[Omega]1)^2 (\[CapitalOmega]+\[Omega]2))-(I E^(-2 I t \[CapitalOmega]-I t \[Omega]1-I t \[Omega]2) p1^4 f[\[Omega]2] J[\[Omega]1])/((\[CapitalOmega]+\[Omega]1)^2 (\[CapitalOmega]+\[Omega]2))+(I E^(2 I t \[CapitalOmega]+I t \[Omega]1+I t \[Omega]2) p1^4 f[\[Omega]2] J[\[Omega]1])/((\[CapitalOmega]+\[Omega]1)^2 (\[CapitalOmega]+\[Omega]2))+(8 I E^(-I t \[Omega]2) p1^2 p3^2 \[CapitalOmega] f[\[Omega]2] J[\[Omega]1])/(\[CapitalOmega]^2 \[Omega]1 \[Omega]2-\[Omega]1^3 \[Omega]2)-(8 I E^(I t \[Omega]2) p1^2 p3^2 \[CapitalOmega] f[\[Omega]2] J[\[Omega]1])/(\[CapitalOmega]^2 \[Omega]1 \[Omega]2-\[Omega]1^3 \[Omega]2)-(8 I E^(-I t \[Omega]1) p1^2 p3^2 \[CapitalOmega] (3 \[CapitalOmega]^2-\[Omega]1^2+\[Omega]2^2) f[\[Omega]2] J[\[Omega]1])/((\[CapitalOmega]-\[Omega]2) (\[CapitalOmega]-\[Omega]1-\[Omega]2) (\[CapitalOmega]+\[Omega]1-\[Omega]2) (\[CapitalOmega]+\[Omega]2) (\[CapitalOmega]-\[Omega]1+\[Omega]2) (\[CapitalOmega]+\[Omega]1+\[Omega]2))+(8 I E^(I t \[Omega]1) p1^2 p3^2 \[CapitalOmega] (3 \[CapitalOmega]^2-\[Omega]1^2+\[Omega]2^2) f[\[Omega]2] J[\[Omega]1])/((\[CapitalOmega]-\[Omega]2) (\[CapitalOmega]-\[Omega]1-\[Omega]2) (\[CapitalOmega]+\[Omega]1-\[Omega]2) (\[CapitalOmega]+\[Omega]2) (\[CapitalOmega]-\[Omega]1+\[Omega]2) (\[CapitalOmega]+\[Omega]1+\[Omega]2))+(2 I E^(I t \[CapitalOmega]-I t \[Omega]1-I t \[Omega]2) p1^2 p3^2 (\[CapitalOmega] (\[CapitalOmega]-\[Omega]1) \[Omega]1+(2 \[CapitalOmega]-\[Omega]1) \[Omega]2^2-\[Omega]2^3) f[\[Omega]2] J[\[Omega]1])/((\[CapitalOmega]-\[Omega]1) \[Omega]1 (\[CapitalOmega]-\[Omega]2) (\[CapitalOmega]-\[Omega]1-\[Omega]2) \[Omega]2^2)-(2 I E^(-I t \[CapitalOmega]+I t \[Omega]1+I t \[Omega]2) p1^2 p3^2 (\[CapitalOmega] (\[CapitalOmega]-\[Omega]1) \[Omega]1+(2 \[CapitalOmega]-\[Omega]1) \[Omega]2^2-\[Omega]2^3) f[\[Omega]2] J[\[Omega]1])/((\[CapitalOmega]-\[Omega]1) \[Omega]1 (\[CapitalOmega]-\[Omega]2) (\[CapitalOmega]-\[Omega]1-\[Omega]2) \[Omega]2^2)+(2 I E^(-I t \[CapitalOmega]-I t \[Omega]1-I t \[Omega]2) p1^2 p3^2 (\[CapitalOmega] \[Omega]1 (\[CapitalOmega]+\[Omega]1)-(2 \[CapitalOmega]+\[Omega]1) \[Omega]2^2-\[Omega]2^3) f[\[Omega]2] J[\[Omega]1])/(\[Omega]1 (\[CapitalOmega]+\[Omega]1) \[Omega]2^2 (\[CapitalOmega]+\[Omega]2) (\[CapitalOmega]+\[Omega]1+\[Omega]2))-(2 I E^(I t \[CapitalOmega]+I t \[Omega]1+I t \[Omega]2) p1^2 p3^2 (\[CapitalOmega] \[Omega]1 (\[CapitalOmega]+\[Omega]1)-(2 \[CapitalOmega]+\[Omega]1) \[Omega]2^2-\[Omega]2^3) f[\[Omega]2] J[\[Omega]1])/(\[Omega]1 (\[CapitalOmega]+\[Omega]1) \[Omega]2^2 (\[CapitalOmega]+\[Omega]2) (\[CapitalOmega]+\[Omega]1+\[Omega]2))-(2 I E^(-I t \[CapitalOmega]+I t \[Omega]1-I t \[Omega]2) p1^2 p3^2 (\[CapitalOmega] (\[CapitalOmega]-\[Omega]1) \[Omega]1+(2 \[CapitalOmega]-\[Omega]1) \[Omega]2^2+\[Omega]2^3) f[\[Omega]2] J[\[Omega]1])/((\[CapitalOmega]-\[Omega]1) \[Omega]1 \[Omega]2^2 (\[CapitalOmega]+\[Omega]2) (\[CapitalOmega]-\[Omega]1+\[Omega]2))+(2 I E^(I t \[CapitalOmega]-I t \[Omega]1+I t \[Omega]2) p1^2 p3^2 (\[CapitalOmega] (\[CapitalOmega]-\[Omega]1) \[Omega]1+(2 \[CapitalOmega]-\[Omega]1) \[Omega]2^2+\[Omega]2^3) f[\[Omega]2] J[\[Omega]1])/((\[CapitalOmega]-\[Omega]1) \[Omega]1 \[Omega]2^2 (\[CapitalOmega]+\[Omega]2) (\[CapitalOmega]-\[Omega]1+\[Omega]2))-(2 I E^(I t \[CapitalOmega]+I t \[Omega]1-I t \[Omega]2) p1^2 p3^2 (\[CapitalOmega] \[Omega]1 (\[CapitalOmega]+\[Omega]1)-(2 \[CapitalOmega]+\[Omega]1) \[Omega]2^2+\[Omega]2^3) f[\[Omega]2] J[\[Omega]1])/(\[Omega]1 (\[CapitalOmega]+\[Omega]1) (\[CapitalOmega]-\[Omega]2) (\[CapitalOmega]+\[Omega]1-\[Omega]2) \[Omega]2^2)+(2 I E^(-I t \[CapitalOmega]-I t \[Omega]1+I t \[Omega]2) p1^2 p3^2 (\[CapitalOmega] \[Omega]1 (\[CapitalOmega]+\[Omega]1)-(2 \[CapitalOmega]+\[Omega]1) \[Omega]2^2+\[Omega]2^3) f[\[Omega]2] J[\[Omega]1])/(\[Omega]1 (\[CapitalOmega]+\[Omega]1) (\[CapitalOmega]-\[Omega]2) (\[CapitalOmega]+\[Omega]1-\[Omega]2) \[Omega]2^2)-(4 I E^(I t \[CapitalOmega]-I t \[Omega]2) p1^2 \[CapitalOmega] ((p3^2 (\[CapitalOmega]-\[Omega]1) (\[CapitalOmega]+\[Omega]1) (\[CapitalOmega]-\[Omega]2) (2 \[CapitalOmega]-\[Omega]2))/((\[CapitalOmega]+\[Omega]1-\[Omega]2) (-\[CapitalOmega]+\[Omega]1+\[Omega]2))+(p1^2 \[Omega]1^2 (\[CapitalOmega]^2+\[Omega]1^2-2 \[Omega]2^2))/((\[Omega]1-\[Omega]2) (\[Omega]1+\[Omega]2))) f[\[Omega]2] J[\[Omega]1])/((\[CapitalOmega]-\[Omega]1)^2 \[Omega]1 (\[CapitalOmega]+\[Omega]1)^2 (-\[CapitalOmega]+\[Omega]2))+(4 I E^(-I t \[CapitalOmega]-I t \[Omega]2) p1^2 \[CapitalOmega] ((p3^2 (\[CapitalOmega]-\[Omega]1) (\[CapitalOmega]+\[Omega]1) (2 \[CapitalOmega]+\[Omega]2))/((\[CapitalOmega]-\[Omega]1+\[Omega]2) (\[CapitalOmega]+\[Omega]1+\[Omega]2))-(p1^2 \[Omega]1^2 (\[CapitalOmega]^2+\[Omega]1^2-2 \[Omega]2^2))/((\[Omega]1-\[Omega]2) (\[CapitalOmega]+\[Omega]2) (\[Omega]1+\[Omega]2))) f[\[Omega]2] J[\[Omega]1])/((\[CapitalOmega]-\[Omega]1)^2 \[Omega]1 (\[CapitalOmega]+\[Omega]1)^2)-(4 I E^(I t \[CapitalOmega]+I t \[Omega]2) p1^2 \[CapitalOmega] ((p3^2 (\[CapitalOmega]-\[Omega]1) (\[CapitalOmega]+\[Omega]1) (2 \[CapitalOmega]+\[Omega]2))/((\[CapitalOmega]-\[Omega]1+\[Omega]2) (\[CapitalOmega]+\[Omega]1+\[Omega]2))-(p1^2 \[Omega]1^2 (\[CapitalOmega]^2+\[Omega]1^2-2 \[Omega]2^2))/((\[Omega]1-\[Omega]2) (\[CapitalOmega]+\[Omega]2) (\[Omega]1+\[Omega]2))) f[\[Omega]2] J[\[Omega]1])/((\[CapitalOmega]-\[Omega]1)^2 \[Omega]1 (\[CapitalOmega]+\[Omega]1)^2)-(4 I E^(-I t \[CapitalOmega]+I t \[Omega]2) p1^2 \[CapitalOmega] ((p3^2 (\[CapitalOmega]-\[Omega]1) (\[CapitalOmega]+\[Omega]1) (2 \[CapitalOmega]-\[Omega]2))/(\[CapitalOmega]+\[Omega]1-\[Omega]2)+(p1^2 \[Omega]1^2 (-\[CapitalOmega]+\[Omega]1+\[Omega]2) (\[CapitalOmega]^2+\[Omega]1^2-2 \[Omega]2^2))/((\[CapitalOmega]-\[Omega]2) (\[Omega]1-\[Omega]2) (\[Omega]1+\[Omega]2))) f[\[Omega]2] J[\[Omega]1])/((\[CapitalOmega]-\[Omega]1)^2 \[Omega]1 (\[CapitalOmega]+\[Omega]1)^2 (-\[CapitalOmega]+\[Omega]1+\[Omega]2))-(2 I E^(I t \[CapitalOmega]-I t \[Omega]1) p1^2 ((2 p3^2 (\[CapitalOmega]-\[Omega]1))/\[Omega]2^2+(p1^2 (\[CapitalOmega]^2 \[Omega]1 (\[Omega]1^2-I t \[CapitalOmega] \[Omega]1^2+\[CapitalOmega]^2 (-2+I t \[Omega]1))+(-I t \[CapitalOmega]^4+\[Omega]1^3+\[CapitalOmega]^2 \[Omega]1 (1-I t \[Omega]1)+\[CapitalOmega] \[Omega]1^2 (-2+I t \[Omega]1)+\[CapitalOmega]^3 (2+I t \[Omega]1)) \[Omega]2^2+I (t \[CapitalOmega] (\[CapitalOmega]-\[Omega]1)+I \[Omega]1) \[Omega]2^4))/((\[CapitalOmega]-\[Omega]2)^2 (\[CapitalOmega]+\[Omega]2)^2 (-\[Omega]1+\[Omega]2) (\[Omega]1+\[Omega]2))) f[\[Omega]2] J[\[Omega]1])/(\[CapitalOmega]-\[Omega]1)^2+(2 I E^(-I t \[CapitalOmega]+I t \[Omega]1) p1^2 ((2 p3^2 (\[CapitalOmega]-\[Omega]1))/\[Omega]2^2+(p1^2 (\[CapitalOmega]^2 \[Omega]1 (\[Omega]1^2+I t \[CapitalOmega] \[Omega]1^2+\[CapitalOmega]^2 (-2-I t \[Omega]1))+(I t \[CapitalOmega]^4+\[Omega]1^3+\[CapitalOmega] \[Omega]1^2 (-2-I t \[Omega]1)+\[CapitalOmega]^3 (2-I t \[Omega]1)+\[CapitalOmega]^2 \[Omega]1 (1+I t \[Omega]1)) \[Omega]2^2+I (I \[Omega]1+t \[CapitalOmega] (-\[CapitalOmega]+\[Omega]1)) \[Omega]2^4))/((\[CapitalOmega]-\[Omega]2)^2 (\[CapitalOmega]+\[Omega]2)^2 (-\[Omega]1+\[Omega]2) (\[Omega]1+\[Omega]2))) f[\[Omega]2] J[\[Omega]1])/(\[CapitalOmega]-\[Omega]1)^2-(2 I E^(-I t \[CapitalOmega]-I t \[Omega]1) p1^2 ((2 p3^2 (\[CapitalOmega]+\[Omega]1))/\[Omega]2^2+(p1^2 (\[CapitalOmega]^2 \[Omega]1 (-\[Omega]1^2-I t \[CapitalOmega] \[Omega]1^2+\[CapitalOmega]^2 (2-I t \[Omega]1))+I (\[CapitalOmega]^3 (-2 I+t \[CapitalOmega])+\[CapitalOmega]^2 (I+t \[CapitalOmega]) \[Omega]1+\[CapitalOmega] (2 I+t \[CapitalOmega]) \[Omega]1^2+(I+t \[CapitalOmega]) \[Omega]1^3) \[Omega]2^2+(\[Omega]1-I t \[CapitalOmega] (\[CapitalOmega]+\[Omega]1)) \[Omega]2^4))/((\[CapitalOmega]-\[Omega]2)^2 (\[CapitalOmega]+\[Omega]2)^2 (-\[Omega]1+\[Omega]2) (\[Omega]1+\[Omega]2))) f[\[Omega]2] J[\[Omega]1])/(\[CapitalOmega]+\[Omega]1)^2+(2 I E^(I t \[CapitalOmega]+I t \[Omega]1) p1^2 ((2 p3^2 (\[CapitalOmega]+\[Omega]1))/\[Omega]2^2+(p1^2 (\[CapitalOmega]^2 \[Omega]1 (-\[Omega]1^2+I t \[CapitalOmega] \[Omega]1^2+\[CapitalOmega]^2 (2+I t \[Omega]1))-I (\[CapitalOmega]^3 (2 I+t \[CapitalOmega])+\[CapitalOmega]^2 (-I+t \[CapitalOmega]) \[Omega]1+\[CapitalOmega] (-2 I+t \[CapitalOmega]) \[Omega]1^2+(-I+t \[CapitalOmega]) \[Omega]1^3) \[Omega]2^2+(\[Omega]1+I t \[CapitalOmega] (\[CapitalOmega]+\[Omega]1)) \[Omega]2^4))/((\[CapitalOmega]-\[Omega]2)^2 (\[CapitalOmega]+\[Omega]2)^2 (-\[Omega]1+\[Omega]2) (\[Omega]1+\[Omega]2))) f[\[Omega]2] J[\[Omega]1])/(\[CapitalOmega]+\[Omega]1)^2-(2 I E^(-I t \[Omega]1-I t \[Omega]2) p1^2 \[CapitalOmega] ((2 p3^2 (\[CapitalOmega]-\[Omega]1) (\[CapitalOmega]+\[Omega]1))/(\[Omega]1 \[Omega]2)+(p1^2 (-3 \[CapitalOmega]^4+\[Omega]1 \[Omega]2 (-2 \[Omega]1^2+\[Omega]1 \[Omega]2-2 \[Omega]2^2)+\[CapitalOmega]^2 (\[Omega]1^2+4 \[Omega]1 \[Omega]2+\[Omega]2^2)))/(\[CapitalOmega]^2-\[Omega]2^2)^2) f[\[Omega]2] J[\[Omega]1])/((\[CapitalOmega]-\[Omega]1)^2 (\[CapitalOmega]+\[Omega]1)^2)+(2 I E^(I t \[Omega]1+I t \[Omega]2) p1^2 \[CapitalOmega] ((2 p3^2 (\[CapitalOmega]-\[Omega]1) (\[CapitalOmega]+\[Omega]1))/(\[Omega]1 \[Omega]2)+(p1^2 (-3 \[CapitalOmega]^4+\[Omega]1 \[Omega]2 (-2 \[Omega]1^2+\[Omega]1 \[Omega]2-2 \[Omega]2^2)+\[CapitalOmega]^2 (\[Omega]1^2+4 \[Omega]1 \[Omega]2+\[Omega]2^2)))/(\[CapitalOmega]^2-\[Omega]2^2)^2) f[\[Omega]2] J[\[Omega]1])/((\[CapitalOmega]-\[Omega]1)^2 (\[CapitalOmega]+\[Omega]1)^2)-(2 I E^(I t \[Omega]1-I t \[Omega]2) p1^2 \[CapitalOmega] ((2 p3^2 (\[CapitalOmega]^2-\[Omega]1^2))/(\[Omega]1 \[Omega]2)-(p1^2 (-3 \[CapitalOmega]^4+\[CapitalOmega]^2 (\[Omega]1^2-4 \[Omega]1 \[Omega]2+\[Omega]2^2)+\[Omega]1 \[Omega]2 (2 \[Omega]1^2+\[Omega]1 \[Omega]2+2 \[Omega]2^2)))/(\[CapitalOmega]^2-\[Omega]2^2)^2) f[\[Omega]2] J[\[Omega]1])/(\[CapitalOmega]^2-\[Omega]1^2)^2+(2 I E^(-I t \[Omega]1+I t \[Omega]2) p1^2 \[CapitalOmega] ((2 p3^2 (\[CapitalOmega]^2-\[Omega]1^2))/(\[Omega]1 \[Omega]2)-(p1^2 (-3 \[CapitalOmega]^4+\[CapitalOmega]^2 (\[Omega]1^2-4 \[Omega]1 \[Omega]2+\[Omega]2^2)+\[Omega]1 \[Omega]2 (2 \[Omega]1^2+\[Omega]1 \[Omega]2+2 \[Omega]2^2)))/(\[CapitalOmega]^2-\[Omega]2^2)^2) f[\[Omega]2] J[\[Omega]1])/(\[CapitalOmega]^2-\[Omega]1^2)^2;


K4BeforeTripleInt = {{0,0,0,0},{4 I E^(-I (t+t1+t2+t3) \[CapitalOmega]) p1 p3 ((-E^(I t1 \[CapitalOmega])+E^(I t2 \[CapitalOmega])) (E^(I (2 t+t1) \[CapitalOmega]) p1^2+E^(I (2 t+t2) \[CapitalOmega]) p1^2+E^(I (t1+2 t3) \[CapitalOmega]) p1^2+E^(I (t2+2 t3) \[CapitalOmega]) p1^2+2 E^(I (2 t+t3) \[CapitalOmega]) p3^2+2 E^(I (t1+t2+t3) \[CapitalOmega]) p3^2) \[Eta][t1-t3] \[Nu][t-t2]+(E^(I t2 \[CapitalOmega])-E^(I t3 \[CapitalOmega])) (E^(I (2 t+t2) \[CapitalOmega]) p1^2+E^(I (2 t1+t2) \[CapitalOmega]) p1^2+E^(I (2 t+t3) \[CapitalOmega]) p1^2+E^(I (2 t1+t3) \[CapitalOmega]) p1^2+2 E^(I (2 t+t1) \[CapitalOmega]) p3^2+2 E^(I (t1+t2+t3) \[CapitalOmega]) p3^2) \[Eta][t-t3] \[Nu][t1-t2]+(-E^(I t1 \[CapitalOmega])+E^(I t2 \[CapitalOmega])) (E^(I (2 t+t1) \[CapitalOmega]) p1^2+E^(I (2 t+t2) \[CapitalOmega]) p1^2+E^(I (t1+2 t3) \[CapitalOmega]) p1^2+E^(I (t2+2 t3) \[CapitalOmega]) p1^2+2 E^(I (2 t+t3) \[CapitalOmega]) p3^2+2 E^(I (t1+t2+t3) \[CapitalOmega]) p3^2) \[Eta][t1-t2] \[Nu][t-t3]),-4 E^(-I (3 t-t1-t2-t3) \[CapitalOmega]) p1^2 p3^2 (E^(I (t-2 (t1+t2+t3)) \[CapitalOmega]) (E^(I t1 \[CapitalOmega])-E^(I t2 \[CapitalOmega])) (E^(2 I t \[CapitalOmega])-E^(I (t1+t2) \[CapitalOmega])) (E^(2 I t \[CapitalOmega])+E^(2 I t3 \[CapitalOmega])) \[Eta][t1-t2] \[Eta][t-t3]+E^(I (t-2 (t1+t2+t3)) \[CapitalOmega]) (E^(I t1 \[CapitalOmega])-E^(I t2 \[CapitalOmega])) (E^(2 I t \[CapitalOmega])-E^(I (t1+t2) \[CapitalOmega])) (E^(2 I t \[CapitalOmega])+E^(2 I t3 \[CapitalOmega])) \[Eta][t-t2] \[Eta][t1-t3]+E^(I (t-t1) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]-E^(I (3 t-t1-2 t2) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]+E^(I (t-t2) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]+E^(I (3 t-2 t1-t2) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]-E^(I (3 t-t1-2 t3) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]+E^(I (3 t-t2-2 t3) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]+E^(I (5 t-2 t1-t2-2 t3) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]-2 E^(I (t-t3) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]-2 E^(I (5 t-2 t1-2 t2-t3) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]+E^(I (5 t-t1-2 (t2+t3)) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]+E^(I (t-t1) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]-E^(I (3 t-t1-2 t2) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]-E^(I (t-t2) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]+E^(I (3 t-2 t1-t2) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]-E^(I (3 t-t1-2 t3) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]+E^(I (3 t-t2-2 t3) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]-E^(I (5 t-2 t1-t2-2 t3) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]+E^(I (5 t-t1-2 (t2+t3)) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]),4 I E^(-I (3 t-t1-t2-t3) \[CapitalOmega]) p1^2 p3^2 (E^(I (t-2 (t1+t2+t3)) \[CapitalOmega]) (E^(I t1 \[CapitalOmega])-E^(I t2 \[CapitalOmega])) (E^(2 I t \[CapitalOmega])-E^(I (t1+t2) \[CapitalOmega])) (E^(2 I t \[CapitalOmega])-E^(2 I t3 \[CapitalOmega])) \[Eta][t1-t2] \[Eta][t-t3]+E^(I (t-2 (t1+t2+t3)) \[CapitalOmega]) (E^(I t1 \[CapitalOmega])-E^(I t2 \[CapitalOmega])) (E^(2 I t \[CapitalOmega])-E^(I (t1+t2) \[CapitalOmega])) (E^(2 I t \[CapitalOmega])-E^(2 I t3 \[CapitalOmega])) \[Eta][t-t2] \[Eta][t1-t3]-E^(I (t-t1) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]+E^(I (3 t-t1-2 t2) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]-E^(I (t-t2) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]+E^(I (3 t-2 t1-t2) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]-E^(I (3 t-t1-2 t3) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]-E^(I (3 t-t2-2 t3) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]+E^(I (5 t-2 t1-t2-2 t3) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]+2 E^(I (t-t3) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]-2 E^(I (5 t-2 t1-2 t2-t3) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]+E^(I (5 t-t1-2 (t2+t3)) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]-E^(I (t-t1) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]+E^(I (3 t-t1-2 t2) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]+E^(I (t-t2) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]-E^(I (3 t-2 t1-t2) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]-E^(I (3 t-t1-2 t3) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]+E^(I (3 t-t2-2 t3) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]-E^(I (5 t-2 t1-t2-2 t3) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]+2 E^(I (3 t-2 t1-t3) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]-2 E^(I (3 t-2 t2-t3) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]+E^(I (5 t-t1-2 (t2+t3)) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]),4 E^(-I (3 t-t1-t2-t3) \[CapitalOmega]) p1 p3 (2 E^(I (2 t-2 t1-2 t2-t3) \[CapitalOmega]) (E^(I t1 \[CapitalOmega])-E^(I t2 \[CapitalOmega])) (E^(2 I t \[CapitalOmega])-E^(I (t1+t2) \[CapitalOmega])) p3^2 \[Eta][t1-t2] \[Eta][t-t3]+E^(-2 I (-4 t+t1+t2+2 t3) \[CapitalOmega]) (2 E^(-3 I (2 t-t3) \[CapitalOmega]) (E^(I t1 \[CapitalOmega])-E^(I t2 \[CapitalOmega])) (E^(2 I t \[CapitalOmega])-E^(I (t1+t2) \[CapitalOmega])) p3^2 \[Eta][t-t2] \[Eta][t1-t3]+E^(-2 I (3 t-t3) \[CapitalOmega]) (-((E^(2 I (t+t1) \[CapitalOmega]) p1^2-E^(2 I (t1+t2) \[CapitalOmega]) p1^2-E^(2 I (t+t3) \[CapitalOmega]) p1^2+E^(2 I (t2+t3) \[CapitalOmega]) p1^2-2 E^(I (2 t+t1+t2) \[CapitalOmega]) p3^2+2 E^(I (2 t+t1+t3) \[CapitalOmega]) p3^2+2 E^(I (t1+2 t2+t3) \[CapitalOmega]) p3^2-2 E^(I (t1+t2+2 t3) \[CapitalOmega]) p3^2) \[Nu][t1-t2] \[Nu][t-t3])-(E^(2 I t1 \[CapitalOmega])-E^(2 I t2 \[CapitalOmega])) (E^(2 I t \[CapitalOmega])-E^(2 I t3 \[CapitalOmega])) p1^2 \[Nu][t-t2] \[Nu][t1-t3])))},{4 E^(-I (t+t1+t2+t3) \[CapitalOmega]) p1 p3 (-((E^(I t1 \[CapitalOmega])-E^(I t2 \[CapitalOmega])) (E^(I (2 t+t1) \[CapitalOmega]) p1^2+E^(I (2 t+t2) \[CapitalOmega]) p1^2-2 E^(I (t+t1+t2) \[CapitalOmega]) p1^2+2 E^(I (t+2 t3) \[CapitalOmega]) p1^2-E^(I (t1+2 t3) \[CapitalOmega]) p1^2-E^(I (t2+2 t3) \[CapitalOmega]) p1^2+2 E^(I (2 t+t3) \[CapitalOmega]) p3^2-2 E^(I (t1+t2+t3) \[CapitalOmega]) p3^2) \[Eta][t1-t3] \[Nu][t-t2])+(E^(I t2 \[CapitalOmega])-E^(I t3 \[CapitalOmega])) (2 E^(I (t+2 t1) \[CapitalOmega]) p1^2+E^(I (2 t+t2) \[CapitalOmega]) p1^2-E^(I (2 t1+t2) \[CapitalOmega]) p1^2+E^(I (2 t+t3) \[CapitalOmega]) p1^2-E^(I (2 t1+t3) \[CapitalOmega]) p1^2-2 E^(I (t+t2+t3) \[CapitalOmega]) p1^2+2 E^(I (2 t+t1) \[CapitalOmega]) p3^2-2 E^(I (t1+t2+t3) \[CapitalOmega]) p3^2) \[Eta][t-t3] \[Nu][t1-t2]-(E^(I t1 \[CapitalOmega])-E^(I t2 \[CapitalOmega])) (E^(I (2 t+t1) \[CapitalOmega]) p1^2+E^(I (2 t+t2) \[CapitalOmega]) p1^2-2 E^(I (t+t1+t2) \[CapitalOmega]) p1^2+2 E^(I (t+2 t3) \[CapitalOmega]) p1^2-E^(I (t1+2 t3) \[CapitalOmega]) p1^2-E^(I (t2+2 t3) \[CapitalOmega]) p1^2+2 E^(I (2 t+t3) \[CapitalOmega]) p3^2-2 E^(I (t1+t2+t3) \[CapitalOmega]) p3^2) \[Eta][t1-t2] \[Nu][t-t3]),4 I E^(-I (3 t-t1-t2-t3) \[CapitalOmega]) p1^2 (E^(I (t-2 (t1+t2+t3)) \[CapitalOmega]) (E^(I t1 \[CapitalOmega])-E^(I t2 \[CapitalOmega])) (E^(2 I t \[CapitalOmega])+E^(2 I t3 \[CapitalOmega])) (E^(I (t+t1) \[CapitalOmega]) p1^2+E^(I (t+t2) \[CapitalOmega]) p1^2+E^(2 I t \[CapitalOmega]) p3^2+E^(I (t1+t2) \[CapitalOmega]) p3^2) \[Eta][t1-t2] \[Eta][t-t3]+E^(I (t-2 (t1+t2+t3)) \[CapitalOmega]) (E^(I t1 \[CapitalOmega])-E^(I t2 \[CapitalOmega])) (E^(2 I t \[CapitalOmega])+E^(2 I t3 \[CapitalOmega])) (E^(I (t+t1) \[CapitalOmega]) p1^2+E^(I (t+t2) \[CapitalOmega]) p1^2+E^(2 I t \[CapitalOmega]) p3^2+E^(I (t1+t2) \[CapitalOmega]) p3^2) \[Eta][t-t2] \[Eta][t1-t3]-E^(2 I (t-t2) \[CapitalOmega]) p1^2 \[Nu][t1-t2] \[Nu][t-t3]-E^(2 I (2 t-t1-t2) \[CapitalOmega]) p1^2 \[Nu][t1-t2] \[Nu][t-t3]+E^(2 I (t-t3) \[CapitalOmega]) p1^2 \[Nu][t1-t2] \[Nu][t-t3]+E^(2 I (2 t-t1-t3) \[CapitalOmega]) p1^2 \[Nu][t1-t2] \[Nu][t-t3]-E^(I (t-t1) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]-E^(I (3 t-t1-2 t2) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]-E^(I (t-t2) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]-E^(I (3 t-2 t1-t2) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]+2 E^(I (2 t-t1-t2) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]+E^(I (3 t-t1-2 t3) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]+E^(I (3 t-t2-2 t3) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]+E^(I (5 t-2 t1-t2-2 t3) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]-2 E^(I (4 t-t1-t2-2 t3) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]+2 E^(I (t-t3) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]-2 E^(I (5 t-2 t1-2 t2-t3) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]-2 E^(I (2 t-t2-t3) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]+2 E^(I (4 t-2 t1-t2-t3) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]+E^(I (5 t-t1-2 (t2+t3)) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]-E^(I (t-t1) \[CapitalOmega]) p3^2 \[Nu][t-t2] \[Nu][t1-t3]-E^(I (3 t-t1-2 t2) \[CapitalOmega]) p3^2 \[Nu][t-t2] \[Nu][t1-t3]+E^(I (t-t2) \[CapitalOmega]) p3^2 \[Nu][t-t2] \[Nu][t1-t3]+E^(I (3 t-2 t1-t2) \[CapitalOmega]) p3^2 \[Nu][t-t2] \[Nu][t1-t3]+E^(I (3 t-t1-2 t3) \[CapitalOmega]) p3^2 \[Nu][t-t2] \[Nu][t1-t3]-E^(I (3 t-t2-2 t3) \[CapitalOmega]) p3^2 \[Nu][t-t2] \[Nu][t1-t3]-E^(I (5 t-2 t1-t2-2 t3) \[CapitalOmega]) p3^2 \[Nu][t-t2] \[Nu][t1-t3]-2 E^(I (3 t-2 t1-t3) \[CapitalOmega]) p3^2 \[Nu][t-t2] \[Nu][t1-t3]+2 E^(I (2 t-t1-t3) \[CapitalOmega]) p3^2 \[Nu][t-t2] \[Nu][t1-t3]+2 E^(I (3 t-2 t2-t3) \[CapitalOmega]) p3^2 \[Nu][t-t2] \[Nu][t1-t3]-2 E^(I (4 t-t1-2 t2-t3) \[CapitalOmega]) p3^2 \[Nu][t-t2] \[Nu][t1-t3]-2 E^(I (2 t-t2-t3) \[CapitalOmega]) p3^2 \[Nu][t-t2] \[Nu][t1-t3]+2 E^(I (4 t-2 t1-t2-t3) \[CapitalOmega]) p3^2 \[Nu][t-t2] \[Nu][t1-t3]+E^(I (5 t-t1-2 (t2+t3)) \[CapitalOmega]) p3^2 \[Nu][t-t2] \[Nu][t1-t3]),4 E^(-I (3 t-t1-t2-t3) \[CapitalOmega]) p1^2 (E^(I (t-2 (t1+t2+t3)) \[CapitalOmega]) (E^(I t1 \[CapitalOmega])-E^(I t2 \[CapitalOmega])) (E^(2 I t \[CapitalOmega])-E^(2 I t3 \[CapitalOmega])) (E^(I (t+t1) \[CapitalOmega]) p1^2+E^(I (t+t2) \[CapitalOmega]) p1^2+E^(2 I t \[CapitalOmega]) p3^2+E^(I (t1+t2) \[CapitalOmega]) p3^2) \[Eta][t1-t2] \[Eta][t-t3]+E^(I (t-2 (t1+t2+t3)) \[CapitalOmega]) (E^(I t1 \[CapitalOmega])-E^(I t2 \[CapitalOmega])) (E^(2 I t \[CapitalOmega])-E^(2 I t3 \[CapitalOmega])) (E^(I (t+t1) \[CapitalOmega]) p1^2+E^(I (t+t2) \[CapitalOmega]) p1^2+E^(2 I t \[CapitalOmega]) p3^2+E^(I (t1+t2) \[CapitalOmega]) p3^2) \[Eta][t-t2] \[Eta][t1-t3]+E^(2 I (t-t2) \[CapitalOmega]) p1^2 \[Nu][t1-t2] \[Nu][t-t3]-E^(2 I (2 t-t1-t2) \[CapitalOmega]) p1^2 \[Nu][t1-t2] \[Nu][t-t3]-E^(2 I (t-t3) \[CapitalOmega]) p1^2 \[Nu][t1-t2] \[Nu][t-t3]+E^(2 I (2 t-t1-t3) \[CapitalOmega]) p1^2 \[Nu][t1-t2] \[Nu][t-t3]+E^(I (t-t1) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]+E^(I (3 t-t1-2 t2) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]+E^(I (t-t2) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]-E^(I (3 t-2 t1-t2) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]-2 E^(I (2 t-t1-t2) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]+E^(I (3 t-t1-2 t3) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]-E^(I (3 t-t2-2 t3) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]+E^(I (5 t-2 t1-t2-2 t3) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]-2 E^(I (4 t-t1-t2-2 t3) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]-2 E^(I (t-t3) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]-2 E^(I (5 t-2 t1-2 t2-t3) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]+2 E^(I (2 t-t2-t3) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]+2 E^(I (4 t-2 t1-t2-t3) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]+E^(I (5 t-t1-2 (t2+t3)) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]+E^(I (t-t1) \[CapitalOmega]) p3^2 \[Nu][t-t2] \[Nu][t1-t3]+E^(I (3 t-t1-2 t2) \[CapitalOmega]) p3^2 \[Nu][t-t2] \[Nu][t1-t3]-E^(I (t-t2) \[CapitalOmega]) p3^2 \[Nu][t-t2] \[Nu][t1-t3]-E^(I (3 t-2 t1-t2) \[CapitalOmega]) p3^2 \[Nu][t-t2] \[Nu][t1-t3]+E^(I (3 t-t1-2 t3) \[CapitalOmega]) p3^2 \[Nu][t-t2] \[Nu][t1-t3]-E^(I (3 t-t2-2 t3) \[CapitalOmega]) p3^2 \[Nu][t-t2] \[Nu][t1-t3]-E^(I (5 t-2 t1-t2-2 t3) \[CapitalOmega]) p3^2 \[Nu][t-t2] \[Nu][t1-t3]-2 E^(I (2 t-t1-t3) \[CapitalOmega]) p3^2 \[Nu][t-t2] \[Nu][t1-t3]-2 E^(I (4 t-t1-2 t2-t3) \[CapitalOmega]) p3^2 \[Nu][t-t2] \[Nu][t1-t3]+2 E^(I (2 t-t2-t3) \[CapitalOmega]) p3^2 \[Nu][t-t2] \[Nu][t1-t3]+2 E^(I (4 t-2 t1-t2-t3) \[CapitalOmega]) p3^2 \[Nu][t-t2] \[Nu][t1-t3]+E^(I (5 t-t1-2 (t2+t3)) \[CapitalOmega]) p3^2 \[Nu][t-t2] \[Nu][t1-t3]),4 I E^(-I (3 t-t1-t2-t3) \[CapitalOmega]) p1 p3 (2 E^(I (2 t-2 t1-2 t2-t3) \[CapitalOmega]) (-E^(I t1 \[CapitalOmega])+E^(I t2 \[CapitalOmega])) (E^(I (t+t1) \[CapitalOmega]) p1^2+E^(I (t+t2) \[CapitalOmega]) p1^2+E^(2 I t \[CapitalOmega]) p3^2+E^(I (t1+t2) \[CapitalOmega]) p3^2) \[Eta][t1-t2] \[Eta][t-t3]+2 E^(I (2 t-2 t1-2 t2-t3) \[CapitalOmega]) (-E^(I t1 \[CapitalOmega])+E^(I t2 \[CapitalOmega])) (E^(I (t+t1) \[CapitalOmega]) p1^2+E^(I (t+t2) \[CapitalOmega]) p1^2+E^(2 I t \[CapitalOmega]) p3^2+E^(I (t1+t2) \[CapitalOmega]) p3^2) \[Eta][t-t2] \[Eta][t1-t3]-E^(2 I (t-t1) \[CapitalOmega]) p1^2 \[Nu][t1-t2] \[Nu][t-t3]+2 E^(I (3 t-t1-2 t2) \[CapitalOmega]) p1^2 \[Nu][t1-t2] \[Nu][t-t3]+2 E^(I (3 t-2 t1-t2) \[CapitalOmega]) p1^2 \[Nu][t1-t2] \[Nu][t-t3]-E^(2 I (2 t-t1-t2) \[CapitalOmega]) p1^2 \[Nu][t1-t2] \[Nu][t-t3]-2 E^(I (3 t-t1-2 t3) \[CapitalOmega]) p1^2 \[Nu][t1-t2] \[Nu][t-t3]-2 E^(I (3 t-t2-2 t3) \[CapitalOmega]) p1^2 \[Nu][t1-t2] \[Nu][t-t3]+E^(2 I (t-t3) \[CapitalOmega]) p1^2 \[Nu][t1-t2] \[Nu][t-t3]+E^(2 I (2 t-t2-t3) \[CapitalOmega]) p1^2 \[Nu][t1-t2] \[Nu][t-t3]+2 E^(I (2 t-t1-t2) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]-2 E^(I (4 t-t1-t2-2 t3) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]-2 E^(I (2 t-t1-t3) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]+2 E^(I (4 t-t1-2 t2-t3) \[CapitalOmega]) p3^2 \[Nu][t1-t2] \[Nu][t-t3]-E^(2 I (t-t1) \[CapitalOmega]) p1^2 \[Nu][t-t2] \[Nu][t1-t3]-2 E^(I (3 t-t1-2 t2) \[CapitalOmega]) p1^2 \[Nu][t-t2] \[Nu][t1-t3]+E^(2 I (t-t2) \[CapitalOmega]) p1^2 \[Nu][t-t2] \[Nu][t1-t3]+2 E^(I (3 t-2 t1-t2) \[CapitalOmega]) p1^2 \[Nu][t-t2] \[Nu][t1-t3]+2 E^(I (3 t-t1-2 t3) \[CapitalOmega]) p1^2 \[Nu][t-t2] \[Nu][t1-t3]-2 E^(I (3 t-t2-2 t3) \[CapitalOmega]) p1^2 \[Nu][t-t2] \[Nu][t1-t3]-E^(2 I (2 t-t1-t3) \[CapitalOmega]) p1^2 \[Nu][t-t2] \[Nu][t1-t3]+E^(2 I (2 t-t2-t3) \[CapitalOmega]) p1^2 \[Nu][t-t2] \[Nu][t1-t3])},{4 I E^(-I (t+t1+t2+t3) \[CapitalOmega]) p1^2 ((-E^(I t1 \[CapitalOmega])+E^(I t2 \[CapitalOmega])) (E^(I (2 t+t1) \[CapitalOmega]) p1^2+E^(I (2 t+t2) \[CapitalOmega]) p1^2+E^(I (t1+2 t3) \[CapitalOmega]) p1^2+E^(I (t2+2 t3) \[CapitalOmega]) p1^2+2 E^(I (2 t+t3) \[CapitalOmega]) p3^2+2 E^(I (t1+t2+t3) \[CapitalOmega]) p3^2) \[Eta][t1-t3] \[Nu][t-t2]+(E^(I t2 \[CapitalOmega])-E^(I t3 \[CapitalOmega])) (E^(I (2 t+t2) \[CapitalOmega]) p1^2+E^(I (2 t1+t2) \[CapitalOmega]) p1^2+E^(I (2 t+t3) \[CapitalOmega]) p1^2+E^(I (2 t1+t3) \[CapitalOmega]) p1^2+2 E^(I (2 t+t1) \[CapitalOmega]) p3^2+2 E^(I (t1+t2+t3) \[CapitalOmega]) p3^2) \[Eta][t-t3] \[Nu][t1-t2]+(-E^(I t1 \[CapitalOmega])+E^(I t2 \[CapitalOmega])) (E^(I (2 t+t1) \[CapitalOmega]) p1^2+E^(I (2 t+t2) \[CapitalOmega]) p1^2+E^(I (t1+2 t3) \[CapitalOmega]) p1^2+E^(I (t2+2 t3) \[CapitalOmega]) p1^2+2 E^(I (2 t+t3) \[CapitalOmega]) p3^2+2 E^(I (t1+t2+t3) \[CapitalOmega]) p3^2) \[Eta][t1-t2] \[Nu][t-t3]),-4 E^(-I (3 t-t1-t2-t3) \[CapitalOmega]) p1^3 p3 (E^(I (t-2 (t1+t2+t3)) \[CapitalOmega]) (E^(I t1 \[CapitalOmega])-E^(I t2 \[CapitalOmega])) (E^(2 I t \[CapitalOmega])-E^(I (t1+t2) \[CapitalOmega])) (E^(2 I t \[CapitalOmega])+E^(2 I t3 \[CapitalOmega])) \[Eta][t1-t2] \[Eta][t-t3]+E^(I (t-2 (t1+t2+t3)) \[CapitalOmega]) (E^(I t1 \[CapitalOmega])-E^(I t2 \[CapitalOmega])) (E^(2 I t \[CapitalOmega])-E^(I (t1+t2) \[CapitalOmega])) (E^(2 I t \[CapitalOmega])+E^(2 I t3 \[CapitalOmega])) \[Eta][t-t2] \[Eta][t1-t3]+E^(I (t-t1) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]-E^(I (3 t-t1-2 t2) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]+E^(I (t-t2) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]+E^(I (3 t-2 t1-t2) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]-E^(I (3 t-t1-2 t3) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]+E^(I (3 t-t2-2 t3) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]+E^(I (5 t-2 t1-t2-2 t3) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]-2 E^(I (t-t3) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]-2 E^(I (5 t-2 t1-2 t2-t3) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]+E^(I (5 t-t1-2 (t2+t3)) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]+E^(I (t-t1) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]-E^(I (3 t-t1-2 t2) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]-E^(I (t-t2) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]+E^(I (3 t-2 t1-t2) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]-E^(I (3 t-t1-2 t3) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]+E^(I (3 t-t2-2 t3) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]-E^(I (5 t-2 t1-t2-2 t3) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]+E^(I (5 t-t1-2 (t2+t3)) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]),4 I E^(-I (3 t-t1-t2-t3) \[CapitalOmega]) p1^3 p3 (E^(I (t-2 (t1+t2+t3)) \[CapitalOmega]) (E^(I t1 \[CapitalOmega])-E^(I t2 \[CapitalOmega])) (E^(2 I t \[CapitalOmega])-E^(I (t1+t2) \[CapitalOmega])) (E^(2 I t \[CapitalOmega])-E^(2 I t3 \[CapitalOmega])) \[Eta][t1-t2] \[Eta][t-t3]+E^(I (t-2 (t1+t2+t3)) \[CapitalOmega]) (E^(I t1 \[CapitalOmega])-E^(I t2 \[CapitalOmega])) (E^(2 I t \[CapitalOmega])-E^(I (t1+t2) \[CapitalOmega])) (E^(2 I t \[CapitalOmega])-E^(2 I t3 \[CapitalOmega])) \[Eta][t-t2] \[Eta][t1-t3]-E^(I (t-t1) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]+E^(I (3 t-t1-2 t2) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]-E^(I (t-t2) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]+E^(I (3 t-2 t1-t2) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]-E^(I (3 t-t1-2 t3) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]-E^(I (3 t-t2-2 t3) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]+E^(I (5 t-2 t1-t2-2 t3) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]+2 E^(I (t-t3) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]-2 E^(I (5 t-2 t1-2 t2-t3) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]+E^(I (5 t-t1-2 (t2+t3)) \[CapitalOmega]) \[Nu][t1-t2] \[Nu][t-t3]-E^(I (t-t1) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]+E^(I (3 t-t1-2 t2) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]+E^(I (t-t2) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]-E^(I (3 t-2 t1-t2) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]-E^(I (3 t-t1-2 t3) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]+E^(I (3 t-t2-2 t3) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]-E^(I (5 t-2 t1-t2-2 t3) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]+2 E^(I (3 t-2 t1-t3) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]-2 E^(I (3 t-2 t2-t3) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]+E^(I (5 t-t1-2 (t2+t3)) \[CapitalOmega]) \[Nu][t-t2] \[Nu][t1-t3]),4 E^(-I (3 t-t1-t2-t3) \[CapitalOmega]) p1^2 (2 E^(I (2 t-2 t1-2 t2-t3) \[CapitalOmega]) (E^(I t1 \[CapitalOmega])-E^(I t2 \[CapitalOmega])) (E^(2 I t \[CapitalOmega])-E^(I (t1+t2) \[CapitalOmega])) p3^2 \[Eta][t1-t2] \[Eta][t-t3]+E^(-2 I (-4 t+t1+t2+2 t3) \[CapitalOmega]) (2 E^(-3 I (2 t-t3) \[CapitalOmega]) (E^(I t1 \[CapitalOmega])-E^(I t2 \[CapitalOmega])) (E^(2 I t \[CapitalOmega])-E^(I (t1+t2) \[CapitalOmega])) p3^2 \[Eta][t-t2] \[Eta][t1-t3]+E^(-2 I (3 t-t3) \[CapitalOmega]) (-((E^(2 I (t+t1) \[CapitalOmega]) p1^2-E^(2 I (t1+t2) \[CapitalOmega]) p1^2-E^(2 I (t+t3) \[CapitalOmega]) p1^2+E^(2 I (t2+t3) \[CapitalOmega]) p1^2-2 E^(I (2 t+t1+t2) \[CapitalOmega]) p3^2+2 E^(I (2 t+t1+t3) \[CapitalOmega]) p3^2+2 E^(I (t1+2 t2+t3) \[CapitalOmega]) p3^2-2 E^(I (t1+t2+2 t3) \[CapitalOmega]) p3^2) \[Nu][t1-t2] \[Nu][t-t3])-(E^(2 I t1 \[CapitalOmega])-E^(2 I t2 \[CapitalOmega])) (E^(2 I t \[CapitalOmega])-E^(2 I t3 \[CapitalOmega])) p1^2 \[Nu][t-t2] \[Nu][t1-t3])))}};


(* ::Chapter:: *)
(*Keep The following section. Important*)


TwoPtTransRule = {
eta[T1_] eta[T2_] -> Hold[EtaIntegrand[\[Omega][3],T1]] Hold[EtaIntegrand[\[Omega][4],T2]], 
nu[T1_] nu[T2_] -> Hold[NuIntegrand[\[Omega][1],T1]] Hold[NuIntegrand[\[Omega][2],T2]],
eta[T1_] nu[T2_] -> Hold[EtaIntegrand[\[Omega][3],T1]] Hold[NuIntegrand[\[Omega][2],T2]]
};

(*
If you use TwoPtTransRule2, which \[Nu] will get \[Omega]1 and which will get \[Omega]2 is unclear.
But this should not make any difference to the result of the integral!
*)

TwoPtTransRule2 = {
\[Eta][T1_] \[Eta][T2_] -> Hold[EtaIntegrand[\[Omega]1,T1]] Hold[EtaIntegrand[\[Omega]2,T2]], 
\[Nu][T1_] \[Nu][T2_] -> Hold[NuIntegrand[\[Omega]1,T1]] Hold[NuIntegrand[\[Omega]2,T2]],
\[Eta][T1_] \[Nu][T2_] -> Hold[EtaIntegrand[\[Omega]1,T1]] Hold[NuIntegrand[\[Omega]2,T2]]
};


(* ::Section:: *)
(*TCL2 expressions*)


azz2Integrand = - Cos[(t-t1) \[CapitalOmega]] Hold[NuIntegrand[\[Omega], t - t1]];
ayx2Integrand = Sin[(t-t1) \[CapitalOmega]] Hold[NuIntegrand[\[Omega], t-t1]];
bz2Integrand = Sin[(t-t1) \[CapitalOmega]] Hold[EtaIntegrand[\[Omega], t-t1]];


Wiki[KapplerAyx2Int]


(* ::Section:: *)
(*TCL4 expressions (for simple Hamiltonian?)*)


Azz4IntTerm1 = Sin[\[CapitalOmega] (t1 - t2)] Sin[\[CapitalOmega] (t - t3)] Hold[NuIntegrand[\[Omega][1], t - t2]] Hold[NuIntegrand[\[Omega][2], t1 - t3]];

Azz4IntTerm2 = Sin[\[CapitalOmega] (t1 - t3)] Sin[\[CapitalOmega] (t - t2)] Hold[NuIntegrand[\[Omega][1], t - t3]] Hold[NuIntegrand[\[Omega][2], t1 - t2]];

Azz4TotalIntegrand = Azz4IntTerm1 + Azz4IntTerm2;

Azz4GenSpecIntResult = -f[\[Omega]1] f[\[Omega]2]((2 \[CapitalOmega]^2 Cos[t \[CapitalOmega]] (t (\[CapitalOmega]^2-\[Omega]1^2) (\[CapitalOmega]^2-\[Omega]2^2) (\[Omega]1^2-\[Omega]2^2) Cos[t \[Omega]1]+2 \[Omega]1 (\[CapitalOmega]^4+\[Omega]1^4-\[Omega]1^2 \[Omega]2^2+\[Omega]2^4-\[CapitalOmega]^2 (\[Omega]1^2+\[Omega]2^2)) Sin[t \[Omega]1]+2 (\[CapitalOmega]^2-\[Omega]1^2) \[Omega]2 (-\[CapitalOmega]^2+\[Omega]2^2) Sin[t \[Omega]2])+(\[CapitalOmega]^2-\[Omega]1^2) (\[Omega]1^2-\[Omega]2^2) Cos[t \[CapitalOmega]]^2 (2 \[Omega]1 (\[CapitalOmega]^2+\[Omega]2^2) Cos[t \[Omega]2] Sin[t \[Omega]1]+4 \[CapitalOmega]^2 \[Omega]2 Cos[t \[Omega]1] Sin[t \[Omega]2])-2 \[CapitalOmega] Cos[t \[Omega]1] ((\[CapitalOmega]^4 (\[Omega]1^2+\[Omega]2^2)+\[Omega]1^2 \[Omega]2^2 (\[Omega]1^2+\[Omega]2^2)+\[CapitalOmega]^2 (\[Omega]1^4-6 \[Omega]1^2 \[Omega]2^2+\[Omega]2^4)) Sin[t \[CapitalOmega]]+(\[CapitalOmega]^2-\[Omega]1^2) (\[Omega]1^2-\[Omega]2^2) (\[CapitalOmega]^2+\[Omega]2^2) Cos[t \[Omega]2] Sin[2 t \[CapitalOmega]]+2 \[CapitalOmega] (\[CapitalOmega]^2-\[Omega]1^2) \[Omega]2 (\[Omega]1^2-\[Omega]2^2) Sin[t \[CapitalOmega]]^2 Sin[t \[Omega]2])+2 (\[CapitalOmega]^2-\[Omega]1^2) (Cos[t \[Omega]2] (2 \[CapitalOmega] \[Omega]2^2 (\[CapitalOmega]^2-\[Omega]2^2) Sin[t \[CapitalOmega]]+\[Omega]1 (\[CapitalOmega]^2-\[Omega]2^2) (\[Omega]1^2-\[Omega]2^2) Sin[t \[Omega]1]-\[Omega]1 (\[Omega]1^2-\[Omega]2^2) (\[CapitalOmega]^2+\[Omega]2^2) Sin[t \[CapitalOmega]]^2 Sin[t \[Omega]1])+\[CapitalOmega] \[Omega]1 (\[Omega]1^2-\[Omega]2^2) Sin[t \[CapitalOmega]] Sin[t \[Omega]1] (t (\[CapitalOmega]^2-\[Omega]2^2)+4 \[Omega]2 Cos[t \[CapitalOmega]] Sin[t \[Omega]2])))/(4 (\[CapitalOmega]-\[Omega]1)^2 (\[CapitalOmega]+\[Omega]1)^2 (\[CapitalOmega]-\[Omega]2)^2 (\[Omega]1-\[Omega]2) (\[CapitalOmega]+\[Omega]2)^2 (\[Omega]1+\[Omega]2)));


Bz4IntTerm1 = Sin[\[CapitalOmega] (t1 - t2)] Cos[\[CapitalOmega] (t - t3)] Hold[NuIntegrand[\[Omega][1], t - t2]] Hold[EtaIntegrand[\[Omega][2],t1 - t3]];
Bz4IntTerm2 = Sin[\[CapitalOmega] (t1 - t3)] Cos[\[CapitalOmega] (t - t2)] Hold[NuIntegrand[\[Omega][1], t - t3]] Hold[EtaIntegrand[\[Omega][2], t1 - t2]];
Bz4IntTerm3 = - Sin[\[CapitalOmega] (t2 - t3)] Cos[\[CapitalOmega] (t - t1)] (Hold[NuIntegrand[\[Omega][1], t - t3]] Hold[EtaIntegrand[\[Omega][2], t1 - t2]] + Hold[EtaIntegrand[\[Omega][2], t - t3]] Hold[NuIntegrand[\[Omega][1], t1 - t2]]); 
(*Is this choice of assigning \[Omega]1 and \[Omega]2 in the 3rd term ok?*)

Bz4TotalIntegrand = Bz4IntTerm1 + Bz4IntTerm2 + Bz4IntTerm3;

Bz4GenSpecIntResult = - 2 (f[\[Omega]1] J[\[Omega]2] (2 \[CapitalOmega]^6 \[Omega]2 Cos[t \[Omega]2] Sin[t \[CapitalOmega]]-5 \[CapitalOmega]^4 \[Omega]1^2 \[Omega]2 Cos[t \[Omega]2] Sin[t \[CapitalOmega]]+\[CapitalOmega]^2 \[Omega]1^4 \[Omega]2 Cos[t \[Omega]2] Sin[t \[CapitalOmega]]+\[CapitalOmega]^4 \[Omega]2^3 Cos[t \[Omega]2] Sin[t \[CapitalOmega]]+2 \[CapitalOmega]^2 \[Omega]1^2 \[Omega]2^3 Cos[t \[Omega]2] Sin[t \[CapitalOmega]]+\[Omega]1^4 \[Omega]2^3 Cos[t \[Omega]2] Sin[t \[CapitalOmega]]-\[CapitalOmega]^2 \[Omega]2^5 Cos[t \[Omega]2] Sin[t \[CapitalOmega]]-\[Omega]1^2 \[Omega]2^5 Cos[t \[Omega]2] Sin[t \[CapitalOmega]]+4 \[CapitalOmega]^3 \[Omega]1^3 \[Omega]2 Cos[t \[Omega]2] Sin[t \[Omega]1]-2 \[CapitalOmega] \[Omega]1^5 \[Omega]2 Cos[t \[Omega]2] Sin[t \[Omega]1]-4 \[CapitalOmega]^3 \[Omega]1 \[Omega]2^3 Cos[t \[Omega]2] Sin[t \[Omega]1]+2 \[CapitalOmega] \[Omega]1 \[Omega]2^5 Cos[t \[Omega]2] Sin[t \[Omega]1]-2 \[CapitalOmega]^3 \[Omega]1^3 \[Omega]2 Cos[t \[Omega]2] Sin[t \[CapitalOmega]]^2 Sin[t \[Omega]1]+2 \[CapitalOmega] \[Omega]1^5 \[Omega]2 Cos[t \[Omega]2] Sin[t \[CapitalOmega]]^2 Sin[t \[Omega]1]+2 \[CapitalOmega]^3 \[Omega]1 \[Omega]2^3 Cos[t \[Omega]2] Sin[t \[CapitalOmega]]^2 Sin[t \[Omega]1]-2 \[CapitalOmega] \[Omega]1^3 \[Omega]2^3 Cos[t \[Omega]2] Sin[t \[CapitalOmega]]^2 Sin[t \[Omega]1]+t \[CapitalOmega]^6 \[Omega]1^2 Sin[t \[CapitalOmega]] Sin[t \[Omega]2]-t \[CapitalOmega]^4 \[Omega]1^4 Sin[t \[CapitalOmega]] Sin[t \[Omega]2]-t \[CapitalOmega]^6 \[Omega]2^2 Sin[t \[CapitalOmega]] Sin[t \[Omega]2]+t \[CapitalOmega]^2 \[Omega]1^4 \[Omega]2^2 Sin[t \[CapitalOmega]] Sin[t \[Omega]2]+t \[CapitalOmega]^4 \[Omega]2^4 Sin[t \[CapitalOmega]] Sin[t \[Omega]2]-t \[CapitalOmega]^2 \[Omega]1^2 \[Omega]2^4 Sin[t \[CapitalOmega]] Sin[t \[Omega]2]+\[CapitalOmega]^4 \[Omega]1^3 Sin[2 t \[CapitalOmega]] Sin[t \[Omega]1] Sin[t \[Omega]2]-\[CapitalOmega]^2 \[Omega]1^5 Sin[2 t \[CapitalOmega]] Sin[t \[Omega]1] Sin[t \[Omega]2]-\[CapitalOmega]^4 \[Omega]1 \[Omega]2^2 Sin[2 t \[CapitalOmega]] Sin[t \[Omega]1] Sin[t \[Omega]2]+2 \[CapitalOmega]^2 \[Omega]1^3 \[Omega]2^2 Sin[2 t \[CapitalOmega]] Sin[t \[Omega]1] Sin[t \[Omega]2]-\[Omega]1^5 \[Omega]2^2 Sin[2 t \[CapitalOmega]] Sin[t \[Omega]1] Sin[t \[Omega]2]-\[CapitalOmega]^2 \[Omega]1 \[Omega]2^4 Sin[2 t \[CapitalOmega]] Sin[t \[Omega]1] Sin[t \[Omega]2]+\[Omega]1^3 \[Omega]2^4 Sin[2 t \[CapitalOmega]] Sin[t \[Omega]1] Sin[t \[Omega]2]+\[CapitalOmega] Cos[t \[CapitalOmega]] (t (\[CapitalOmega]^2-\[Omega]1^2) \[Omega]2 (\[CapitalOmega]^2-\[Omega]2^2) (\[Omega]1^2-\[Omega]2^2) Cos[t \[Omega]2]+2 \[Omega]1 (-\[CapitalOmega]^2+\[Omega]1^2) \[Omega]2 (-\[CapitalOmega]^2+2 \[Omega]1^2-\[Omega]2^2) Sin[t \[Omega]1]+2 (-\[Omega]1^4 \[Omega]2^2+\[CapitalOmega]^4 (\[Omega]1^2-2 \[Omega]2^2)+\[CapitalOmega]^2 \[Omega]2^2 (\[Omega]1^2+\[Omega]2^2)) Sin[t \[Omega]2])+1/2 \[CapitalOmega] (\[CapitalOmega]^2-\[Omega]1^2) (\[Omega]1^2-\[Omega]2^2) Cos[t \[CapitalOmega]]^2 (4 \[Omega]1 \[Omega]2 Cos[t \[Omega]2] Sin[t \[Omega]1]+2 (\[CapitalOmega]^2+\[Omega]2^2) Cos[t \[Omega]1] Sin[t \[Omega]2])+\[CapitalOmega] Cos[t \[Omega]1] (-2 \[CapitalOmega] (\[CapitalOmega]^2-\[Omega]1^2) \[Omega]2 (\[CapitalOmega]^2-2 \[Omega]1^2+\[Omega]2^2) Sin[t \[CapitalOmega]]+(\[CapitalOmega]^2-\[Omega]1^2) (\[CapitalOmega]^2+\[Omega]2^2) (-\[Omega]1^2+\[Omega]2^2) Sin[t \[CapitalOmega]]^2 Sin[t \[Omega]2]+(\[Omega]1^2-\[Omega]2^2) (-2 \[CapitalOmega] (\[CapitalOmega]^2-\[Omega]1^2) \[Omega]2 Cos[t \[Omega]2] Sin[2 t \[CapitalOmega]]+(-3 \[CapitalOmega]^4+\[Omega]1^2 \[Omega]2^2+\[CapitalOmega]^2 (\[Omega]1^2+\[Omega]2^2)) Sin[t \[Omega]2]))))/(4 (\[CapitalOmega]-\[Omega]1)^2 (\[CapitalOmega]+\[Omega]1)^2 (\[CapitalOmega]-\[Omega]2)^2 (\[Omega]1-\[Omega]2) (\[CapitalOmega]+\[Omega]2)^2 (\[Omega]1+\[Omega]2));


(* ::Chapter:: *)
(*Functions*)


(* ::Section:: *)
(*TCL4 Functions*)


GenHAssump = {t>0, \[Theta] \[Element] Reals, \[Omega]2>0, \[Omega]1>0, \[Omega]>0, \[CapitalOmega] > 0, p1 \[Element] Reals, p3 \[Element] Reals}; 

TCL4TripleInt[TotalIntegrand_, Assump_:GenHAssump] := 
  Assuming[Assump,
    Integrate[
    Integrate[Integrate[TotalIntegrand, {t3, 0, t2}], {t2, 0, t1}], {t1, 0, 
     t}]]


(*
Both the following functions assume that the argument of the exponent is not of the form I t (\[Omega]1 + \[Omega]2)
because in that case it will not be able to identify and isolate the individual parts.
*)

PowIter[term_,var_,m_] := Coefficient[ExprNiceForm[term], Exp[ Simplify[ I var t] ], m]  Exp[ Simplify[ m I var t] ];

PowIter2[term_,var_,m_] := Coefficient[term, Exp[ Simplify[ I var t] ], m]  Exp[ Simplify[ m I var t] ];


FreeArgument[testExpr_, f_] := ReleaseHold[testExpr //. f[x_] -> Hold[f[Expand[x]]]] (*For f in testExpr, expands the arguments of f*)

ExprNiceForm[e0_]:= Module[{e1, e2, e3, e4, e5, e6},
	(*e1 = FreeArgument[e0, Exp];*)
	(*e2 = Coefficient[e1, Exp[I t \[CapitalOmega]]] ExpToTrig[Exp[I t \[CapitalOmega]]] + Coefficient[e1,Exp[-I t \[CapitalOmega]]] ExpToTrig[Exp[-I t \[CapitalOmega]]] + Coefficient[e1, Exp[I t \[CapitalOmega]], 0];*)
	e2 = TrigToExp[e0];(*Gets rid of trignometric terms.*)
	e3 = Expand[e2]; (*Writes all terms as sum of terms so that there is just one exp term.*)
	e4 = collectCoefficientT[e3, Exp]; (*Factors t out of the exponents*)
	e5 = CollectTerms[e4, Exp]; (*Collects terms with same exp argument together*)
	FreeArgument[e5, Exp]
	
]


(*
The z input to FindDenSolutions is to exclude the pole at x = z.

This part is there because in the \[Omega]2 integration loop, we are assuming that
z means a pole in the bulk region. So, if the solution for \[Omega]2 contains
z, then that solution can't exist on the number line.

Notes: Inbuilt FunctionPoles function of mathematica fails for me.
*)

FindDenSolutions[expr_, x_, z_:Unique["noMatch"]] := Module[
  {terms, denominators, relevantDenominators, solutions, PoleList, SimpSol},
  
  If[expr ===0, Return[{}]];
  
  (* Expand the expression and extract all individual terms *)
  terms = List @@ Expand[expr];
  
  (* Extract the denominator of each term *)
  denominators = Denominator /@ terms;
  
  (* Select denominators that contain the variable x *)
  relevantDenominators = Select[denominators, ! FreeQ[#, x] &];
  
  (* Solve each relevant denominator equal to zero for x and collect solutions *)
  solutions = Flatten[
    Table[
      x /. Solve[denominator == 0, x],
      {denominator, relevantDenominators}
    ]
  ];
  
  SimpSol = Simplify[solutions, Assumptions -> GenHAssump];
  
  
  PoleList = deleteTermsWithZ[Union[SimpSol], z]; (* What is this code? Read the comment above. *)
  
  Print["PoleList for ", x, " was ",  PoleList];
  (* Return only unique solutions *)
  PoleList
]


(* Define the function *)
deleteTermsWithZ[list_List, z_] := Select[list, FreeQ[#, z] &]


(* Helper function: Checks one omega symbol in an expression based on the name string *)
checkOmega[expr_, name_String, sym_Symbol] := Module[{integrated},
  (* We assume that the occurrence of "om1" or "om2" in the name means it was integrated, i.e. it should not appear. *)
  integrated = StringContainsQ[name, ToString[sym]];
  If[integrated && !FreeQ[expr, sym],
    Print["WARNING: In expression '", name, "', variable ", sym, " was integrated (name contains '", ToString[sym],
          "') but appears in the expression. See the expression: ", expr];
  ];
  If[!integrated && FreeQ[expr, sym],
    Print["WARNING: In expression '", name, "', variable ", sym, " was not integrated (name does not contain '", ToString[sym],
          "') but it is missing from the expression. See the expression: ", expr];
  ];
]

(* Main function: Checks both om1 and om2 in an expression given its name *)
CheckOmegaIntegration[expr_, name_String] := Module[{},
  If[FullSimplify[expr] === 0, Return[0]];
  checkOmega[expr, name, \[Omega]1];
  checkOmega[expr, name, \[Omega]2];
];


(*SkipZero is now legacy. I don't trigger it from AllResCareful function*)

OneResExtractSeg[expr_, x_, i_, PoleListX_, z_:z1] := Module[
{ResX, ExtractedTerm, XBulk, TotalXInregral},
	
	ZeroSelect[m_]:= If[m==0, 1, 0];
	
	(*Fix this. R and R2 functions are unnecessary and ugly.*)
	
	BulkPole[term_, var_, m_]:= ZeroSelect[m] term //. {var -> z};
	
	PoleSign[m_] := If[m == 0, 1, Sign[m]];
	
	ExtractedTerm = PowIter2[expr, x, i]; (*extract the exp(I m x t) term from expression*);
	
	MyPrint["For variable = ", x, ", m = ",i];
	MyPrint["Osci Expression = ", ExtractedTerm];
	
	ResX = \[Pi] I PoleSign[i] Rsidu[ExtractedTerm, x, PoleListX]; (*Do re4sidue calcs on this term and add a sign. For no exponential term, positive sign.*)
	
	XBulk = BulkPole[ExtractedTerm, x, i ]; 
	
	{ResX, XBulk}
	];


(*This resloop code is better because it can help debug*)
ResLoopSeg[expr_, xval_, PoleListz_, z_:z1]:= Module[{PoleSum, i, PartSum, ExpandedExpr, Bulk, PartBulk, Divergence, PoleListX, PoleListX2},
    If[expr ===0, Return[{0,0,0}]];
    
    ExpandedExpr = ExprNiceForm[expr];
    PoleListX2 = FindDenSolutions[expr, xval, PoleListz];
    
    PoleListX = FindDenSolutions[ExpandedExpr, xval, PoleListz];
    
    If[PoleListX != PoleListX2, Print["Warning. Pole calculation may have issue."]];
	
	MyPrint["Nice form Expr = ", ExpandedExpr];
	PoleSum = 0;
	Bulk = 0;
	For[i = -1, i < 2, i++, {PartSum, PartBulk} = OneResExtractSeg[ExpandedExpr, xval, i, PoleListX, z];
	PoleSum = PoleSum + PartSum;
	Bulk = Bulk + PartBulk;	
	];
	
	Divergence = O2Res[expr, xval, PoleListX];
	{PoleSum, Bulk, Divergence}
]



AllResSeg[expr_]:= Module[{PoleListX, PoleListY1,PoleListY2, om1Poles, om2Poles, startTime, elapsedTime, endTime, Divergence1,Divergence2, Divergence3, Result, 
Om1Bulk, Om2Bulk, Only\[Omega]1Done, \[Omega]1\[Omega]2Done, Only\[Omega]2Done, NoneDone, om1Done, Om1NotDone, SingleInt},
	startTime = AbsoluteTime[];
	(*It is important that this z be supplied to the subsequent functions because
	it will be later used here for further calcs.*)
	
	If[expr === 0, Return[{0,0, 0, 0}]];
	
	 (*Specifying the z here does not make sense. But it does not hurt either because expr is expected to be devoid of z at this stage*)

	{om1Done, Om1NotDone, Divergence1} = ResLoopSeg[expr, \[Omega]1, z, z];
	

	{\[Omega]1\[Omega]2Done, Only\[Omega]1Done, Divergence2} = ResLoopSeg[om1Done, \[Omega]2, z, z1] //. {z -> \[Omega]1, z1 -> \[Omega]2};
	
	(*PoleListY2 = FindDenSolutions[Om1NotDone, \[Omega]2, z]*);
	{Only\[Omega]2Done, NoneDone, Divergence3} = ResLoopSeg[Om1NotDone, \[Omega]2, z, z1] //. {z -> \[Omega]1, z1 -> \[Omega]2};


	(* Check the expressions *)
	CheckOmegaIntegration[Only\[Omega]1Done, "Only\[Omega]1Done"];
	CheckOmegaIntegration[Only\[Omega]2Done, "Only\[Omega]2Done"];
	CheckOmegaIntegration[\[Omega]1\[Omega]2Done, "\[Omega]1\[Omega]2Done"];
	CheckOmegaIntegration[NoneDone, "NoneDone"];
	
	
	SingleInt = (Only\[Omega]1Done + Only\[Omega]2Done) //. { \[Omega]1 -> \[Omega], \[Omega]2 -> \[Omega]};
    
    endTime = AbsoluteTime[];
    elapsedTime = endTime - startTime;
  
    Print["Time taken to run AllRes: ", NumberForm[elapsedTime, {5, 3}], " seconds"];
    
    {\[Omega]1\[Omega]2Done, SingleInt, NoneDone, Divergence1 + Divergence2 + Divergence3}
]


(* ::Section:: *)
(*Other Functions*)


(*The problem with this code is, once you have a complicated function, it is difficult to extricate parts of it that are independent of some variable.

This is because you could have (x-y)/(ax - ay). if you expand, it becomes (x/(ax-ay)) - (y/(ax-ay)). Mathematica would fail it pick that it is independent of both x and y.

But the following function is needed in drude cutoff evaluation to separate the nu and eta terms in the begining itself. So, this code should be in drude cutoff function actually.
*)

decomposeExpr[expr_, x_, y_] := Module[{A, terms, expr2},
  (* Initialize a 2\[Times]2 list. 
     Note: Mathematica lists start at index 1.
     Here we use:
         A[[1,1]] for terms free of both x and y,
         A[[2,1]] for terms containing x but not y,
         A[[1,2]] for terms containing y but not x,
         A[[2,2]] for terms containing both x and y. *)
  A = {{0, 0}, {0, 0}};
  
  expr;
  
  (* Expand the expression and decompose it into addends *)
  terms = TermsAsList[Expand[expr]];
  
  (* Process each term *)
  Do[
    Module[{i, j},
      (* Set i = 1 if term is free of x, else 2 *)
      i = If[FreeQ[Simplify[term], x], 1, 2];
      (* Set j = 1 if term is free of y, else 2 *)
      j = If[FreeQ[Simplify[term], y], 1, 2];
      (* Add the term into the appropriate slot *)
      
      A[[i, j]] += term;
    ],
    {term, terms}
  ];
  
  A
]


Min[\[Pi], 5]


(* ::Section:: *)
(*Numerical Integral Eval*)


Coth[\[Beta] \[Omega]/2] /. \[Omega] -> 2 \[Pi] I /\[Beta]


ContourIntegral[expr_, poleList_:{}, prec_:16] := Module[{e2, \[Delta]},
  
  \[Delta]  = 0.5*Min[ poleList ];
  
  Print["\[Delta] = ",\[Delta]];
  
  e2 = expr /. {\[Omega] -> \[Omega] + I\[ThinSpace]\[Delta]};
  NIntegrate[e2, {\[Omega], -\[Infinity], \[Infinity]}, WorkingPrecision -> prec]
];

DoubleContourIntegral[expr_, poleList_:{}, prec_:16] := Module[{e2, \[Delta]},
  \[Delta]  = 0.5*Min[poleList];
  e2 = expr /. {
    \[Omega]1 -> \[Omega]1 + I\[ThinSpace]\[Delta],
    \[Omega]2 -> \[Omega]2 + 0.5\[ThinSpace]I\[ThinSpace]\[Delta]
  };
  NIntegrate[
    e2,
    {\[Omega]2, -\[Infinity], \[Infinity]},
    {\[Omega]1, -\[Infinity], \[Infinity]}, WorkingPrecision -> prec
  ]
];


(* ::Text:: *)
(*I can set the WorkingPrecision -> 20 above because for some values, *)


ClearAll[getDefaultWPandPG]
getDefaultWPandPG[f_] := Module[{wp, pg},
  wp = WorkingPrecision /. Options[f, WorkingPrecision];
  pg = PrecisionGoal    /. Options[f, PrecisionGoal];
  {If[wp === MachinePrecision, $MachinePrecision, wp],
   Replace[pg, Automatic :> Floor[If[wp === MachinePrecision,
                                      $MachinePrecision, wp]/2]]}
]

(* Example calls *)
getDefaultWPandPG[NIntegrate]
getDefaultWPandPG[NDSolve]


ContourIntegralBetter[expr_, LambdaNumVal_, BetaNumVal_, prec_:16] := Module[{e2, \[Delta]},
  \[Delta]  = 0.5*Min[LambdaNumVal, 2 \[Pi]/BetaNumVal];
  e2 = expr /. {\[Omega] -> \[Omega] + I\[ThinSpace]\[Delta]};
  NIntegrate[e2, {\[Omega], -\[Infinity], \[Infinity]}, WorkingPrecision -> prec]
];


DoubleContourIntegralBetter[expr_, LambdaNumVal_, BetaNumVal_, prec_:16] := Module[{e2, \[Delta]},
  \[Delta]  = 0.5*Min[LambdaNumVal,  2 \[Pi] /BetaNumVal];
  e2 = expr /. {
    \[Omega]1 -> \[Omega]1 + I\[ThinSpace]\[Delta],
    \[Omega]2 -> \[Omega]2 + 0.5\[ThinSpace]I\[ThinSpace]\[Delta]
  };
  NIntegrate[
    e2,
    {\[Omega]2, -\[Infinity], \[Infinity]},
    {\[Omega]1, -\[Infinity], \[Infinity]},
    WorkingPrecision -> prec
  ]
];


modifyMatrix[A_?MatrixQ, c_] := Module[{B = A},
  B[[2, All]] = c\[ThinSpace]B[[4, All]];
  B
];


(* \[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash] *)
(* 2.  ProcessGen with auto\[Hyphen]double detection                  *)
(* \[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash] *)
ClearAll[ProcessGen];
ProcessGen[t : {_List ..}, poleList_:{}, p3Bp1_:0, prec_:16] := Module[
  { nTerms, cons, single, dbl, total},
  (* 2) inspect how many slots each entry has *)
  nTerms = Length@t[[1, 1]];
  (* 3) constant and single always present *)
  cons   = t[[All, All, 1]];
  
  single = ParallelMap[ContourIntegral[#, poleList, prec] &, t[[All, All, 2]], {2}];
  (* 4) only do double if there are at least 4 elements *)
  If[nTerms >= 4,
    dbl   = ParallelMap[DoubleContourIntegral[#, poleList, prec] &, t[[All, All, 3]], {2}];
    total = modifyMatrix[cons + single + dbl,p3Bp1],
    (* else skip double *)
    dbl   = ConstantArray[0, Dimensions[cons]];
    total = cons + single
  ];
  (* 5) return an Association *)
  <|
    "Const"  -> cons,
    "Single" -> single,
    "Double" -> dbl,
    "Total"  -> total
  |>
];

TCLNumericalEval[t : {_List ..}, poleList_:{}, p3Bp1_:0, prec_:16]:= ProcessGen[t , poleList, p3Bp1, prec]["Total"]


(* ::Text:: *)
(*The argument, poleList to the function TCLNumericalEval is an optional list of smallest imaginary pole location that the spectral density may have which you do not want to overshoot while integrating by going about the real axis. For example, in the case of Drude cutoff, there is a pole at I \[CapitalLambda] above the real axis.*)


CleanInt[expr_] := Module[{e1, e2, Result},
	e1 = ReleaseHold[expr //. TwoPtTransRule2];
	e2 = ExprNiceForm[e1];
	Print["Integrand = ", e2];
	Result = TCL4TripleInt[e2, GenHAssump];
	Print["final result = ", Result];
	Result
]

TCL2NuEtaReplace = {\[Eta][T1_] -> Hold[EtaIntegrand[\[Omega],T1]], \[Nu][T1_]-> Hold[NuIntegrand[\[Omega],T1]]};

TCL2TimeIntegral[expr_]:= Module[{e1},
e1 = ReleaseHold[expr //. TCL2NuEtaReplace];
Integrate[e1, {t1, 0,t}, Assumptions -> GenHAssump]
]


(* ::Section:: *)
(*Nom Markovianity Functions*)


(* Density matrix definition*)
den[{x_,y_,z_}]:=0.5 (IdentityMatrix[2]+x PauliMatrix[1]+y PauliMatrix[2]+z PauliMatrix[3])


(*trace distance definition*)
trdist[{x1_,y1_,z1_},{x2_,y2_,z2_}]:=0.5*Re[Sqrt[(x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2]]


(* ::Text:: *)
(*2 v -1 ranges from -1 to 1. This means v ranges from 0 to 1. This means ArcCos[2 v-1] ranges from 0 to Pi*)
(**)
(*u ranges from 0 to 1 and hence 2Pi u ranges from 0 to 2 Pi*)
(**)
(*z = 2 v - 1. Ranges from -1 to 1.*)
(**)
(*Think of 2v - 1 as a *)


(* parametrising the arbitrary state on the bloch sphere*)
xcomp[u_,v_]:=Cos[2Pi u]Sin[ArcCos[2 v-1]]
ycomp[u_,v_]:=Sin[2Pi u]Sin[ArcCos[2v-1]]
zcomp[v_]:= 2v-1


vComp[z_]:= (z + 1)/2
uComp[x_, z_]:= 1/(2Pi) ArcCos[x/Sin[ArcCos[2 vComp[z] - 1]]]


(*Antipodal state for given state *)
xcompa[u_,v_]:=Cos[Pi+2Pi u]Sin[Pi-ArcCos[-1+2 v]]
ycompa[u_,v_]:=Sin[Pi+2Pi u]Sin[Pi-ArcCos[-1+2 v]]
zcompa[v_]:=Cos[Pi-ArcCos[-1+2 v]]


(* \[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash] *)
(* 2.  NumericalIntegralNM with auto\[Hyphen]double detection                  *)
(* \[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash] *)
ClearAll[NumericalIntegralNM];
NumericalIntegralNM[gen : {_List ..}, lambdaVal_, TVal_, prec_:16] := Module[
  {t, nTerms, cons, single, dbl, total},
  (* 1) substitute numerical values *)
  t      = Map[NumReplace, gen, {2}] //. {\[CapitalLambda] -> lambdaVal, \[Beta] -> 1/TVal};
  
  (* 2) inspect how many slots each entry has *)
  nTerms = Length@t[[1, 1]];
  (* 3) constant and single always present *)
  cons   = t[[All, All, 1]];
  single = Map[ContourIntegralBetter[#, lambdaVal, 1/TVal, prec] & , t[[All, All, 2]], {2}];
  (* 4) only do double if there are at least 4 elements *)
  If[nTerms >= 4,
    dbl   = Map[DoubleContourIntegralBetter[#, lambdaVal, 1/TVal, prec] &, t[[All, All, 3]], {2}];
    total = cons + single + dbl,
    (* else skip double *)
    dbl   = ConstantArray[0, Dimensions[cons]];
    total = cons + single
  ];
  (* 5) return an Association *)
	total
];


SSCalcUV[Ltotalff_]:= Module[{zeroEigenstates,u,v},
zeroEigenstates = Re[ NullSpace[Ltotalff][[1]]/NullSpace[Ltotalff][[1,1]]];
v = vComp[zeroEigenstates[[4]]];
u = uComp[zeroEigenstates[[2]], zeroEigenstates[[4]]];
{{u},{v},{u + 1/2},{1 - 2 v}}
]


TrDist[x_,y_]:= Re[0.5 Sqrt[Total[(x - y)^2]]]

StateEvolve[t_, vin_, L2_]:= Re[MatrixExp[L2 t] . vin]

FindThresholdTime[nul2_, vin_, L2_] := Module[{t = 0},
  While[TrDist[nul2, StateEvolve[t, vin, L2]] >= 0.001,
   t = t + 20
   ];
   
  Return[t]
]

TMaxCalc[L2_]:= Module[{nul1, nul2, antinul},
nul1 = NullSpace[L2][[1]];

nul2 = Re[nul1/nul1[[1]]];
antinul = {1,0,0, - Sign[nul2[[4]]]};

FindThresholdTime[nul2, antinul, L2]
]


NonMarkovianityEval[lambdaVal_,TVal_, numPoint_:2, dt_:0.1, L0gen_:L0gen, p3Bp1_: p3NumVal/p1NumVal]:=Module[{uvList, raw4, TCL4GenFinal, TCL2Final, LtotalModified, Ltotalff, T1, T2, \[CapitalDelta]T,L2gen, tmax},
(*Parameter setup*)

(*Process generator matrix*)

TCL2Final     = NumericalIntegralNM[TCL2Generator, lambdaVal, TVal];

L2gen = L0gen + TCL2Final;

Print["TCL2 generator calculated"];

tmax = TMaxCalc[L2gen];

Print["TMax = ", tmax];

raw4 = NumericalIntegralNM[TCL4Generator, lambdaVal, TVal];

TCL4GenFinal=modifyMatrix[raw4, p3Bp1];

Print["TCL4 generator calculated"];

(*Combine matrices*)
LtotalModified = L0gen + TCL2Final + TCL4GenFinal; (*L0gen here is hard coded*)
Ltotalff = N[LtotalModified];




uvList = Table[RandomReal[{0, 1}, numPoint], {2}];

{NonMarkCalc[Ltotalff, uvList, tmax, dt], NonMarkCalc[L2gen, uvList, tmax, dt]}
]


(*I have commented out parts of this code that did the calcs for general ensemble.
For the sake of speed, I now made restricted it to antipodal.
*)

NonMarkCalc[Ltotalff_, uvList_, tmax_, dt_:0.1]:=Module[{map,vec,vec1,allTrace,backflowIntegrals,maxBackflow,tmin,
tvals,evolved,evolved1,x,y,z,x1,y1,z1, traceSeries, derivatives, positives, uList, vList, u1List,
v1List, T1, T2, \[CapitalDelta]T, currentBackflowIntegral, u, v, u1, v1, maxIndices},



T1 = AbsoluteTime[];
(*Time evolution setup*)
tmin=0;

tvals=Range[tmin,tmax,dt];
map[t_]:=MatrixExp[Ltotalff*t];
(*State vectors*)
vec={{1},{xcomp[u,v]},{ycomp[u,v]},{zcomp[v]}};
vec1={{1},{xcompa[u,v]},{ycompa[u,v]},{zcompa[v]}} (*/. {u -> u1, v -> v1}*);
(*Trace distance calculation*)

(*For all the u, v and t values speciefied, create a list of evolved trace distance.
The earlier calcs followed with ; will not be considered while creating the table. *)

{uList, vList
(*, u1List, v1List*)
} = uvList;


(*Print["uv values = ", {uList, vList(*, u1List, v1List*)}];*)


allTrace=Table[evolved=map[t] . vec;
evolved1=map[t] . vec1;
(*extract the x, y and z components from the evolved vector.*)
{x,y,z}=Flatten[evolved[[2;;4]]]; 
{x1,y1,z1}=Flatten[evolved1[[2;;4]]];
trdist[{x,y,z},{x1,y1,z1}],{u, uList},{v,vList},
(*{u1,u1List},{v1,v1List},*)
{t,tvals}];

(*Backflow analysis*)




backflowIntegrals=Table[traceSeries=allTrace[[i,j(*, i1, j1*)]];
derivatives=Differences[traceSeries];
derivatives=Join[derivatives,{0}]; (*Just add a 0 entry to the end*)
positives=Select[derivatives,#>0&];
Total[positives],{i,Length[allTrace]},{j,Length[allTrace[[1]]]}
(*, {i1,Length[allTrace]},{j1,Length[allTrace[[1]]]}*)
];

maxBackflow=Max[backflowIntegrals];

(* Find the {i, j} position(s) of that maximum value *)
maxIndices = Position[backflowIntegrals, maxBackflow];

T2 = AbsoluteTime[];
\[CapitalDelta]T = (T2 - T1)/60;

Print["non markovianity calc Done"];

{maxBackflow, allTrace[[maxIndices[[1,1]],maxIndices[[1,2]]]]}
]


(* ::Chapter:: *)
(*Integral results*)


(* ::Section:: *)
(*General J(\[Omega]), special hamiltonian*)


Wiki[KapplerAzz2Int]
Wiki[Azz4GeneralSpecOmFirst]
Wiki[Bz4GeneralSpecOmFirst]


(* ::Text:: *)
(*The extra factor of 2 in the following because later, the vectorization was changed from expectation value of sigma/2 to sigma. Hence, both the vector and generator have to be scaled by a factor of 2.*)


Azz2GenJ = - \[Pi] f[\[CapitalOmega]]/2; (*The factor of half should be here because there is No integral here!*)

Azz4GenJ = 2 ((\[Pi] \[CapitalOmega]^2 f[\[Omega]] f[\[CapitalOmega]])/(4 (-\[Omega]^2+\[CapitalOmega]^2)^2)-(\[Pi] \[CapitalOmega] f[\[Omega]] Derivative[1][f][\[CapitalOmega]])/(8 (-\[Omega]^2+\[CapitalOmega]^2)));

Bz4GenJ = 2 (-((\[Pi] \[Omega] \[CapitalOmega] f[\[CapitalOmega]] J[\[Omega]])/(4 (-\[Omega]+\[CapitalOmega])^2 (\[Omega]+\[CapitalOmega])^2))+(\[Pi] (\[Omega]^2+3 \[CapitalOmega]^2) f[\[Omega]] J[\[CapitalOmega]])/(8 (\[Omega]-\[CapitalOmega])^2 (\[Omega]+\[CapitalOmega])^2)+(\[Pi] \[CapitalOmega] f[\[Omega]] Derivative[1][J][\[CapitalOmega]])/(8 \[Omega]^2-8 \[CapitalOmega]^2));


(* ::Text:: *)
(*Integrals over \[Omega] is assumed while avoiding all poles at the real line by moving above the pole.*)
(**)
(*Also, these integrals have to be done from 0 to + Infinity*)


(* ::Text:: *)
(*Verification part*)


Wiki[Azz4General2Drude]
Wiki[Bz4GeneralVsDrude]


GenHJCoh = (J[\[Omega]] Sin[2 \[Theta]] (-\[CapitalOmega]+\[Omega] Coth[(\[Beta] \[Omega])/2] Tanh[(\[Beta] \[CapitalOmega])/2]))/(2 \[Omega]^3-2 \[Omega] \[CapitalOmega]^2)


SimpHDiagPop = -((\[Omega] \[CapitalOmega] J[\[Omega]])/((\[Omega]-\[CapitalOmega])^2 (\[Omega]+\[CapitalOmega])^2))+(\[Beta] \[CapitalOmega] Coth[(\[Beta] \[Omega])/2] J[\[Omega]] Sech[(\[Beta] \[CapitalOmega])/2]^2)/(4 \[Omega]^2-4 \[CapitalOmega]^2)+((\[Omega]^2+\[CapitalOmega]^2) Coth[(\[Beta] \[Omega])/2] J[\[Omega]] Tanh[(\[Beta] \[CapitalOmega])/2])/(2 (\[Omega]-\[CapitalOmega])^2 (\[Omega]+\[CapitalOmega])^2)


SigmaXSpinBosonMFGS = 2*\[Lambda]^2*Sin[2*\[Theta]]/\[Omega]q * 
  (
    Tanh[\[Beta]*\[Omega]q/2] * J[\[Omega]] * Coth[\[Beta]*\[Omega]/2] * \[Omega]q/(\[Omega]^2 - \[Omega]q^2) 
    - J[\[Omega]] * \[Omega]/(\[Omega]^2 - \[Omega]q^2) + J[\[Omega]]/\[Omega]
  )
  
(*eqn c4 of Cresser anders*)


Wiki[GenHCoherenceCalcs]


O2GenHJGenerator = {{0,0,0,0},{-(1/2) \[Pi] \[Lambda]^2 Cos[\[Theta]] J[\[CapitalOmega]] Sin[\[Theta]],-\[Pi] \[Lambda]^2 Cos[\[Theta]]^2 f[0],-\[CapitalOmega],-(1/2) \[Pi] \[Lambda]^2 Cos[\[Theta]] f[\[CapitalOmega]] Sin[\[Theta]]},{1/4 \[Lambda]^2 ((2 \[CapitalOmega]^2 J[\[Omega]])/(\[Omega]^3-\[Omega] \[CapitalOmega]^2)+I \[Pi] J[\[CapitalOmega]]) Sin[2 \[Theta]],\[CapitalOmega]+1/2 \[Lambda]^2 (-I \[Pi] f[\[CapitalOmega]]+(2 \[CapitalOmega] Coth[(\[Beta] \[Omega])/2] J[\[Omega]])/(-\[Omega]^2+\[CapitalOmega]^2)) Sin[\[Theta]]^2,-(1/2) \[Pi] \[Lambda]^2 (2 Cos[\[Theta]]^2 f[0]+f[\[CapitalOmega]] Sin[\[Theta]]^2),1/4 \[Lambda]^2 (I \[Pi] f[\[CapitalOmega]]+(2 \[CapitalOmega] Coth[(\[Beta] \[Omega])/2] J[\[Omega]])/(\[Omega]^2-\[CapitalOmega]^2)) Sin[2 \[Theta]]},{-(1/2) \[Pi] \[Lambda]^2 J[\[CapitalOmega]] Sin[\[Theta]]^2,-\[Pi] \[Lambda]^2 Cos[\[Theta]] f[0] Sin[\[Theta]],0,-(1/2) \[Pi] \[Lambda]^2 f[\[CapitalOmega]] Sin[\[Theta]]^2}};
O2GenHJGenerator // MatrixForm;


(* ::Chapter:: *)
(*Cresser Anders MFGS results*)


(*eqn c5 of Cresser anders*)
SigmaZSpinBosonMFGS = 2*\[Lambda]^2*Sin[\[Theta]]^2 *
  (
    Tanh[\[Beta]*\[Omega]q/2] * J[\[Omega]] * Coth[\[Beta]*\[Omega]/2] * (\[Omega]^2 + \[Omega]q^2)/((\[Omega]^2 - \[Omega]q^2)^2)
    - J[\[Omega]] * (2*\[Omega]*\[Omega]q)/((\[Omega]^2 - \[Omega]q^2)^2)
    + (\[Beta]/2) * Sech[\[Beta]*\[Omega]q/2]^2 * J[\[Omega]] * Coth[\[Beta]*\[Omega]/2] * \[Omega]q/(\[Omega]^2 - \[Omega]q^2)
  )


(*eqn c4 of Cresser anders*)
SigmaXSpinBosonMFGS = 2 Sin[2 \[Theta]]/\[CapitalOmega] J[\[Omega]](Tanh[\[Beta] \[CapitalOmega]/2] Coth[\[Beta] \[Omega]/2] \[CapitalOmega]/(\[Omega]^2 - \[CapitalOmega]^2) - \[Omega]/(\[Omega]^2 - \[CapitalOmega]^2) + 1/\[Omega])


(* ::Chapter:: *)
(*Rough*)


(*Why are these functions here and not in the DrudeCutoffFunctions?
I need to move it there! Later!
*)


(*Too many simplify operations in this code may be slowing it down. But I am afraid to remove them.*)


(*The following code, LargeTLimit, was not required! Simply taking the limit would have worked and would also had been robust!!!!! *)

OneExprOp[term_, Assump_] := Module[{exponents, exponent, exponentSimplified, coeffT, realPart},
  (* Extract exponents from all exponential terms. These assumes that the exponents are only in the numerator? *)
  exponents = Cases[term, E^(arg_) :> arg, {0, Infinity}];
  
  (* If no exponentials are found, return the term as is *)
  If[exponents === {},
    Return[term],
    
    (* Simplify the total exponent *)
    (
    exponent = Total[exponents];
    exponentSimplified = Simplify[exponent, Assumptions -> Assump];
    
    (* Extract the coefficient of t *)
    coeffT = Simplify[Coefficient[exponentSimplified, t]];
    
    (* Compute the real part of the coefficient *)
    realPart = Simplify[Re[coeffT], Assumptions -> Assump];
    
    (* Check if the real part is negative under the assumptions *)
    If[Simplify[realPart < 0, Assumptions -> Assump],
      (* If negative, set the term to zero *)
      (Print["Setting ",Exp[exponentSimplified], " to zero"];
      Return[0];);,
      (Print["did not set ", Exp[exponentSimplified], " to zero"];
      Return[term];);,
      Return[term];
      (*Simplify[term]*)
      ]);,
      Return[term];
    ];
];

LargeTLimit[expr_,Assump_] := Module[{TermList, processedTerms},
	TermList = TermsAsList[ExprNiceForm[expr]]; (*ExprNiceForm gets rid of trignometric terms, which will destroy this function.*)
	processedTerms = OneExprOp[#, Assump] & /@ TermList;
	Total[processedTerms]
]


(*SkipZero is now legacy. I don't trigger it from AllResCareful function*)

OneResExtract[expr_, x_, i_, PoleListX_, QtrPoleList_:{}, z_:z1, SkipZero_:False] := Module[
{ResX, ExtractedTerm, XBulk, TotalXInregral},
	
	ZeroSelect[m_]:= If[m==0, 1, 0];
	
	(*Fix this. R and R2 functions are unnecessary and ugly.*)
	
	BulkPole[term_, var_, m_]:= ZeroSelect[m] term //. {var -> z};
	
	PoleSign[m_] := If[m == 0, 1, Sign[m]];
	
	ExtractedTerm = PowIter2[expr, x, i]; (*extract the exp(I m x t) term from expression*);
	
	MyPrint["For variable = ", x, ", m = ",i];
	MyPrint["Osci Expression = ", ExtractedTerm];
	
	If[And[SkipZero, i ==0], Return[ExtractedTerm], (ResX = \[Pi] I PoleSign[i] Rsidu[ExtractedTerm, x, PoleListX]; (*Do re4sidue calcs on this term and add a sign. For no exponential term, positive sign.*)
	
	If[Length[QtrPoleList] >0, ResX = ResX + \[Pi]/2 I PoleSign[i] Rsidu[ExtractedTerm, x, QtrPoleList], ResX = ResX];
	
	XBulk = BulkPole[ExtractedTerm, x, i ]; 
	(*Bulk pole code only gets activated when i = 0. i=0 means no exponential oscillation for that variable.
	The idea is that when there is oscillating exponential, there are no bulk poles because of decaying effect.
	Note that for negative i, we complete the contour from below.
	In either case, only the poles at the real line contribute.
	For i=0, bulk contri can't be ignored.
	
	For i=0 case, this code replaces \[Omega]1 -> z and \[Omega]2 -> z1.
	
	
	Also, bulk pole lives in upper plane.*)
	
	TotalXInregral = ResX + XBulk;
	
	Return[TotalXInregral])
	];
];



(*This resloop code is better because it can help debug*)
ResLoop[expr_, xval_, PoleListX_, QtrPoleList_:{}, z_:z1, SkipZero_:False]:= Module[{PoleSum, i, PartSum, ExpandedExpr},
	ExpandedExpr = ExprNiceForm[expr];
	MyPrint["Nice form Expr = ", ExpandedExpr];
	PoleSum = 0;
	For[i = -1, i < 2, i++, PartSum = OneResExtract[ExpandedExpr, xval, i, PoleListX, QtrPoleList, z, SkipZero];
	PoleSum = PoleSum + PartSum; 
	];
	Divergence = O2Res[expr, xval, PoleListX];
	{PoleSum, Divergence}
]



(*Legacy code. Right now using AllResCareful*)
AllRes[expr_]:= Module[{PoleListX, PoleListY, om1Poles, om2Poles, startTime, result, elapsedTime, endTime, Divergence1,Divergence2, Result},
	startTime = AbsoluteTime[];
	(*It is important that this z be supplied to the subsequent functions because
	it will be later used here for further calcs.*)
	PoleListX = FindDenSolutions[expr, \[Omega]1, z]; (*Specifying the z here does not make sense.*)
	
	{om1Poles, Divergence1} = ResLoop[expr, \[Omega]1, PoleListX, {}, z];
	
	PoleListY = FindDenSolutions[om1Poles, \[Omega]2, z];
	
	{om2Poles, Divergence2} = ResLoop[om1Poles, \[Omega]2, PoleListY, {}, z];
	MyPrint["Final Result = ", om2Poles];
	
	(*
	RealPoleList = FindDenSolutions[om2Poles, z];
	
	result = ConditionalOp[Expand[om2Poles], {z}, RealResSubtract, RealPoleList];
	*)
	
	MyPrint["As a Real integral = ", result];
	
	Print["Total Divergence = ", Divergence1 + Divergence2];	

    Result = FullSimplify[om2Poles //. z-> \[Omega]];
    
    endTime = AbsoluteTime[];
    elapsedTime = endTime - startTime;
  
    Print["Time taken to run AllRes: ", NumberForm[elapsedTime, {5, 3}], " seconds"];
    
    Result
]



(* ::Section:: *)
(*Debugging code*)


(* ::Text:: *)
(*See*)


Wiki[FunctionWithinFunctionMathematica]
