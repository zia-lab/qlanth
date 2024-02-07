BeginPackage["qonstants`"];

(* Physical Constants*)
bohrRadius = 5.29177210903 * 10^-9;
ee = 1.602176634 * 10^-19;

(* Spectroscopic niceties*)
theLanthanides = {"Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb"};
theActinides   = {"Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"};
theTrivalents  = {"Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm"};
specAlphabet   = "SPDFGHIKLMNOQRTUV"

(* SI *)
\[Mu]0  = 4 \[Pi]*10^-7; (* magnetic permeability in vacuum in SI *)
hPlanck = 6.62607*10^-34;  (* Planck's constant in SI *)
\[Mu]B  = 9.27401*10^-24;  (* Bohr magneton in SI *)
me      = 9.10938*10^-31;  (* electron mass in SI *)
cLight  = 2.99792*10^8;    (* speed of light in SI *)
eCharge = 1.60218*10^-19;  (* elementary charge in SI *)
\[Alpha]Fine = 1/137.036; (* fine structure constant in SI *)
bohrRadius   = 5.29177*10^-11; (* Bohr radius in SI *)

(* Hartree atomic units *)
hPlanckHartree = 2 \[Pi]; (* Planck's constant in Hartree *)
meHartree      = 1;       (* electron mass in Hartree *)
cLightHartree  = 137.036; (* speed of light in Hartree *)
eChargeHartree = 1;       (* elementary charge in Hartree *)
\[Mu]0Hartree  = \[Alpha]Fine^2; (* magnetic permeability in vacuum in Hartree *)

EndPackage[];