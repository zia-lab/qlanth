BeginPackage["qonstants`"];

(* Physical Constants*)
bohrRadius = 5.29177210903 * 10^-9;
ee = 1.602176634 * 10^-19;

(* Spectroscopic niceties*)
theLanthanides = {"Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb"};
theActinides   = {"Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"};
theTrivalents  = {"Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm"};
specAlphabet   = "SPDFGHIKLMNOQRTUV"

\[Mu]0  = 4 \[Pi]*10^-7; (* magnetic permeability in vacuum in SI *)
hPlanck = 6.626*10^-34;  (* Planck's constant in SI *)
\[Mu]B  = 9.274*10^-24;  (* Bohr magneton in SI *)
me      = 9.109*10^-31;  (* electron mass in SI *)
cLight  = 2.998*10^8;    (* speed of light in SI *)
eCharge = 1.602*10^-19;  (* elementary charge in SI *)
\[Alpha]Fine = 1/137.036; (* fine structure constant in SI *)
bohrRadius   = 5.29177210903*10^-11; (* Bohr radius in SI *)

EndPackage[];