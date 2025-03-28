(* -----------------------------------------------------------------
+------------------------------------------------------------------+
|                                                                  |
|                      __                     __     __            |
|            ____ _   / /  ____ _   ____     / /_   / /_           |
|           / __ `/  / /  / __ `/  / __ \   / __/  / __ \          |
|          / /_/ /  / /  / /_/ /  / / / /  / /_   / / / /          |
|          \__, /  /_/   \__, _/  /_/ /_/   \__/  /_/ /_/          |
|            /_/                                                   |
|                                                                  |
|                                                                  |
+------------------------------------------------------------------+
This code was initially authored by Christopher M. Dodson and Rashid
Zia, and then rewritten and expanded by Juan David Lizarazo Ferro in
the  years  2022-2024  under  the advisory of Dr. Rashid Zia. It has
also  benefited  from  the  discussions  with  Tharnier  Puel at the
University of Iowa.

It   grew   out  of  a  collaboration  sponsored  by  the  NSF  (NSF
DMR-1922025)   between  the  groups  of  Dr.  Rashid  Zia  at  Brown
University,  the Quantum Engineering Laboratory at the University of
Pennsylvania  led  by  Dr. Lee Bassett, and the group of Dr. Michael
Flatté at the University of Iowa.

It   uses  an  effective  Hamiltonian  to  describe  the  electronic
structure of lanthanide ions in crystals. This effective Hamiltonian
includes  terms representing the following interactions/relativistic
corrections: spin-orbit, electrostatic repulsion, spin-spin, crystal
field, and spin-other-orbit.

The  Hilbert  space used in this effective Hamiltonian is limited to
single   f^n   configurations.   The   inaccuracy   of  this  single
configuration  description is partially compensated by the inclusion
of  configuration  interaction  terms as parametrized by the Casimir
operators  of  SO(3),  G(2),  and SO(7), and by three-body effective
operators ti.

The  parameters  included  in  this  model  are listed in the string
paramAtlas.

The  notebook  "qlanth.nb"  contains  a  gallery  with  many  of the
functions included in this module with some simple use cases.

The notebook "The Lanthanides in LaF3.nb" is an example in which the
results from this code are compared against the published results by
Carnall  et.  al for the energy levels of lanthinde ions in crystals
of lanthanum trifluoride.

VERSION: MARCH 2025

REFERENCES:

+ Condon, E U, and G Shortley. "The Theory of Atomic Spectra." 1935.

+ Racah,  Giulio.   "Theory of Complex Spectra. II." Physical Review
62,      no.      9-10      (November      1,     1942):     438-62.
https://doi.org/10.1103/PhysRev.62.438.

+ Racah,   Giulio. "Theory of Complex Spectra. III." Physical Review
63,       no.       9-10       (May      1,      1943):      367-82.
https://doi.org/10.1103/PhysRev.63.367.

+ Judd,  B. R.  "Optical Absorption Intensities of Rare-Earth Ions."
Physical   Review   127,   no.   3   (August   1,   1962):   750-61.
https://doi.org/10.1103/PhysRev.127.750.

+ Ofelt,  GS.   "Intensities of Crystal Spectra of Rare-Earth Ions."
The Journal of Chemical Physics 37, no. 3 (1962): 511-20.

+ Rajnak,  K, and BG Wybourne. "Configuration Interaction Effects in
l^N  Configurations."  Physical  Review  132,  no.  1  (1963):  280.
https://doi.org/10.1103/PhysRev.132.280.

+ Nielson,   C. W., and George F Koster. "Spectroscopic Coefficients
for the p^n, d^n, and f^n Configurations", 1963.

+ Wybourne, Brian. "Spectroscopic Properties of Rare Earths." 1965.

+ Carnall,  W To, PR Fields, and BG Wybourne. "Spectral Intensities
of  the  Trivalent  Lanthanides  and Actinides in Solution. I. Pr3+,
Nd3+, Er3+, Tm3+, and Yb3+." The Journal of Chemical Physics 42, no.
11 (1965): 3797-3806. 

+ Judd,   BR.  "Three-Particle  Operators for Equivalent Electrons."
Physical       Review      141,      no.      1      (1966):      4.
https://doi.org/10.1103/PhysRev.141.4.

+ Judd, BR, HM  Crosswhite,  and  Hannah  Crosswhite.  "Intra-Atomic
Magnetic  Interactions  for f Electrons." Physical Review 169, no. 1
(1968): 130. https://doi.org/10.1103/PhysRev.169.130.

+ (TASS)  Cowan,  Robert  Duane. "The Theory of Atomic Structure and
Spectra." Los  Alamos  Series  in  Basic  and  Applied  Sciences  3.
Berkeley: University of California Press, 1981.

+ Judd,   BR,  and  MA  Suskin.  "Complete  Set of Orthogonal Scalar
Operators  for  the  Configuration  f^3."  JOSA  B  1, no. 2 (1984):
261-65. https://doi.org/10.1364/JOSAB.1.000261.

+ Carnall,   W.  T.,  G.  L.  Goodman, K. Rajnak, and R. S. Rana. "A
Systematic  Analysis  of  the  Spectra of the Lanthanides Doped into
Single  Crystal  LaF3."  The  Journal  of Chemical Physics 90, no. 7
(1989): 3443-57. https://doi.org/10.1063/1.455853.

+ Thorne, Anne, Ulf Litzen, and Sveneric Johansson. "Spectrophysics:
Principles  and  Applications." Springer  Science  & Business Media,
1999.

+ Hansen,   JE,  BR Judd, and Hannah Crosswhite. "Matrix Elements of
Scalar Three-Electron Operators for the Atomic f-Shell." Atomic Data
and    Nuclear    Data    Tables    62,    no.   1   (1996):   1-49.
https://doi.org/10.1006/adnd.1996.0001.

+ Velkov,   Dobromir.  "Multi-Electron  Coefficients  of  Fractional
Parentage  for  the  p,  d,  and f Shells." John Hopkins University,
2000. The B1F_ALL.TXT file is from this thesis.

+  Dodson,  Christopher  M.,  and  Rashid  Zia. "Magnetic Dipole and
Electric  Quadrupole Transitions in the Trivalent Lanthanide Series:
Calculated Emission Rates and Oscillator Strengths." Physical Review
B     86,     no.     12     (September     5,     2012):    125102.
https://doi.org/10.1103/PhysRevB.86.125102.

+  Hehlen,  Markus  P,  Mikhail  G  Brik,  and  Karl W Kramer. "50th
Anniversary  of  the Judd-Ofelt Theory: An Experimentalist's View of
the  Formalism  and  Its  Application."  Journal of Luminescence 136
(2013): 221-39.

+ Rudzikas, Zenonas. Theoretical Atomic Spectroscopy, 2007.

+ Benelli, Cristiano, and Dante Gatteschi. Introduction to Molecular
Magnetism: From Transition Metals to Lanthanides. John Wiley & Sons,
2015.

+  Newman, DJ, and G Balasubramanian. "Parametrization of Rare-Earth
Ion  Transition  Intensities."  Journal  of  Physics  C: Solid State
Physics 8, no. 1 (1975): 37.

----------------------------------------------------------------- *)

BeginPackage["qlanth`"];
Get[FileNameJoin[{DirectoryName[$InputFileName], "qonstants.m"}]]
Get[FileNameJoin[{DirectoryName[$InputFileName], "misc.m"}]]

paramAtlas = "
E0: linear combination of F_k, see eqn. (2-80) in Wybourne 1965
E1: linear combination of F_k, see eqn. (2-80) in Wybourne 1965
E2: linear combination of F_k, see eqn. (2-80) in Wybourne 1965
E3: linear combination of F_k, see eqn. (2-80) in Wybourne 1965

\[Zeta]: spin-orbit strength parameter.

F0: Direct Slater integral F^0, produces an overall shift of all energy levels.
F2: Direct Slater integral F^2
F4: Direct Slater integral F^4, possibly constrained by ratio to F^2
F6: Direct Slater integral F^6, possibly constrained by ratio to F^2

M0: 0th Marvin integral
M2: 2nd Marvin integral
M4: 4th Marvin integral
\[Sigma]SS: spin-spin override, if 0 spin-spin is omitted, if 1 then spin-spin is included

T2:  three-body effective operator parameter T^2  (non-orthogonal)
T2p: three-body effective operator parameter T^2' (orthogonalized T2)
T3:  three-body effective operator parameter T^3
T4:  three-body effective operator parameter T^4
T6:  three-body effective operator parameter T^6
T7:  three-body effective operator parameter T^7
T8:  three-body effective operator parameter T^8

T11p: three-body effective operator parameter T^11' (orthogonalized T11)
T12:  three-body effective operator parameter T^12
T14:  three-body effective operator parameter T^14
T15:  three-body effective operator parameter T^15
T16:  three-body effective operator parameter T^16
T17:  three-body effective operator parameter T^17
T18:  three-body effective operator parameter T^18
T19:  three-body effective operator parameter T^19

P0: pseudo-magnetic parameter P^0
P2: pseudo-magnetic parameter P^2
P4: pseudo-magnetic parameter P^4
P6: pseudo-magnetic parameter P^6

gs: electronic gyromagnetic ratio

\[Alpha]: Trees' parameter \[Alpha] describing configuration interaction via the Casimir operator of SO(3)
\[Beta]: Trees' parameter \[Beta] describing configuration interaction via the Casimir operator of G(2)
\[Gamma]: Trees' parameter \[Gamma] describing configuration interaction via the Casimir operator of SO(7)

B02: crystal field parameter B_0^2 (real)
B04: crystal field parameter B_0^4 (real)
B06: crystal field parameter B_0^6 (real)
B12: crystal field parameter B_1^2 (real)
B14: crystal field parameter B_1^4 (real)

B16: crystal field parameter B_1^6 (real)
B22: crystal field parameter B_2^2 (real)
B24: crystal field parameter B_2^4 (real)
B26: crystal field parameter B_2^6 (real)
B34: crystal field parameter B_3^4 (real)

B36: crystal field parameter B_3^6 (real)
B44: crystal field parameter B_4^4 (real)
B46: crystal field parameter B_4^6 (real)
B56: crystal field parameter B_5^6 (real)
B66: crystal field parameter B_6^6 (real)

S12: crystal field parameter S_1^2 (real)
S14: crystal field parameter S_1^4 (real)
S16: crystal field parameter S_1^6 (real)
S22: crystal field parameter S_2^2 (real)

S24: crystal field parameter S_2^4 (real)
S26: crystal field parameter S_2^6 (real)
S34: crystal field parameter S_3^4 (real)
S36: crystal field parameter S_3^6 (real)

S44: crystal field parameter S_4^4 (real)
S46: crystal field parameter S_4^6 (real)
S56: crystal field parameter S_5^6 (real)
S66: crystal field parameter S_6^6 (real)

\[Epsilon]: ground level baseline shift
t2Switch: controls the usage of the t2 operator beyond f7 (1 for f7 or below, 0 for f8 or above)
wChErrA: If 1 then the type-A errors in Chen are used, if 0 then not.
wChErrB: If 1 then the type-B errors in Chen are used, if 0 then not.

Bx: x component of external magnetic field (in T)
By: y component of external magnetic field (in T)
Bz: z component of external magnetic field (in T)

\[CapitalOmega]2: Judd-Ofelt intensity parameter k=2 (in cm^2)
\[CapitalOmega]4: Judd-Ofelt intensity parameter k=4 (in cm^2)
\[CapitalOmega]6: Judd-Ofelt intensity parameter k=6 (in cm^2)

nE: number of electrons in a configuration

E0p: orthogonalized E0
E1p: orthogonalized E1
E2p: orthogonalized E2
E3p: orthogonalized E3
\[Alpha]p: orthogonalized \[Alpha]
\[Beta]p: orthogonalized \[Beta]
\[Gamma]p: orthogonalized \[Gamma]
";
paramSymbols   = StringSplit[paramAtlas, "\n"];
paramSymbols   = Select[paramSymbols, # != ""& ];
paramSymbols   = ToExpression[StringSplit[#, ":"][[1]]] & /@ paramSymbols;
Protect /@ paramSymbols;

(* Parameter families *)
Unprotect[racahSymbols, chenSymbols, slaterSymbols, controlSymbols, cfSymbols, TSymbols, pseudoMagneticSymbols, marvinSymbols, casimirSymbols, magneticSymbols, juddOfeltIntensitySymbols,mostlyOrthogonalSymbols];
racahSymbols   = {E0, E1, E2, E3};
chenSymbols    = {wChErrA, wChErrB};
slaterSymbols  = {F0, F2, F4, F6};
controlSymbols = {t2Switch, \[Sigma]SS};
cfSymbols      = {B02, B04, B06, B12, B14, B16, B22, B24, B26, B34, B36, 
                  B44, B46, B56, B66, 
                  S12, S14, S16, S22, S24, S26, S34, S36, S44, S46, S56, S66};
TSymbols       = {T2, T2p, T3, T4, T6, T7, T8, T11p, T12, T14, T15, T16, T17, T18, T19};
pseudoMagneticSymbols = {P0, P2, P4, P6};
marvinSymbols         = {M0, M2, M4};
magneticSymbols       = {Bx, By, Bz, gs, \[Zeta]};
casimirSymbols        = {\[Alpha], \[Beta], \[Gamma]};
juddOfeltIntensitySymbols = {\[CapitalOmega]2, \[CapitalOmega]4, \[CapitalOmega]6};
mostlyOrthogonalSymbols = Join[{E0p, E1p, E2p, E3p, \[Alpha]p, \[Beta]p, \[Gamma]p}, 
  {T2p, T3, T4, T6, T7, T8, T11p, T12, T14, T15, T16, T17, T18, T19},
  cfSymbols];

paramFamilies         = Hold[{racahSymbols, chenSymbols, slaterSymbols, controlSymbols, cfSymbols, TSymbols, pseudoMagneticSymbols, marvinSymbols, casimirSymbols, magneticSymbols, juddOfeltIntensitySymbols, mostlyOrthogonalSymbols}];
ReleaseHold[Protect /@ paramFamilies];
crystalGroups = {"C1","Ci","C2","Cs","C2h","D2","C2v","D2h","C4","S4","C4h","D4","C4v","D3d","D4h","C3","C3i","D3","C3v","D3d","C6","C3h","C6h","D6","C6v","D3h","D6h","T","Th","O","Td","Oh"};

(* Parameter usage *)
paramLines = Select[StringSplit[paramAtlas, "\n"], # != "" &];
usageTemplate = StringTemplate["`paramSymbol`::usage=\"`paramSymbol` : `paramUsage`\";"];
Do[(
  {paramString, paramUsage} = StringSplit[paramLine, ":"];
  paramUsage                = StringTrim[paramUsage];
  expressionString          = usageTemplate[<|"paramSymbol" -> paramString, "paramUsage" -> paramUsage|>];
  ToExpression[usageTemplate[<|"paramSymbol" -> paramString, "paramUsage" -> paramUsage|>]]
), 
{paramLine, paramLines}
];

AllowedJ;
AllowedMforJ;
AllowedNKSLJMforJMTerms;
AllowedNKSLJMforJTerms;
AllowedNKSLJTerms;

AllowedNKSLTerms;
AllowedNKSLforJTerms;
AllowedSLJMTerms;
AllowedSLJTerms;
AllowedSLTerms;

AngularMomentumMatrices;
BasisLSJ;
BasisLSJMJ;
Bqk;
CFP;

CFPAssoc;
CFPTable;
CFPTerms;
Carnall;
CasimirG2;

CasimirSO3;
CasimirSO7;
Ck;
Cqk;
CrystalField;
CrystalFieldForm;

Dk;
EigenLever;
EffectiveHamiltonian;
EffectiveElectricDipole;
Electrostatic;
ElectrostaticConfigInteraction;

ElectrostaticTable;
EnergyLevelDiagram;
EnergyStates;
EtoF;
ExportMZip;

FindNKLSTerm;
FindSL;
FreeHam;
FreeIonTable;
ElectricDipLineStrength;
FromArrayToTable;
FtoE;

GG2U;
GSO7W;
GenerateCFP;
GenerateCFPAssoc;
GenerateCFPTable;

GenerateCrystalFieldTable;
GenerateElectrostaticTable;
GenerateFreeIonTable;
GenerateReducedUkTable;
GenerateReducedV1kTable;
GenerateSOOandECSOLSTable;

GenerateSOOandECSOTable;
GenerateSpinOrbitTable;
GenerateSpinSpinTable;
GenerateT22Table;
GenerateThreeBodyTables;

GroundMagDipoleOscillatorStrength;
HamiltonianForm;

HamiltonianMatrixPlot;
HoleElectronConjugation;
ImportMZip;
IonSolver;
JJBlockMagDip;

JJBlockMatrix;
JJBlockMatrixFileName;
JJBlockMatrixTable;
JuddCFPPhase;
JuddOfeltUkSquared;
LabeledGrid;
LandeFactor;

LevelElecDipoleOscillatorStrength;
LevelJJBlockMagDipole;
LevelMagDipoleLineStrength;
LevelMagDipoleMatrixAssembly;
LevelMagDipoleOscillatorStrength;

LevelMagDipoleSpontaneousDecayRates;
LevelSimplerEffectiveHamiltonian;
LevelSolver;
LoadAll;

LoadCFP;
LoadCarnall;
LoadChenDeltas;
LoadElectrostatic;
LoadFreeIon;

LoadLaF3Parameters;
LoadLiYF4Parameters;
LoadSOOandECSO;
LoadSOOandECSOLS;
LoadSpinOrbit;
LoadSpinSpin;

LoadT22;
LoadTermLabels;
LoadThreeBody;

LoadUk;
LoadV1k;
MagDipLineStrength;
MagDipoleMatrixAssembly;
MagDipoleRates;

MagneticInteractions;
MapToSparseArray;
MaxJ;
MinJ;
NKCFPPhase;

ParamPad;
ParseBenelli2015;
ParseStates;
ParseStatesByNumBasisVecs;
ParseStatesByProbabilitySum;

ParseTermLabels;
Phaser;
PrettySaundersSL;
PrettySaundersSLJ;
PrettySaundersSLJmJ;

PrintL;
PrintSLJ;
PrintSLJM;
ReducedSOOandECSOinf2;
ReducedSOOandECSOinfn;

ReducedT11inf2;
ReducedT22inf2;
ReducedT22infn;
ReducedUk;
ReducedUkTable;

ReducedV1k;
ReducedV1kTable;
Reducedt11inf2;
ReplaceInSparseArray;
ScalarLSJMFromLS;
SOOandECSO;
SOOandECSOLSTable;

SOOandECSOTable;
Seniority;
ShiftedLevels;
SimplerEffectiveHamiltonian;

SixJay;
SpinOrbit;
SpinOrbitTable;
SpinSpin;
SpinSpinTable;

Sqk;
SquarePrimeToNormal;
TPO;
T22Table;
TabulateJJBlockMagDipTable;
TabulateJJBlockMatrixTable;

TabulateManyJJBlockMagDipTables;
TabulateManyJJBlockMatrixTables;
ThreeBodyTable;
ThreeBodyTables;
ThreeJay;
UkOperator;

chenDeltas;
fnTermLabels;
fsubk;

fsupk;
moduleDir;
symbolicHamiltonians;

(* this selects the function that is applied to calculated matrix elements which helps keep down the complexity of the resulting algebraic expressions *)
SimplifyFun = Expand;

Begin["`Private`"]

  ListRepeater;
  TotalCFIters;
  
  moduleDir = ParentDirectory[DirectoryName[$InputFileName]];
  frontEndAvailable = (Head[$FrontEnd] === FrontEndObject);

  (* ########################################################### *)
  (* ########################## MISC ########################### *)

  TPO::usage = "TPO[x, y, ...] gives the product of 2x+1, 2y+1, ...";
  TPO[args__] := Times @@ ((2*# + 1) & /@ {args});

  Phaser::usage = "Phaser[x] gives (-1)^x.";
  Phaser[exponent_] := ((-1)^exponent);

  TriangleCondition::usage = "TriangleCondition[a, b, c] evaluates the triangle condition on a, b, and c.";
  TriangleCondition[a_, b_, c_] := (Abs[b - c] <= a <= (b + c));

  TriangleAndSumCondition::usage = "TriangleAndSumCondition[a, b, c] evaluates the joint satisfaction of the triangle and sum conditions.";
  TriangleAndSumCondition[a_, b_, c_] := (
    And[
      Abs[b - c] <= a <= (b + c),
      IntegerQ[a + b + c]
      ]
    );

  SquarePrimeToNormal::usage = "SquarePrimeToNormal[squarePrime] evaluates the standard representation of a number from the squared prime representation given in the list squarePrime. For squarePrime of the form {c0, c1, c2, c3, ...} this function returns the number c0 * Sqrt[p1^c1 * p2^c2 * p3^c3 * ...] where pi is the ith prime number. Exceptionally some of the ci might be letters in which case they have to be one of \"A\", \"B\", \"C\", \"D\" with them corresponding to 10, 11, 12, and 13, respectively.";
  SquarePrimeToNormal[squarePrime_] :=
  (
    radical = Product[Prime[idx1 - 1] ^ Part[squarePrime, idx1], {idx1, 2, Length[squarePrime]}];
    radical = radical /. {"A" -> 10, "B" -> 11, "C" -> 12, "D" -> 13};
    val = squarePrime[[1]] * Sqrt[radical];
    Return[val];
  );

  ParamPad::usage = "ParamPad[params] takes an association params whose keys are a subset of paramSymbols. The function returns a new association where all the keys not present in paramSymbols, will now be included in the returned association with their values set to zero.
  The function additionally takes an option \"Print\" that if set to True, will print the symbols that were not present in the given association. The default is True.";
  Options[ParamPad] = {"PrintFun" -> PrintTemporary};
  ParamPad[params_, OptionsPattern[]] := (
    notPresentSymbols = Complement[paramSymbols, Keys[params]];
    PrintFun = OptionValue["PrintFun"];
    PrintFun["Following symbols were not given and are being set to 0: ",
      notPresentSymbols];
    newParams = Transpose[{paramSymbols, ConstantArray[0, Length[paramSymbols]]}];
    newParams = (#[[1]] -> #[[2]]) & /@ newParams;
    newParams = Association[newParams];
    newParams = Join[newParams, params];
    Return[newParams];
    )

  (* ########################################################### *)
  (* #################### Angular Momentum ##################### *)
  
  AngularMomentumMatrices::usage = "AngularMomentumMatrices[j] gives the matrix representation for the angular momentum operators Jx, Jy, and Jz for a given angular momentum j in the basis of eigenvectors of jz. j may be a half-integer or an integer.
  The options are \"Sparse\" which defaults to False and \"Order\" which defaults to \"HighToLow\". 
  The option \"Order\" can be set to \"LowToHigh\" to get the matrices in the order from -jay to jay otherwise they are returned in the order jay to -jay.
  The function returns a list {JxMatrix, JyMatrix, JzMatrix} with the matrix representations for the cartesian components of the angular momentum operator.";
  Options[AngularMomentumMatrices] = {
    "Sparse" -> False,
    "Order" -> "HighToLow"};
  AngularMomentumMatrices[jay_, OptionsPattern[]] := Module[
    {
      JxMatrix, JyMatrix, JzMatrix, 
      JPlusMatrix, JMinusMatrix, 
      m1, m2,
      ArrayInverter
    }, 
    (
      ArrayInverter = #[[-1 ;; 1 ;; -1, -1 ;; 1 ;; -1]] &;
      JPlusMatrix   = Table[
        If[m2 == m1 + 1, 
          Sqrt[(jay - m1) (jay + m1 + 1)], 
          0
          ], 
        {m1, jay, -jay, -1}, 
        {m2, jay, -jay, -1}];
      JMinusMatrix = Table[
        If[m2 == m1 - 1, 
          Sqrt[(jay + m1) (jay - m1 + 1)], 
          0
          ], 
        {m1, jay, -jay, -1}, 
        {m2, jay, -jay, -1}];
      JxMatrix = (JPlusMatrix + JMinusMatrix)/2;
      JyMatrix = (JMinusMatrix - JPlusMatrix)/(2 I);
      JzMatrix = DiagonalMatrix[Table[m, {m, jay, -jay, -1}]];
      If[OptionValue["Sparse"],
        {JxMatrix, JyMatrix, JzMatrix} = SparseArray /@ {JxMatrix, JyMatrix, JzMatrix}
      ];
      If[OptionValue["Order"] == "LowToHigh",
        {JxMatrix, JyMatrix, JzMatrix} = ArrayInverter /@ {JxMatrix, JyMatrix, JzMatrix};
      ];
      Return[{JxMatrix, JyMatrix, JzMatrix}];
    )
  ];
 
  LandeFactor::usage="LandeFactor[J, L, S] gives the Lande factor for a given total angular momentum J, orbital angular momentum L, and spin S.";
  LandeFactor[J_, L_, S_] := (3/2) + (S*(S + 1) - L*(L + 1))/(2*J*(J + 1));
  
  (* #################### Angular Momentum ##################### *)
  (* ########################################################### *)


  (* ########################################################### *)
  (* ###################### Racah Algebra ###################### *)

  ReducedUk::usage = "ReducedUk[n, l, SL, SpLp, k] gives the reduced matrix element of the symmetric unit tensor operator U^(k). See equation 11.53 in TASS.";
  ReducedUk[numE_, l_, SL_, SpLp_, k_] := Module[
    {spin, orbital, Uk, S, L,
    Sp, Lp, Sb, Lb, parentSL,
    cfpSL, cfpSpLp, Ukval,
    SLparents, SLpparents,
    commonParents, phase},
      {spin, orbital} = {1/2, 3};
      {S, L}          = FindSL[SL];
      {Sp, Lp}        = FindSL[SpLp];
      If[Not[S == Sp],
        Return[0]
      ];
      cfpSL      = CFP[{numE, SL}];
      cfpSpLp    = CFP[{numE, SpLp}];
      SLparents  = First /@ Rest[cfpSL];
      SLpparents = First /@ Rest[cfpSpLp];
      commonParents = Intersection[SLparents, SLpparents];
      Uk = Sum[(
        {Sb, Lb} = FindSL[\[Psi]b];
        Phaser[Lb] *
          CFPAssoc[{numE, SL, \[Psi]b}] *
          CFPAssoc[{numE, SpLp, \[Psi]b}] *
          SixJay[{orbital, k, orbital}, {L, Lb, Lp}]
      ),
      {\[Psi]b, commonParents}
      ];
      phase     = Phaser[orbital + L + k];
      prefactor = numE * phase * Sqrt[TPO[L,Lp]];
      Ukval     = prefactor*Uk;
      Return[Ukval];
  ]

  Ck::usage = "Ck[orbital, k] gives the diagonal reduced matrix element <l||C^(k)\[VerticalSeparator]\[VerticalSeparator]l> where the Subscript[C, q]^(k) are renormalized spherical harmonics. See equation 11.23 in TASS with l=l'.";
  Ck[orbital_, k_] := (-1)^orbital * TPO[orbital] * ThreeJay[{orbital, 0}, {k, 0}, {orbital, 0}];

  SixJay::usage = "SixJay[{j1, j2, j3}, {j4, j5, j6}] provides the value for SixJSymbol[{j1, j2, j3}, {j4, j5, j6}] with memorization of computed values and short-circuiting values based on triangle conditions.";
  SixJay[{j1_, j2_, j3_}, {j4_, j5_, j6_}] := (
    sixJayval = Which[
      Not[TriangleAndSumCondition[j1, j2, j3]],
      0,
      Not[TriangleAndSumCondition[j1, j5, j6]],
      0,
      Not[TriangleAndSumCondition[j4, j2, j6]],
      0,
      Not[TriangleAndSumCondition[j4, j5, j3]],
      0,
      True,
      SixJSymbol[{j1, j2, j3}, {j4, j5, j6}]];
    SixJay[{j1, j2, j3}, {j4, j5, j6}] = sixJayval);

  ThreeJay::usage = "ThreeJay[{j1, m1}, {j2, m2}, {j3, m3}] gives the value of the Wigner 3j-symbol and memorizes the computed value.";
  ThreeJay[{j1_, m1_}, {j2_, m2_}, {j3_, m3_}] := (
   threejval = Which[
     Not[(m1 + m2 + m3) == 0],
     0,
     Not[TriangleCondition[j1,j2,j3]],
     0,
     True,
     ThreeJSymbol[{j1, m1}, {j2, m2}, {j3, m3}]
     ];
   ThreeJay[{j1, m1}, {j2, m2}, {j3, m3}] = threejval);

  ReducedV1k::usage = "ReducedV1k[n, SL, SpLp, k] gives the reduced matrix element of the spherical tensor operator V^(1k). See equation 2-101 in Wybourne 1965.";
  ReducedV1k[numE_, SL_, SpLp_, k_] := Module[
    {Vk1, S, L, Sp, Lp,
    Sb, Lb, spin, orbital,
    cfpSL, cfpSpLp, 
    SLparents, SpLpparents,
    commonParents, prefactor},
    (
      {spin, orbital} = {1/2, 3};
      {S, L}        = FindSL[SL];
      {Sp, Lp}      = FindSL[SpLp];
      cfpSL         = CFP[{numE, SL}];
      cfpSpLp       = CFP[{numE, SpLp}];
      SLparents     = First /@ Rest[cfpSL];
      SpLpparents   = First /@ Rest[cfpSpLp];
      commonParents = Intersection[SLparents, SpLpparents];
      Vk1 = Sum[(
          {Sb, Lb} = FindSL[\[Psi]b];
          Phaser[(Sb + Lb + S + L + orbital + k - spin)] *
          CFPAssoc[{numE, SL, \[Psi]b}] *
          CFPAssoc[{numE, SpLp, \[Psi]b}] *
          SixJay[{S, Sp, 1}, {spin, spin, Sb}] *
          SixJay[{L, Lp, k}, {orbital, orbital, Lb}]
        ), 
      {\[Psi]b, commonParents}
      ];
      prefactor = numE * Sqrt[spin * (spin + 1) * TPO[spin, S, L, Sp, Lp] ];
      Return[prefactor * Vk1];
    )
  ];

  GenerateReducedUkTable::usage = "GenerateReducedUkTable[numEmax] can be used to generate the association of reduced matrix elements for the unit tensor operators Uk from f^1 up to f^numEmax. If the option \"Export\" is set to True then the resulting data is saved to ./data/ReducedUkTable.m.";
  Options[GenerateReducedUkTable] = {"Export" -> True, "Progress" -> True};
  GenerateReducedUkTable[numEmax_Integer:7, OptionsPattern[]] := (
    numValues = Total[Length[AllowedNKSLTerms[#]]*Length[AllowedNKSLTerms[#]]&/@Range[1, numEmax]] * 4;
    Echo["Calculating " <> ToString[numValues] <> " values for Uk k=0,2,4,6."];
    counter = 1;
    If[And[OptionValue["Progress"], frontEndAvailable],
    progBar = PrintTemporary[
        Dynamic[Row[{ProgressIndicator[counter, {0, numValues}], " ", 
          counter}]]]
      ];
    ReducedUkTable = Table[
      (
        counter = counter+1;
        {numE, 3, SL, SpLp, k} -> SimplifyFun[ReducedUk[numE, 3, SL, SpLp, k]]
      ),
      {numE, 1, numEmax},
      {SL,   AllowedNKSLTerms[numE]}, 
      {SpLp, AllowedNKSLTerms[numE]}, 
      {k, {0, 2, 4, 6}}
    ];
    ReducedUkTable = Association[Flatten[ReducedUkTable]];
    ReducedUkTableFname = FileNameJoin[{moduleDir, "data", "ReducedUkTable.m"}];
    If[And[OptionValue["Progress"], frontEndAvailable],
      NotebookDelete[progBar]
    ];
    If[OptionValue["Export"],
      (
        Echo["Exporting to file " <> ToString[ReducedUkTableFname]];
        Export[ReducedUkTableFname, ReducedUkTable];
      )
    ];
    Return[ReducedUkTable];
  )

  GenerateReducedV1kTable::usage = "GenerateReducedV1kTable[nmax] calculates values for Vk1 and returns an association where the keys are lists of the form {n, SL, SpLp, 1}. If the option \"Export\" is set to True then the resulting data is saved to ./data/ReducedV1kTable.m.";
  Options[GenerateReducedV1kTable] = {"Export" -> True, "Progress" -> True};
  GenerateReducedV1kTable[numEmax_Integer:7, OptionsPattern[]] := (
    numValues = Total[Length[AllowedNKSLTerms[#]]*Length[AllowedNKSLTerms[#]]&/@Range[1, numEmax]];
    Echo["Calculating " <> ToString[numValues] <> " values for Vk1."];
    counter = 1;
    If[And[OptionValue["Progress"], frontEndAvailable],
    progBar = PrintTemporary[
        Dynamic[Row[{ProgressIndicator[counter, {0, numValues}], " ", 
          counter}]]]
      ];
    ReducedV1kTable = Table[
      (
        counter = counter+1;
        {n, SL, SpLp, 1} -> SimplifyFun[ReducedV1k[n, SL, SpLp, 1]]
      ),
      {n, 1, numEmax},
      {SL, AllowedNKSLTerms[n]}, 
      {SpLp, AllowedNKSLTerms[n]}
    ];
    ReducedV1kTable = Association[ReducedV1kTable];
    If[And[OptionValue["Progress"], frontEndAvailable],
      NotebookDelete[progBar]
    ];
    exportFname = FileNameJoin[{moduleDir, "data", "ReducedV1kTable.m"}];
    If[OptionValue["Export"],
      (
        Echo["Exporting to file " <> ToString[exportFname]];
        Export[exportFname, ReducedV1kTable];
      )
    ];
    Return[ReducedV1kTable];
  )

  (* ###################### Racah Algebra ###################### *)
  (* ########################################################### *)

  (* ########################################################### *)
  (* ###################### Electrostatic ###################### *)

  fsubk::usage = "fsubk[numE, orbital, SL, SLp, k] gives the Slater integral f_k for the given configuration and pair of SL terms. See equation 12.17 in TASS.";
  fsubk[numE_, orbital_, NKSL_, NKSLp_, k_] := Module[
    {terms, S, L, Sp, Lp,
    termsWithSameSpin, SL,
    fsubkVal, spinMultiplicity,
    prefactor, summand1, summand2},
    (
      {S, L}   = FindSL[NKSL];
      {Sp, Lp} = FindSL[NKSLp];
      terms    = AllowedNKSLTerms[numE];
      (* sum for summand1 is over terms with same spin *)
      spinMultiplicity  = 2*S + 1;
      termsWithSameSpin = StringCases[terms, ToString[spinMultiplicity] ~~ __];
      termsWithSameSpin = Flatten[termsWithSameSpin];
      If[Not[{S, L} == {Sp, Lp}],
        Return[0]
      ];
      prefactor = 1/2 * Abs[Ck[orbital, k]]^2;
      summand1 = Sum[(
          ReducedUkTable[{numE, orbital, SL, NKSL,  k}] *
          ReducedUkTable[{numE, orbital, SL, NKSLp, k}]
          ),
        {SL, termsWithSameSpin}
      ];
      summand1 = 1 / TPO[L] * summand1;
      summand2 = (
        KroneckerDelta[NKSL, NKSLp] *
          (numE *(4*orbital + 2 - numE)) /
          ((2*orbital + 1) * (4*orbital + 1))
        );
      fsubkVal = prefactor*(summand1 - summand2);
      Return[fsubkVal];
    )
  ];

  fsupk::usage = "fsupk[numE, orbital, SL, SLp, k] gives the superscripted Slater integral f^k = Subscript[f, k] * Subscript[D, k].";
  fsupk[numE_, orbital_, NKSL_, NKSLp_ ,k_] := (
    Dk[k] * fsubk[numE, orbital, NKSL, NKSLp, k]
    )

  Dk::usage = "D[k] gives the ratio between the super-script and sub-scripted Slater integrals (F^k / F_k). k must be even. See table 6-3 in TASS, and also section 2-7 of Wybourne (1965). See also equation 6.41 in TASS.";
  Dk[k_] := {1, 225, 1089, 184041/25}[[k/2+1]];

  FtoE::usage = "FtoE[F0, F2, F4, F6] calculates the Racah parameters {E0, E1, E2, E3} corresponding to the given Slater integrals.
  See eqn. 2-80 in Wybourne.
  Note that in that equation the subscripted Slater integrals are used but since this function assumes the the input values are superscripted Slater integrals, it is necessary to convert them using Dk.";
  FtoE[F0_, F2_, F4_, F6_] := Module[
    {E0, E1, E2, E3},
    (
      E0 = (F0 - 10*F2/Dk[2] - 33*F4/Dk[4] - 286*F6/Dk[6]);
      E1 = (70*F2/Dk[2] + 231*F4/Dk[4] + 2002*F6/Dk[6])/9;
      E2 = (F2/Dk[2] - 3*F4/Dk[4] + 7*F6/Dk[6])/9;
      E3 = (5*F2/Dk[2] + 6*F4/Dk[4] - 91*F6/Dk[6])/3;
      Return[{E0, E1, E2, E3}];
    )
  ];

  EtoF::usage = "EtoF[E0, E1, E2, E3] calculates the Slater integral parameters {F0, F2, F4, F6} corresponding to the given Racah parameters {E0, E1, E2, E3}. This is the inverse of the FtoE function.";
  EtoF[E0_, E1_, E2_, E3_] := Module[
    {F0, F2, F4, F6},
    (
      F0 = 1/7      (7 E0 + 9 E1);
      F2 = 75/14    (E1 + 143 E2 + 11 E3);
      F4 = 99/7     (E1 - 130 E2 + 4 E3);
      F6 = 5577/350 (E1 + 35 E2 - 7 E3);
      Return[{F0, F2, F4, F6}];
    )
  ];

  Electrostatic::usage = "Electrostatic[{numE, NKSL, NKSLp}] returns the LS reduced matrix element for repulsion matrix element for equivalent electrons. See equation 2-79 in Wybourne (1965). The option \"Coefficients\" can be set to \"Slater\" or \"Racah\". If set to \"Racah\" then E_k parameters and e^k operators are assumed, otherwise the Slater integrals F^k and operators f_k. The default is \"Slater\".";
  Options[Electrostatic] = {"Coefficients" -> "Slater"};
  Electrostatic[{numE_, NKSL_, NKSLp_}, OptionsPattern[]] := Module[
    {fsub0, fsub2, fsub4, fsub6, 
     esub0, esub1, esub2, esub3, 
     fsup0, fsup2, fsup4, fsup6, 
     eMatrixVal, orbital},
    (
      orbital = 3;
      Which[
        OptionValue["Coefficients"] == "Slater",
        (
          fsub0 = fsubk[numE, orbital, NKSL, NKSLp, 0];
          fsub2 = fsubk[numE, orbital, NKSL, NKSLp, 2];
          fsub4 = fsubk[numE, orbital, NKSL, NKSLp, 4];
          fsub6 = fsubk[numE, orbital, NKSL, NKSLp, 6];
          eMatrixVal = fsub0*F0 + fsub2*F2 + fsub4*F4 + fsub6*F6;
        ),
        OptionValue["Coefficients"] == "Racah",
        (
          fsup0 = fsupk[numE, orbital, NKSL, NKSLp, 0];
          fsup2 = fsupk[numE, orbital, NKSL, NKSLp, 2];
          fsup4 = fsupk[numE, orbital, NKSL, NKSLp, 4];
          fsup6 = fsupk[numE, orbital, NKSL, NKSLp, 6];
          esub0 = fsup0;
          esub1 = 9/7*fsup0 +   1/42*fsup2 +   1/77*fsup4 +  1/462*fsup6;
          esub2 =             143/42*fsup2 - 130/77*fsup4 + 35/462*fsup6;
          esub3 =              11/42*fsup2 + 4/77*fsup4   -  7/462*fsup6;
          eMatrixVal = esub0*E0 + esub1*E1 + esub2*E2 + esub3*E3;
        )
      ];
      Return[eMatrixVal];
    )
  ];

  GenerateElectrostaticTable::usage = "GenerateElectrostaticTable[numEmax] can be used to generate the table for the electrostatic interaction from f^1 to f^numEmax. If the option \"Export\" is set to True then the resulting data is saved to ./data/ElectrostaticTable.m.";
  Options[GenerateElectrostaticTable] = {"Export" -> True, "Coefficients" -> "Slater"};
  GenerateElectrostaticTable[numEmax_Integer:7, OptionsPattern[]] := (
    ElectrostaticTable = Table[
      {numE, SL, SpLp} -> SimplifyFun[Electrostatic[{numE, SL, SpLp}, "Coefficients" -> OptionValue["Coefficients"]]], 
      {numE, 1, numEmax},
      {SL, AllowedNKSLTerms[numE]}, 
      {SpLp, AllowedNKSLTerms[numE]}
    ];
    ElectrostaticTable = Association[Flatten[ElectrostaticTable]];
    If[OptionValue["Export"],
      Export[FileNameJoin[{moduleDir, "data", "ElectrostaticTable.m"}],
      ElectrostaticTable];
    ];
    Return[ElectrostaticTable];
  );

  (* ###################### Electrostatic ###################### *)
  (* ########################################################### *)

  (* ########################################################### *)
  (* ########################## Bases ########################## *)

  BasisTableGenerator::usage = "BasisTableGenerator[numE] gives an association whose keys are lists of the form {numE, J} and whose values are lists with elements of the form {LS, J, MJ} representing the elements of the LSJM coupled basis.";
  BasisTableGenerator[numE_] := Module[
    {energyStatesTable, allowedJ, J, Jp}, 
    (
      energyStatesTable = <||>;
      allowedJ = AllowedJ[numE];
      Do[
      (
        energyStatesTable[{numE, J}] = EnergyStates[numE, J];
      ),
      {Jp, allowedJ},
      {J,  allowedJ}];
      Return[energyStatesTable]
    )
  ];

  BasisLSJMJ::usage = "BasisLSJMJ[numE] returns the ordered basis in L-S-J-MJ with the total orbital angular momentum L and total spin angular momentum S coupled together to form J. The function returns a list with each element representing the quantum numbers for each basis vector. Each element is of the form {SL (string in spectroscopic notation),J, MJ}.
  The option \"AsAssociation\" can be set to True to return the basis as an association with the keys corresponding to values of J and the values lists with the corresponding {L, S, J, MJ} list. The default of this option is False.
  ";
  Options[BasisLSJMJ] = {"AsAssociation" -> False};
  BasisLSJMJ[numE_, OptionsPattern[]] := Module[
    {energyStatesTable, basis, idx1},
    (
      energyStatesTable = BasisTableGenerator[numE];
      basis = Table[
        energyStatesTable[{numE, AllowedJ[numE][[idx1]]}],
        {idx1, 1, Length[AllowedJ[numE]]}];
      basis = Flatten[basis, 1];
      If[OptionValue["AsAssociation"],
        (
          Js    = AllowedJ[numE];
          basis = Table[(J -> Select[basis, #[[2]] == J &]), {J, Js}];
          basis = Association[basis];
        )
      ];
      Return[basis];
    )
  ];

  BasisLSJ::usage = "BasisLSJ[numE] returns the level basis LSJ. The function returns a list with each element representing the quantum numbers for each basis vector. Each element is of the form {SL (string in spectroscopic notation), J}.
  The option \"AsAssociation\" can be set to True to return the basis as an association with the keys being the allowed J values. The default is False.
    ";
  Options[BasisLSJ]={"AsAssociation"->False};
  BasisLSJ[numE_,OptionsPattern[]] := Module[
    {Js, basis},
    (
      Js    = AllowedJ[numE];
      basis = BasisLSJMJ[numE,"AsAssociation"->False];
      basis = DeleteDuplicates[{#[[1]],#[[2]]}& /@ basis];
      If[OptionValue["AsAssociation"],
        (
          basis = Association @ Table[(J->Select[basis,#[[2]]==J&]),{J,Js}]
        )
      ];
      Return[basis];
    )
  ];

  (* ########################## Bases ########################## *)
  (* ########################################################### *)

  (* ########################################################### *)
  (* ########### Coefficients of Fracional Parentage ########### *)

  GenerateCFP::usage = "GenerateCFP[] generates the association for the coefficients of fractional parentage. Result is exported to the file ./data/CFP.m. The coefficients of fractional parentage are taken beyond the half-filled shell using the phase convention determined by the option \"PhaseFunction\". The default is \"NK\" which corresponds to the phase convention of Nielson and Koster. The other option is \"Judd\" which corresponds to the phase convention of Judd.";
  Options[GenerateCFP] = {"Export" -> True, "PhaseFunction"-> "NK"};
  GenerateCFP[OptionsPattern[]] := (
    CFP = Table[
      {numE, NKSL} -> First[CFPTerms[numE, NKSL]],
      {numE, 1, 7},
      {NKSL, AllowedNKSLTerms[numE]}];
    CFP = Association[CFP];
    (* Go all the way to f14 *)
    CFP = CFPExpander["Export" -> False, "PhaseFunction"-> OptionValue["PhaseFunction"]];
    If[OptionValue["Export"],
      Export[FileNameJoin[{moduleDir, "data", "CFPs.m"}], CFP];
    ];
    Return[CFP];
  );

  JuddCFPPhase::usage = "Phase between conjugate coefficients of fractional parentage according to Velkov's thesis, page 40.";
  JuddCFPPhase[parent_, parentS_, parentL_, daughterS_, daughterL_, parentSeniority_, daughterSeniority_] := Module[
    {spin, orbital, expo, phase},
    (
      {spin, orbital} = {1/2, 3};
      expo = (
          (parentS + parentL + daughterS + daughterL) - 
          (orbital + spin) + 
          1/2 * (parentSeniority + daughterSeniority - 1)
      );
      phase = Phaser[-expo];
      Return[phase];
    )
  ];

  NKCFPPhase::usage = "Phase between conjugate coefficients of fractional parentage according to Nielson and Koster page viii. Note that there is a typo on there the expression for zeta should be (-1)^((v-1)/2) instead of (-1)^(v - 1/2).";
  NKCFPPhase[parent_, parentS_, parentL_, daughterS_, daughterL_, parentSeniority_, daughterSeniority_] := Module[
    {spin, orbital, expo, phase},
    (
      {spin, orbital} = {1/2, 3};
      expo = (
          (parentS + parentL + daughterS + daughterL) - 
          (orbital + spin)
      );
      phase = Phaser[-expo];
      If[parent == 2*orbital,
          phase = phase * Phaser[(daughterSeniority-1)/2]];
      Return[phase];
    )
  ];

  Options[CFPExpander] = {"Export" -> True, "PhaseFunction" -> "NK"};
  CFPExpander::usage = "Using the coefficients of fractional parentage up to f7 this function calculates them up to f14. 
  The coefficients of fractional parentage are taken beyond the half-filled shell using the phase convention determined by the option \"PhaseFunction\". The default is \"NK\" which corresponds to the phase convention of Nielson and Koster. The other option is \"Judd\" which corresponds to the phase convention of Judd. The result is exported to the file ./data/CFPs_extended.m.";
  CFPExpander[OptionsPattern[]] := Module[
    {orbital, halfFilled, fullShell, parentMax, PhaseFun, 
    complementaryCFPs, daughter, conjugateDaughter, 
    conjugateParent, parentTerms, daughterTerms, 
    parentCFPs, daughterSeniority, daughterS, daughterL, 
    parentCFP, parentTerm, parentCFPval, 
    parentS, parentL, parentSeniority, phase, prefactor, 
    newCFPval, key, extendedCFPs, exportFname},
    (
      orbital    = 3;
      halfFilled = 2 * orbital + 1;
      fullShell  = 2 * halfFilled;
      parentMax  = 2 * orbital;
      
      PhaseFun  = <|
          "Judd" -> JuddCFPPhase,
          "NK" -> NKCFPPhase|>[OptionValue["PhaseFunction"]];
      PrintTemporary["Calculating CFPs using the phase system from ", PhaseFun];
      (* Initialize everything with lists to be filled in the next Do*)
      complementaryCFPs = 
          Table[
          ({numE, term} -> {term}),
          {numE, halfFilled + 1, fullShell - 1, 1},
          {term, AllowedNKSLTerms[numE]
          }];
      complementaryCFPs = Association[Flatten[complementaryCFPs]];
      Do[(
          daughter          = parent + 1;
          conjugateDaughter = fullShell - parent;
          conjugateParent   = conjugateDaughter - 1;
          parentTerms       = AllowedNKSLTerms[parent];
          daughterTerms     = AllowedNKSLTerms[daughter];
          Do[
          (
              parentCFPs             = Rest[CFP[{daughter, daughterTerm}]];
              daughterSeniority      = Seniority[daughterTerm];
              {daughterS, daughterL} = FindSL[daughterTerm];
              Do[
              (
                  {parentTerm, parentCFPval} = parentCFP;
                  {parentS, parentL}         = FindSL[parentTerm];
                  parentSeniority            = Seniority[parentTerm];
                  phase = PhaseFun[parent, parentS, parentL, 
                                  daughterS, daughterL,
                                  parentSeniority, daughterSeniority];
                  prefactor = (daughter * TPO[daughterS, daughterL]) / 
                              (conjugateDaughter * TPO[parentS, parentL]);
                  prefactor = Sqrt[prefactor];
                  newCFPval = phase * prefactor * parentCFPval;
                  key = {conjugateDaughter, parentTerm};
                  complementaryCFPs[key] = Append[complementaryCFPs[key], {daughterTerm, newCFPval}]
              ),
              {parentCFP, parentCFPs}
              ]
          ),
          {daughterTerm, daughterTerms}
          ]
          ),
      {parent, 1, parentMax}
      ];
      
      complementaryCFPs[{14, "1S"}] = {"1S", {"2F",1}};
      extendedCFPs         = Join[CFP, complementaryCFPs];
      If[OptionValue["Export"];,
      (
          exportFname = FileNameJoin[{moduleDir, "data", "CFPs_extended.m"}];
          Echo["Exporting to " <> exportFname];
          Export[exportFname, extendedCFPs];
      )
      ];
      Return[extendedCFPs];
  )
  ];

  GenerateCFPTable::usage = "GenerateCFPTable[] generates the table for the coefficients of fractional parentage. If the optional parameter \"Export\" is set to True then the resulting data is saved to ./data/CFPTable.m.
  The data being parsed here is the file attachment B1F_ALL.TXT which comes from Velkov's thesis.";
  Options[GenerateCFPTable] = {"Export" -> True};
  GenerateCFPTable[OptionsPattern[]] := Module[
    {rawText, rawLines, leadChar, configIndex, line, daughter,
    lineParts, numberCode, parsedNumber, toAppend, CFPTablefname},
    (
      CleanWhitespace[string_]     := StringReplace[string,RegularExpression["\\s+"]->" "];
      AddSpaceBeforeMinus[string_] := StringReplace[string,RegularExpression["(?<!\\s)-"]->" -"];
      ToIntegerOrString[list_]     := Map[If[StringMatchQ[#, NumberString], ToExpression[#], #] &, list];
      CFPTable      = ConstantArray[{},7];
      CFPTable[[1]] = {{"2F",{"1S",1}}};
      
      (* Cleaning before processing is useful *)
      rawText  = Import[FileNameJoin[{moduleDir, "data", "B1F_ALL.TXT"}]];
      rawLines = StringTrim/@StringSplit[rawText,"\n"];
      rawLines = Select[rawLines,#!=""&];
      rawLines = CleanWhitespace/@rawLines;
      rawLines = AddSpaceBeforeMinus/@rawLines;
      
      Do[(
        (* the first character can be used to identify the start of a block *)
        leadChar=StringTake[line,{1}];
        (* ..FN, N is at position 50 in that line *)
        If[leadChar=="[",
        (
          configIndex=ToExpression[StringTake[line,{50}]];
          Continue[];
        )
        ];
        (* Identify which daughter term is being listed *)
        If[StringContainsQ[line,"[DAUGHTER TERM]"],
          daughter=StringSplit[line,"["][[1]];
          CFPTable[[configIndex]]=Append[CFPTable[[configIndex]],{daughter}];
          Continue[];
        ];
        (* Once we get here we are already parsing a row with coefficient data *)
        lineParts    = StringSplit[line," "];
        parent       = lineParts[[1]];
        numberCode   = ToIntegerOrString[lineParts[[3;;]]];
        parsedNumber = SquarePrimeToNormal[numberCode];
        toAppend     = {parent,parsedNumber};
        CFPTable[[configIndex]][[-1]] = Append[CFPTable[[configIndex]][[-1]], toAppend]
      ),
      {line,rawLines}];
      If[OptionValue["Export"],
        (
        CFPTablefname = FileNameJoin[{moduleDir, "data", "CFPTable.m"}];
        Export[CFPTablefname, CFPTable];
        )
      ];
      Return[CFPTable];
    )
  ];

  GenerateCFPAssoc::usage = "GenerateCFPAssoc[] converts the coefficients of fractional parentage into an association in which zero values are explicit. If the option \"Export\" is set to True, the association is exported to the file /data/CFPAssoc.m. This function requires that the association CFP be defined.";
  Options[GenerateCFPAssoc] = {"Export" -> True};
  GenerateCFPAssoc[OptionsPattern[]] := (
    CFPAssoc = Association[];
    Do[
      (daughterTerms = AllowedNKSLTerms[numE];
      parentTerms    = AllowedNKSLTerms[numE - 1];
      Do[
        (
        cfps = CFP[{numE, daughter}];
        cfps = cfps[[2 ;;]];
        parents = First /@ cfps;
        Do[
          (
          key = {numE, daughter, parent};
          cfp = If[
            MemberQ[parents, parent],
            (
              idx = Position[parents, parent][[1, 1]];
              cfps[[idx]][[2]]
            ),
            0
            ];
          CFPAssoc[key] = cfp;
          ),
          {parent, parentTerms}
          ]
        ),
        {daughter, daughterTerms}
        ]
      ),
      {numE, 1, 14}
      ];
    If[OptionValue["Export"],
      (
      CFPAssocfname = FileNameJoin[{moduleDir, "data", "CFPAssoc.m"}];
      Export[CFPAssocfname, CFPAssoc];
      )
    ];
    Return[CFPAssoc];
  );

  CFPTerms::usage = "CFPTerms[numE] gives all the daughter and parent terms, together with the corresponding coefficients of fractional parentage, that correspond to the the f^n configuration.
  CFPTerms[numE, SL] gives all the daughter and parent terms, together with the corresponding coefficients of fractional parentage, that are compatible with the given string SL in the f^n configuration.
  CFPTerms[numE, L, S] gives all the daughter and parent terms, together with the corresponding coefficients of fractional parentage, that correspond to the given total orbital angular momentum L and total spin S n the f^n configuration. L being an integer, and S being integer or half-integer.
  In all cases the output is in the shape of a list with enclosed lists having the format {daughter_term, {parent_term_1, CFP_1}, {parent_term_2, CFP_2}, ...}.
  Only the one-body coefficients for f-electrons are provided.
  In all cases it must be that 1 <= n <= 7.
  These are according to the tables from Nielson & Koster.
  ";
  CFPTerms[numE_] := Part[CFPTable, numE]
  CFPTerms[numE_, SL_] := Module[
    {NKterms, CFPconfig},
    (
      NKterms   = {{}};
      CFPconfig = CFPTable[[numE]];
      Map[
        If[StringFreeQ[First[#], SL],
          Null, 
          NKterms = Join[NKterms, {#}, 1]
        ] &,
      CFPconfig
      ];
      NKterms = DeleteCases[NKterms, {}]
    )
  ];
  CFPTerms[numE_, L_, S_] := Module[
  {NKterms, SL, CFPconfig},
  (
    SL = StringJoin[ToString[2 S + 1], PrintL[L]];
    NKterms = {{}};
    CFPconfig = Part[CFPTable, numE];
    Map[
      If[StringFreeQ[First[#], SL],
        Null, 
        NKterms = Join[NKterms, {#}, 1]
      ]&,
    CFPconfig
    ];
    NKterms = DeleteCases[NKterms, {}]
  )
  ];

  (* ########### Coefficients of Fracional Parentage ########### *)
  (* ########################################################### *)

  (* ########################################################### *)
  (* ###################### Spin Orbit ######################### *)

  SpinOrbit::usage = "SpinOrbit[numE, SL, SpLp, J] returns the LSJ reduced matrix element \[Zeta] <SL, J|L.S|SpLp, J>. These are given as a function of \[Zeta]. This function requires that the association ReducedV1kTable be defined. 
  See equations 2-106 and 2-109 in Wybourne (1965). Equivalently see eqn. 12.43 in TASS.";
  SpinOrbit[numE_, SL_, SpLp_, J_] := Module[
    {S, L, Sp, Lp, orbital, sign, prefactor, val},
    (
      orbital   = 3;
      {S, L}    = FindSL[SL];
      {Sp, Lp}  = FindSL[SpLp];
      prefactor = Sqrt[orbital * (orbital+1) * (2*orbital+1)] *
                  SixJay[{L, Lp, 1}, {Sp, S, J}];
      sign      = Phaser[J + L + Sp];
      val       = sign * prefactor * \[Zeta] * ReducedV1kTable[{numE, SL, SpLp, 1}];
      Return[val];
    )
  ];

  SpinOrbitTable::usage="An association containing the matrix elements for the spin-orbit interaction for f^n configurations. The keys are lists of the form {n, SL, SpLp, J}.";

  GenerateSpinOrbitTable::usage = "GenerateSpinOrbitTable[nmax] computes the matrix elements for the spin-orbit interaction for f^n configurations up to n = nmax. The function returns an association whose keys are lists of the form {n, SL, SpLp, J}. If \"Export\" is set to True, then the result is exported to the data folder. It requires ReducedV1kTable to be defined.";
  Options[GenerateSpinOrbitTable] = {"Export" -> True};
  GenerateSpinOrbitTable[nmax_Integer:7, OptionsPattern[]] := Module[
    {numE, J, SL, SpLp, exportFname},
    (
      SpinOrbitTable = 
        Table[
          {numE, SL, SpLp, J} -> SpinOrbit[numE, SL, SpLp, J], 
        {numE, 1, nmax},
        {J, MinJ[numE], MaxJ[numE]},
        {SL,   Map[First, AllowedNKSLforJTerms[numE, J]]},
        {SpLp, Map[First, AllowedNKSLforJTerms[numE, J]]}
        ];
      SpinOrbitTable = Association[SpinOrbitTable];
      
      exportFname = FileNameJoin[{moduleDir, "data", "SpinOrbitTable.m"}];
      If[OptionValue["Export"],
        (
          Echo["Exporting to file " <> ToString[exportFname]];
          Export[exportFname, SpinOrbitTable];
        )
      ];
      Return[SpinOrbitTable];
    )
  ];

  (* ###################### Spin Orbit ######################### *)
  (* ########################################################### *)

  (* ########################################################### *)
  (* ################# Three Body Operators #################### *)

  ParseJudd1984::usage = "This function parses the data from tables 1 and 2 of Judd from Judd, BR, and MA Suskin. \"Complete Set of Orthogonal Scalar Operators for the Configuration f^3\". JOSA B 1, no. 2 (1984): 261-65.";
  Options[ParseJudd1984] = {"Export" -> False};
  ParseJudd1984[OptionsPattern[]] := (
    ParseJuddTab1[str_] := (
      strR = ToString[str];
      strR = StringReplace[strR, ".5" -> "^(1/2)"];
      num = ToExpression[strR];
      sign = Sign[num];
      num = sign*Simplify[Sqrt[num^2]];
      If[Round[num] == num, num = Round[num]];
      Return[num]);
    
    (* Parse table 1 from Judd 1984 *)
    judd1984Fname1 = FileNameJoin[{moduleDir, "data", "Judd1984-1.csv"}];
    data = Import[judd1984Fname1, "CSV", "Numeric" -> False];
    headers = data[[1]];
    data = data[[2 ;;]];
    data = Transpose[data];
    \[Psi] = Select[data[[1]], # != "" &];
    \[Psi]p = Select[data[[2]], # != "" &];
    matrixKeys = Transpose[{\[Psi], \[Psi]p}];
    data = data[[3 ;;]];
    cols = Table[ParseJuddTab1 /@ Select[col, # != "" &], {col, data}];
    cols = Select[cols, Length[#] == 21 &];
    tab1 = Prepend[Prepend[cols, \[Psi]p], \[Psi]];
    tab1 = Transpose[Prepend[Transpose[tab1], headers]];
    
    (* Parse table 2 from Judd 1984 *)
    judd1984Fname2 = FileNameJoin[{moduleDir, "data", "Judd1984-2.csv"}];
    data = Import[judd1984Fname2, "CSV", "Numeric" -> False];
    headers = data[[1]];
    data = data[[2 ;;]];
    data = Transpose[data];
    {operatorLabels, WUlabels, multiFactorSymbols, multiFactorValues} = data[[;; 4]];
    multiFactorValues = ParseJuddTab1 /@ multiFactorValues;
    multiFactorValues = AssociationThread[multiFactorSymbols -> multiFactorValues];
    
    (*scale values of table 1 given the values in table 2*)
    oppyS = {};
    normalTable = 
      Table[header = col[[1]];
        If[StringContainsQ[header, " "],
          (
            multiplierSymbol = StringSplit[header, " "][[1]];
            multiplierValue = multiFactorValues[multiplierSymbol];
            operatorSymbol = StringSplit[header, " "][[2]];
            oppyS = Append[oppyS, operatorSymbol];
          ), 
          (
            multiplierValue = 1;
            operatorSymbol = header;
          )
        ];
        normalValues = 1/multiplierValue*col[[2 ;;]];
        Join[{operatorSymbol}, normalValues], {col, tab1[[3 ;;]]}
      ];
    
    (*Create an association for the reduced matrix elements in the f^3 config*)
    juddOperators = Association[];
    Do[(
      col      = normalTable[[colIndex]];
      opLabel  = col[[1]];
      opValues = col[[2 ;;]];
      opMatrix = AssociationThread[matrixKeys -> opValues];
      Do[(
        opMatrix[Reverse[mKey]] = opMatrix[mKey]
        ), 
      {mKey, matrixKeys}
      ];
      juddOperators[{3, opLabel}] = opMatrix),
      {colIndex, 1, Length[normalTable]}
    ];
    
    (* special case of t2 in f3 *)
    (* this is the same as getting the reduced matrix elements from Judd 1966 *)
    numE = 3;
    e3Op      = juddOperators[{3, "e_{3}"}];
    t2prime   = juddOperators[{3, "t_{2}^{'}"}];
    prefactor = 1/(70 Sqrt[2]);
    t2Op = (# -> (t2prime[#] + prefactor*e3Op[#])) & /@ Keys[t2prime];
    t2Op = Association[t2Op];
    juddOperators[{3, "t_{2}"}] = t2Op;
    
    (*Special case of t11 in f3*)
    t11 = juddOperators[{3, "t_{11}"}];
    e\[Beta]primeOp = juddOperators[{3, "e_{\\beta}^{'}"}];
    t11primeOp = (# -> (t11[#] + Sqrt[3/385] e\[Beta]primeOp[#])) & /@ Keys[t11];
    t11primeOp = Association[t11primeOp];
    juddOperators[{3, "t_{11}^{'}"}] = t11primeOp;
    If[OptionValue["Export"],
      (
        (*export them*)
        PrintTemporary["Exporting ..."];
        exportFname = FileNameJoin[{moduleDir, "data", "juddOperators.m"}];
        Export[exportFname, juddOperators];
      )
    ];
    Return[juddOperators];
  );
  
  ThreeBodyTable::usage="ThreeBodyTable is an association containing the LS-reduced matrix elements for the three-body operators for f^n configurations. The keys are lists of the form {n, SL, SpLp}.";

  ThreeBodyTables::usage="ThreeBodyTables is an association whose keys are integers n from 1 to 7 and whose values are associations whose keys are symbols for the different three-body operators, and whose keys are of the form {LS, LpSp} where LS and LpSp are strings for LS-terms in f^n.";

  GenerateThreeBodyTables::usage = "This function generates the reduced matrix elements for the three body operators using the coefficients of fractional parentage, including those beyond f^7.";
  Options[GenerateThreeBodyTables] = {"Export" -> False};
  GenerateThreeBodyTables[OptionsPattern[]] := (
    tiKeys        = (StringReplace[ToString[#], {"T" -> "t_{", "p" -> "}^{'"}] <> "}") & /@ TSymbols;
    TSymbolsAssoc = AssociationThread[tiKeys -> TSymbols];
    juddOperators = ParseJudd1984[];
    (* op3MatrixElement[SL, SpLp, opSymbol] returns the value for the reduced matrix element of the operator opSymbol for the terms {SL, SpLp} in the f^3 configuration. *)
    op3MatrixElement[SL_, SpLp_, opSymbol_] := (
      jOP = juddOperators[{3, opSymbol}];
      key = {SL, SpLp};
      val = If[MemberQ[Keys[jOP], key],
        jOP[key],
        0];
      Return[val];
    );
    (* ti: This is the implementation of formula (2) in Judd & Suskin 1984. It computes the reduced matrix elements of ti in f^n by using the reduced matrix elements in f^3 and the coefficients of fractional parentage. If the option \"Fast\" is set to True then the values for n>7 are simply computed as the negatives of the values in the complementary configuration; this except for t2 and t11 which are treated as special cases. *)
    Options[ti] = {"Fast" -> True};
    ti[nE_, SL_, SpLp_, tiKey_, opOrder_ : 3, OptionsPattern[]] := Module[
      {nn, S, L, Sp, Lp, 
      cfpSL, cfpSpLp,
      parentSL, parentSpLp,
      tnk, tnks},
      (
        {S, L}   = FindSL[SL];
        {Sp, Lp} = FindSL[SpLp];
        fast     = OptionValue["Fast"];
        numH     = 14 - nE;
        If[fast && Not[MemberQ[{"t_{2}","t_{11}"},tiKey]] && nE > 7,
          Return[-tktable[{numH, SL, SpLp, tiKey}]]
        ];
        If[(S == Sp && L == Lp),
        (
          cfpSL   = CFP[{nE, SL}];
          cfpSpLp = CFP[{nE, SpLp}];
          tnks = Table[(
              parentSL   = cfpSL[[nn, 1]];
              parentSpLp = cfpSpLp[[mm, 1]];
              cfpSL[[nn, 2]] * cfpSpLp[[mm, 2]] *
              tktable[{nE - 1, parentSL, parentSpLp, tiKey}]
              ),
              {nn, 2, Length[cfpSL]},
              {mm, 2, Length[cfpSpLp]}
              ];
          tnk = Total[Flatten[tnks]];
          ),
        tnk = 0;
        ];
        Return[ nE / (nE - opOrder) * tnk];
      )
    ];
    (* Calculate the reduced matrix elements of t^i for n up to 14 *)
    tktable = <||>;
    Do[(
      Do[(
        tkValue = Which[numE <= 2,
          (*Initialize n=1,2 with zeros*)
          0,
          numE == 3,
          (* Grab matrix elem in f^3 from Judd 1984 *)
          SimplifyFun[op3MatrixElement[SL, SpLp, opKey]],
          True,
          SimplifyFun[ti[numE, SL, SpLp, opKey, If[opKey == "e_{3}", 2, 3]]]
          ];
        tktable[{numE, SL, SpLp, opKey}] = tkValue;
        ),
      {SL, AllowedNKSLTerms[numE]},
      {SpLp, AllowedNKSLTerms[numE]},
      {opKey, Append[tiKeys, "e_{3}"]}
      ];
      PrintTemporary[StringJoin["\[ScriptF]", ToString[numE], " configuration complete"]];
      ),
      {numE, 1, 14}
    ];
    
    (* Now use those reduced matrix elements to determine their sum as weighted by their corresponding strengths Ti *)
    ThreeBodyTable = <||>;
    Do[
      Do[
      (
        ThreeBodyTable[{numE, SL, SpLp}] = (
          Sum[(
            If[tiKey == "t_{2}", t2Switch, 1] * 
            tktable[{numE, SL, SpLp, tiKey}] * 
            TSymbolsAssoc[tiKey] +
            If[tiKey == "t_{2}", 1 - t2Switch, 0] * 
            (-tktable[{14 - numE, SL, SpLp, tiKey}]) * 
            TSymbolsAssoc[tiKey]
            ),
          {tiKey, tiKeys}
          ]
        );
      ),
      {SL, AllowedNKSLTerms[numE]},
      {SpLp, AllowedNKSLTerms[numE]}
      ];
    PrintTemporary[StringJoin["\[ScriptF]", ToString[numE], " matrix complete"]];,
    {numE, 1, 7}
    ];
    
    ThreeBodyTables = Table[(
      terms = AllowedNKSLTerms[numE];
      singleThreeBodyTable = 
        Table[
          {SL, SLp} -> ThreeBodyTable[{numE, SL, SLp}],
          {SL, terms},
          {SLp, terms}
        ];
      singleThreeBodyTable  = Flatten[singleThreeBodyTable];
      singleThreeBodyTables = Table[(
          notNullPosition = Position[TSymbols, notNullSymbol][[1, 1]];
          reps = ConstantArray[0, Length[TSymbols]];
          reps[[notNullPosition]] = 1;
          rep = AssociationThread[TSymbols -> reps];
          notNullSymbol -> Association[(singleThreeBodyTable /. rep)]
          ),
        {notNullSymbol, TSymbols}
        ];
      singleThreeBodyTables = Association[singleThreeBodyTables];
      numE -> singleThreeBodyTables),
      {numE, 1, 7}
  ];
  
  ThreeBodyTables = Association[ThreeBodyTables];
  If[OptionValue["Export"],
    (
      threeBodyTablefname = FileNameJoin[{moduleDir, "data", "ThreeBodyTable.m"}];
      Export[threeBodyTablefname, ThreeBodyTable];
      threeBodyTablesfname = FileNameJoin[{moduleDir, "data", "ThreeBodyTables.m"}];
      Export[threeBodyTablesfname, ThreeBodyTables];
    )
   ];
  Return[{ThreeBodyTable, ThreeBodyTables}];
  );


  (* ################# Three Body Operators #################### *)
  (* ########################################################### *)

  (* ########################################################### *)
  (* ################## Reduced SOO and ECSO ################### *)

  ReducedT11inf2::usage = "ReducedT11inf2[SL, SpLp] returns the reduced matrix element of the scalar component of the double tensor T11 for the given SL terms SL, SpLp.
  Data used here for m0, m2, m4 is from Table II of Judd, BR, HM Crosswhite, and Hannah Crosswhite. Intra-Atomic Magnetic Interactions for f Electrons. Physical Review 169, no. 1 (1968): 130.
  ";
  ReducedT11inf2[SL_, SpLp_] := Module[
    {T11inf2},
    (
      T11inf2 = <|
        {"1S", "3P"} -> 6 M0 + 2 M2 + 10/11 M4,
        {"3P", "3P"} -> -36 M0 - 72 M2 - 900/11 M4,
        {"3P", "1D"} -> -Sqrt[(2/15)] (27 M0 + 14 M2 + 115/11 M4),
        {"1D", "3F"} -> Sqrt[2/5] (23 M0 + 6 M2 - 195/11 M4),
        {"3F", "3F"} -> 2 Sqrt[14] (-15 M0 - M2 + 10/11 M4),
        {"3F", "1G"} -> Sqrt[11] (-6 M0 + 64/33 M2 - 1240/363 M4),
        {"1G", "3H"} -> Sqrt[2/5] (39 M0 - 728/33 M2 - 3175/363 M4),
        {"3H", "3H"} -> 8/Sqrt[55] (-132 M0 + 23 M2 + 130/11 M4),
        {"3H", "1I"} -> Sqrt[26] (-5 M0 - 30/11 M2 - 375/1573 M4)
        |>;
      Which[
        MemberQ[Keys[T11inf2],{SL,SpLp}],
          Return[T11inf2[{SL,SpLp}]],
        MemberQ[Keys[T11inf2],{SpLp,SL}],
          Return[T11inf2[{SpLp,SL}]],
        True,
          Return[0]
        ]
    )
  ];

  Reducedt11inf2::usage = "Reducedt11inf2[SL, SpLp] returns the reduced matrix element in f^2 of the double tensor operator t11 for the corresponding given terms {SL, SpLp}.
  Values given here are those from Table VII of \"Judd, BR, HM Crosswhite, and Hannah Crosswhite. \"Intra-Atomic Magnetic Interactions for f Electrons.\" Physical Review 169, no. 1 (1968): 130.\"
  ";
  Reducedt11inf2[SL_, SpLp_] := Module[
    {t11inf2},
    (
      t11inf2 = <|
        {"1S", "3P"} -> -2 P0 - 105 P2 - 231 P4 - 429 P6,
        {"3P", "3P"} -> -P0 - 45 P2 - 33 P4 + 1287 P6,
        {"3P", "1D"} -> Sqrt[15/2] (P0 + 32 P2 - 33 P4 - 286 P6),
        {"1D", "3F"} -> Sqrt[10] (-P0 - 9/2 P2 + 66 P4 - 429/2 P6),
        {"3F", "3F"} -> Sqrt[14] (-P0 + 10 P2 + 33 P4 + 286 P6),
        {"3F", "1G"} -> Sqrt[11] (P0 - 20 P2 + 32 P4 - 104 P6),
        {"1G", "3H"} -> Sqrt[10] (-P0 + 55/2 P2 - 23 P4 - 65/2 P6),
        {"3H", "3H"} -> Sqrt[55] (-P0 + 25 P2 + 51 P4 + 13 P6),
        {"3H", "1I"} -> Sqrt[13/2] (P0 - 21 P4 - 6 P6)
        |>;
      Which[
        MemberQ[Keys[t11inf2],{SL,SpLp}],
          Return[t11inf2[{SL,SpLp}]],
        MemberQ[Keys[t11inf2],{SpLp,SL}],
          Return[t11inf2[{SpLp,SL}]],
        True,
          Return[0]
      ]
    )
  ];

  ReducedSOOandECSOinf2::usage = "ReducedSOOandECSOinf2[SL, SpLp] returns the reduced matrix element corresponding to the operator (T11 + t11 - a13 * z13 / 6) for the terms {SL, SpLp}. This combination of operators corresponds to the spin-other-orbit plus ECSO interaction.
  The T11 operator corresponds to the spin-other-orbit interaction, and the t11 operator (associated with electrostatically-correlated spin-orbit) originates from configuration interaction analysis. To their sum a factor proportional to the operator z13 is subtracted since its effect is redundant to the spin-orbit interaction. The factor of 1/6 is not on Judd's 1968 paper, but it is on \"Chen, Xueyuan, Guokui Liu, Jean Margerie, and Michael F Reid. \"A Few Mistakes in Widely Used Data Files for Fn Configurations Calculations.\" Journal of Luminescence 128, no. 3 (2008): 421-27\".
  The values for the reduced matrix elements of z13 are obtained from Table IX of the same paper. The value for a13 is from table VIII.
  Rigorously speaking the Pk parameters here are subscripted. The conversion to superscripted parameters is performed elsewhere with the Prescaling replacement rules.
  ";
  ReducedSOOandECSOinf2[SL_, SpLp_]  := Module[
    {a13, z13, z13inf2, matElement, redSOOandECSOinf2}, 
    (
      a13 = (-33 M0 + 3 M2 + 15/11 M4 -
            6 P0 + 3/2 (35 P2 + 77 P4 + 143 P6));    
      z13inf2 = <|
          {"1S","3P"} -> 2,
          {"3P","3P"} -> 1,
          {"3P","1D"} -> -Sqrt[(15/2)],
          {"1D","3F"} -> Sqrt[10],
          {"3F","3F"} -> Sqrt[14],
          {"3F","1G"} -> -Sqrt[11],
          {"1G","3H"} -> Sqrt[10],
          {"3H","3H"} -> Sqrt[55],
          {"3H","1I"} -> -Sqrt[(13/2)]
          |>;
      matElement = Which[
        MemberQ[Keys[z13inf2],{SL,SpLp}],
          z13inf2[{SL,SpLp}],
        MemberQ[Keys[z13inf2],{SpLp,SL}],
          z13inf2[{SpLp,SL}],
        True,
          0
      ];
      redSOOandECSOinf2 = (
          ReducedT11inf2[SL, SpLp] +
          Reducedt11inf2[SL, SpLp] -
          a13 / 6 * matElement
      );
      redSOOandECSOinf2 = SimplifyFun[redSOOandECSOinf2];
      Return[redSOOandECSOinf2];
    )
  ];

  ReducedSOOandECSOinfn::usage = "ReducedSOOandECSOinfn[numE, SL, SpLp] calculates the reduced matrix elements of the (spin-other-orbit + ECSO) operator for the f^numE configuration corresponding to the terms SL and SpLp. This is done recursively, starting from tabulated values for f^2 from \"Judd, BR, HM Crosswhite, and Hannah Crosswhite. \"Intra-Atomic Magnetic Interactions for f Electrons.\" Physical Review 169, no. 1 (1968): 130.\", and by using equation (4) of that same paper. 
  ";
  ReducedSOOandECSOinfn[numE_, SL_, SpLp_] := Module[
    {spin, orbital, t, S, L, Sp, Lp, 
    idx1, idx2, cfpSL, cfpSpLp, parentSL,
    Sb, Lb, Sbp, Lbp, parentSpLp, funval},
    (
      {spin, orbital} = {1/2, 3};
      {S, L}   = FindSL[SL];
      {Sp, Lp} = FindSL[SpLp];
      t = 1;
      cfpSL    = CFP[{numE, SL}];
      cfpSpLp  = CFP[{numE, SpLp}];
      funval   = Sum[
        (
          parentSL = cfpSL[[idx2, 1]];
          parentSpLp = cfpSpLp[[idx1, 1]];
          {Sb, Lb}   = FindSL[parentSL];
          {Sbp, Lbp} = FindSL[parentSpLp];
          phase = Phaser[Sb + Lb + spin + orbital + Sp + Lp];
          (
            phase *
            cfpSpLp[[idx1, 2]] * cfpSL[[idx2, 2]] *
            SixJay[{S, t, Sp}, {Sbp, spin, Sb}] *
            SixJay[{L, t, Lp}, {Lbp, orbital, Lb}] *
            SOOandECSOLSTable[{numE - 1, parentSL, parentSpLp}]
          )
        ),
      {idx1, 2, Length[cfpSpLp]},
      {idx2, 2, Length[cfpSL]}
      ];
      funval *= numE / (numE - 2) * Sqrt[TPO[S, Sp, L, Lp]];
      Return[funval];
    )
  ];

  GenerateSOOandECSOLSTable::usage = "GenerateSOOandECSOLSTable[nmax] generates the LS reduced matrix elements of the spin-other-orbit + ECSO for the f^n configurations up to n=nmax. The values for n=1 and n=2 are taken from \"Judd, BR, HM Crosswhite, and Hannah Crosswhite. \"Intra-Atomic Magnetic Interactions for f Electrons.\" Physical Review 169, no. 1 (1968): 130.\", and the values for n>2 are calculated recursively using equation (4) of that same paper. The values are then exported to a file \"ReducedSOOandECSOLSTable.m\" in the data folder of this module. The values are also returned as an association.";
  Options[GenerateSOOandECSOLSTable] = {"Progress" -> True, "Export" -> True};
  GenerateSOOandECSOLSTable[nmax_Integer, OptionsPattern[]] := (
    If[And[OptionValue["Progress"], frontEndAvailable],
      (
        numItersai = Association[Table[numE->Length[AllowedNKSLTerms[numE]]^2, {numE, 1, nmax}]];
        counters   = Association[Table[numE->0, {numE, 1, nmax}]];
        totalIters = Total[Values[numItersai[[1;;nmax]]]];
        template1  = StringTemplate["Iteration `numiter` of `totaliter`"];
        template2  = StringTemplate["`remtime` min remaining"];template3 = StringTemplate["Iteration speed = `speed` ms/it"];
        template4  = StringTemplate["Time elapsed = `runtime` min"];
        progBar = PrintTemporary[
          Dynamic[
            Pane[
              Grid[{
                    {Superscript["f", numE]},
                    {template1[<|"numiter"->numiter, "totaliter"->totalIters|>]},
                    {template4[<|"runtime"->Round[QuantityMagnitude[UnitConvert[(Now-startTime), "min"]], 0.1]|>]},
                    {template2[<|"remtime"->Round[QuantityMagnitude[UnitConvert[(Now-startTime)/(numiter)*(totalIters-numiter), "min"]], 0.1]|>]},
                    {template3[<|"speed"->Round[QuantityMagnitude[Now-startTime, "ms"]/(numiter), 0.01]|>]}, {ProgressIndicator[Dynamic[numiter], {1, totalIters}]}
                  },
                  Frame->All
              ],
              Full,
              Alignment->Center
            ]
          ]
        ];
      )
    ];
    SOOandECSOLSTable = <||>;
    numiter   = 1;
    startTime = Now;
    Do[
      (
        numiter+= 1;
        SOOandECSOLSTable[{numE, SL, SpLp}] = Which[
          numE==1,
          0,
          numE==2,
          SimplifyFun[ReducedSOOandECSOinf2[SL, SpLp]],
          True,
          SimplifyFun[ReducedSOOandECSOinfn[numE,  SL, SpLp]]
        ];
      ),
    {numE, 1, nmax},
    {SL, AllowedNKSLTerms[numE]}, 
    {SpLp, AllowedNKSLTerms[numE]}
    ];
    If[And[OptionValue["Progress"], frontEndAvailable],
      NotebookDelete[progBar]];
    If[OptionValue["Export"],
      (fname = FileNameJoin[{moduleDir, "data", "ReducedSOOandECSOLSTable.m"}];
      Export[fname, SOOandECSOLSTable];
      )
    ];
    Return[SOOandECSOLSTable];
  );

  (* ################## Reduced SOO and ECSO ################### *)
  (* ########################################################### *)

  (* ########################################################### *)
  (* ###################### Spin-Spin ########################## *)

  ReducedT22inf2::usage = "ReducedT22inf2[SL, SpLp] returns the reduced matrix element of the scalar component of the double tensor T22 for the terms SL, SpLp in f^2.
  Data used here for m0, m2, m4 is from Table I of Judd, BR, HM Crosswhite, and Hannah Crosswhite. Intra-Atomic Magnetic Interactions for f Electrons. Physical Review 169, no. 1 (1968): 130.
  ";
  ReducedT22inf2[SL_, SpLp_] := Module[
    {statePosition, PsiPsipStates, m0, m2, m4, Tkk2m},
    (
      T22inf2 = <|
      {"3P", "3P"} -> -12 M0 - 24 M2 - 300/11 M4,
      {"3P", "3F"} -> 8/Sqrt[3] (3 M0 + M2 - 100/11 M4),
      {"3F", "3F"} -> 4/3 Sqrt[14] (-M0 + 8 M2 - 200/11 M4),
      {"3F", "3H"} -> 8/3 Sqrt[11/2] (2 M0 - 23/11 M2 - 325/121 M4),
      {"3H", "3H"} -> 4/3 Sqrt[143] (M0 - 34/11 M2 - (1325/1573) M4)
      |>;
      Which[
        MemberQ[Keys[T22inf2],{SL,SpLp}],
          Return[T22inf2[{SL,SpLp}]],
        MemberQ[Keys[T22inf2],{SpLp,SL}],
          Return[T22inf2[{SpLp,SL}]],
        True,
          Return[0]
      ]
    )
  ];

  ReducedT22infn::usage = "ReducedT22infn[n, SL, SpLp] calculates the reduced matrix element of the T22 operator for the f^n configuration corresponding to the terms SL and SpLp.
  This is done by using equation (4) of \"Judd, BR, HM Crosswhite, and Hannah Crosswhite. \"Intra-Atomic Magnetic Interactions for f Electrons.\" Physical Review 169, no. 1 (1968): 130.\"
  ";
  ReducedT22infn[numE_, SL_, SpLp_] := Module[
    {spin, orbital, t, idx1, idx2, S, L,
    Sp, Lp, cfpSL, cfpSpLp, parentSL,
    parentSpLp, Sb, Lb, Tnkk, phase, Sbp, Lbp},
    (
      {spin, orbital} = {1/2, 3};
      {S, L}   = FindSL[SL];
      {Sp, Lp} = FindSL[SpLp];
      t = 2;
      cfpSL    = CFP[{numE, SL}];
      cfpSpLp  = CFP[{numE, SpLp}];
      Tnkk = Sum[(
        parentSL   = cfpSL[[idx2, 1]]; 
        parentSpLp = cfpSpLp[[idx1, 1]];
        {Sb, Lb}   = FindSL[parentSL];
        {Sbp, Lbp} = FindSL[parentSpLp];
        phase = Phaser[Sb + Lb + spin + orbital + Sp + Lp];
        (
          phase *
          cfpSpLp[[idx1, 2]] * cfpSL[[idx2, 2]] *
          SixJay[{S, t, Sp}, {Sbp, spin, Sb}] *
          SixJay[{L, t, Lp}, {Lbp, orbital, Lb}] *
          T22Table[{numE - 1, parentSL, parentSpLp}]
        )
      ),
      {idx1, 2, Length[cfpSpLp]},
      {idx2, 2, Length[cfpSL]}
      ];
      Tnkk *= numE / (numE - 2) * Sqrt[TPO[S,Sp,L,Lp]];
      Return[Tnkk];
    )
  ];

  GenerateT22Table::usage = "GenerateT22Table[nmax] generates the LS reduced matrix elements for the double tensor operator T22 in f^n up to n=nmax. If the option \"Export\" is set to true then the resulting association is saved to the data folder. The values for n=1 and n=2 are taken from \"Judd, BR, HM Crosswhite, and Hannah Crosswhite. \"Intra-Atomic Magnetic Interactions for f Electrons.\" Physical Review 169, no. 1 (1968): 130.\", and the values for n>2 are calculated recursively using equation (4) of that same paper.
  This is an intermediate step to the calculation of the reduced matrix elements of the spin-spin operator.";
  Options[GenerateT22Table] = {"Export" -> True, "Progress" -> True};
  GenerateT22Table[nmax_Integer, OptionsPattern[]] := (
    If[And[OptionValue["Progress"], frontEndAvailable],
      (
        numItersai = Association[Table[numE->Length[AllowedNKSLTerms[numE]]^2, {numE, 1, nmax}]];
        counters = Association[Table[numE->0, {numE, 1, nmax}]];
        totalIters = Total[Values[numItersai[[1;;nmax]]]];
        template1 = StringTemplate["Iteration `numiter` of `totaliter`"];
        template2 = StringTemplate["`remtime` min remaining"];template3 = StringTemplate["Iteration speed = `speed` ms/it"];
        template4 = StringTemplate["Time elapsed = `runtime` min"];
        progBar = PrintTemporary[
          Dynamic[
            Pane[
              Grid[{{Superscript["f", numE]},
                    {template1[<|"numiter"->numiter, "totaliter"->totalIters|>]},
                    {template4[<|"runtime"->Round[QuantityMagnitude[UnitConvert[(Now-startTime), "min"]], 0.1]|>]},
                    {template2[<|"remtime"->Round[QuantityMagnitude[UnitConvert[(Now-startTime)/(numiter)*(totalIters-numiter), "min"]], 0.1]|>]},
                    {template3[<|"speed"->Round[QuantityMagnitude[Now-startTime, "ms"]/(numiter), 0.01]|>]},
                    {ProgressIndicator[Dynamic[numiter], {1, totalIters}]}},
                    Frame->All],
                    Full,
                    Alignment->Center]
                  ]
                ];
      )
    ];
    T22Table = <||>;
    startTime = Now;
    numiter = 1;
    Do[
      (
        numiter+= 1;
        T22Table[{numE, SL, SpLp}] = Which[
          numE==1,
          0,
          numE==2,
          SimplifyFun[ReducedT22inf2[SL, SpLp]],
          True,
          SimplifyFun[ReducedT22infn[numE,  SL, SpLp]]
        ];
      ),
    {numE, 1, nmax},
    {SL,   AllowedNKSLTerms[numE]}, 
    {SpLp, AllowedNKSLTerms[numE]}
    ];
    If[And[OptionValue["Progress"],frontEndAvailable],
      NotebookDelete[progBar]
    ];
    If[OptionValue["Export"],
      (
        fname = FileNameJoin[{moduleDir, "data", "ReducedT22Table.m"}];
        Export[fname, T22Table];
      )
    ];
    Return[T22Table];
  );

  SpinSpin::usage = "SpinSpin[n, SL, SpLp, J] returns the matrix element <|SL,J|spin-spin|SpLp,J|> for the spin-spin operator within the configuration f^n. This matrix element is independent of MJ. This is obtained by querying the relevant reduced matrix element from the association T22Table, putting in the adequate phase, and 6-j symbol.
  This is calculated according to equation (3) in \"Judd, BR, HM Crosswhite, and Hannah Crosswhite. \"Intra-Atomic Magnetic Interactions for f Electrons.\" Physical Review 169, no. 1 (1968): 130.\"
  \".
  ";
  SpinSpin[numE_, SL_, SpLp_, J_] := Module[
    {S, L, Sp, Lp, \[Alpha], val},
    (
      \[Alpha] = 2;
      {S, L} = FindSL[SL];
      {Sp, Lp} = FindSL[SpLp];
      val = (
              Phaser[Sp + L + J] * 
              SixJay[{Sp, Lp, J}, {L, S, \[Alpha]}] *
              T22Table[{numE, SL, SpLp}]
            );
      Return[val]
    )
  ];

  GenerateSpinSpinTable::usage = "GenerateSpinSpinTable[nmax] generates the reduced matrix elements in the |LSJ> basis for the spin-spin operator. It returns an association where the keys are of the form {numE, SL, SpLp, J}. If the option \"Export\" is set to True then the resulting object is saved to the data folder. Since this is a scalar operator, there is no MJ dependence. This dependence only comes into play when the crystal field contribution is taken into account.";
  Options[GenerateSpinSpinTable] = {"Export"->False};
  GenerateSpinSpinTable[nmax_, OptionsPattern[]] :=
    (
      SpinSpinTable = <||>;
      PrintTemporary[Dynamic[numE]];
      Do[
        SpinSpinTable[{numE, SL, SpLp, J}] = (SpinSpin[numE, SL, SpLp, J]);,
      {numE, 1, nmax},
      {J, MinJ[numE], MaxJ[numE]},
      {SL,   First /@ AllowedNKSLforJTerms[numE, J]},
      {SpLp, First /@ AllowedNKSLforJTerms[numE, J]}
      ];
      If[OptionValue["Export"],
      (fname = FileNameJoin[{moduleDir, "data", "SpinSpinTable.m"}];
        Export[fname, SpinSpinTable];
        )
      ];
      Return[SpinSpinTable];
    );

  (* ###################### Spin-Spin ########################## *)
  (* ########################################################### *)

  (* ################################################################## *)
  (* ## Spin-Other-Orbit and Electrostatically-Correlated-Spin-Orbit ## *)

  SOOandECSO::usage = "SOOandECSO[n, SL, SpLp, J] returns the matrix element <|SL,J|spin-spin|SpLp,J|> for the combined effects of the spin-other-orbit interaction and the electrostatically-correlated-spin-orbit (which originates from configuration interaction effects) within the configuration f^n. This matrix element is independent of MJ. This is obtained by querying the relevant reduced matrix element by querying the association SOOandECSOLSTable and putting in the adequate phase and 6-j symbol. The SOOandECSOLSTable puts together the reduced matrix elements from three operators.
  This is calculated according to equation (3) in \"Judd, BR, HM Crosswhite, and Hannah Crosswhite. \"Intra-Atomic Magnetic Interactions for f Electrons.\" Physical Review 169, no. 1 (1968): 130.\".
  ";
  SOOandECSO[numE_, SL_, SpLp_, J_] := Module[
    {S, Sp, L, Lp, \[Alpha], val},
    (
      \[Alpha] = 1;
      {S, L}   = FindSL[SL];
      {Sp, Lp} = FindSL[SpLp];
      val = (
              Phaser[Sp + L + J] *
              SixJay[{Sp, Lp, J}, {L, S, \[Alpha]}] *
              SOOandECSOLSTable[{numE, SL, SpLp}]
            );
      Return[val];
    )
  ];

  Prescaling = {P2 -> P2/225, P4 -> P4/1089, P6 -> 25 * P6 / 184041};

  GenerateSOOandECSOTable::usage = "GenerateSOOandECSOTable[nmax] generates the reduced matrix elements in the |LSJ> basis for the (spin-other-orbit + electrostatically-correlated-spin-orbit) operator. It returns an association where the keys are of the form {n, SL, SpLp, J}. If the option \"Export\" is set to True then the resulting object is saved to the data folder. Since this is a scalar operator, there is no MJ dependence. This dependence only comes into play when the crystal field contribution is taken into account.";
  Options[GenerateSOOandECSOTable] = {"Export"->False};
  GenerateSOOandECSOTable[nmax_, OptionsPattern[]] := (
    SOOandECSOTable = <||>;
    Do[
      SOOandECSOTable[{numE, SL, SpLp, J}] = (SOOandECSO[numE, SL, SpLp, J] /. Prescaling);,
    {numE, 1, nmax},
    {J, MinJ[numE], MaxJ[numE]},
    {SL,   First /@ AllowedNKSLforJTerms[numE, J]},
    {SpLp, First /@ AllowedNKSLforJTerms[numE, J]}
    ];
    If[OptionValue["Export"],
    (
      fname = FileNameJoin[{moduleDir, "data", "SOOandECSOTable.m"}];
      Export[fname, SOOandECSOTable];
    )
    ];
    Return[SOOandECSOTable];
  );

  (* ## Spin-Other-Orbit and Electrostatically-Correlated-Spin-Orbit ## *)
  (* ################################################################## *)

  (* ########################################################### *)
  (* ################## Magnetic Interactions ################## *)

  MagneticInteractions::usage = "MagneticInteractions[{numE, SL, SLp, J}] returns the matrix element of the magnetic interaction between the terms SL and SLp in the f^numE configuration for the given value of J. The interaction is given by the sum of the spin-spin, the spin-other-orbit,  and the electrostatically-correlated-spin-orbit interactions.
  The part corresponding to the spin-spin interaction is provided by SpinSpin[{numE, SL, SLp, J}].
  The part corresponding to SOO and ECSO is provided by the function SOOandECSO[{numE, SL, SLp, J}].
  The option \"ChenDeltas\" can be used to include or exclude the Chen deltas from the calculation. The default is to exclude them. If this option is used, then the chenDeltas association needs to be loaded into the session with LoadChenDeltas[].";
  Options[MagneticInteractions] = {"ChenDeltas" -> False};
  MagneticInteractions[{numE_, SL_, SLp_, J_}, OptionsPattern[]] := Module[
    {key, ss, sooandecso, total,
    S, L, Sp, Lp, phase, sixjay, 
    M0v, M2v, M4v, 
    P2v, P4v, P6v}, 
    (
      key        = {numE, SL, SLp, J};
      ss         = \[Sigma]SS * SpinSpinTable[key];
      sooandecso = SOOandECSOTable[key];
      total = ss + sooandecso;
      total = SimplifyFun[total];
      If[
        Not[OptionValue["ChenDeltas"]],
        Return[total]
      ];
      (* In the type A errors the wrong values are different *)
      If[MemberQ[Keys[chenDeltas["A"]], {numE, SL, SLp}],
        (
          {S, L}   = FindSL[SL];
          {Sp, Lp} = FindSL[SLp];
          phase    = Phaser[Sp + L + J];
          sixjay   = SixJay[{Sp, Lp, J}, {L, S, 1}];
          {M0v, M2v, M4v, P2v, P4v, P6v} = chenDeltas["A"][{numE, SL, SLp}]["wrong"];
          total    = (
                        phase * sixjay * 
                        (
                          M0v*M0 + M2v*M2 + M4v*M4 +
                          P2v*P2 + P4v*P4 + P6v*P6
                        )
                    );
          total    = wChErrA * total + (1 - wChErrA) * (ss + sooandecso)
        )
      ];
      (* In the type B errors the wrong values are zeros all around *)
      If[MemberQ[chenDeltas["B"], {numE, SL, SLp}],
        (
          total  = (1 - wChErrB) * (ss + sooandecso)
        )
      ];
      Return[total];
    )
  ];

  (* ################## Magnetic Interactions ################## *)
  (* ########################################################### *)

  (* ########################################################### *)
  (* #################### Free-Ion Energies #################### *)
  
  GenerateFreeIonTable::usage="GenerateFreeIonTable[] generates an association for free-ion energies in terms of Slater integrals Fk and spin-orbit parameter \[Zeta]. It returns an association where the keys are of the form {nE, SL, SpLp}. If the option \"Export\" is set to True then the resulting object is saved to the data folder. The free-ion Hamiltonian is the sum of the electrostatic and spin-orbit interactions. The electrostatic interaction is given by the function Electrostatic[{numE, SL, SpLp}] and the spin-orbit interaction is given by the function SpinOrbitTable[{numE, SL, SpLp}]. The values for the electrostatic interaction are taken from the data file ElectrostaticTable.m and the values for the spin-orbit interaction are taken from the data file SpinOrbitTable.m. The values for the free-ion Hamiltonian are then exported to a file \"FreeIonTable.m\" in the data folder of this module. The values are also returned as an association.";
  Options[GenerateFreeIonTable] = {"Export" -> False};
  GenerateFreeIonTable[OptionsPattern[]] := Module[
    {terms, numEH, zetaSign, fname, FreeIonTable},
    (
      If[Not[ValueQ[ElectrostaticTable]],
        LoadElectrostatic[]
        ];    
      If[Not[ValueQ[SpinOrbitTable]],
        LoadSpinOrbit[]
      ];
      If[Not[ValueQ[ReducedUkTable]],
        LoadUk[]
      ];
      FreeIonTable = <||>;
      Do[
        (
          terms = AllowedNKSLJTerms[nE];
          numEH = Min[nE, 14 - nE];
          zetaSign = If[nE > 7, -1, 1];
          Do[
            FreeIonTable[{nE,term[[1]],term[[2]]}] = (
              Electrostatic[{numEH, term[[1]], term[[1]]}] + 
                zetaSign * SpinOrbitTable[{numEH, term[[1]], term[[1]], term[[2]]}]
                ),
            {term, terms}];
        ),
        {nE, 1, 14}
      ];
      If[OptionValue["Export"],
        (
          fname = FileNameJoin[{moduleDir, "data", "FreeIonTable.m"}];
          Export[fname, FreeIonTable];
        )
      ];
      Return[FreeIonTable];
    )
  ];

  LoadFreeIon::usage = "LoadFreeIon[] loads the free-ion energies from the data folder. The values are stored in the association FreeIonTable.";
  LoadFreeIon[] := (
    If[ValueQ[FreeIonTable],
      Return[]
    ];
    PrintTemporary["Loading the association of free-ion energies ..."];
    FreeIonTableFname = FileNameJoin[{moduleDir, "data", "FreeIonTable.m"}];
    FreeIonTable = If[!FileExistsQ[FreeIonTableFname],
      (
        PrintTemporary[">> FreeIonTable.m not found, generating ..."];
        GenerateFreeIonTable["Export" -> True]
      ),
        Import[FreeIonTableFname]
    ];
  );
 
  (* #################### Free-Ion Energies #################### *)
  (* ########################################################### *)

  (* ########################################################### *)
  (* ##################### Crystal Field ####################### *)

  Cqk::usage = "Cqk[numE, q, k, NKSL, J, M, NKSLp, Jp, Mp]. In Wybourne (1965) see equations 6-3, 6-4, and 6-5. Also in TASS see equation 11.53.";
  Cqk[numE_, q_, k_, NKSL_, J_, M_, NKSLp_, Jp_, Mp_] := Module[
    {S, Sp, L, Lp, orbital, val},
    (
      orbital = 3;
      {S, L}   = FindSL[NKSL];
      {Sp, Lp} = FindSL[NKSLp];
      f1  = ThreeJay[{J, -M}, {k, q}, {Jp, Mp}];
      val = 
        If[f1==0,
          0,
          (
            f2 = SixJay[{L, J, S}, {Jp, Lp, k}] ;
            If[f2==0,
              0,
              (
                f3 = ReducedUkTable[{numE, orbital, NKSL, NKSLp, k}];
                If[f3==0,
                  0,
                  (
                    ( 
                      Phaser[J - M + S + Lp + J + k] *
                      Sqrt[TPO[J, Jp]] *
                      f1 *
                      f2 *
                      f3 *
                      Ck[orbital, k]
                    )
                  )
                ]
              )
            ]
          )
        ];
      Return[val];
    )
  ];

  Bqk::usage = "Real part of the Bqk coefficients.";
  Bqk[q_, 2] := {B02/2, B12, B22}[[q + 1]];
  Bqk[q_, 4] := {B04/2, B14, B24, B34, B44}[[q + 1]];
  Bqk[q_, 6] := {B06/2, B16, B26, B36, B46, B56, B66}[[q + 1]];

  Sqk::usage = "Imaginary part of the Bqk coefficients.";
  Sqk[q_, 2] := {0, S12, S22}[[q + 1]];
  Sqk[q_, 4] := {0, S14, S24, S34, S44}[[q + 1]];
  Sqk[q_, 6] := {0, S16, S26, S36, S46, S56, S66}[[q + 1]];

  CrystalField::usage = "CrystalField[n, NKSL, J, M, NKSLp, Jp, Mp] calculates the matrix element of the crystal field in terms of Bqk and Sqk parameters for configuration f^numE. It is calculated as an association with keys of the form {n, NKSL, J, M, NKSLp, Jp, Mp}.";
  CrystalField[numE_, NKSL_, J_, M_, NKSLp_, Jp_, Mp_] := (
    Sum[
      (
        cqk  = Cqk[numE,  q, k, NKSL, J, M, NKSLp, Jp, Mp];
        cmqk = Cqk[numE, -q, k, NKSL, J, M, NKSLp, Jp, Mp];
        Bqk[q, k]   * (cqk + (-1)^q * cmqk) +
        I*Sqk[q, k] * (cqk - (-1)^q * cmqk)
        ),
    {k, {2, 4, 6}},
    {q, 0, k}
    ]
  )

  TotalCFIters::usage = "TotalIters[i, j] returns total number of function evaluations for calculating all the matrix elements for the \!\(\*SuperscriptBox[\(f\), \(i\)]\) to the \!\(\*SuperscriptBox[\(f\), \(j\)]\) configurations.";
  TotalCFIters[i_, j_] := (
    numIters = {196, 8281, 132496, 1002001, 4008004, 9018009, 11778624};
    Return[Total[numIters[[i ;; j]]]];
    )

  GenerateCrystalFieldTable::usage = "GenerateCrystalFieldTable[{numEs]] computes the matrix values for the crystal field interaction for f^n configurations the given list of numE in  numEs. The function calculates the association CrystalFieldTable with keys of the form {numE, NKSL, J, M, NKSLp, Jp, Mp}. If the option \"Export\" is set to True, then the result is exported to the data subfolder for the folder in which this package is in. If the option \"Progress\" is set to True then an interactive progress indicator is shown. If \"Compress\" is set to true the exported values are compressed when exporting.";
  Options[GenerateCrystalFieldTable] = {"Export" -> False, "Progress" -> True, "Compress" -> True, "Overwrite" -> False};
  GenerateCrystalFieldTable[numEs_List:{1,2,3,4,5,6,7}, OptionsPattern[]] := (
    ExportFun = 
    If[OptionValue["Compress"],
      ExportMZip,
      Export
    ];
    numiter = 1;
    template1 = StringTemplate["Iteration `numiter` of `totaliter`"];
    template2 = StringTemplate["`remtime` min remaining"];
    template3 = StringTemplate["Iteration speed = `speed` ms/it"];
    template4 = StringTemplate["Time elapsed = `runtime` min"];
    totalIter = Total[TotalCFIters[#, #] & /@ numEs];
    freebies = 0;
    startTime = Now;
    If[And[OptionValue["Progress"], frontEndAvailable],
      progBar = PrintTemporary[
        Dynamic[
          Pane[
            Grid[
              {
                {Superscript["f", numE]},
                {template1[<|"numiter" -> numiter, "totaliter" -> totalIter|>]},
                {template4[<|"runtime" -> Round[QuantityMagnitude[UnitConvert[(Now - startTime), "min"]], 0.1]|>]},
                {template2[<|"remtime" -> Round[QuantityMagnitude[UnitConvert[(Now - startTime)/(numiter - freebies) * (totalIter - numiter), "min"]], 0.1]|>]},
                {template3[<|"speed" -> Round[QuantityMagnitude[Now - startTime, "ms"]/(numiter-freebies), 0.01]|>]},
                {ProgressIndicator[Dynamic[numiter], {1, totalIter}]}
              },
            Frame -> All
            ],
          Full,
          Alignment -> Center
          ]
        ]
      ];
    ];
    Do[
      (
        exportFname = FileNameJoin[{moduleDir, "data", "CrystalFieldTable_f"<>ToString[numE]<>".m"}];
        If[And[FileExistsQ[exportFname],Not[OptionValue["Overwrite"]]],
          Echo["File exists, skipping ..."];
          numiter+=  TotalCFIters[numE, numE];
          freebies+= TotalCFIters[numE, numE];
          Continue[];
        ];
        CrystalFieldTable = <||>;
        Do[
          (
            numiter+=1;
            {S,L}   = FindSL[NKSL];
            {Sp,Lp} = FindSL[NKSLp];
            CrystalFieldTable[{numE,NKSL,J,M,NKSLp,Jp,Mp}] = 
              Which[
                Abs[M-Mp]>6,
                0,
                Abs[S-Sp]!=0,
                0,
                True,
                CrystalField[numE,NKSL,J,M,NKSLp,Jp,Mp]
              ];
          ),
        {J,  MinJ[numE], MaxJ[numE]},
        {Jp, MinJ[numE], MaxJ[numE]},
        {M,  AllowedMforJ[J]},
        {Mp, AllowedMforJ[Jp]},
        {NKSL, First/@AllowedNKSLforJTerms[numE,J]},
        {NKSLp, First/@AllowedNKSLforJTerms[numE,Jp]}
        ];
        If[And[OptionValue["Progress"],frontEndAvailable],
          NotebookDelete[progBar]
        ];
        If[OptionValue["Export"],
          (
            Echo["Exporting to file " <> ToString[exportFname]];
            ExportFun[exportFname, CrystalFieldTable];
          )
        ];
      ),
    {numE, numEs}
    ]
  )

  ParseBenelli2015::usage="ParseBenelli2015[] parses the data from file /data/benelli_and_gatteschi_table3p3.csv which corresponds to Table 3.3 of Benelli and Gatteshi, Introduction to Molecular. This data provides the form that the crystal field has under different point group symmetries. This function parses that data into an association with keys equal to strings representing any of the 32 crystallographic point groups.";
  Options[ParseBenelli2015] = {"Export"->False};
  ParseBenelli2015[OptionsPattern[]] := Module[
    {fname, crystalSym,
    crystalSymmetries, parseFun,
    chars, qk, groupName, family,
    groupNum, params},
  (
    fname      = FileNameJoin[{moduleDir,"data","benelli_and_gatteschi_table3p3.csv"}];
    crystalSym = Import[fname][[2;;33]];
    crystalSymmetries = <||>;
    parseFun[txt_] := (
      chars = Characters[txt];
      qk    = chars[[-2;;]];
      If[chars[[1]]=="R",
      (
        Return[{ToExpression@StringJoin[{"B",qk[[1]],qk[[2]]}]}]
      ),
      (
        If[qk[[1]]=="0",
          Return[{ToExpression@StringJoin[{"B",qk[[1]],qk[[2]]}]}]
        ];
        Return[{
          ToExpression@StringJoin[{"B",qk[[1]],qk[[2]]}],
          ToExpression@StringJoin[{"S",qk[[1]],qk[[2]]}]
        }]
      )]
    );
    Do[
    (
      groupNum  = Round@ToExpression@row[[1]];
      groupName = row[[2]];
      family    = row[[3]];
      params    = Select[row[[4;;]],And[FreeQ[#,Missing],#=!=""]&];
      params    = parseFun/@params;
      params    = <|"BqkSqk"->Sort@Flatten[params],
        "aliases"->{groupNum},
        "constraints"->{}|>;
      If[MemberQ[{"T","Th","O","Td","Oh"},groupName],
        params["constraints"] = {
            {B44 ->  Sqrt[5/14] B04, B46 -> -Sqrt[7/2] B06},
            {B44 -> -Sqrt[5/14] B04, B46 ->  Sqrt[7/2] B06}
            }
      ];
      If[StringContainsQ[groupName,","],
      (
        {alias1, alias2}          = StringSplit[groupName,","];
        crystalSymmetries[alias1] = params;
        crystalSymmetries[alias1]["aliases"] = {groupNum, alias2};
        crystalSymmetries[alias2] = params;
        crystalSymmetries[alias2]["aliases"] = {groupNum, alias1};
      ),
      (
      crystalSymmetries[groupName] = params;
      )
      ]
    ),
    {row,crystalSym}];
    crystalSymmetries["source"] = "Benelli and Gatteschi, 2015, Introduction to Molecular Magnetism, table 3.3.";
    If[OptionValue["Export"],
      Export[FileNameJoin[{moduleDir,"data","crystalFieldFunctionalForms.m"}],crystalSymmetries];
    ];
    Return[crystalSymmetries];
  )
  ]; 

  CrystalFieldForm::usage = "CrystalFieldForm[symmetryGroup] returns an association that describes the crystal field parameters that are necessary to describe a crystal field for the given symmetry group.
  
  The symmetry group must be given as a string in Schoenflies notation and must be one of C1, Ci, S2, Cs, C1h, C2, C2h, C2v, D2, D2h, S4, C4, C4h, D2d, C4v, D4, D4h, C3, S6, C3h, C3v, D3, D3d, D3h, C6, C6h, C6v, D6, D6h, T, Th, Td, O, Oh.
   
  The returned association has three keys: 
    \"BqkSqk\" whose values is a list with the nonzero Bqk and Sqk parameters;
    \"constraints\" whose value is either an empty list, or a lists of replacements rules that are constraints on the Bqk and Sqk parameters;
    \"simplifier\" whose value is an association that can be used to set to zero the crystal field parameters that are zero for the given symmetry group;
    \"aliases\" whose value is a list with the integer by which the point group is also known for and an alternate Schoenflies symbol if it exists.
  
  This uses data from table 3.3 in Benelli and Gatteschi, 2015.";
  CrystalFieldForm[symmetryGroupString_] := (
    If[Not@ValueQ[crystalFieldFunctionalForms],
      crystalFieldFunctionalForms = Import[FileNameJoin[{moduleDir, "data", "crystalFieldFunctionalForms.m"}]];
    ];
    cfForm = crystalFieldFunctionalForms[symmetryGroupString];
    simplifier = Association[(# -> 0) &/@ Complement[cfSymbols, cfForm["BqkSqk"]]];
    Return[Join[cfForm, <|"simplifier"->simplifier|>]];
  )

  (* ##################### Crystal Field ####################### *)
  (* ########################################################### *)

  (* ########################################################### *)
  (* #### Configuration-Interaction via Casimir Operators ###### *)

  CasimirSO3::usage = "CasimirSO3[{SL, SpLp}] returns LS reduced matrix element of the configuration interaction term corresponding to the Casimir operator of R3.";
  CasimirSO3[{SL_, SpLp_}] := (
    {S, L} = FindSL[SL];
    If[SL == SpLp,
      \[Alpha] * L * (L + 1),
      0
    ]
  )

  GG2U::usage = "GG2U is an association whose keys are labels for the irreducible representations of group G2 and whose values are the eigenvalues of the corresponding Casimir operator.
  Reference: Wybourne, \"Spectroscopic Properties of Rare Earths\", table 2-6.";
  GG2U = Association[{
      "00" -> 0,
      "10" -> 6/12 ,
      "11" -> 12/12,
      "20" -> 14/12,
      "21" -> 21/12,
      "22" -> 30/12,
      "30" -> 24/12,
      "31" -> 32/12,
      "40" -> 36/12}
    ];

  CasimirG2::usage = "CasimirG2[{SL, SpLp}] returns LS reduced matrix element of the configuration interaction term corresponding to the Casimir operator of G2.";
  CasimirG2[{SL_, SpLp_}] := (
    Ulabel = FindNKLSTerm[SL][[1]][[4]];
    If[SL==SpLp,
      \[Beta] * GG2U[Ulabel],
      0
    ]
  )

  GSO7W::usage = "GSO7W is an association whose keys are labels for the irreducible representations of group R7 and whose values are the eigenvalues of the corresponding Casimir operator.
  Reference: Wybourne, \"Spectroscopic Properties of Rare Earths\", table 2-7.";
  GSO7W := Association[
    {
      "000" -> 0,
      "100" -> 3/5,
      "110" -> 5/5,
      "111" -> 6/5,
      "200" -> 7/5,
      "210" -> 9/5,
      "211" -> 10/5,
      "220" -> 12/5,
      "221" -> 13/5,
      "222" -> 15/5
    }
  ];

  CasimirSO7::usage = "CasimirSO7[{SL, SpLp}] returns the LS reduced matrix element of the configuration interaction term corresponding to the Casimir operator of R7.";
  CasimirSO7[{SL_, SpLp_}] := (
    Wlabel = FindNKLSTerm[SL][[1]][[3]];
    If[SL==SpLp,
      \[Gamma] * GSO7W[Wlabel],
      0
    ]
  )

  ElectrostaticConfigInteraction::usage = "ElectrostaticConfigInteraction[numE, {SL, SpLp}] returns the matrix element for configuration interaction as approximated by the Casimir operators of the groups R3, G2, and R7. SL and SpLp are strings that represent terms under LS coupling.";
  ElectrostaticConfigInteraction[numE_, {SL_, SpLp_}] := Module[
    {S, L, val},
    (
      If[
        Or[numE == 1, numE==13],
        Return[0];
      ];
      {S, L} = FindSL[SL];
      val = (
        If[SL == SpLp,
          CasimirSO3[{SL, SL}] +
          CasimirSO7[{SL, SL}] +
          CasimirG2[{SL, SL}],
          0
        ]
        );
      ElectrostaticConfigInteraction[numE, {S, L}] = val;
      Return[val];
    )
  ];

  (* #### Configuration-Interaction via Casimir Operators ###### *)
  (* ########################################################### *)

  (* ########################################################### *)
  (* ##################### Block assembly ###################### *)

  JJBlockMatrix::usage = "For given J, J' in the f^n configuration JJBlockMatrix[numE, J, J'] determines all the SL S'L' terms that may contribute to them and using those it provides the matrix elements <J, LS | H | J', LS'>. H having contributions from the following interactions: Coulomb, spin-orbit, spin-other-orbit, electrostatically-correlated-spin-orbit, spin-spin, three-body interactions, and crystal-field.";
  Options[JJBlockMatrix] = {"Sparse"->True, "ChenDeltas"->False};
  JJBlockMatrix[numE_, J_, Jp_, CFTable_, OptionsPattern[]] := Module[
    {NKSLJMs, NKSLJMps, NKSLJM, NKSLJMp,
    SLterm, SpLpterm,
    MJ, MJp,
    subKron, matValue, eMatrix},
    (
      NKSLJMs  = AllowedNKSLJMforJTerms[numE, J];
      NKSLJMps = AllowedNKSLJMforJTerms[numE, Jp];
      eMatrix = 
        Table[
          (*Condition for a scalar matrix op*)
          SLterm   = NKSLJM[[1]];
          SpLpterm = NKSLJMp[[1]];
          MJ       = NKSLJM[[3]];
          MJp      = NKSLJMp[[3]];
          subKron  = (
              KroneckerDelta[J, Jp] *
              KroneckerDelta[MJ, MJp]
                    );
          matValue = 
            If[subKron==0,
              0,
                (
                  ElectrostaticTable[{numE, SLterm, SpLpterm}] +
                  ElectrostaticConfigInteraction[numE, {SLterm, SpLpterm}] +
                  SpinOrbitTable[{numE, SLterm, SpLpterm, J}] +
                  MagneticInteractions[{numE, SLterm, SpLpterm, J}, 
                    "ChenDeltas" -> OptionValue["ChenDeltas"]] +
                  ThreeBodyTable[{numE, SLterm, SpLpterm}]
                )
            ];
          matValue += CFTable[{numE, SLterm, J, MJ, SpLpterm, Jp, MJp}];
          matValue,
        {NKSLJMp, NKSLJMps},
        {NKSLJM , NKSLJMs}
        ];
    If[OptionValue["Sparse"],
      eMatrix = SparseArray[eMatrix]
    ];
    Return[eMatrix]
  )
  ];

  EnergyStates::usage = "Alias for AllowedNKSLJMforJTerms. At some point may be used to redefine states used in basis.";
  EnergyStates[numE_, J_]:= AllowedNKSLJMforJTerms[numE, J];

  JJBlockMatrixFileName::usage = "JJBlockMatrixFileName[numE] gives the filename for the energy matrix table for an atom with numE f-electrons. The function admits an optional parameter \"FilenameAppendix\" which can be used to modify the filename.";
  Options[JJBlockMatrixFileName] = {"FilenameAppendix" -> ""};
  JJBlockMatrixFileName[numE_Integer, OptionsPattern[]] := (
    fileApp = OptionValue["FilenameAppendix"];
    fname = FileNameJoin[{moduleDir,
        "hams",
        StringJoin[{"f", ToString[numE], "_JJBlockMatrixTable", fileApp ,".m"}]}];
    Return[fname];
  );

  TabulateJJBlockMatrixTable::usage = "TabulateJJBlockMatrixTable[numE, CFTable] returns an association JJBlockMatrixTable with keys equal to lists of the form {numE, J, Jp} and values equal to JJ blocks of the semi-empirical Hamiltonian. The function admits an optional parameter \"Sparse\" which can be used to control whether the matrix is returned as a SparseArray or not. The default is True. Another admitted option is \"ChenDeltas\" which can be used to include or exclude the Chen deltas from the calculation. The default is to exclude them. If this option is used, then the chenDeltas association needs to be loaded into the session with LoadChenDeltas[].";
  Options[TabulateJJBlockMatrixTable] = {"Sparse"->True, "ChenDeltas"->False};
  TabulateJJBlockMatrixTable[numE_, CFTable_, OptionsPattern[]] := (
    JJBlockMatrixTable = <||>;
    totalIterations = Length[AllowedJ[numE]]^2;
    template1 = StringTemplate["Iteration `numiter` of `totaliter`"];
    template2 = StringTemplate["`remtime` min remaining"];
    template4 = StringTemplate["Time elapsed = `runtime` min"];
    numiter   = 0;
    startTime = Now;
    If[$FrontEnd =!= Null,
      (
        temp = PrintTemporary[
          Dynamic[
            Grid[
              {
                {template1[<|"numiter"->numiter, "totaliter"->totalIterations|>]},
                {template2[<|"remtime"->Round[QuantityMagnitude[UnitConvert[(Now-startTime)/(Max[1,numiter])*(totalIterations-numiter), "min"]], 0.1]|>]},
                {template4[<|"runtime"->Round[QuantityMagnitude[UnitConvert[(Now-startTime), "min"]], 0.1]|>]},
                {ProgressIndicator[numiter, {1, totalIterations}]}
              }
            ]
        ]
        ];
      )
    ];
    Do[
      (
        JJBlockMatrixTable[{numE, J, Jp}] = JJBlockMatrix[numE, J, Jp, CFTable, "Sparse"->OptionValue["Sparse"], "ChenDeltas" -> OptionValue["ChenDeltas"]];
        numiter += 1;
      ),
    {Jp, AllowedJ[numE]},
    {J, AllowedJ[numE]}
    ];
    If[$FrontEnd =!= Null,
      NotebookDelete[temp]
    ];
    Return[JJBlockMatrixTable];
  );

  TabulateManyJJBlockMatrixTables::usage = "TabulateManyJJBlockMatrixTables[{n1, n2, ...}] calculates the tables of matrix elements for the requested f^n_i configurations. The function does not return the matrices themselves. It instead returns an association whose keys are numE and whose values are the filenames where the output of TabulateJJBlockMatrixTables was saved to.The output consists of an association whose keys are of the form {n, J, Jp} and whose values are rectangular arrays given the values of <|LSJMJa|H|L'S'J'MJ'a'|>.";
  Options[TabulateManyJJBlockMatrixTables] = {"Overwrite"->False, "Sparse"->True, "ChenDeltas"->False, "FilenameAppendix"-> "", "Compressed" -> False};
  TabulateManyJJBlockMatrixTables[ns_, OptionsPattern[]] := (
    overwrite = OptionValue["Overwrite"];
    fNames = <||>;
    fileApp = OptionValue["FilenameAppendix"];
    ExportFun = If[OptionValue["Compressed"], ExportMZip, Export];
    Do[
      (
        CFdataFilename = FileNameJoin[{moduleDir, "data", "CrystalFieldTable_f"<>ToString[numE]<>".zip"}];
        PrintTemporary["Importing CrystalFieldTable from ", CFdataFilename, " ..."];
        CrystalFieldTable = ImportMZip[CFdataFilename];
        
        PrintTemporary["#------- numE = ", numE, " -------#"];
        exportFname = JJBlockMatrixFileName[numE, "FilenameAppendix" -> fileApp];
        fNames[numE] = exportFname;
        If[FileExistsQ[exportFname] && Not[overwrite],
          Continue[]
        ];
        JJBlockMatrixTable = TabulateJJBlockMatrixTable[numE, CrystalFieldTable, "Sparse"->OptionValue["Sparse"], "ChenDeltas" -> OptionValue["ChenDeltas"]];
        If[FileExistsQ[exportFname] && overwrite,
          DeleteFile[exportFname]
        ];
        ExportFun[exportFname, JJBlockMatrixTable];
        
        ClearAll[CrystalFieldTable];
      ),
    {numE, ns}
    ];
    Return[fNames];
  );

  EffectiveHamiltonian::usage = "EffectiveHamiltonian[numE] returns the Hamiltonian matrix for the f^numE configuration. The matrix is returned as a SparseArray. 
  The function admits an optional parameter \"FilenameAppendix\" which can be used to control which variant of the JJBlocks is used to assemble the matrix.
  It also admits an optional parameter \"IncludeZeeman\" which can be used to include the Zeeman interaction. The default is False.
  The option \"Set t2Switch\" can be used to toggle on or off setting the t2 selector automatically or not, the default is True, which replaces the parameter according to numE.
  The option \"ReturnInBlocks\" can be used to return the matrix in block or flattened form. The default is to return it in flattened form.";
  Options[EffectiveHamiltonian] = {
        "FilenameAppendix" -> "",
        "IncludeZeeman"    -> False,
        "Set t2Switch"     -> True,
        "ReturnInBlocks"   -> False,
        "OperatorBasis"    -> "Legacy"};
  EffectiveHamiltonian[nf_, OptionsPattern[]] := Module[
    {numE, ii, jj, howManyJs, Js, blockHam, opBasis},
    (
      (*#####################################*)
      ImportFun = ImportMZip;
      opBasis = OptionValue["OperatorBasis"];
      If[Not[MemberQ[{"Legacy", "MostlyOrthogonal", "Orthogonal"}, opBasis]],
        Echo["Operator basis " <> opBasis <> " not recognized, using \"Legacy\" basis."];
        opBasis = "Legacy";
      ];
      If[opBasis == "Orthogonal",
        Echo["Operator basis \"Orthogonal\" not implemented yet, aborting ..."];
        Return[Null];
      ];
      (*#####################################*)
      If[opBasis == "MostlyOrthogonal",
        (
          blockHam = EffectiveHamiltonian[nf, 
          "OperatorBasis" -> "Legacy", 
          "FilenameAppendix" -> OptionValue["FilenameAppendix"], 
          "IncludeZeeman" -> OptionValue["IncludeZeeman"], 
          "Set t2Switch" -> OptionValue["Set t2Switch"], 
          "ReturnInBlocks" -> OptionValue["ReturnInBlocks"]];
          paramChanger = Which[
            nf < 7,
            <|
              F0 -> 1/91 (54 E1p+91 E0p+78 \[Gamma]p),
              F2 -> (15/392 * 
                (
                  140   E1p  +
                  20020 E2p +
                  1540  E3p  +
                  770 \[Alpha]p -
                  70 \[Gamma]p  +
                  22 Sqrt[2] T2p -
                  11 Sqrt[2] nf T2p t2Switch -
                  11 Sqrt[2] (14-nf) T2p (1 - t2Switch)
                )
              ),
              F4 -> (99/490 *
                (
                  70   E1p - 
                  9100 E2p +
                  280  E3p +
                  140 \[Alpha]p - 
                  35 \[Gamma]p  + 
                  4 Sqrt[2] T2p -
                  2 Sqrt[2] nf T2p t2Switch -
                  2 Sqrt[2] (14-nf) T2p (1-t2Switch)
                )
                ),
              F6 -> (5577/7000 * 
                (
                  20  E1p +
                  700 E2p -
                  140 E3p -
                  70 \[Alpha]p -
                  10 \[Gamma]p -
                  2 Sqrt[2] T2p +
                  Sqrt[2] nf T2p t2Switch +
                  Sqrt[2] (14-nf) T2p (1-t2Switch)
                )
                ),
              \[Zeta]  -> \[Zeta],
              \[Alpha] -> (5 \[Alpha]p)/4,
              \[Beta]  -> -6  (5 \[Alpha]p + \[Beta]p),
              \[Gamma] -> 5/2 (2 \[Beta]p + 5 \[Gamma]p),
              T2 -> 0
            |>,
          nf >= 7,
            <|
              F0 -> 1/91 (54 E1p+91 E0p+78 \[Gamma]p),
              F2 -> (15/392 * 
                (
                  140   E1p  +
                  20020 E2p +
                  1540  E3p  +
                  770 \[Alpha]p -
                  70 \[Gamma]p  +
                  22 Sqrt[2] T2p -
                  11 Sqrt[2] nf T2p
                )
              ),
              F4 -> (99/490 *
                (
                  70   E1p - 
                  9100 E2p +
                  280  E3p +
                  140 \[Alpha]p - 
                  35 \[Gamma]p  + 
                  4 Sqrt[2] T2p -
                  2 Sqrt[2] nf T2p
                )
                ),
              F6 -> (5577/7000 * 
                (
                  20  E1p +
                  700 E2p -
                  140 E3p -
                  70 \[Alpha]p -
                  10 \[Gamma]p -
                  2 Sqrt[2] T2p +
                  Sqrt[2] nf T2p
                )
                ),
              \[Zeta]  -> \[Zeta],
              \[Alpha] -> (5 \[Alpha]p)/4,
              \[Beta]  -> -6  (5 \[Alpha]p + \[Beta]p),
              \[Gamma] -> 5/2 (2 \[Beta]p + 5 \[Gamma]p),
              T2 -> 0
            |>     
          ];
          blockHamMO = Which[
            OptionValue["ReturnInBlocks"] == False,
            ReplaceInSparseArray[blockHam, paramChanger],
            OptionValue["ReturnInBlocks"] == True,
            Map[ReplaceInSparseArray[#, paramChanger]&, blockHam, {2}]
          ];
          Return[blockHamMO];
        )
      ];
      (*#####################################*)
      
      (*hole-particle equivalence enforcement*)
      numE = nf;
      allVars = {E0, E1, E2, E3, \[Zeta], F0, F2, F4, F6, M0, M2, M4, T2, T2p,
        T3, T4, T6, T7, T8, P0, P2, P4, P6, gs, 
        \[Alpha], \[Beta], \[Gamma], B02, B04, B06, B12, B14, B16, 
        B22, B24, B26, B34, B36, B44, B46, B56, B66, S12, S14, S16, S22, 
        S24, S26, S34, S36, S44, S46, S56, S66, T11p, T12, T14, T15, T16, 
        T17, T18, T19, Bx, By, Bz};
      params0 = AssociationThread[allVars, allVars];
      If[nf > 7, 
        (
          numE = 14 - nf;
          params = HoleElectronConjugation[params0];
          If[OptionValue["Set t2Switch"], params[t2Switch] = 0];
        ),
        params = params0;
        If[OptionValue["Set t2Switch"], params[t2Switch] = 1];
      ];
      (* Load symbolic expressions for LS,J,J' energy sub-matrices. *)
      emFname = JJBlockMatrixFileName[numE, "FilenameAppendix" -> OptionValue["FilenameAppendix"]];
      JJBlockMatrixTable = ImportFun[emFname];
      (*Patch together the entire matrix representation using J,J' blocks.*)
      PrintTemporary["Patching JJ blocks ..."];
      Js        = AllowedJ[numE];
      howManyJs = Length[Js];
      blockHam  = ConstantArray[0, {howManyJs, howManyJs}];
      Do[
        blockHam[[jj, ii]] = JJBlockMatrixTable[{numE, Js[[ii]], Js[[jj]]}];,
      {ii, 1, howManyJs},
      {jj, 1, howManyJs}
      ];
      (* Once the block form is created flatten it *)
      If[Not[OptionValue["ReturnInBlocks"]],
        (blockHam = ArrayFlatten[blockHam];
        blockHam  = ReplaceInSparseArray[blockHam, params];
        ),
        (blockHam  = Map[ReplaceInSparseArray[#, params]&,blockHam,{2}];)
      ];
      
      If[OptionValue["IncludeZeeman"],
        (
          PrintTemporary["Including Zeeman terms ..."];
          {magx, magy, magz} = MagDipoleMatrixAssembly[numE, "ReturnInBlocks" -> OptionValue["ReturnInBlocks"]];
          blockHam += - teslaToKayser * (Bx * magx + By * magy + Bz * magz);
        )
      ];
      Return[blockHam];
    )
  ];

  SimplerEffectiveHamiltonian::usage = "SimplerEffectiveHamiltonian[numE, simplifier] is a simple addition to EffectiveHamiltonian that applies a given simplification to the full Hamiltonian. simplifier is a list of replacement rules. 
  If the option \"Export\" is set to True, then the function also exports the resulting sparse array to the ./hams/ folder.
  The option \"PrependToFilename\" can be used to prepend a string to the filename to which the function may export to.
  The option \"Return\" can be used to choose whether the function returns the matrix or not.
  The option \"Overwrite\" can be used to overwrite the file if it already exists, if this options is set to False then this function simply reloads a file that it assumed to be present already in the ./hams folder.
  The option \"IncludeZeeman\" can be used to toggle the inclusion of the Zeeman interaction with an external magnetic field.
  The option \"OperatorBasis\" can be used to choose the basis in which the operator is expressed. The default is the \"Legacy\" basis. Order alternatives being: \"MostlyOrthogonal\" and \"Orthogonal\". In the \"Legacy\" alternative the operators used are the same as in Carnall's work. In the \"MostlyOrthogonal\"  all operators are orthogonal except those corresponding to the Mk and Pk parameters. In the \"Orthogonal\" basis all operators are orthogonal, with the operators corresponding to the Mk and Pk parameters replaced by zi operators and accompanying ai coefficients. The \"Orthogonal\" option has not been implemented yet.";
  Options[SimplerEffectiveHamiltonian] = {
    "Export"->True, 
    "PrependToFilename"->"", 
    "EorF"->"F",
    "Overwrite" -> False,
    "Return" -> True,
    "Set t2Switch" -> False,
    "IncludeZeeman" -> False,
    "OperatorBasis" -> "Legacy"
    };
  SimplerEffectiveHamiltonian[numE_, simplifier_, OptionsPattern[]] := Module[
    {thisHam,fname, fnamemx},
    (
      If[Not[ValueQ[ElectrostaticTable]],
        LoadElectrostatic[]
        ];
      If[Not[ValueQ[SOOandECSOTable]],
        LoadSOOandECSO[]
      ];
      If[Not[ValueQ[SpinOrbitTable]],
        LoadSpinOrbit[]
      ];
      If[Not[ValueQ[SpinSpinTable]],
        LoadSpinSpin[]
      ];
      If[Not[ValueQ[ThreeBodyTable]],
        LoadThreeBody[]
      ];
      
      opBasis = OptionValue["OperatorBasis"];
      If[Not[MemberQ[{"Legacy", "MostlyOrthogonal", "Orthogonal"}, opBasis]],
        Echo["Operator basis " <> opBasis <> " not recognized, using \"Legacy\" basis."];
        opBasis = "Legacy";
      ];
      If[opBasis == "Orthogonal",
        Echo["Operator basis \"Orthogonal\" not implemented yet, aborting ..."];
        Return[Null];
      ];
      
      fnamePrefix = Which[
        opBasis == "Legacy",
        "SymbolicMatrix-f",
        opBasis == "MostlyOrthogonal",
        "SymbolicMatrix-mostly-orthogonal-f",
        opBasis == "Orthogonal",
        "SymbolicMatrix-orthogonal-f"
      ];
      fname   = FileNameJoin[{moduleDir, "hams",
        OptionValue["PrependToFilename"] <> fnamePrefix <> ToString[numE] <> ".m"}];
      fnamemx = FileNameJoin[{moduleDir, "hams", 
        OptionValue["PrependToFilename"] <> fnamePrefix <> ToString[numE] <> ".mx"}];
      fnamezip = FileNameJoin[{moduleDir, "hams", 
        OptionValue["PrependToFilename"] <> fnamePrefix <> ToString[numE] <> ".zip"}];
      If[Or[FileExistsQ[fname], 
            FileExistsQ[fnamemx],
            FileExistsQ[fnamezip]] 
        && Not[OptionValue["Overwrite"]],
        (
          If[OptionValue["Return"],
            (
              Which[
                FileExistsQ[fnamezip],
                (
                  Echo["File " <> fnamezip <> " already exists, and option \"Overwrite\" is set to False, loading file ..."];
                  thisHam = ImportMZip[fnamezip];
                  Return[thisHam];
                ),
                FileExistsQ[fnamemx],
                (
                  Echo["File " <> fnamemx <> " already exists, and option \"Overwrite\" is set to False, loading file ..."];
                  thisHam = Import[fnamemx];
                  Return[thisHam];
                ),
                FileExistsQ[fname],
                (
                  Echo["File " <> fname <> " already exists, and option \"Overwrite\" is set to False, loading file ..."];
                  thisHam = Import[fname];
                  Echo["Exporting to file " <> fnamemx <> " for quicker loading."];
                  Export[fnamemx, thisHam];
                  Return[thisHam];
                )
              ]
            ),
            (
              Echo["File " <> fname <> " already exists, skipping ..."];
              Return[Null];
            )
          ]
        )
      ];
      
      thisHam = EffectiveHamiltonian[numE, "Set t2Switch" -> OptionValue["Set t2Switch"], "IncludeZeeman" -> OptionValue["IncludeZeeman"], "OperatorBasis" -> opBasis];
      If[Length[simplifier] > 0,
        thisHam = ReplaceInSparseArray[thisHam, simplifier];
      ];
      (* This removes zero entries from being included in the sparse array *)
      thisHam = SparseArray[thisHam];
      If[OptionValue["Export"],
      (
        If[Not@FileExistsQ[fname],
        (
          Echo["Exporting to file " <> fname];
          Export[fname,   thisHam];
        )
        ];
        If[Not@FileExistsQ[fnamemx],
        (
          Echo["Exporting to file " <> fnamemx];
          Export[fnamemx, thisHam];
        )
        ];
      )
      ];
      If[OptionValue["Return"],
        Return[thisHam],
        Return[Null]
      ];
    )
  ];

  ScalarLSJMFromLS::usage = "ScalarLSJMFromLS[numE, LSReducedMatrixElements]. Given the LS-reduced matrix elements LSReducedMatrixElements of a scalar operator, this function returns the corresponding LSJM representation. This is returned as a SparseArray.";
  ScalarLSJMFromLS[numE_, LSReducedMatrixElements_] := Module[
    {jjBlocktable, NKSLJMs, NKSLJMps, J, Jp, eMatrix, SLterm, SpLpterm, 
    MJ, MJp, subKron, matValue, Js, howManyJs, blockHam, ii, jj},
    (
    SparseDiagonalArray[diagonalElements_] := SparseArray[
      Table[{i, i} -> diagonalElements[[i]],
        {i, 1, Length[diagonalElements]}]
      ];
    SparseZeroArray[width_, height_] := (
      SparseArray[
        Join[
        Table[{1, i} -> 0, {i, 1, width}],
        Table[{i, 1} -> 0, {i, 1, height}]
        ]
        ]
      );
    jjBlockTable = <||>;
    Do[
      NKSLJMs  = AllowedNKSLJMforJTerms[numE, J];
      NKSLJMps = AllowedNKSLJMforJTerms[numE, Jp];
      If[J != Jp,
        jjBlockTable[{numE, J, Jp}] = SparseZeroArray[Length[NKSLJMs], Length[NKSLJMps]];
      Continue[];
      ];
      eMatrix = Table[
        (* Condition for a scalar matrix op *)
        SLterm = NKSLJM[[1]];
        SpLpterm = NKSLJMp[[1]];
        MJ = NKSLJM[[3]];
        MJp = NKSLJMp[[3]];
        subKron = (KroneckerDelta[MJ, MJp]);
        matValue = If[subKron == 0, 
        0,
          (
            Which[MemberQ[Keys[LSReducedMatrixElements], {numE, SLterm, SpLpterm}],
              LSReducedMatrixElements[{numE, SLterm, SpLpterm}],
              MemberQ[Keys[LSReducedMatrixElements], {numE, SpLpterm, SLterm}],
              LSReducedMatrixElements[{numE, SpLpterm, SLterm}],
              True, 
              0
            ]
          )
        ];
        matValue,
        {NKSLJMp, NKSLJMps},
        {NKSLJM, NKSLJMs}
        ];
      jjBlockTable[{numE, J, Jp}] = SparseArray[eMatrix],
      {J, AllowedJ[numE]},
      {Jp, AllowedJ[numE]}
      ];
    
    Js = AllowedJ[numE];
    howManyJs = Length[Js];
    blockHam = ConstantArray[0, {howManyJs, howManyJs}];
    Do[blockHam[[jj, ii]] = jjBlockTable[{numE, Js[[ii]], Js[[jj]]}];,
      {ii, 1, howManyJs},
      {jj, 1, howManyJs}];
    blockHam = ArrayFlatten[blockHam];
    blockHam = SparseArray[blockHam];
    Return[blockHam];
    )
    ];

  (* ##################### Block assembly ###################### *)
  (* ########################################################### *)

  (* ########################################################### *)
  (* #################### Level Description #################### *)
  
  FreeHam::usage = "FreeHam[JJBlocks, numE] given the JJ blocks of the Hamiltonian for f^n, this function returns a list with all the scalar-simplified versions of the blocks.";
  FreeHam[JJBlocks_List, numE_Integer] := Module[
    {Js, basisJ, pivot, freeHam, idx, J,
    thisJbasis, shrunkBasisPositions, theBlock},
    (
      Js      = AllowedJ[numE];
      basisJ  = BasisLSJMJ[numE, "AsAssociation" -> True];
      pivot   = If[OddQ[numE], 1/2, 0];
      freeHam = Table[(
        J          = Js[[idx]];
        theBlock   = JJBlocks[[idx]];
        thisJbasis = basisJ[J];
        (* find the basis vectors that end with pivot *)
        shrunkBasisPositions = Flatten[Position[thisJbasis, {_ ..., pivot}]];
        (* take only those rows and columns *)
        theBlock[[shrunkBasisPositions, shrunkBasisPositions]]
      ),
      {idx, 1, Length[Js]}
      ];
      Return[freeHam];
    )
  ];
  
  ListRepeater::usage = "ListRepeater[list, reps] repeats each element of list reps times.";
  ListRepeater[list_List, repeats_Integer] := (
    Flatten[ConstantArray[#, repeats] & /@ list]
  );

  ListLever::usage = "ListLever[vecs, multiplicity] takes a list of vectors and returns all interleaved shifted versions of them.";
  ListLever[vecs_, multiplicity_] := Module[
  {uppityVecs, uppityVec},
  (
    uppityVecs = Table[(
      uppityVec = PadRight[{#}, multiplicity] & /@ vec;
      uppityVec = Permutations /@ uppityVec;
      uppityVec = Transpose[uppityVec];
      uppityVec = Flatten /@ uppityVec
      ),
    {vec, vecs}
    ];
    Return[Flatten[uppityVecs, 1]];
  )
  ];

  EigenLever::usage = "EigenLever[eigenSys, multiplicity] takes a list eigenSys of the form {eigenvalues, eigenvectors} and returns the eigenvalues repeated multiplicity times and the eigenvectors interleaved and shifted accordingly.";
  EigenLever[eigenSys_, multiplicity_] := Module[
    {eigenVals, eigenVecs,
    leveledEigenVecs, leveledEigenVals},
    (
      {eigenVals, eigenVecs} = eigenSys;
      leveledEigenVals       = ListRepeater[eigenVals, multiplicity];
      leveledEigenVecs       = ListLever[eigenVecs, multiplicity];
      Return[{Flatten[leveledEigenVals], leveledEigenVecs}]
    )
  ];

 
  LevelSimplerEffectiveHamiltonian::usage = "LevelSimplerEffectiveHamiltonian[numE] is a variation of EffectiveHamiltonian that returns the diagonal JJ Hamiltonian blocks applying a simplifier and with simplifications adequate for the level description. The keys of the given association correspond to the different values of J that are possible for f^numE, the values are sparse array that are meant to be interpreted in the basis provided by BasisLSJ.
  The option \"Simplifier\" is a list of symbols that are set to zero. At a minimum this has to include the crystal field parameters. By default this includes everything except the Slater parameters Fk and the spin orbit coupling \[Zeta].
  The option \"Export\" controls whether the resulting association is saved to disk, the default is True and the resulting file is saved to the ./hams/ folder. A hash is appended to the filename that corresponds to the simplifier used in the resulting expression. If the option \"Overwrite\" is set to False then these files may be used to quickly retrieve a previously computed case. The file is saved both in .m and .mx format.
  The option \"PrependToFilename\" can be used to append a string to the filename to which the function may export to.
  The option \"Return\" can be used to choose whether the function returns the matrix or not.
  The option \"Overwrite\" can be used to overwrite the file if it already exists.";
  Options[LevelSimplerEffectiveHamiltonian] = {
    "Export" -> True,
    "PrependToFilename" -> "",
    "Overwrite" -> False,
    "Return" -> True,
    "Simplifier" -> Join[
      {F0,\[Sigma]SS},
      cfSymbols,
      TSymbols,
      casimirSymbols,
      pseudoMagneticSymbols,
      marvinSymbols,
      DeleteCases[magneticSymbols,\[Zeta]]
    ]
  };
  LevelSimplerEffectiveHamiltonian[numE_Integer, OptionsPattern[]] := Module[
    {thisHamAssoc, Js, fname,
    fnamemx, hash, simplifier},
    (
      simplifier = (#->0)&/@Sort[OptionValue["Simplifier"]];
      hash       = Hash[simplifier];
      If[Not[ValueQ[ElectrostaticTable]],LoadElectrostatic[]];
      If[Not[ValueQ[SOOandECSOTable]],LoadSOOandECSO[]];
      If[Not[ValueQ[SpinOrbitTable]],LoadSpinOrbit[]];
      If[Not[ValueQ[SpinSpinTable]],LoadSpinSpin[]];
      If[Not[ValueQ[ThreeBodyTable]],LoadThreeBody[]];
      fname   = FileNameJoin[{moduleDir,"hams",OptionValue["PrependToFilename"]<>"Level-SymbolicMatrix-f"<>ToString[numE]<>"-"<>ToString[hash]<>".m"}];
      fnamemx = FileNameJoin[{moduleDir,"hams",OptionValue["PrependToFilename"]<>"Level-SymbolicMatrix-f"<>ToString[numE]<>"-"<>ToString[hash]<>".mx"}];
      If[Or[FileExistsQ[fname],FileExistsQ[fnamemx]]&&Not[OptionValue["Overwrite"]],
      (
        If[OptionValue["Return"],
        (
          Which[FileExistsQ[fnamemx],
          (
            Echo["File " <> fnamemx <> " already exists, and option \"Overwrite\" is set to False, loading file ..."];
            thisHamAssoc=Import[fnamemx];
            Return[thisHamAssoc];
          ),
          FileExistsQ[fname],
          (
            Echo["File " <> fname <> " already exists, and option \"Overwrite\" is set to False, loading file ..."];
            thisHamAssoc=Import[fname];
            Echo["Exporting to file " <> fnamemx <> " for quicker loading."];
            Export[fnamemx,thisHamAssoc];
            Return[thisHamAssoc];
          )
          ]
        ),
        (
          Echo["File " <> fname <> " already exists, skipping ..."];
          Return[Null];
        )
        ]
      )
      ];
      Js           = AllowedJ[numE];
      thisHamAssoc = EffectiveHamiltonian[numE,
        "Set t2Switch"->True,
        "IncludeZeeman"->False,
        "ReturnInBlocks"->True
      ];
      thisHamAssoc = Diagonal[thisHamAssoc];
      thisHamAssoc = Map[SparseArray[ReplaceInSparseArray[#,simplifier]]&,thisHamAssoc,{1}];
      thisHamAssoc = FreeHam[thisHamAssoc,numE];
      thisHamAssoc = AssociationThread[Js->thisHamAssoc];
      If[OptionValue["Export"],
        (
          Echo["Exporting to file " <> fname <> " and to " <> fnamemx];
          Export[fname,thisHamAssoc];
          Export[fnamemx,thisHamAssoc];
        )
      ];
      If[OptionValue["Return"],
        Return[thisHamAssoc],
        Return[Null]
      ];
    )
  ];

  LevelSolver::usage = "LevelSolver[numE, params] puts together (or retrieves from disk) the symbolic level Hamiltonian for the f^numE configuration and solves it for the given params returning the resultant energies and eigenstates.
  If the option \"Return as states\" is set to False, then the function returns an association whose keys are values for J in f^numE, and whose values are lists with two elements. The first element being equal to the ordered basis for the corresponding subspace, given as a list of lists of the form {LS string, J}. The second element being another list of two elements, the first element being equal to the energies and the second being equal to the corresponding normalized eigenvectors. The energies given have been subtracted the energy of the ground state.
  If the option \"Return as states\" is set to True, then the function returns a list with three elements. The first element is the global level basis for the f^numE configuration, given as a list of lists of the form {LS string, J}. The second element are the mayor LSJ components in the returned eigenstates. The third element is a list of lists with three elements, in each list the first element being equal to the energy, the second being equal to the value of J, and the third being equal to the corresponding normalized eigenvector (given as a row). The energies given have been subtracted the energy of the ground state, and the states have been sorted in order of increasing energy.
  The following options are admitted:
    - \"Overwrite Hamiltonian\", if set to True the function will overwrite the symbolic Hamiltonian. Default is False.
    - \"Return as states\", see description above. Default is True.
    - \"Simplifier\", this is a list with symbols that are set to zero for defining the parameters kept in the level  description.
  ";
  Options[LevelSolver] = {
    "Overwrite Hamiltonian" -> False,
    "Return as states" -> True,
    "Simplifier" -> Join[
      cfSymbols,
      TSymbols,
      casimirSymbols,
      pseudoMagneticSymbols,
      marvinSymbols,
      DeleteCases[magneticSymbols,\[Zeta]]
    ],
    "PrintFun" -> PrintTemporary
  };
  LevelSolver[numE_Integer, params0_Association, OptionsPattern[]] := Module[
    {ln, simplifier, simpleHam, basis,
    numHam, eigensys, startTime, endTime,
    diagonalTime, params=params0, globalBasis,
    eigenVectors, eigenEnergies, eigenJs,
    states, groundEnergy, allEnergies, PrintFun},
    (
      ln         = theLanthanides[[numE]];
      basis      = BasisLSJ[numE, "AsAssociation"->True];
      simplifier = OptionValue["Simplifier"];
      PrintFun   = OptionValue["PrintFun"];
      PrintFun["> LevelSolver for ",ln," with ",numE," f-electrons."];
      PrintFun["> Loading the symbolic level Hamiltonian ..."];
      simpleHam  = LevelSimplerEffectiveHamiltonian[numE,
        "Simplifier" -> simplifier,
        "Overwrite"  -> OptionValue["Overwrite Hamiltonian"]
      ];
      (* Everything that is not given is set to zero *)
      PrintFun["> Setting to zero every parameter not given ..."];
      params    = ParamPad[params, "PrintFun" -> PrintFun];
      PrintFun[params];
      (* Create the numeric hamiltonian *)
      PrintFun["> Replacing parameters in the J-blocks of the Hamiltonian to produce numeric arrays ..."];
      numHam    = N /@ Map[ReplaceInSparseArray[#, params]&, simpleHam];
      Clear[simpleHam];
      (* Eigensolver *)
      PrintFun["> Diagonalizing the numerical Hamiltonian within each separate J-subspace ..."];
      startTime    = Now;
      eigensys     = Eigensystem /@ numHam;
      endTime      = Now;
      diagonalTime = QuantityMagnitude[endTime-startTime,"Seconds"];
      allEnergies  = Flatten[First/@Values[eigensys]];
      groundEnergy = Min[allEnergies];
      eigensys     = Map[Chop[{#[[1]]-groundEnergy,#[[2]]}]&,eigensys];
      eigensys     = Association@KeyValueMap[#1->{basis[#1],#2} &,eigensys];
      PrintFun[">> Diagonalization took ",diagonalTime," seconds."];
      If[OptionValue["Return as states"],
        (
          PrintFun["> Padding the eigenvectors to correspond to the level basis ..."];
          eigenVectors  = SparseArray @ BlockDiagonalMatrix[Values[#[[2, 2]] & /@ eigensys]];
          globalBasis   = Flatten[Values[basis], 1];
          eigenEnergies = Flatten[Values[#[[2, 1]] & /@ eigensys]];
          eigenJs       = Flatten[KeyValueMap[ConstantArray[#1, Length[#2[[2, 2]]]] &, eigensys]];
          states        = Transpose[{eigenEnergies, eigenJs, eigenVectors}];
          states        = SortBy[states, First];
          eigenVectors  = Last /@ states;
          LSJmultiplets = (RemoveTrailingDigits[#[[1]]] <> ToString[InputForm[#[[2]]]]) & /@ globalBasis;
          majorComponentIndices = Ordering[Abs[#]][[-1]] & /@ eigenVectors;
          levelLabels           = LSJmultiplets[[majorComponentIndices]];
          Return[{globalBasis, levelLabels, states}];
        ),
        Return[{basis, eigensys}]
      ];
    )
  ]; 


  (* #################### Level Description #################### *)
  (* ########################################################### *)

  (* ########################################################### *)
  (* #################### Optical Operators #################### *)

  magOp = <||>;

  JJBlockMagDip::usage = "JJBlockMagDip[numE, J, Jp] returns an array for the LSJM matrix elements of the magnetic dipole operator between states with given J and Jp. The option \"Sparse\" can be used to return a sparse matrix. The default is to return a sparse matrix.
  See eqn 15.7 in TASS.
  Here it is provided in atomic units in which the Bohr magneton is 1/2.
  \[Mu] = -(1/2) (L + gs S)
  We are using the Racah convention for the reduced matrix elements in the Wigner-Eckart theorem. See TASS eqn 11.15.
  ";
  Options[JJBlockMagDip]={"Sparse"->True};
  JJBlockMagDip[numE_, braJ_, ketJ_, OptionsPattern[]] := Module[
    {braSLJs, ketSLJs,
    braSLJ,   ketSLJ,
    braSL,    ketSL,
    braS,     braL,
    ketS,     ketL,
    braMJ,    ketMJ,
    matValue, magMatrix,
    summand1, summand2,
    threejays},
    (
      braSLJs  = AllowedNKSLJMforJTerms[numE,braJ];
      ketSLJs  = AllowedNKSLJMforJTerms[numE,ketJ];
      magMatrix  = Table[
        braSL     = braSLJ[[1]];
        ketSL     = ketSLJ[[1]];
        {braS, braL}   = FindSL[braSL];
        {ketS, ketL}   = FindSL[ketSL];
        braMJ       = braSLJ[[3]];
        ketMJ       = ketSLJ[[3]];
        summand1    = If[Or[braJ  != ketJ,
                            braSL != ketSL],
          0,
          Sqrt[braJ*(braJ+1)*TPO[braJ]]
        ];
        (* looking at the string includes checking L=L', S=S', and \alpha=\alpha'*)
        summand2 = If[braSL!= ketSL,
          0,
          (gs-1) *
            Phaser[braS+braL+ketJ+1] *
              Sqrt[TPO[braJ]*TPO[ketJ]] *
                SixJay[{braJ,1,ketJ},{braS,braL,braS}] *
                  Sqrt[braS(braS+1)TPO[braS]]
        ];
        matValue = summand1 + summand2;
        (* We are using the Racah convention for red matrix elements in Wigner-Eckart *)
        threejays = (ThreeJay[{braJ, -braMJ}, {1, #}, {ketJ, ketMJ}] &) /@ {-1,0,1};
        threejays *= Phaser[braJ-braMJ];
        matValue  = - 1/2 * threejays * matValue;
        matValue,
      {braSLJ, braSLJs},
      {ketSLJ, ketSLJs}
      ];
      If[OptionValue["Sparse"],
        magMatrix = SparseArray[magMatrix]
      ];
      Return[magMatrix];
    )
  ]; 

  Options[TabulateJJBlockMagDipTable]={"Sparse"->True};
  TabulateJJBlockMagDipTable[numE_,OptionsPattern[]] := (
    JJBlockMagDipTable=<||>;
    Js=AllowedJ[numE];
    Do[
      (
        JJBlockMagDipTable[{numE,braJ,ketJ}]=
          JJBlockMagDip[numE,braJ,ketJ,"Sparse"->OptionValue["Sparse"]]
      ),
    {braJ, Js},
    {ketJ, Js}
    ];
    Return[JJBlockMagDipTable];
  );

  TabulateManyJJBlockMagDipTables::usage = "TabulateManyJJBlockMagDipTables[{n1, n2, ...}] calculates the tables of matrix elements for the requested f^n_i configurations. The function does not return the matrices themselves. It instead returns an association whose keys are numE and whose values are the filenames where the output of TabulateManyJJBlockMagDipTables was saved to. The output consists of an association whose keys are of the form {n, J, Jp} and whose values are rectangular arrays given the values of <|LSJMJa|H_dip|L'S'J'MJ'a'|>.";
  Options[TabulateManyJJBlockMagDipTables]={"FilenameAppendix"->"","Overwrite"->False,"Compressed"->True};
  TabulateManyJJBlockMagDipTables[ns_,OptionsPattern[]] := (
    fnames=<||>;
    Do[
    (
      ExportFun=If[OptionValue["Compressed"],ExportMZip,Export];
      PrintTemporary["#------- numE = ",numE," -------#"];
      appendTo     = (OptionValue["FilenameAppendix"]<>"-magDip");
      exportFname  = JJBlockMatrixFileName[numE,"FilenameAppendix"->appendTo];
      fnames[numE] = exportFname;
      If[FileExistsQ[exportFname]&&Not[OptionValue["Overwrite"]],
        Continue[]
      ];
      JJBlockMatrixTable = TabulateJJBlockMagDipTable[numE];
      If[FileExistsQ[exportFname] && OptionValue["Overwrite"],
        DeleteFile[exportFname]
      ];
      ExportFun[exportFname,  JJBlockMatrixTable];
    ),
    {numE,ns}
    ];
    Return[fnames];
  );

  MagDipoleMatrixAssembly::usage = "MagDipoleMatrixAssembly[numE] returns the matrix representation of the operator - 1/2 (L + gs S) in the f^numE configuration. The function returns a list with three elements corresponding to the x,y,z components of this operator. The option \"FilenameAppendix\" can be used to append a string to the filename from which the function imports from in order to patch together the array. For numE beyond 7 the function returns the same as for the complementary configuration. The option \"ReturnInBlocks\" can be used to return the matrices in blocks. The default is to return the matrices in flattened form and as sparse array.";
  Options[MagDipoleMatrixAssembly]={
    "FilenameAppendix"->"",
    "ReturnInBlocks"->False};
  MagDipoleMatrixAssembly[nf_Integer, OptionsPattern[]] := Module[
    {ImportFun, numE, appendTo,
    emFname, JJBlockMagDipTable, 
    Js, howManyJs, blockOp,
    rowIdx, colIdx},
    (
      ImportFun = ImportMZip;
      numE      = nf;
      numH      = 14 - numE;
      numE      = Min[numE,numH];
      
      appendTo  = (OptionValue["FilenameAppendix"]<>"-magDip");
      emFname   = JJBlockMatrixFileName[numE,"FilenameAppendix"->appendTo];
      JJBlockMagDipTable = ImportFun[emFname];
      
      Js        = AllowedJ[numE];
      howManyJs = Length[Js];
      blockOp   = ConstantArray[0,{howManyJs,howManyJs}];
      Do[
        blockOp[[rowIdx,colIdx]] = JJBlockMagDipTable[{numE,Js[[rowIdx]],Js[[colIdx]]}],
        {rowIdx,1,howManyJs},
        {colIdx,1,howManyJs}
      ];
      If[OptionValue["ReturnInBlocks"],
        (
          opMinus = Map[#[[1]]&, blockOp, {4}];
          opZero  = Map[#[[2]]&, blockOp, {4}];
          opPlus  = Map[#[[3]]&, blockOp, {4}];
          opX =   (opMinus - opPlus)/Sqrt[2];
          opY = I (opPlus + opMinus)/Sqrt[2];
          opZ = opZero;  
        ),
          blockOp = ArrayFlatten[blockOp];
          opMinus = blockOp[[;; , ;; , 1]];
          opZero  = blockOp[[;; , ;; , 2]];
          opPlus  = blockOp[[;; , ;; , 3]];
          opX =   (opMinus - opPlus)/Sqrt[2];
          opY = I (opPlus + opMinus)/Sqrt[2];
          opZ = opZero;  
      ];
      Return[{opX, opY, opZ}];
    )
  ];

  MagDipLineStrength::usage = "MagDipLineStrength[theEigensys, numE] takes the eigensystem of an ion and the number numE of f-electrons that correspond to it and calculates the line strength array Stot.
  The option \"Units\" can be set to either \"SI\" (so that the units of the returned array are (A m^2)^2) or to \"Hartree\".
  The option \"States\" can be used to limit the states for which the line strength is calculated. The default, All, calculates the line strength for all states. A second option for this is to provide an index labelling a specific state, in which case only the line strengths between that state and all the others are computed.
  The returned array should be interpreted in the eigenbasis of the Hamiltonian. As such the element Stot[[i,i]] corresponds to the line strength states between states |i> and |j>.";
  Options[MagDipLineStrength]={"Reload MagOp" -> False, "Units"->"SI", "States" -> All};
  MagDipLineStrength[theEigensys_List, numE0_Integer, OptionsPattern[]] := Module[
    {numE, allEigenvecs, Sx, Sy, Sz, Stot, factor},
    (
      numE = Min[14-numE0, numE0];
      (*If not loaded then load it, *)
      If[Or[
          Not[MemberQ[Keys[magOp], numE]],
          OptionValue["Reload MagOp"]],
      (
        magOp[numE] = ReplaceInSparseArray[#, {gs->2}]& /@ MagDipoleMatrixAssembly[numE];
      )
      ];
      allEigenvecs = Transpose[Last /@ theEigensys];
      Which[OptionValue["States"] === All,
        (
          {Sx, Sy, Sz} = (ConjugateTranspose[allEigenvecs].#.allEigenvecs) & /@ magOp[numE];
          Stot         = Abs[Sx]^2+Abs[Sy]^2+Abs[Sz]^2;
        ),
        IntegerQ[OptionValue["States"]],
        (
          singleState  = theEigensys[[OptionValue["States"],2]];
          {Sx, Sy, Sz} = (ConjugateTranspose[allEigenvecs].#.singleState) & /@ magOp[numE];
          Stot         = Abs[Sx]^2+Abs[Sy]^2+Abs[Sz]^2;
        )
      ];
      Which[
        OptionValue["Units"] == "SI",
          Return[4 \[Mu]B^2 * Stot],
        OptionValue["Units"] == "Hartree",
          Return[Stot],
        True,
        (
          Echo["Invalid option for \"Units\". Options are \"SI\" and \"Hartree\"."];
          Abort[];
        )
      ];
    )
  ];

  MagDipoleRates::usage = "MagDipoleRates[eigenSys, numE] calculates the magnetic dipole transition rate array for the provided eigensystem. The option \"Units\" can be set to \"SI\" or to \"Hartree\". If the option \"Natural Radiative Lifetimes\" is set to true then the reciprocal of the rate is returned instead. eigenSys is a list of lists with two elements, in each list the first element is the energy and the second one the corresponding eigenvector.
  Based on table 7.3 of Thorne 1999, using g2=1.
  The energy unit assumed in eigenSys is kayser.
  The returned array should be interpreted in the eigenbasis of the Hamiltonian. As such the element AMD[[i,i]] corresponds to the transition rate (or the radiative lifetime, depending on options) between eigenstates |i> and |j>. 
  By default this assumes that the refractive index is unity, this may be changed by setting the option \"RefractiveIndex\" to the desired value.
  The option \"Lifetime\" can be used to return the reciprocal of the transition rates. The default is to return the transition rates.";
  Options[MagDipoleRates]={"Units"->"SI", "Lifetime"->False, "RefractiveIndex"->1};
  MagDipoleRates[eigenSys_List, numE0_Integer,OptionsPattern[]] := Module[
    {AMD, Stot, eigenEnergies,
    transitionWaveLengthsInMeters, nRefractive},
    (
      nRefractive   = OptionValue["RefractiveIndex"];
      numE          = Min[14-numE0, numE0];
      Stot          = MagDipLineStrength[eigenSys, numE, "Units" -> OptionValue["Units"]];
      eigenEnergies = Chop[First/@eigenSys];
      energyDiffs   = Outer[Subtract, eigenEnergies, eigenEnergies];
      energyDiffs   = ReplaceDiagonal[energyDiffs, Indeterminate];
      (* Energies assumed in kayser.*)
      transitionWaveLengthsInMeters = 0.01/energyDiffs;
      
      unitFactor    = Which[
      OptionValue["Units"]=="Hartree",
      (
        (* The bohrRadius factor in SI needed to convert the wavelengths which are assumed in m*)
        16 \[Pi]^3 (\[Mu]0Hartree /(3 hPlanckHartree)) * bohrRadius^3
      ),
      OptionValue["Units"]=="SI",
      (
        16 \[Pi]^3 \[Mu]0/(3 hPlanck)
      ),
      True,
      (
        Echo["Invalid option for \"Units\". Options are \"SI\" and \"Hartree\"."];
        Abort[];
      )
      ];
      AMD = unitFactor / transitionWaveLengthsInMeters^3 * Stot * nRefractive^3;
      Which[OptionValue["Lifetime"],
        Return[1/AMD],
        True,
        Return[AMD]
      ]
    )
  ];

  GroundMagDipoleOscillatorStrength::usage = "GroundMagDipoleOscillatorStrength[eigenSys, numE] calculates the magnetic dipole oscillator strengths between the ground state and the excited states as given by eigenSys.
  Based on equation 8 of Carnall 1965, removing the 2J+1 factor since this degeneracy has been removed by the crystal field. 
  eigenSys is a list of lists with two elements, in each list the first element is the energy and the second one the corresponding eigenvector.
  The energy unit assumed in eigenSys is Kayser.
  The oscillator strengths are dimensionless.
  The returned array should be interpreted in the eigenbasis of the Hamiltonian. As such the element fMDGS[[i]] corresponds to the oscillator strength between ground state and eigenstate |i>.
  By default this assumes that the refractive index is unity, this may be changed by setting the option \"RefractiveIndex\" to the desired value.";
  Options[GroundMagDipoleOscillatorStrength]={"RefractiveIndex"->1};
  GroundMagDipoleOscillatorStrength[eigenSys_List, numE_Integer, OptionsPattern[]] := Module[
    {eigenEnergies, SMDGS, GSEnergy, energyDiffs,
    transitionWaveLengthsInMeters, unitFactor, nRefractive},
    (
      eigenEnergies    = First/@eigenSys;
      nRefractive      = OptionValue["RefractiveIndex"];
      SMDGS            = MagDipLineStrength[eigenSys,numE, "Units"->"SI", "States"->1];
      GSEnergy         = eigenSys[[1,1]];
      energyDiffs      = eigenEnergies-GSEnergy;
      energyDiffs[[1]] = Indeterminate;
      transitionWaveLengthsInMeters = 0.01/energyDiffs;
      unitFactor       = (8\[Pi]^2 me)/(3 hPlanck eCharge^2 cLight);
      fMDGS            = unitFactor / transitionWaveLengthsInMeters * SMDGS * nRefractive;
      Return[fMDGS];
    )
  ];

  UkOperator::usage = "UkOperator[numE, k] provides the matrix representation of the symmetric unit tensor operator U^k in the LSJMJ basis. The function returns an association with keys of the form {k,q} where q ranges from -k to k and  values representing the different components of the tensor operator. This is computed by using the doubly reduced matrix elements provided by ReducedUkTable with adequate coupling coefficients.";
  UkOperator[numE0_, k0_] := Module[
    {
      k = k0,
      basis,
      numE = numE0,
      opComponents,
      basisSize,
      val,
      coupling1,
      coupling2,
      phase,
      braLS,
      braJ,
      braM,
      ketLS,
      ketJ,
      ketM,
      braS,
      braL,
      ketS,
      ketL,
      q
      },
    (
      If[Not@ValueQ[ReducedUkTable],
      LoadUk[]
      ];
      basis = BasisLSJMJ[numE];
      basisSize = Length[basis];
      opComponents = <||>;
      numE = Min[14-numE, numE];
      Do[
      (
        nonZ = Reap[
          Do[(
            {braLS, braJ, braM} = basis[[rowIdx]];
            {ketLS, ketJ, ketM} = basis[[colIdx]];
            
            redUk = ReducedUkTable[{numE, 3, braLS, ketLS, k}];
            If[redUk == 0, Continue[]];
            
            coupling1 = ThreeJay[{ketJ, -ketM}, {k, q}, {braJ, braM}];
            If[coupling1 == 0, Continue[]];
            
            {braS, braL} = FindSL[braLS];
            {ketS, ketL} = FindSL[ketLS];
            
            coupling2 = Sqrt[TPO[braJ]*TPO[ketJ]] * SixJay[{ketL, ketJ, ketS}, {braJ, braL, k}];
            If[coupling2 == 0, Continue[]];
            
            phase = Phaser[braS + k + braJ + ketL + braJ - braM];
            val = N[phase * coupling1 * coupling2 * redUk];
            Sow[({rowIdx, colIdx} -> val)]
            ),
            {rowIdx, 1, basisSize},
            {colIdx, 1, basisSize}]
          ][[2, 1]];
        opComponents[{k, q}] = SparseArray[nonZ, {basisSize, basisSize}];
        ), 
      {q, -k, k}];
      Return[opComponents]
      )
    ];

  EffectiveElectricDipole::usage = "EffectiveElectricDipole[numE, Aparams] calculates the effective dipole operator in configuration f^numE for the given intensity parameters Aparams (given as an Association). The function returns an Association with three keys {-1,0,1} representing the three components of the operator and with values equal to arrays for the corresponding matrix representations in the standard LSJMJ basis. The units being equal to those of the units assumed for Aparams.";
  EffectiveElectricDipole[numE0_, Aparams0_] := Module[
    {
      Aparams = Aparams0,
      numE = numE0,
      Deff,
      keys,
      key,
      altKey,
      Aparam,
      qComponents
      },
    (
      keys = Keys[Aparams];
      Deff = <||>;
      qComponents = {-1, 0, 1};
      If[Not[ValueQ[Uks]],
        Uks = <||>
      ];
      Do[
      (
        sumTerms = Reap[
          Do[
            (
            key    = {\[Lambda], t, p};
            altKey = {\[Lambda], t, -p};
            (* the way the iterator is carried out, some CG are non-physical, leave them out early *)
            If[
              Or[Abs[p] > t,
              Abs[p + q] > \[Lambda]
              ],
            Continue[]
            ];
            (* enforce relation between conjugate values *)
            Aparam = Which[
              MemberQ[keys, key],
              Aparams[key],
              MemberQ[keys, altKey],
              Conjugate[Phaser[t + p + 1] * Aparams[altKey]],
              True,
              Continue[]
              ];
            (* calculate the Uk operator if it has not been calculated *)
            If[Not@MemberQ[Keys[Uks], \[Lambda]],
              Uks[\[Lambda]] = UkOperator[numE, \[Lambda]];
              ];
            
            clebscG = ClebschGordan[
              {\[Lambda], p + q}, 
              {1, -q}, 
              {t, p}];
            If[clebscG === 0, Continue[]];
            
            Sow[Phaser[q]*Aparam*clebscG*Uks[\[Lambda]][{\[Lambda], p + q}]];
            ),
            {\[Lambda], {2, 4, 6}},
            {t, {\[Lambda] - 1, \[Lambda], \[Lambda] + 1}},
            {p, -t, t, 1}
            ]
          ][[2, 1]];
        Deff[q] = Total[N@sumTerms];
        ),
      {q, qComponents}
      ];
      Return[Deff]
      )
    ];

  ElectricDipLineStrength::usage = "ElectricDipLineStrength[theEigensys, numE, Aparams] takes the eigensystem of an ion with configuration f^numE together with the intensity parameters A that define the effective electric dipole operator.
  The option \"Units\" can be set to either \"SI\" (so that the units of the returned array are (units of A)^2 * C^2) or to \"Hartree-partial\" in which case the units of the returned array are determined by the units of Aparams and equal to (units of Aparams)^2.
  The option \"States\" can be used to limit the states for which the line strength is calculated. The default, All, calculates the line strength for all states. A second option for this is to provide an index labelling a specific state, in which case only the line strengths between that state and all the others are computed.
  The returned array should be interpreted in the eigenbasis of the Hamiltonian. As such the element Sed[[i,i]] corresponds to the line strength states between states |i> and |j>.
  In what is returned no sum is made over degenerate states.
  ";
  Options[ElectricDipLineStrength] = {"Units" -> "SI","States" -> All}; 
  ElectricDipLineStrength[theEigensys_List, numE0_Integer, Aparams_Association, OptionsPattern[]] := Module[
    {
      numE = numE0,
      Deffcalc,
      allEigenvecs,
      Sed
    },
    (
      numE = Min[14-numE, numE];
      Deff = EffectiveElectricDipole[numE, Aparams];
      allEigenvecs = Transpose[Last /@ theEigensys];
      
      Which[OptionValue["States"] === All,
        (
          Deffcalc = (ConjugateTranspose[allEigenvecs] . # . allEigenvecs) & /@ Deff;
        ),
        IntegerQ[OptionValue["States"]],
        (
          singleState  = theEigensys[[OptionValue["States"],2]];
          Deffcalc = (ConjugateTranspose[allEigenvecs] . # . singleState) & /@ Deff; 
        )
      ];
      Sed  = (Abs[Deffcalc[-1]^2] + Abs[Deffcalc[0]^2] + Abs[Deffcalc[1]^2]);
      Which[
        OptionValue["Units"] == "SI",
        Return[eCharge^2 * Sed],
        OptionValue["Units"] == "Hartree-partial",
        Return[Sed]
      ];
    )
  ];

  (* #################### Optical Operators #################### *)
  (* ########################################################### *)

  (* ########################################################### *)
  (* ################### Printers and Labels ################### *)

  PrintL::usage = "PrintL[L] give the string representation of a given angular momentum.";
  PrintL[L_] := If[StringQ[L], L, StringTake[specAlphabet, {L + 1}]]

  FindSL::usage = "FindSL[LS] gives the spin and orbital angular momentum that corresponds to the provided string LS.";
  FindSL[SL_] := (
    FindSL[SL] = 
    If[StringQ[SL],
      {
        (ToExpression[StringTake[SL, 1]]-1)/2,
        StringPosition[specAlphabet, StringTake[SL, {2}]][[1, 1]]-1
      },
      SL
    ]
  );

  PrintSLJ::usage = "Given a list with three elements {S, L, J} this function returns a symbol where the spin multiplicity is presented as a superscript, the orbital angular momentum as its corresponding spectroscopic letter, and J as a subscript. Function does not check to see if the given J is compatible with the given S and L.";
  PrintSLJ[SLJ_] := (
    RowBox[{
      SuperscriptBox[" ", 2 SLJ[[1]] + 1], 
      SubscriptBox[PrintL[SLJ[[2]]], SLJ[[3]]]
      }
      ] // DisplayForm
    );

  PrintSLJM::usage = "Given a list with four elements {S, L, J, MJ} this function returns a symbol where the spin multiplicity is presented as a superscript, the orbital angular momentum as its corresponding spectroscopic letter, and {J, MJ} as a subscript. No attempt is made to guarantee that the given input is consistent.";
  PrintSLJM[SLJM_] := (
    RowBox[{
      SuperscriptBox[" ", 2 SLJM[[1]] + 1], 
      SubscriptBox[PrintL[SLJM[[2]]], {SLJM[[3]], SLJM[[4]]}]
      }
      ] // DisplayForm
    );

  (* ################### Printers and Labels ################### *)
  (* ########################################################### *)

  (* ########################################################### *)
  (* ##################### Term management ##################### *)

  AllowedSLTerms::usage = "AllowedSLTerms[numE] returns a list with the allowed terms in the f^numE configuration, the terms are given as lists in the format {S, L}. This list may have redundancies which are compatible with the degeneracies that might correspond to the given case.";
  AllowedSLTerms[numE_] := Map[FindSL[First[#]] &, CFPTerms[Min[numE, 14-numE]]];

  AllowedNKSLTerms::usage = "AllowedNKSLTerms[numE] returns a list with the allowed terms in the f^numE configuration, the terms are given as strings in spectroscopic notation. The integers in the last positions are used to distinguish cases with degeneracy.";
  AllowedNKSLTerms[numE_] := Map[First, CFPTerms[Min[numE, 14-numE]]];
  AllowedNKSLTerms[0]  = {"1S"};
  AllowedNKSLTerms[14] = {"1S"};

  MaxJ::usage = "MaxJ[numE] gives the maximum J = S+L that corresponds to the configuration f^numE.";
  MaxJ[numE_] := Max[Map[Total, AllowedSLTerms[Min[numE, 14-numE]]]];

  MinJ::usage = "MinJ[numE] gives the minimum J = S+L that corresponds to the configuration f^numE.";
  MinJ[numE_] := Min[Map[Abs[Part[#, 1] - Part[#, 2]] &, AllowedSLTerms[Min[numE, 14-numE]]]]

  AllowedSLJTerms::usage = "AllowedSLJTerms[numE] returns a list with the allowed {S, L, J} terms in the f^n configuration, the terms are given as lists in the format {S, L, J}. This list may have repeated elements which account for possible degeneracies of the related term.";
  AllowedSLJTerms[numE_] := Module[
    {idx1, allowedSL, allowedSLJ},
    (
      allowedSL = AllowedSLTerms[numE];
      allowedSLJ = {}; 
      For[
        idx1 = 1,
        idx1 <= Length[allowedSL],
        termSL = allowedSL[[idx1]];
        termsSLJ = 
          Table[
            {termSL[[1]], termSL[[2]], J},
          {J,Abs[termSL[[1]] - termSL[[2]]], Total[termSL]}
          ];
        allowedSLJ = Join[allowedSLJ, termsSLJ];
        idx1++
      ];
      SortBy[allowedSLJ, Last]
    )
  ];

  AllowedNKSLJTerms::usage = "AllowedNKSLJTerms[numE] returns a list with the allowed {SL, J} terms in the f^n configuration, the terms are given as lists in the format {SL, J} where SL is a string in spectroscopic notation.";
  AllowedNKSLJTerms[numE_] := Module[
    {allowedSL, allowedNKSL, allowedSLJ, nn},
    (
      allowedNKSL = AllowedNKSLTerms[numE];
      allowedSL = AllowedSLTerms[numE];
      allowedSLJ = {}; 
      For[
        nn = 1,
        nn <= Length[allowedSL],
        (
          termSL = allowedSL[[nn]];
          termNKSL = allowedNKSL[[nn]];
          termsSLJ = 
            Table[{termNKSL, J}, 
            {J, Abs[termSL[[1]] - termSL[[2]]], Total[termSL]}
            ];
          allowedSLJ = Join[allowedSLJ, termsSLJ];
          nn++
        )
      ];
      SortBy[allowedSLJ, Last]
    )
  ];

  AllowedNKSLforJTerms::usage = "AllowedNKSLforJTerms[numE, J] gives the terms that correspond to the given total angular momentum J in the f^n configuration. The result is a list whose elements are lists of length 2, the first element being the SL term in spectroscopic notation, and the second element being J.";
  AllowedNKSLforJTerms[numE_, J_] := Module[
    {allowedSL, allowedNKSL, allowedSLJ,
    nn, termSL, termNKSL, termsSLJ},
    (
      allowedNKSL = AllowedNKSLTerms[numE];
      allowedSL = AllowedSLTerms[numE];
      allowedSLJ = {}; 
      For[
        nn = 1,
        nn <= Length[allowedSL],
        (
          termSL = allowedSL[[nn]];
          termNKSL = allowedNKSL[[nn]];
          termsSLJ = If[Abs[termSL[[1]] - termSL[[2]]] <= J <= Total[termSL], 
            {{termNKSL, J}},
            {}
            ];
          allowedSLJ = Join[allowedSLJ, termsSLJ];
          nn++
        )
      ];
      Return[allowedSLJ]
    )
  ];

  AllowedSLJMTerms::usage = "AllowedSLJMTerms[numE] returns a list with all the states that correspond to the configuration f^n. A list is returned whose elements are lists of the form {S, L, J, MJ}.";
  AllowedSLJMTerms[numE_] := Module[
    {allowedSLJ, allowedSLJM,
    termSLJ, termsSLJM, nn},
    (
    allowedSLJ = AllowedSLJTerms[numE];
    allowedSLJM = {}; 
    For[
      nn = 1,
      nn <= Length[allowedSLJ],
      nn++,
      (
        termSLJ = allowedSLJ[[nn]];
        termsSLJM = 
          Table[{termSLJ[[1]], termSLJ[[2]], termSLJ[[3]], M},
          {M, - termSLJ[[3]], termSLJ[[3]]}
          ];
        allowedSLJM = Join[allowedSLJM, termsSLJM];
      )
    ];
    Return[SortBy[allowedSLJM, Last]];
    )
  ];

  AllowedNKSLJMforJMTerms::usage = "AllowedNKSLJMforJMTerms[numE, J, MJ] returns a list with all the terms that contain states of the f^n configuration that have a total angular momentum J, and a projection along the z-axis MJ. The returned list has elements of the form {SL (string in spectroscopic notation), J, MJ}.";
  AllowedNKSLJMforJMTerms[numE_, J_, MJ_] := Module[
    {allowedSL, allowedNKSL,
    allowedSLJM, nn},
    (
      allowedNKSL = AllowedNKSLTerms[numE];
      allowedSL   = AllowedSLTerms[numE];
      allowedSLJM = {}; 
      For[
        nn = 1,
        nn <= Length[allowedSL],
        termSL = allowedSL[[nn]];
        termNKSL = allowedNKSL[[nn]];
        termsSLJ = If[(Abs[termSL[[1]] - termSL[[2]]] 
                      <= J 
                      <= Total[termSL] 
                      && (Abs[MJ] <= J)
                      ),
                      {{termNKSL, J, MJ}},
                      {}];
        allowedSLJM = Join[allowedSLJM, termsSLJ];
        nn++
      ];
      Return[allowedSLJM];
    )
  ];

  AllowedNKSLJMforJTerms::usage = "AllowedNKSLJMforJTerms[numE, J] returns a list with all the states that have a total angular momentum J. The returned list has elements of the form {SL (string in spectroscopic notation), J, MJ}.";
  AllowedNKSLJMforJTerms[numE_, J_] := Module[
    {MJs, labelsAndMomenta, termsWithJ},
    (
      MJs = AllowedMforJ[J];
      (* Pair LS labels and their {S,L} momenta *)
      labelsAndMomenta = ({#, FindSL[#]}) & /@ AllowedNKSLTerms[numE];
      (* A given term will contain J if |L-S|<=J<=L+S *)
      ContainsJ[{SL_String, {S_, L_}}] := (Abs[S - L] <= J <= (S + L));
      (* Keep just the terms that satisfy this condition *)
      termsWithJ = Select[labelsAndMomenta, ContainsJ];
      (* We don't want to keep the {S,L} *)
      termsWithJ = {#[[1]], J} & /@ termsWithJ;
      (* This is just a quick way of including up all the MJ values *)
      Return[Flatten /@ Tuples[{termsWithJ, MJs}]]
    )
  ];

  AllowedMforJ::usage = "AllowedMforJ[J] is shorthand for Range[-J, J, 1].";
  AllowedMforJ[J_] := Range[-J, J, 1];
  
  AllowedJ::usage = "AllowedJ[numE] returns the total angular momenta J that appear in the f^numE configuration.";
  AllowedJ[numE_] := Table[J, {J, MinJ[numE], MaxJ[numE]}];

  Seniority::usage = "Seniority[LS] returns the seniority of the given term.";
  Seniority[LS_] := FindNKLSTerm[LS][[1, 2]];

  FindNKLSTerm::usage = "Given the string LS FindNKLSTerm[SL] returns all the terms that are compatible with it. This is only for f^n configurations. The provided terms might belong to more than one configuration. The function returns a list with elements of the form {LS, seniority, W, U}.";
  FindNKLSTerm[SL_] := Module[
    {NKterms, n},
    (
      n = 7;
      NKterms = {{}};
      Map[
        If[! StringFreeQ[First[#], SL], 
          If[ToExpression[Part[#, 2]] <= n, 
            NKterms = Join[NKterms, {#}, 1]
          ]
        ] &,
      fnTermLabels
      ];
      NKterms = DeleteCases[NKterms, {}];
      NKterms
    )
  ];

  ParseTermLabels::usage = "ParseTermLabels[] parses the labels for the terms in the f^n configurations based on the labels for the f6 and f7 configurations. The function returns a list whose elements are of the form {LS, seniority, W, U}.";
  Options[ParseTermLabels] = {"Export" -> True};
  ParseTermLabels[OptionsPattern[]] := Module[
    {labelsTextData, fNtextLabels, nielsonKosterLabels,
    seniorities, RacahW, RacahU},
    (
      labelsTextData = FileNameJoin[{moduleDir, "data", "NielsonKosterLabels_f6_f7.txt"}];
      fNtextLabels   = Import[labelsTextData];
      nielsonKosterLabels = Partition[StringSplit[fNtextLabels], 3];
      termLabels = Map[Part[#, {1}] &, nielsonKosterLabels];
      seniorities = Map[ToExpression[Part[# , {2}]] &, nielsonKosterLabels];
      racahW = 
        Map[
          StringTake[
            Flatten[StringCases[Part[# , {3}], 
              "(" ~~ DigitCharacter ~~ DigitCharacter ~~ DigitCharacter ~~ ")"]],
            {2, 4}
          ] &,
        nielsonKosterLabels];
      racahU = 
        Map[
          StringTake[
            Flatten[StringCases[Part[# , {3}], 
              "(" ~~ DigitCharacter ~~ DigitCharacter ~~ ")"]], 
            {2, 3}
          ] &,
        nielsonKosterLabels];
      fnTermLabels = Join[termLabels, seniorities, racahW, racahU, 2];
      fnTermLabels = Sort[fnTermLabels];
      If[OptionValue["Export"],
        (
          broadFname = FileNameJoin[{moduleDir,"data","fnTerms.m"}];
          Export[broadFname, fnTermLabels];
        )
      ];
      Return[fnTermLabels];
    )
  ];

  (* ##################### Term management ##################### *)
  (* ########################################################### *)

  LoadLaF3Parameters::usage="LoadLaF4Parameters[ln] gives the models for the semi-empirical Hamiltonian for the trivalent lanthanide ion with symbol ln.";
  Options[LoadLaF3Parameters] = {
    "Vintage"  -> "Standard",
    "With Uncertainties"->False
  };
  LoadLaF3Parameters[Ln_String, OptionsPattern[]]:= Module[
    {paramData, params},
    (
      paramData = Which[
        OptionValue["Vintage"] == "Carnall",
        Import[FileNameJoin[{moduleDir, "data", "carnallParams.m"}]],
        OptionValue["Vintage"] == "Standard",
        Import[FileNameJoin[{moduleDir, "data", "standardParams.m"}]] 
      ];
    params = If[OptionValue["With Uncertainties"],
      paramData[Ln],
      If[Head[#] === Around, #["Value"], #] & /@ paramData[Ln]
    ];
    Return[params];
    )
  ];  

  HoleElectronConjugation::usage = "HoleElectronConjugation[params] takes the parameters (as an association) that define a configuration and converts them so that they may be interpreted as corresponding to a complementary hole configuration. Some of this can be simply done by changing the sign of the model parameters. In the case of the effective three body interaction of T2 the relationship is more complex and is controlled by the value of the t2Switch variable.";
  HoleElectronConjugation[params_] := Module[
    {newparams = params},
    (
      flipSignsOf = Join[{\[Zeta]} , cfSymbols, TSymbols];
      flipped = Table[
        (
          flipper -> - newparams[flipper]
        ), 
        {flipper, flipSignsOf}
        ];
      nonflipped = Table[
        (
          flipper -> newparams[flipper]
        ),
        {flipper, Complement[Keys[newparams], flipSignsOf]}
        ];
      flippedParams = Association[Join[nonflipped, flipped]];
      flippedParams = Select[flippedParams, FreeQ[#, Missing]&];
      Return[flippedParams];
    )
  ];

  IonSolver::usage = "IonSolver[numE, params, host] puts together (or retrieves from disk) the symbolic Hamiltonian for the f^numE configuration and solves it for the given params. 
  params is an Association with keys equal to parameter symbols and values their numerical values. The function will replace the symbols in the symbolic Hamiltonian with their numerical values and then diagonalize the resulting matrix. Any parameter that is not defined in the params Association is assumed to be zero.
  host is an optional string that may be used to prepend the filename of the symbolic Hamiltonian that is saved to disk. The default is \"Ln\".
  The function returns the eigensystem as a list of lists where in each list the first element is the energy and the second element the corresponding eigenvector.
  The ordered basis in which these eigenvectors are to be interpreted is the one corresponding to BasisLSJMJ[numE].
  The function admits the following options:
  \"Include Spin-Spin\" (bool) : If True then the spin-spin interaction is included as a contribution to the m_k operators. The default is True.
  \"Overwrite Hamiltonian\" (bool) : If True then the function will overwrite the symbolic Hamiltonian that is saved to disk to expedite calculations. The default is False. The symbolic Hamiltonian is saved to disk to the ./hams/ folder preceded by the string host.
  \"Zeroes\" (list) : A list with symbols assumed to be zero.
  "; 
  Options[IonSolver] = {
    "Include Spin-Spin" -> True,
    "Overwrite Hamiltonian" -> False,
    "Zeroes" -> {}
  };
  IonSolver[numE_Integer, params0_Association, host_String:"Ln", OptionsPattern[]] := Module[
    {ln, simplifier, simpleHam, numHam, eigensys,
    startTime, endTime, diagonalTime,
    params=params0, zeroSymbols},
    (
      ln = theLanthanides[[numE]];
      
      (* This could be done when replacing values, but this produces smaller saved arrays. *)
      simplifier = (#-> 0) & /@ OptionValue["Zeroes"];
      simpleHam = SimplerEffectiveHamiltonian[numE, 
        simplifier,
        "PrependToFilename" -> host,
        "Overwrite" -> OptionValue["Overwrite Hamiltonian"]
      ];
      
      (* Note that we don't have to flip signs of parameters for fn beyond f7 since the matrix produced
      by SimplerEffectiveHamiltonian has already accounted for this. *)
      
      (* Everything that is not given is set to zero *)
      params = ParamPad[params];
      PrintFun[params];
      
      (* Enforce the override to the spin-spin contribution to the magnetic interactions *)
      params[\[Sigma]SS] = If[OptionValue["Include Spin-Spin"], 1, 0];
      
      (* Create the numeric hamiltonian *)
      numHam = ReplaceInSparseArray[simpleHam, params];
      Clear[simpleHam];
      
      (* Eigensolver *)
      PrintFun["> Diagonalizing the numerical Hamiltonian ..."];
      startTime = Now;
      eigensys  = Eigensystem[numHam];
      endTime   = Now;
      diagonalTime = QuantityMagnitude[endTime - startTime, "Seconds"];
      PrintFun[">> Diagonalization took ", diagonalTime, " seconds."];
      eigensys = Chop[eigensys];
      eigensys = Transpose[eigensys];
      
      (* Shift the baseline energy *)
      eigensys = ShiftedLevels[eigensys];
      (* Sort according to energy *)
      eigensys = SortBy[eigensys, First];
      Return[eigensys];
    )
  ];

  ShiftedLevels::usage = "ShiftedLevels[eigenSys] takes a list of levels of the form
  {{energy_1, coeff_vector_1}, {energy_2, coeff_vector_2},...}} and returns the same input except that now to every energy the minimum of all of them has been subtracted.";
  ShiftedLevels[originalLevels_] := Module[
    {groundEnergy, shifted},
    (
      groundEnergy = Sort[originalLevels][[1,1]];
      shifted      = Map[{#[[1]] - groundEnergy, #[[2]]} &, originalLevels];
      Return[shifted];
    )
  ];


  (* ########################################################### *)
  (* ############# Optical Transitions for Levels ############## *)

  JuddOfeltUkSquared::usage = "JuddOfeltUkSquared[numE, params] calculates the matrix elements of the Uk operator in the level basis. These are calculated according to equation (7) in Carnall 1965. 
  The function returns a list with the following elements:
    - basis : A list with the allowed {SL, J} terms in the f^numE configuration. Equal to BasisLSJ[numE].
    - eigenSys : A list with the eigensystem of the Hamiltonian for the f^n configuration.
    - levelLabels : A list with the labels of the major components of the level eigenstates.
    - LevelUkSquared : An association with the squared matrix elements of the Uk operators in the level eigenbasis. The keys being {2, 4, 6} corresponding to the rank of the Uk operator. The basis in which the matrix elements are given is the one corresponding to the level eigenstates given in eigenSys and whose major SLJ components are given in levelLabels. The matrix is symmetric and given as a SymmetrizedArray.
    The function admits the following options:
    \"PrintFun\" : A function that will be used to print the progress of the calculations. The default is PrintTemporary.";
  Options[JuddOfeltUkSquared] = {"PrintFun" -> PrintTemporary};
  JuddOfeltUkSquared[numE_, params_, OptionsPattern[]] := Module[
    {eigenChanger, numEH, basis, eigenSys,
    Js, Ukmat, LevelUkSquared, kRank,
    S, L, Sp, Lp, J, Jp, phase,
    braTerm, ketTerm, levelLabels,
    eigenVecs, majorComponentIndices},
    (
      If[Not[ValueQ[ReducedUkTable]],
        LoadUk[]
      ];
      numEH = Min[numE, 14-numE];
      PrintFun = OptionValue["PrintFun"];
      PrintFun["> Calculating the levels for the given parameters ..."];
      {basis, majorComponents, eigenSys} = LevelSolver[numE, params];
      (* The change of basis matrix to the eigenstate basis *)
      eigenChanger = Transpose[Last /@ eigenSys];
      PrintFun["Calculating the matrix elements of Uk in the physical coupling basis ..."];
      LevelUkSquared = <||>;
      Do[(
        Ukmat = Table[(
            {S, L}   = FindSL[braTerm[[1]]];
            J        = braTerm[[2]];
            Jp       = ketTerm[[2]];
            {Sp, Lp} = FindSL[ketTerm[[1]]];
            phase    = Phaser[S + Lp + J + kRank];
            Simplify @ (
              phase *
              Sqrt[TPO[J]*TPO[Jp]] *
              SixJay[{J, Jp, kRank}, {Lp, L, S}] *
              ReducedUkTable[{numEH, 3, braTerm[[1]], ketTerm[[1]], kRank}]
              )
            ),
          {braTerm, basis},
          {ketTerm, basis}
        ];
        Ukmat = (Transpose[eigenChanger] . Ukmat . eigenChanger)^2;
        Ukmat = Chop@Ukmat;
        LevelUkSquared[kRank] = SymmetrizedArray[Ukmat, Dimensions[eigenChanger], Symmetric[{1, 2}]];
        ),
        {kRank, {2, 4, 6}}
      ];
      LSJmultiplets = (RemoveTrailingDigits[#[[1]]] <> ToString[InputForm[#[[2]]]]) & /@ basis;
      eigenVecs = Last /@ eigenSys;
      majorComponentIndices = Ordering[Abs[#]][[-1]] & /@ eigenVecs;
      levelLabels = LSJmultiplets[[majorComponentIndices]];
      Return[{basis, eigenSys, levelLabels, LevelUkSquared}];
    )
  ];

  LevelElecDipoleOscillatorStrength::usage = "LevelElecDipoleOscillatorStrength[numE, levelParams, juddOfeltParams] uses Judd-Ofelt theory to estimate the forced electric dipole oscillator strengths ions whose level description is determined by levelParams. 
  The third parameter juddOfeltParams is an association with keys equal to the three Judd-Ofelt intensity parameters {\[CapitalOmega]2, \[CapitalOmega]4, \[CapitalOmega]6} and corresponding values in cm^2.
  The local field correction implemented here corresponds to the one given by the virtual cavity model of Lorentz.
  The function returns a list with the following elements:
    - basis : A list with the allowed {SL, J} terms in the f^numE configuration. Equal to BasisLSJ[numE].
    - eigenSys : A list with the eigensystem of the Hamiltonian for the f^n configuration in the level description.
    - levelLabels : A list with the labels of the major components of the calculated levels.
    - oStrengthArray : A square array whose elements represent the oscillator strengths between levels such that the element oStrengthArray[[i,j]] is the oscillator strength between the levels |Subscript[\[Psi], i]> and |Subscript[\[Psi], j]>. In this array, the elements below the diagonal represent emission oscillator strengths, and elements above the diagonal represent absorption oscillator strengths.
  The function admits the following three options:
    \"PrintFun\" : A function that will be used to print the progress of the calculations. The default is PrintTemporary.
    \"RefractiveIndex\" : The refractive index of the medium where the transitions are taking place. This may be a number or a function. If a number then the oscillator strengths are calculated for assuming a wavelength-independent refractive index. If a function then the refractive indices are calculated accordingly to the wavelength of each transition (the function must admit a single argument equal to the wavelength in nm). The default is 1.
    \"LocalFieldCorrection\" : The local field correction to be used. The default is \"VirtualCavity\". The options are: \"VirtualCavity\" and \"EmptyCavity\".
  The equation implemented here is the one given in eqn. 29 from the review article of Hehlen (2013). See that same article for a discussion on the local field correction.
  ";
  Options[LevelElecDipoleOscillatorStrength]={
    "PrintFun"        -> PrintTemporary,
    "RefractiveIndex" -> 1,
    "LocalFieldCorrection" -> "VirtualCavity"
  };
  LevelElecDipoleOscillatorStrength[numE_, levelParams_Association, juddOfeltParams_Association, OptionsPattern[]] := Module[
    {PrintFun, basis, eigenSys, levelLabels,
    LevelUkSquared, eigenEnergies, energyDiffs,
    oStrengthArray, nRef, \[Chi], nRefs,
    \[Chi]OverN, groundLevel, const,
    transitionFrequencies, wavelengthsInNM,
    fieldCorrectionType},
    (
      PrintFun = OptionValue["PrintFun"];
      nRef     = OptionValue["RefractiveIndex"];
      PrintFun["Calculating the Uk^2 matrix elements for the given parameters ..."];
      {basis, eigenSys, levelLabels, LevelUkSquared} = JuddOfeltUkSquared[numE, levelParams, "PrintFun"->PrintFun];
      eigenEnergies = First/@eigenSys;
      (* converted to cm^-2 *)
      const         = (8\[Pi]^2)/3 me/hPlanck * 10^(-4);
      energyDiffs   = Transpose@Outer[Subtract,eigenEnergies,eigenEnergies];
      (* since energies are assumed in Kayser, speed of light needs to be in cm/s, so that the frequencies are in 1/s *)
      transitionFrequencies = energyDiffs*cLight*100;
      (* grab the J for each level *)
      levelJs        = #[[2]] & /@ eigenSys;
      oStrengthArray = (
        juddOfeltParams[\[CapitalOmega]2]*LevelUkSquared[2]+
        juddOfeltParams[\[CapitalOmega]4]*LevelUkSquared[4]+
        juddOfeltParams[\[CapitalOmega]6]*LevelUkSquared[6]
      );
      oStrengthArray = Abs@(const * transitionFrequencies * oStrengthArray);
      (* it is necessary to divide each oscillator strength by the degeneracy of the initial level *)
      oStrengthArray = MapIndexed[1/(2 levelJs[[#2[[1]]]]+1) #1 &, oStrengthArray,{2}];
      (* including the effects of the refractive index *)
      fieldCorrectionType = OptionValue["LocalFieldCorrection"];
      Which[
        nRef === 1,
          True,
        NumberQ[nRef],
        (
          \[Chi] = Which[
            fieldCorrectionType == "VirtualCavity",
            (
              ( (nRef^2 + 2) / 3 )^2
            ),
            fieldCorrectionType == "EmptyCavity",
            (
              ( 3 * nRef^2 / ( 2 * nRef^2 + 1 ) )^2
            )
          ];
          \[Chi]OverN    = \[Chi] / nRef;
          oStrengthArray = \[Chi]OverN * oStrengthArray;
          (* the refractive index participates differently in absorption and in emission *)
          aFunction      = If[#2[[1]]>#2[[2]], #1 * nRef^2, #1]&;
          oStrengthArray = MapIndexed[aFunction, oStrengthArray, {2}];
        ),
        True,
        (
          wavelengthsInNM = Abs[1 / energyDiffs] * 10^7;
          nRefs           = Map[nRef, wavelengthsInNM];
          Echo["Calculating the oscillator strengths for the given refractive index ..."];
          \[Chi] = Which[
            fieldCorrectionType == "VirtualCavity",
            (
              ( (nRefs^2 + 2) / 3 )^2
            ),
            fieldCorrectionType == "EmptyCavity",
            (
              ( 3 * nRefs^2 / ( 2*nRefs^2 + 1 ) )^2
            )
          ];
          \[Chi]OverN     = \[Chi] / nRefs;
          oStrengthArray  = \[Chi]OverN * oStrengthArray
        )
      ];
      Return[{basis, eigenSys, levelLabels, oStrengthArray}];
    )
  ];

  LevelJJBlockMagDipole::usage = "LevelJJBlockMagDipole[numE, J, Jp] returns an array of the LSJ reduced matrix elements of the magnetic dipole operator between states with given J and Jp. The option \"Sparse\" can be used to return a sparse matrix. The default is to return a sparse matrix.";
  Options[LevelJJBlockMagDipole] = {"Sparse"->True};
  LevelJJBlockMagDipole[numE_, braJ_, ketJ_, OptionsPattern[]] := Module[
    {
      braSLJs, ketSLJs,
      braSLJ, ketSLJ,
      braSL, ketSL,
      braS, braL, 
      ketS, ketL,
      matValue, magMatrix,
      summand1, summand2
    },
  (
    braSLJs     = AllowedNKSLforJTerms[numE,braJ];
    ketSLJs     = AllowedNKSLforJTerms[numE,ketJ];
    magMatrix   = Table[
      (
        braSL       = braSLJ[[1]];
        ketSL       = ketSLJ[[1]];
        {braS,braL} = FindSL[braSL];
        {ketS,ketL} = FindSL[ketSL];
        summand1    = If[Or[braJ!=ketJ,braSL!=ketSL],
          0,
          Sqrt[braJ*(braJ+1)*TPO[braJ]
          ]
        ];
        (*looking at the string includes checking L=L',S=S',and\alpha=\alpha'*)
        summand2    = If[braSL!=ketSL,
          0,
          (gs-1)*
          Phaser[braS+braL+ketJ+1]*
          Sqrt[TPO[braJ]*TPO[ketJ]]*
          SixJay[{braJ,1,ketJ},{braS,braL,braS}]*
          Sqrt[braS(braS+1)TPO[braS]]
        ];
        matValue = summand1 + summand2;
        matValue = -1/2 * matValue;
        matValue
      ),
      {braSLJ,braSLJs},
      {ketSLJ,ketSLJs}
    ];
    If[OptionValue["Sparse"],
      magMatrix = SparseArray[magMatrix]];
    Return[magMatrix];
  )
  ];

  LevelMagDipoleMatrixAssembly::usage = "LevelMagDipoleMatrixAssembly[numE] puts together an array with the reduced matrix elements of the magnetic dipole operator in the level basis for the f^numE configuration. The function admits the two following options:
    \"Flattened\": If True then the returned matrix is flattened. The default is True.
    \"gs\": The electronic gyromagnetic ratio. The default is 2.";
  Options[LevelMagDipoleMatrixAssembly] = {
    "Flattened" -> True,
    gs -> 2
  };
  LevelMagDipoleMatrixAssembly[numE_, OptionsPattern[]] := Module[
    {Js, magDip, braJ, ketJ},
    (
      Js      = AllowedJ[numE];
      magDip  = Table[
        ReplaceInSparseArray[LevelJJBlockMagDipole[numE, braJ, ketJ], {gs -> OptionValue[gs]}],
        {braJ,Js},
        {ketJ,Js}
      ];
      If[OptionValue["Flattened"],
        magDip = ArrayFlatten[magDip];
      ];
      Return[magDip];
    )
  ];

  LevelMagDipoleLineStrength::usage = "LevelMagDipoleLineStrength[eigenSys, numE] calculates the magnetic dipole line strengths for an ion whose level description is determined by levelParams. The function returns a square array whose elements represent the magnetic dipole line strengths between the levels given in eigenSys such that the element magDipoleLineStrength[[i,j]] is the line strength between the levels |Subscript[\[Psi], i]> and |Subscript[\[Psi], j]>. Eigensys must be such that it consists of a lists of lists where in each list the last element corresponds to the eigenvector of a level (given as a row) in the standard basis for levels of the f^numE configuration.
  The function admits the following options:
    \"Units\" : The units in which the line strengths are given. The default is \"SI\". The options are \"SI\" and \"Hartree\". If \"SI\" then the unit of the line strength is (A m^2)^2 = (J/T)^2. If \"Hartree\" then the line strength is given in units of 2 \[Mu]B.";
  Options[LevelMagDipoleLineStrength] = {
    "Units"->"SI"
    };
  LevelMagDipoleLineStrength[theEigensys_List, numE0_Integer, OptionsPattern[]] := Module[
    {numE, levelMagOp, allEigenvecs, magDipoleLineStrength, units},
    (
      numE         = Min[14-numE0,numE0];
      levelMagOp   = LevelMagDipoleMatrixAssembly[numE];
      allEigenvecs = Transpose[Last/@theEigensys];
      units        = OptionValue["Units"];
      magDipoleLineStrength         = Transpose[allEigenvecs].levelMagOp.allEigenvecs;
      magDipoleLineStrength         = Abs[magDipoleLineStrength]^2;
      Which[
        units=="SI",
          Return[4 \[Mu]B^2 * magDipoleLineStrength],
        units=="Hartree",
          Return[magDipoleLineStrength]
      ];
    )
  ];

  LevelMagDipoleOscillatorStrength::usage = "LevelMagDipoleOscillatorStrength[eigenSys, numE] calculates the magnetic dipole oscillator strengths for an ion whose levels are described by eigenSys in configuration f^numE. The refractive index of the medium is relevant, but here it is assumed to be 1, this can be changed through the option \"RefractiveIndex\". eigenSys must consist of a lists of lists with three elements: the first element being the energy of the level, the second element being the J of the level, and the third element being the eigenvector of the level.
    The function returns a list with the following elements:
      - basis : A list with the allowed {SL, J} terms in the f^numE configuration. Equal to BasisLSJ[numE].
      - eigenSys : A list with the eigensystem of the Hamiltonian for the f^n configuration in the level description.
      - magDipoleOstrength : A square array whose elements represent the magnetic dipole oscillator strengths between the levels given in eigenSys such that the element magDipoleOstrength[[i,j]] is the oscillator strength between the levels |Subscript[\[Psi], i]> and |Subscript[\[Psi], j]>. In this array the elements below the diagonal represent emission oscillator strengths, and elements above the diagonal represent absorption oscillator strengths. The emission oscillator strengths are negative. The oscillator strength is a dimensionless quantity.
    The function admits the following option:
      \"RefractiveIndex\" : The refractive index of the medium where the transitions are taking place. This may be a number or a function. If a number then the oscillator strengths are calculated assuming a wavelength-independent refractive index as given. If a function then the refractive indices are calculated accordingly to the vaccum wavelength of each transition (the function must admit a single argument equal to the wavelength in nm). The default is 1.
    For reference see equation (27.8) in Rudzikas (2007). The expression for the line strenght is the simplest when using atomic units, (27.8) is missing a factor of \[Alpha]^2.";
  Options[LevelMagDipoleOscillatorStrength]={
    "RefractiveIndex" -> 1
  };
  LevelMagDipoleOscillatorStrength[eigenSys_, numE_, OptionsPattern[]] := Module[
    {eigenEnergies, eigenVecs, levelJs,
    energyDiffs, magDipoleOstrength, nRef,
    wavelengthsInNM, nRefs, degenDivisor},
    (
      basis         = BasisLSJ[numE];
      eigenEnergies = First/@eigenSys;
      nRef          = OptionValue["RefractiveIndex"];
      eigenVecs     = Last/@eigenSys;
      levelJs       = #[[2]]&/@eigenSys;
      energyDiffs   = -Outer[Subtract,eigenEnergies,eigenEnergies];
      energyDiffs  *= kayserToHartree;
      magDipoleOstrength = LevelMagDipoleLineStrength[eigenSys, numE, "Units"->"Hartree"];
      magDipoleOstrength = 2/3 * \[Alpha]Fine^2 * energyDiffs * magDipoleOstrength;
      degenDivisor       = #1 / ( 2 * levelJs[[#2[[1]]]] + 1 ) &;
      magDipoleOstrength = MapIndexed[degenDivisor, magDipoleOstrength, {2}];
      Which[nRef===1,
        True,
        NumberQ[nRef],
          (
            magDipoleOstrength = nRef * magDipoleOstrength;
          ),
        True,
          (
            wavelengthsInNM    = Abs[kayserToHartree / energyDiffs] * 10^7;
            nRefs              = Map[nRef, wavelengthsInNM];
            magDipoleOstrength = nRefs * magDipoleOstrength;
          )
      ];
      Return[{basis, eigenSys, magDipoleOstrength}];
    )
  ];

  LevelMagDipoleSpontaneousDecayRates::usage = "LevelMagDipoleSpontaneousDecayRates[eigenSys, numE] calculates the spontaneous emission rates for the magnetic dipole transitions between the levels given in eigenSys. The function returns a square array whose elements represent the spontaneous emission rates between the levels given in eigenSys such that the element [[i,j]] of the returned array is the rate of spontaneous emission from the level |Subscript[\[Psi], i]> to the level |Subscript[\[Psi], j]>. In this array the elements below the diagonal represent emission rates, and elements above the diagonal are identically zero.
  The function admits two optional arguments:
    + \"Units\" : The units in which the rates are given. The default is \"SI\". The options are \"SI\" and \"Hartree\". If \"SI\" then the rates are given in s^-1. If \"Hartree\" then the rates are given in the atomic unit of frequency.
    + \"RefractiveIndex\" : The refractive index of the medium where the transitions are taking place. This may be a number or a function. If a number then the rates are calculated assuming a wavelength-independent refractive index as given. If a function then the refractive indices are calculated accordingly to the vaccum wavelength of each transition (the function must admit a single argument equal to the wavelength in nm). The default is 1.";
  Options[LevelMagDipoleSpontaneousDecayRates] = {
    "Units" -> "SI",
    "RefractiveIndex" -> 1};
  LevelMagDipoleSpontaneousDecayRates[eigenSys_List, numE_Integer, OptionsPattern[]] := Module[
    {levMDlineStrength, eigenEnergies, energyDiffs,
    levelJs, spontaneousRatesInHartree, spontaneousRatesInSI,
    degenDivisor, units, nRef, nRefs, wavelengthsInNM},
    (
      nRef              = OptionValue["RefractiveIndex"];
      units             = OptionValue["Units"];
      levMDlineStrength = LowerTriangularize@LevelMagDipoleLineStrength[eigenSys,numE,"Units"->"Hartree"];
      levMDlineStrength = SparseArray[levMDlineStrength];
      eigenEnergies     = First /@ eigenSys;
      energyDiffs       = Outer[Subtract, eigenEnergies, eigenEnergies];
      energyDiffs       = kayserToHartree * energyDiffs;
      energyDiffs       = SparseArray[LowerTriangularize[energyDiffs]];
      levelJs           = #[[2]]& /@ eigenSys;
      spontaneousRatesInHartree = 4/3 \[Alpha]Fine^5 * energyDiffs^3 * levMDlineStrength;
      degenDivisor              = #1 / (2*levelJs[[#2[[1]]]] + 1)&;
      spontaneousRatesInHartree = MapIndexed[degenDivisor, spontaneousRatesInHartree, {2}];
      Which[nRef===1,
        True,
        NumberQ[nRef],
          (
            spontaneousRatesInHartree = nRef^3 * spontaneousRatesInHartree;
          ),
        True,
          (
            wavelengthsInNM    = Abs[kayserToHartree / energyDiffs] * 10^7;
            nRefs              = Map[nRef, wavelengthsInNM];
            spontaneousRatesInHartree = nRefs^3 * spontaneousRatesInHartree;
          )
      ];
      If[units == "SI",
        (
          spontaneousRatesInSI = 1/hartreeTime * spontaneousRatesInHartree;
          Return[SparseArray@spontaneousRatesInSI];
        ),
        Return[SparseArray@spontaneousRatesInHartree];
      ];
    )
  ];

  (* ############# Optical Transitions for Levels ############## *)
  (* ########################################################### *)

  (* ########################################################### *)
  (* ################## Eigensystem analysis ################### *)

  PrettySaundersSL::usage = "PrettySaundersSL[SL] produces a human-readeable symbol for the spectroscopic term SL. SL can be either a string (in RS notation for the term) or a list of two numbers {S, L}. The option \"Representation\" can be used to specify whether the output is given as a symbol or as a ket. The default is \"Ket\".";
  Options[PrettySaundersSL] = {"Representation" -> "Ket"};
  PrettySaundersSL[SL_, OptionsPattern[]] := (
    If[StringQ[SL],
    (
      {S,L} = FindSL[SL];
      L     = StringTake[SL,{2,-1}];
    ),
      {S,L}=SL
    ];
    pretty = RowBox[{
        AdjustmentBox[Style[2*S+1,Smaller], BoxBaselineShift->-1, BoxMargins->0],
        AdjustmentBox[PrintL[L]]
      }
    ];
    pretty = DisplayForm[pretty];
    pretty = Which[
      OptionValue["Representation"]=="Ket",
        Ket[{pretty}],
      OptionValue["Representation"]=="Symbol",
        pretty
    ];
    Return[pretty];
  );

  PrettySaundersSLJmJ::usage = "PrettySaundersSLJmJ[{SL, J, mJ}] produces a formatted symbol for the given basis vector {SL, J, mJ}.";
  Options[PrettySaundersSLJmJ] = {"Representation" -> "Ket"};
  PrettySaundersSLJmJ[{SL_, J_, mJ_}, OptionsPattern[]] := (If[
    StringQ[SL], 
    ({S, L} = FindSL[SL];
      L = StringTake[SL, {2, -1}];
      ), 
    {S, L} = SL];
    pretty = RowBox[{AdjustmentBox[Style[2*S + 1, Smaller], 
        BoxBaselineShift -> -1, BoxMargins -> 0], 
        AdjustmentBox[PrintL[L], BoxMargins -> -0.2], 
        AdjustmentBox[
        Style[Row[{InputForm[J], ",", mJ}], Small], 
        BoxBaselineShift -> 1, 
        BoxMargins -> {{0.7, 0}, {0.4, 0.4}}]}];
    pretty = DisplayForm[pretty];
    If[OptionValue["Representation"] == "Ket",
      pretty = Ket[{pretty}]
    ];
    Return[pretty];
  );

  PrettySaundersSLJ::usage = "PrettySaundersSLJ[{SL, J}] produces a formatted symbol for the given level {SL, J}. SL can be either a list of two numbers representing S and L or a string representing the spin multiplicity in spectroscopic notation. J is the total angular momentum to which S and L have coupled to. The option \"Representation\" can be used to specify whether the output is given as a symbol or as a ket. The default is \"Ket\".";
  Options[PrettySaundersSLJ] = {"Representation"->"Ket"};
  PrettySaundersSLJ[{SL_,J_},OptionsPattern[]] := (
    If[StringQ[SL],
    (
      {S,L}=FindSL[SL];
      L=StringTake[SL,{2,-1}];
    ),
      {S,L}=SL
    ];
    pretty = RowBox[{
        AdjustmentBox[Style[2*S+1,Smaller],BoxBaselineShift->-1,BoxMargins->0],
        AdjustmentBox[PrintL[L],BoxMargins->-0.2],
        AdjustmentBox[Style[InputForm[J],Small,FontTracking->"Narrow"],BoxBaselineShift->1,BoxMargins->{{0.7,0},{0.4,0.4}}]
      }
    ];
    pretty = DisplayForm[pretty];
    pretty = Which[
      OptionValue["Representation"]=="Ket",
        Ket[{pretty}],
      OptionValue["Representation"]=="Symbol",
        pretty
    ];
    Return[pretty];
  );

  BasisVecInRusselSaunders::usage = "BasisVecInRusselSaunders[basisVec] takes a basis vector in the format {LSstring, Jval, mJval} and returns a formatted symbol for the corresponding Russel-Saunders term.";
  BasisVecInRusselSaunders[basisVec_] := (
    {LSstring, Jval, mJval} = basisVec;
    Ket[PrettySaundersSLJmJ[basisVec]]
  );

  LSJMJTemplate = 
    StringTemplate[
    "\!\(\*TemplateBox[{\nRowBox[{\"`LS`\", \",\", \nRowBox[{\"J\", \
  \"=\", \"`J`\"}], \",\", \nRowBox[{\"mJ\", \"=\", \"`mJ`\"}]}]},\n\
  \"Ket\"]\)"];

  BasisVecInLSJMJ::usage = "BasisVecInLSJMJ[basisVec] takes a basis vector in the format {{{LSstring, Jval}, mJval}, nucSpin} and returns a formatted symbol for the corresponding LSJMJ term in the form |LS, J=..., mJ=...>.";
  BasisVecInLSJMJ[basisVec_] := (
    {LSstring, Jval, mJval} = basisVec;
    LSJMJTemplate[<|
      "LS" -> LSstring,
      "J" -> ToString[Jval, InputForm], 
      "mJ" -> ToString[mJval, InputForm]|>]
  );

  ParseStates::usage = "ParseStates[eigenSys, basis] takes a list of eigenstates in terms of their coefficients in the given basis and returns a list of the same states in terms of their energy, LSJMJ symbol, J, mJ, S, L, LSJ symbol, and LS symbol. eigenSys is a list of lists with two elements, in each list the first element is the energy and the second one the corresponding eigenvector. The LS symbol returned corresponds to the term with the largest coefficient in the given basis.";
  ParseStates[states_, basis_, OptionsPattern[]] := Module[
    {parsedStates},
    (
      parsedStates = Table[(
        {energy, eigenVec}      = state;
        maxTermIndex            = Ordering[Abs[eigenVec]][[-1]];
        {LSstring, Jval, mJval} = basis[[maxTermIndex]];
        LSJsymbol               = Subscript[LSstring, {Jval, mJval}];
        LSJMJsymbol             = LSstring <> ToString[Jval, InputForm];
        {S, L}                  = FindSL[LSstring];
        {energy, LSstring, Jval, mJval, S, L, LSJsymbol, LSJMJsymbol}
      ),
      {state, states}
      ];
      Return[parsedStates];
    )
  ];

  ParseStatesByNumBasisVecs::usage = "ParseStatesByNumBasisVecs[eigenSys, basis, numBasisVecs, roundTo] takes a list of eigenstates (given in eigenSys) in terms of their coefficients in the given basis and returns a list of the same states in terms of their energy and the coefficients at most numBasisVecs basis vectors. By default roundTo is 0.01 and this is the value used to round the amplitude coefficients. eigenSys is a list of lists with two elements, in each list the first element is the energy and the second one the corresponding eigenvector.
  The option \"Coefficients\" can be used to specify whether the coefficients are given as \"Amplitudes\" or \"Probabilities\". The default is \"Amplitudes\".
  ";
  Options[ParseStatesByNumBasisVecs] = {
    "Coefficients" -> "Amplitudes",
    "Representation" -> "Ket",
    "ReturnAs" -> "Dot"
    };
  ParseStatesByNumBasisVecs[eigensys_List, basis_List, numBasisVecs_Integer, roundTo_Real : 0.01, OptionsPattern[]] := Module[
    {parsedStates, energy, eigenVec, 
    probs, amplitudes, ordering, 
    returnAs, 
    chosenIndices, majorComponents, 
    majorAmplitudes, majorRep},
    (
      returnAs     = OptionValue["ReturnAs"];
      parsedStates = Table[(
        {energy, eigenVec} = state;
        energy             = Chop[energy];
        probs              = Round[Abs[eigenVec^2], roundTo];
        amplitudes         = Round[eigenVec, roundTo];
        ordering           = Ordering[probs];
        chosenIndices      = ordering[[-numBasisVecs ;;]];
        majorComponents    = basis[[chosenIndices]];
        majorThings = If[OptionValue["Coefficients"] == "Probabilities",
          (
            probs[[chosenIndices]]
          ),
          (
            amplitudes[[chosenIndices]]
          )
        ];
        majorComponents    = PrettySaundersSLJmJ[#, "Representation" -> OptionValue["Representation"]] & /@ majorComponents;
        nonZ               = (# != 0.) & /@ majorThings;
        majorThings        = Pick[majorThings, nonZ];
        majorComponents    = Pick[majorComponents, nonZ];
        If[OptionValue["Coefficients"] == "Probabilities",
          (
            majorThings = majorThings * 100 * "%";
          )
        ];
        majorRep  = Which[
                    returnAs == "Dot",
                      majorThings . majorComponents,
                    returnAs == "List",
                      Transpose[{Reverse@majorThings, Reverse@majorComponents}]
                    ];
        {energy, majorRep}
        ),
        {state, eigensys}];
      Return[parsedStates]
    )
  ];

  FindThresholdPosition::usage = "FindThresholdPosition[list, threshold] returns the position of the first element in list that is greater than or equal to threshold. If no such element exists, it returns the length of list. The elements of the given list must be in ascending order.";
  FindThresholdPosition[list_, threshold_] := Module[
    {position},
    (
    position = Position[list, _?(# >= threshold &), 1, 1];
    thrPos = If[Length[position] > 0,
      position[[1, 1]],
      Length[list]];
    If[thrPos == 0,
      Return[1],
      Return[thrPos]]
    )
  ]; 

  ParseStatesByProbabilitySum::usage = "ParseStatesByProbabilitySum[eigensys, basis, probSum] takes a list of eigenstates in terms of their coefficients in the given basis and returns a list of the same states in terms of their energy and the coefficients of the basis vectors that sum to at least probSum.";
  ParseStatesByProbabilitySum[eigensys_, basis_, probSum_, roundTo_ : 0.01, maxParts_: 20] := Module[
    {parsedByProb, numStates, state, energy,
    eigenVec, amplitudes, probs, ordering,
    orderedProbs, truncationIndex, accProb,
    thresholdIndex, chosenIndices, majorComponents,
    majorAmplitudes, absMajorAmplitudes, notnullAmplitudes, majorRep},
    (
      numStates    = Length[eigensys];
      parsedByProb = Table[(
        state              = eigensys[[idx]];
        {energy, eigenVec} = state;
        (*Round them up*)
        amplitudes         = Round[eigenVec, roundTo];
        probs              = Round[Abs[eigenVec^2], roundTo];
        ordering           = Reverse[Ordering[probs]];
        (*Order the probabilities from high to low*)
        orderedProbs       = probs[[ordering]];
        (*To speed up Accumulate, assume that only as much as maxParts will be needed*)
        truncationIndex    = Min[maxParts, Length[orderedProbs]];
        orderedProbs       = orderedProbs[[;;truncationIndex]];
        (*Accumulate the probabilities*)
        accProb            = Accumulate[orderedProbs];
        (*Find the index of the first element in accProb that is greater than probSum*)
        thresholdIndex     = Min[Length[accProb], FindThresholdPosition[accProb, probSum]];
        (*Grab all the indicees up till that one*)
        chosenIndices      = ordering[[;; thresholdIndex]];
        (*Select the corresponding elements from the basis*)
        majorComponents    = basis[[chosenIndices]];
        (*Select the corresponding amplitudes*)
        majorAmplitudes    = amplitudes[[chosenIndices]];
        (*Take their absolute value*)
        absMajorAmplitudes = Abs[majorAmplitudes];
        (*Make sure that there are no effectively zero contributions*)
        notnullAmplitudes  = Flatten[Position[absMajorAmplitudes, x_ /; x != 0]];
        (* majorComponents    = PrettySaundersSLJmJ[{#[[1]],#[[2]],#[[3]]}] & /@ majorComponents; *)
        majorComponents    = PrettySaundersSLJmJ /@ majorComponents;
        majorAmplitudes    = majorAmplitudes[[notnullAmplitudes]];
        majorComponents    = majorComponents[[notnullAmplitudes]];
        (*Multiply and add to build the final Ket*)
        majorRep           = majorAmplitudes . majorComponents;
        {energy, majorRep}
        ), {idx, numStates}];
      Return[parsedByProb];
    )
  ];

  (* ################## Eigensystem analysis ################### *)
  (* ########################################################### *)

  (* ########################################################### *)
  (* ########################## Misc ########################### *)

  SymbToNum::usage = "SymbToNum[expr, numAssociation] takes an expression expr and returns what results after making the replacements defined in the given replacementAssociation. If replacementAssociation doesn't define values for expected keys, they are taken to be zero.";
  SymbToNum[expr_, replacementAssociation_] := (
    includedKeys = Keys[replacementAssociation];
    (*If a key is not defined, make its value zero.*)
    fullAssociation = Table[(
      If[MemberQ[includedKeys, key],
        ToExpression[key]->replacementAssociation[key],
        ToExpression[key]->0
      ]
    ),
    {key, paramSymbols}];
    Return[expr/.fullAssociation];
  );

  SimpleConjugate::usage = "SimpleConjugate[expr] takes an expression and applies a simplified version of the conjugate in that all it does is that it replaces the imaginary unit I with -I. It assumes that every other symbol is real so that it remains the same under complex conjugation. Among other expressions it is valid for any rational or polynomial expression with complex coefficients and real variables.";
  SimpleConjugate[expr_] := expr /. Complex[a_, b_] :> a - I b;

  ExportMZip::usage = "ExportMZip[\"dest.[zip,m]\", expr] saves a compressed version of expr to the given destination.";
  ExportMZip[filename_, expr_] := Module[
    {baseName, exportName, mImportName, zipImportName},
    (
      baseName    = FileBaseName[filename];
      exportName  = StringReplace[filename,".m"->".zip"];
      mImportName = StringReplace[exportName,".zip"->".m"];
      If[FileExistsQ[mImportName],
      (
        PrintTemporary[mImportName<>" exists already, deleting"];
        DeleteFile[mImportName];
        Pause[2];
      )
      ];
      Export[exportName, (baseName<>".m") -> expr];
    )
  ];

  ImportMZip::usage = "ImportMZip[filename] imports a .m file inside a .zip file with corresponding filename. If the Option \"Leave Uncompressed\" is set to True (the default) then this function also leaves an umcompressed version of the object in the same folder of filename";
  Options[ImportMZip]={"Leave Uncompressed" -> True};
  ImportMZip[filename_String, OptionsPattern[]] := Module[
    {baseName, importKey, zipImportName, mImportName, mxImportName, imported},
    (
      baseName      = FileBaseName[filename];
      (*Function allows for the filename to be .m or .zip*)
      importKey     = baseName <> ".m";
      zipImportName = StringReplace[filename, ".m"->".zip"];
      mImportName   = StringReplace[zipImportName, ".zip"->".m"];
      mxImportName  = StringReplace[zipImportName, ".zip"->".mx"];
      Which[
        FileExistsQ[mxImportName],
      (
        PrintTemporary[".mx version exists already, importing that instead ..."];
        Return[Import[mxImportName]];
      ),
        FileExistsQ[mImportName],
      (
        PrintTemporary[".m version exists already, importing that instead ..."];
        imported = Import[mImportName];
        If[OptionValue["Leave Uncompressed"],
          (
            Export[mxImportName, imported];
          )
        ];
        Return[Import[mImportName]];
      )
      ];
      imported = Import[zipImportName, importKey];
      If[OptionValue["Leave Uncompressed"],
        (
          Export[mImportName, imported];
          Export[mxImportName, imported];
        )
      ];
      Return[imported];
    )
  ];

  ReplaceInSparseArray::usage = "ReplaceInSparseArray[sparseArray, rules] takes a sparse array that may contain symbolic quantities and returns a sparse array in which the given rules have been used on every element.";
  ReplaceInSparseArray[sparseA_SparseArray, rules_] := (
    SparseArray[Automatic,
      Dimensions[sparseA],
      sparseA["Background"] /. rules,
      {
        1,
        {sparseA["RowPointers"], sparseA["ColumnIndices"]},
        sparseA["NonzeroValues"] /. rules
      }
      ]
  );

  MapToSparseArray::usage = "MapToSparseArray[sparseArray, function] takes a sparse array and returns a sparse array after the function has been applied to it.";
  MapToSparseArray[sparseA_SparseArray, func_] := Module[
    {nonZ, backg, mapped},
    (
      nonZ   = func /@ sparseA["NonzeroValues"];
      backg  = func[sparseA["Background"]];
      mapped = SparseArray[Automatic, 
        Dimensions[sparseA], 
        backg, 
        {
          1,
          {sparseA["RowPointers"], sparseA["ColumnIndices"]},
          nonZ
        }
      ];
      Return[mapped];
    )
  ];

  ParseTeXLikeSymbol::usage = "ParseTeXLikeSymbol[string] parses a string for a symbol given in LaTeX notation and returns a corresponding mathematica symbol. The string may have expressions for several symbols, they need to be separated by single spaces. In addition the _ and ^ symbols used in LaTeX notation need to have arguments that are enclosed in parenthesis, for example \"x_2\" is invalid, instead \"x_{2}\" should have been given.";
  Options[ParseTeXLikeSymbol] = {"Form" -> "List"};
  ParseTeXLikeSymbol[bigString_, OptionsPattern[]] := Module[
    {form, mainSymbol, symbols},
    (
      form = OptionValue["Form"];
      (* parse greek *)
      symbols = Table[(
          str = StringReplace[string, {"\\alpha" -> "\[Alpha]",
            "\\beta" -> "\[Beta]",
            "\\gamma" -> "\[Gamma]",
            "\\psi" -> "\[Psi]"}];
          symbol = Which[
            StringContainsQ[str, "_"] && Not[StringContainsQ[str, "^"]],
            (
              (*yes sub no sup*)
              mainSymbol = StringSplit[str, "_"][[1]];
              mainSymbol = ToExpression[mainSymbol];
              
              subPart = 
                StringCases[str, 
                  RegularExpression@"\\{(.*?)\\}" -> "$1"][[1]];
              Subscript[mainSymbol, subPart]
            ),
            Not[StringContainsQ[str, "_"]] && StringContainsQ[str, "^"],
            (
              (*no sub yes sup*)
              mainSymbol = StringSplit[str, "^"][[1]];
              mainSymbol = ToExpression[mainSymbol];
              
              supPart = 
                StringCases[str, 
                  RegularExpression@"\\{(.*?)\\}" -> "$1"][[1]];
              Superscript[mainSymbol, supPart]
            ),
            StringContainsQ[str, "_"] && StringContainsQ[str, "^"],
            (
              (*yes sub yes sup*)
              mainSymbol = StringSplit[str, "_"][[1]];
              mainSymbol = ToExpression[mainSymbol];
              {subPart, supPart} = 
                StringCases[str, RegularExpression@"\\{(.*?)\\}" -> "$1"];
              Subsuperscript[mainSymbol, subPart, supPart]
            ),
            True,
            (
              (*no sup or sub*)
              str
            )
            ];
          symbol
          ),
        {string, StringSplit[bigString, " "]}
      ];
      Which[
        form == "Row",
        Return[Row[symbols]],
        form == "List",
        Return[symbols]
      ]
    )
  ];

  FromArrayToTable::usage = "FromArrayToTable[array, labels, energies] takes a square array of values and returns a table with the labels of the rows and columns, the energies of the initial and final levels, the level energies, the vacuum wavelength of the transition, and the value of the array. The array must be square and the labels and energies must be compatible with the order implied by the array. The array must be a square array of values. The function returns a list of lists with the following elements:
    - Initial level index
    - Final level index
    - Initial level label
    - Final level label
    - Initial level energy
    - Final level energy
    - Vacuum wavelength
    - Value of the array element.
    - The reciprocal of the value of the array element.
  Elements in which the array is zero are not included in the return of this function.";
  FromArrayToTable[array_,labels_,energies_] := Module[
    {tableFun, atl},
    (
      tableFun = {
        #2[[1]],
        #2[[2]],
        labels[[#2[[1]]]],
        labels[[#2[[2]]]],
        energies[[#2[[1]]]],
        energies[[#2[[2]]]],
        If[#2[[1]]==#2[[2]],"--",10^7/(energies[[#2[[1]]]]-energies[[#2[[2]]]])],
        #1
      }&;
      atl = Select[Flatten[MapIndexed[tableFun,array,{2}],1],#[[-1]]=!=0.&];
      atl = Append[#,1/#[[-1]]]&/@atl;
      Return[atl]
    )
  ];
  (* ########################## Misc ########################### *)
  (* ########################################################### *)

  (* ########################################################### *)
  (* ################# Some Plotting Routines ################## *)

  EnergyLevelDiagram::usage = "EnergyLevelDiagram[states] takes states and produces a visualization of its energy spectrum.
  The resultant visualization can be navigated by clicking and dragging to zoom in on a region, or by clicking and dragging horizontally while pressing Ctrl. Double-click to reset the view.";
  Options[EnergyLevelDiagram] = {
    "Title"->"", 
    "ImageSize"->1000, 
    "AspectRatio" -> 1/8, 
    "Background"->"Automatic",
    "Epilog"->{},
    "Explorer"->True,
    "Energy Unit" -> "cm^-1"
    };
  EnergyLevelDiagram[states_, OptionsPattern[]] := Module[
    {energies, epi, explora},
    (
      energies = If[Length[Dimensions[states]]==1,
        states,
        First/@states
      ];
      epi = OptionValue["Epilog"];
      explora = If[OptionValue["Explorer"],
        ExploreGraphics,
        Identity
      ];
      frameLabel = "E (" <> OptionValue["Energy Unit"] <> ")";
      plotTips = Which[
        OptionValue["Energy Unit"] == "cm^-1",
        Tooltip[{{#, 0}, {#, 1}}, {Quantity[#/8065.54429, "eV"], Quantity[#, 1/"Centimeters"]}] &/@ energies,
        OptionValue["Energy Unit"] == "eV",
        Tooltip[{{#, 0}, {#, 1}}, {Quantity[# * 8065.54429, 1/"Centimeters"], Quantity[#, "eV"]}] &/@ energies,
        True,
        Tooltip[{{#, 0}, {#, 1}}, Quantity[# , OptionValue["Energy Unit"]]] &/@ energies
      ];
      explora@ListPlot[plotTips,
        Joined       -> True,
        PlotStyle    -> Black,
        AspectRatio  -> OptionValue["AspectRatio"],
        ImageSize    -> OptionValue["ImageSize"],
        Frame        -> True,
        PlotRange    -> {All, {0, 1}},
        FrameTicks   -> {{None, None}, {Automatic, Automatic}},
        FrameStyle   -> Directive[15, Dashed, Thin],
        PlotLabel    -> Style[OptionValue["Title"], 15, Bold],
        Background   -> OptionValue["Background"],
        FrameLabel   -> {frameLabel},
        Epilog       -> epi]
    )
  ];

  ExploreGraphics::usage = "Pass a Graphics object to explore it. Zoom by clicking and dragging a rectangle. Pan by clicking and dragging while pressing Ctrl. Click twice to reset view. 
  Based on ZeitPolizei @ https://mathematica.stackexchange.com/questions/7142/how-to-manipulate-2d-plots.
  The option \"OptAxesRedraw\" can be used to specify whether the axes should be redrawn. The default is False.";
  Options[ExploreGraphics] = {OptAxesRedraw -> False};
  ExploreGraphics[graph_Graphics, opts : OptionsPattern[]] := With[
    {
      gr  = First[graph],
      opt = DeleteCases[Options[graph], 
            PlotRange -> PlotRange | AspectRatio | AxesOrigin -> _],
      plr = PlotRange /. AbsoluteOptions[graph, PlotRange],
      ar  = AspectRatio /. AbsoluteOptions[graph, AspectRatio],
      ao  = AbsoluteOptions[AxesOrigin],
      rectangle = {Dashing[Small], 
        Line[{#1, 
              {First[#2], Last[#1]},
              #2, 
              {First[#1], Last[#2]},
              #1}]} &,
      optAxesRedraw = OptionValue[OptAxesRedraw]
    },
    DynamicModule[
        {dragging=False, first, second, rx1, rx2, ry1, ry2,
        range = plr}, 
        {{rx1, rx2}, {ry1, ry2}} = plr;
      Panel@
      EventHandler[
        Dynamic@Graphics[
          If[dragging, {gr, rectangle[first, second]}, gr],
          PlotRange -> Dynamic@range,
          AspectRatio -> ar,
          AxesOrigin -> If[optAxesRedraw,
            Dynamic@Mean[range\[Transpose]], ao],
          Sequence @@ opt],
        {{"MouseDown", 1} :> (
          first = MousePosition["Graphics"]
          ),
        {"MouseDragged", 1} :> (
          dragging = True;
          second = MousePosition["Graphics"]
          ),
        "MouseClicked" :> (
          If[CurrentValue@"MouseClickCount"==2,
            range = plr];
          ),
        {"MouseUp", 1} :> If[dragging,
          dragging = False;
          
          range = {{rx1, rx2}, {ry1, ry2}} = 
            Transpose@{first, second};
          range[[2]] = {0, 1}],
        {"MouseDown", 2} :> (
          first = {sx1, sy1} = MousePosition["Graphics"]
          ),
        {"MouseDragged", 2} :> (
          second = {sx2, sy2} = MousePosition["Graphics"];
          rx1 = rx1 - (sx2 - sx1);
          rx2 = rx2 - (sx2 - sx1);
          ry1 = ry1 - (sy2 - sy1);
          ry2 = ry2 - (sy2 - sy1);
          range = {{rx1, rx2}, {ry1, ry2}};
          range[[2]] = {0, 1};
          )}]]];

  LabeledGrid::usage = "LabeledGrid[data, rowHeaders, columnHeaders] provides a grid of given data interpreted as a matrix of values whose rows are labeled by rowHeaders and whose columns are labeled by columnHeaders. When hovering with the mouse over the grid elements, the row and column labels are displayed with the given separator between them.";
  Options[LabeledGrid]={
      ItemSize->Automatic,
      Alignment->Center,
      Frame->All,
      "Separator"->",",
      "Pivot"->""
  };
  LabeledGrid[data_,rowHeaders_,columnHeaders_,OptionsPattern[]] := Module[
    {gridList=data, rowHeads=rowHeaders, colHeads=columnHeaders},
    (
      separator=OptionValue["Separator"];
      pivot=OptionValue["Pivot"];
      gridList=Table[
              Tooltip[
                data[[rowIdx,colIdx]],
                DisplayForm[
                  RowBox[{rowHeads[[rowIdx]],
                          separator,
                          colHeads[[colIdx]]}
                        ]
                        ]
                      ],
          {rowIdx,Dimensions[data][[1]]},
          {colIdx,Dimensions[data][[2]]}];
      gridList=Transpose[Prepend[gridList,colHeads]];
      rowHeads=Prepend[rowHeads,pivot];
      gridList=Prepend[gridList,rowHeads]//Transpose;
      Grid[gridList,
          Frame->OptionValue[Frame],
          Alignment->OptionValue[Alignment],
          Frame->OptionValue[Frame],
          ItemSize->OptionValue[ItemSize]
          ]
    )
  ];

  HamiltonianForm::usage = "HamiltonianForm[hamMatrix, basisLabels] takes the matrix representation of a hamiltonian together with a set of symbols representing the ordered basis in which the operator is represented. With this it creates a displayed form that has adequately labeled row and columns together with informative values when hovering over the matrix elements using the mouse cursor.";
  Options[HamiltonianForm]={"Separator"->"","Pivot"->""};
  HamiltonianForm[hamMatrix_, basisLabels_List, OptionsPattern[]] := (
      braLabels=DisplayForm[RowBox[{"\[LeftAngleBracket]",#,"\[RightBracketingBar]"}]]& /@ basisLabels;
      ketLabels=DisplayForm[RowBox[{"\[LeftBracketingBar]",#,"\[RightAngleBracket]"}]]& /@ basisLabels;
      LabeledGrid[hamMatrix,braLabels,ketLabels,"Separator"->OptionValue["Separator"],"Pivot"->OptionValue["Pivot"]]
      )

  HamiltonianMatrixPlot::usage = "HamiltonianMatrixPlot[hamMatrix, basisLabels] creates a matrix plot of the given hamiltonian matrix with the given basis labels. The matrix elements can be hovered over to display the corresponding row and column labels together with the value of the matrix element. The option \"OverlayValues\" can be used to specify whether the matrix elements should be displayed on top of the matrix plot.";
  Options[HamiltonianMatrixPlot] = Join[Options[MatrixPlot], {"Hover" -> True, "OverlayValues" -> True}];
  HamiltonianMatrixPlot[hamMatrix_, basisLabels_, opts : OptionsPattern[]] := (
    braLabels = DisplayForm[RowBox[{"\[LeftAngleBracket]", #, "\[RightBracketingBar]"}]] & /@ basisLabels;
    ketLabels = DisplayForm[Rotate[RowBox[{"\[LeftBracketingBar]", #, "\[RightAngleBracket]"}],\[Pi]/2]] & /@ basisLabels;
    ketLabelsUpright = DisplayForm[RowBox[{"\[LeftBracketingBar]", #, "\[RightAngleBracket]"}]] & /@ basisLabels;
    numRows = Length[hamMatrix];
    numCols = Length[hamMatrix[[1]]];
    epiThings = Which[
      And[OptionValue["Hover"], Not[OptionValue["OverlayValues"]]],
      Flatten[
        Table[
          Tooltip[
            {
              Transparent,
              Rectangle[
              {j - 1, numRows - i},
              {j - 1, numRows - i} + {1, 1}
              ]
            },
          Row[{braLabels[[i]],ketLabelsUpright[[j]],"=",hamMatrix[[i, j]]}]
          ],
        {i, 1, numRows},
        {j, 1, numCols}
        ]
      ],
      And[OptionValue["Hover"], OptionValue["OverlayValues"]],
      Flatten[
        Table[
          Tooltip[
            {
              Transparent,
              Rectangle[
              {j - 1, numRows - i},
              {j - 1, numRows - i} + {1, 1}
              ]
            },
          DisplayForm[RowBox[{"\[LeftAngleBracket]", basisLabels[[i]], "\[LeftBracketingBar]", basisLabels[[j]], "\[RightAngleBracket]"}]]
          ],
        {i, numRows},
        {j, numCols}
        ]
      ],
      True,
      {}
      ];
    textOverlay = If[OptionValue["OverlayValues"],
      (
        Flatten[
        Table[
          Text[hamMatrix[[i, j]],
            {j - 1/2, numRows - i + 1/2}
          ],
        {i, 1, numRows},
        {j, 1, numCols}
        ]
        ]
      ),
      {}
      ];
    epiThings = Join[epiThings, textOverlay];
    MatrixPlot[hamMatrix,
      FrameTicks -> {
        {Transpose[{Range[Length[braLabels]], braLabels}], None}, 
        {None, Transpose[{Range[Length[ketLabels]], ketLabels}]}
        },
      Evaluate[FilterRules[{opts}, Options[MatrixPlot]]],
      Epilog -> epiThings
    ]
  );

  (* ################# Some Plotting Routines ################## *)
  (* ########################################################### *)

  (* ########################################################### *)
  (* ##################### Load Functions ###################### *)
  
  LoadAll::usage = "LoadAll[] executes most Load* functions.";
  LoadAll[] := (
    LoadTermLabels[];
    LoadCFP[];
    LoadUk[];
    LoadV1k[];
    LoadT22[];
    LoadSOOandECSOLS[];
    
    LoadElectrostatic[];
    LoadSpinOrbit[];
    LoadSOOandECSO[];
    LoadSpinSpin[];
    LoadThreeBody[];
    LoadChenDeltas[];
    LoadCarnall[];
  );

  fnTermLabels::usage = "This list contains the labels of f^n configurations. Each element of the list has four elements {LS, seniority, W, U}. At first sight this seems to only include the labels for the f^6 and f^7 configuration, however, all is included in these two.";

  LoadTermLabels::usage = "LoadTermLabels[] loads into the session the labels for the terms in the f^n configurations.";
  LoadTermLabels[] := (
    If[ValueQ[fnTermLabels], Return[]];
    PrintTemporary["Loading data for state labels in the f^n configurations..."];
    fnTermsFname = FileNameJoin[{moduleDir, "data", "fnTerms.m"}];
    
    If[!FileExistsQ[fnTermsFname],
      (PrintTemporary[">> fnTerms.m not found, generating ..."];
        fnTermLabels = ParseTermLabels["Export"->True];
      ),
      fnTermLabels = Import[fnTermsFname];
    ];
  );

  Carnall::usage = "Association of data from Carnall et al (1989) with the following keys: {data, annotations, paramSymbols, elementNames, rawData, rawAnnotations, annnotatedData, appendix:Pr:Association, appendix:Pr:Calculated, appendix:Pr:RawTable, appendix:Headings}";

  LoadCarnall::usage = "LoadCarnall[] loads data for trivalent lanthanides in LaF3 using the data from Bill Carnall's 1989 paper.";
  LoadCarnall[] := (
    If[ValueQ[Carnall], Return[]];
    carnallFname = FileNameJoin[{moduleDir, "data", "Carnall.m"}];
    If[!FileExistsQ[carnallFname],
      (PrintTemporary[">> Carnall.m not found, generating ..."];
        Carnall = ParseCarnall[];
      ),
      Carnall = Import[carnallFname];
    ];
  );

  LoadChenDeltas::usage = "LoadChenDeltas[] loads the differences noted by Chen.";
  LoadChenDeltas[] := (
    If[ValueQ[chenDeltas], Return[]];
    PrintTemporary["Loading the association of discrepancies found by Chen ..."];
    chenDeltasFname = FileNameJoin[{moduleDir, "data", "chenDeltas.m"}];
    If[!FileExistsQ[chenDeltasFname],
      (PrintTemporary[">> chenDeltas.m not found, generating ..."];
        chenDeltas = ParseChenDeltas[];
      ),
      chenDeltas = Import[chenDeltasFname];
    ];
  );

  ParseChenDeltas::usage = "ParseChenDeltas[] parses the data found in ./data/the-chen-deltas-A.csv and ./data/the-chen-deltas-B.csv. If the option \"Export\" is set to True (True is the default), then the parsed data is saved to ./data/chenDeltas.m";
  Options[ParseChenDeltas] = {"Export" -> True};
  ParseChenDeltas[OptionsPattern[]] := (
    chenDeltasRaw = Import[FileNameJoin[{moduleDir, "data", "the-chen-deltas-A.csv"}]];
    chenDeltasRaw = chenDeltasRaw[[2 ;;]];
    chenDeltas = <||>;
    chenDeltasA = <||>;
    Off[Power::infy];
    Do[
      ({right, wrong} = {chenDeltasRaw[[row]][[4 ;;]], 
        chenDeltasRaw[[row + 1]][[4 ;;]]};
      key = chenDeltasRaw[[row]][[1 ;; 3]];
      repRule = (#[[1]] -> #[[2]]*#[[1]]) & /@ 
        Transpose[{{M0, M2, M4, P2, P4, P6}, right/wrong}];
      chenDeltasA[key] = <|"right" -> right, "wrong" -> wrong, 
        "repRule" -> repRule|>;
      chenDeltasA[{key[[1]], key[[3]], key[[2]]}] = <|"right" -> right, 
        "wrong" -> wrong, "repRule" -> repRule|>;
      ),
      {row, 1, Length[chenDeltasRaw], 2}];
    chenDeltas["A"] = chenDeltasA;
    
    chenDeltasRawB = Import[FileNameJoin[{moduleDir, "data", "the-chen-deltas-B.csv"}], "Text"];
    chenDeltasB = StringSplit[chenDeltasRawB, "\n"];
    chenDeltasB = StringSplit[#, ","] & /@ chenDeltasB;
    chenDeltasB = {ToExpression[StringTake[#[[1]], {2}]], #[[2]], #[[3]]} & /@ chenDeltasB;
    chenDeltas["B"] = chenDeltasB;
    On[Power::infy];
    If[OptionValue["Export"],
      (chenDeltasFname = FileNameJoin[{moduleDir, "data", "chenDeltas.m"}];
      Export[chenDeltasFname, chenDeltas];
      )
      ];
    Return[chenDeltas];
  );

  ParseCarnall::usage = "ParseCarnall[] parses the data found in ./data/Carnall.xls. If the option \"Export\" is set to True (True is the default), then the parsed data is saved to ./data/Carnall. This data is from the tables and appendices of Carnall et al (1989).";
  Options[ParseCarnall] = {"Export" -> True};
  ParseCarnall[OptionsPattern[]] := (
    ions         = {"Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb"};
    templates    = StringTemplate/@StringSplit["appendix:`ion`:Association appendix:`ion`:Calculated appendix:`ion`:RawTable appendix:`ion`:Headings"," "];
    
    (* How many unique eigenvalues, after removing Kramer's degeneracy *)
    fullSizes    = AssociationThread[ions, {7, 91, 182, 1001, 1001, 3003, 1716, 3003, 1001, 1001, 182, 91, 7}];
    carnall      = Import[FileNameJoin[{moduleDir,"data","Carnall.xls"}]][[2]];
    carnallErr   = Import[FileNameJoin[{moduleDir,"data","Carnall.xls"}]][[3]];
    
    elementNames = carnall[[1]][[2;;]];
    carnall      = carnall[[2;;]];
    carnallErr   = carnallErr[[2;;]];
    carnall      = Transpose[carnall];
    carnallErr   = Transpose[carnallErr];
    paramNames   = ToExpression/@carnall[[1]][[1;;]];
    carnall      = carnall[[2;;]];
    carnallErr   = carnallErr[[2;;]];
    carnallData  = Table[(
                    data         = carnall[[i]];
                    data         = (#[[1]]->#[[2]])&/@Select[Transpose[{paramNames,data}],#[[2]]!=""&];
                    elementNames[[i]]->data
                    ),
                    {i,1,13}
                    ];
    carnallData  = Association[carnallData];
    carnallNotes = Table[(
                    data        = carnallErr[[i]];
                    elementName = elementNames[[i]];
                    dataFun     = (
                        #[[1]] -> If[#[[2]]=="[]",
                        "Not allowed to vary in fitting.",
                        If[#[[2]]=="[R]",
                            "Ratio constrained by: " <> <|"Eu"->"F4/F2=0.713; F6/F2=0.512",
                                "Gd"->"F4/F2=0.710]",
                                "Tb"->"F4/F2=0.707"|>[elementName],
                            If[#[[2]]=="i",
                                "Interpolated",
                                #[[2]]
                            ]
                        ]
                        ]) &;
                    data = dataFun /@ Select[Transpose[{paramNames,data}],#[[2]]!=""&];
                    elementName->data
                      ),
                  {i,1,13}
                  ];
    carnallNotes  = Association[carnallNotes];
    
    annotatedData = Table[
                    If[NumberQ[#[[1]]],Tooltip[#[[1]],#[[2]]],""] & /@ Transpose[{paramNames/.carnallData[element],
                      paramNames/.carnallNotes[element]
                        }],
                    {element,elementNames}
                    ];
    annotatedData = Transpose[annotatedData];
    
    Carnall = <|"data"      -> carnallData,
        "annotations"       -> carnallNotes,
        "paramSymbols"      -> paramNames,
        "elementNames"      -> elementNames,
        "rawData"           -> carnall,
        "rawAnnotations"    -> carnallErr,
        "includedTableIons" -> ions,
        "annnotatedData"    -> annotatedData
    |>;
    
    Do[(
        carnallData   = Import[FileNameJoin[{moduleDir,"data","Carnall.xls"}]][[sheetIdx]];
        headers       = carnallData[[1]];
        calcIndex     = Position[headers,"Calc (1/cm)"][[1,1]];
        headers       = headers[[2;;]];
        carnallLabels = carnallData[[1]];
        carnallData   = carnallData[[2;;]];
        carnallTerms  = DeleteDuplicates[First/@carnallData];
        parsedData    = Table[(
                          rows = Select[carnallData,#[[1]]==term&];
                          rows = #[[2;;]]&/@rows;
                          rows = Transpose[rows];
                          rows = Transpose[{headers,rows}];
                          rows = Association[(#[[1]]->#[[2]])&/@rows];
                          term->rows
                          ),
                      {term,carnallTerms}
                      ];
        carnallAssoc        = Association[parsedData];
        carnallCalcEnergies = #[[calcIndex]]&/@carnallData;
        carnallCalcEnergies = If[NumberQ[#],#,Missing[]]&/@carnallCalcEnergies;
        ion                 = ions[[sheetIdx-3]];
        carnallCalcEnergies = PadRight[carnallCalcEnergies, fullSizes[ion], Missing[]];
        keys                = #[<|"ion"->ion|>]&/@templates;
        Carnall[keys[[1]]]  = carnallAssoc;
        Carnall[keys[[2]]]  = carnallCalcEnergies;
        Carnall[keys[[3]]]  = carnallData;
        Carnall[keys[[4]]]  = headers;
        ),
    {sheetIdx,4,16}
    ];
    
    goodions = Select[ions,#!="Pm"&];
    expData  = Select[Transpose[Carnall["appendix:"<>#<>":RawTable"]][[1+Position[Carnall["appendix:"<>#<>":Headings"],"Exp (1/cm)"][[1,1]]]],NumberQ]&/@goodions;
    Carnall["All Experimental Data"]=AssociationThread[goodions,expData];
    If[OptionValue["Export"],
      (
        carnallFname = FileNameJoin[{moduleDir, "data", "Carnall.m"}];
        Echo["Exporting to " <> carnallFname];
        Export[carnallFname, Carnall];
      )
      ];
    Return[Carnall];
  );

  CFP::usage = "CFP[{n, NKSL}] provides a list whose first element echoes NKSL and whose other elements are lists with two elements the first one being the symbol of a parent term and the second being the corresponding coefficient of fractional parentage. n must satisfy 1 <= n <= 7.
  These are according to the tables from Nielson & Koster.";

  CFPAssoc::usage = "CFPAssoc is an association where keys are of lists of the form {num_electrons, daughterTerm, parentTerm} and values are the corresponding coefficients of fractional parentage. The terms given in string-spectroscopic notation. If a certain daughter term does not have a parent term, the value is 0. Loaded using LoadCFP[].
  These are according to the tables from Nielson & Koster.";

  LoadCFP::usage = "LoadCFP[] loads CFP, CFPAssoc, and CFPTable into the session.";
  LoadCFP[] := (
    If[And[ValueQ[CFP], ValueQ[CFPTable], ValueQ[CFPAssoc]],Return[]];
    
    PrintTemporary["Loading CFPtable ..."];
    CFPTablefname = FileNameJoin[{moduleDir, "data", "CFPTable.m"}];
    If[!FileExistsQ[CFPTablefname],
      (PrintTemporary[">> CFPTable.m not found, generating ..."];
        CFPTable = GenerateCFPTable["Export"->True];
      ),
      CFPTable = Import[CFPTablefname];
    ];
    
    PrintTemporary["Loading CFPs.m ..."];
    CFPfname = FileNameJoin[{moduleDir, "data", "CFPs.m"}];
    If[!FileExistsQ[CFPfname],
      (PrintTemporary[">> CFPs.m not found, generating ..."];
        CFP = GenerateCFP["Export"->True];
      ),
      CFP = Import[CFPfname];
    ];
    
    PrintTemporary["Loading CFPAssoc.m ..."];
    CFPAfname = FileNameJoin[{moduleDir, "data", "CFPAssoc.m"}];
    If[!FileExistsQ[CFPAfname],
      (PrintTemporary[">> CFPAssoc.m not found, generating ..."];
        CFPAssoc = GenerateCFPAssoc["Export"->True];
      ),
      CFPAssoc = Import[CFPAfname];
    ];
  );

  ReducedUkTable::usage = "ReducedUkTable[{n, l = 3, SL, SpLp, k}] provides reduced matrix elements of the unit spherical tensor operator Uk. See TASS section 11-9 \"Unit Tensor Operators\". Loaded using LoadUk[].";

  LoadUk::usage = "LoadUk[] loads into session the reduced matrix elements for unit tensor operators.";
  LoadUk[] := (
    If[ValueQ[ReducedUkTable], Return[]];
    PrintTemporary["Loading the association of reduced matrix elements for unit tensor operators ..."];
    ReducedUkTableFname = FileNameJoin[{moduleDir, "data", "ReducedUkTable.m"}];
    If[!FileExistsQ[ReducedUkTableFname],
      (PrintTemporary[">> ReducedUkTable.m not found, generating ..."];
        ReducedUkTable = GenerateReducedUkTable[7];
      ),
      ReducedUkTable = Import[ReducedUkTableFname];
    ];
  );

  ElectrostaticTable::usage = "ElectrostaticTable[{n, SL, SpLp}] provides the calculated result of Electrostatic[{n, SL, SpLp}]. Load using LoadElectrostatic[].";

  LoadElectrostatic::usage = "LoadElectrostatic[] loads the reduced matrix elements for the electrostatic interaction.";
  LoadElectrostatic[] := (
    If[ValueQ[ElectrostaticTable], Return[]];
    PrintTemporary["Loading the association of matrix elements for the electrostatic interaction ..."];
    ElectrostaticTablefname = FileNameJoin[{moduleDir, "data", "ElectrostaticTable.m"}];
    If[!FileExistsQ[ElectrostaticTablefname],
      (PrintTemporary[">> ElectrostaticTable.m not found, generating ..."];
        ElectrostaticTable = GenerateElectrostaticTable[7];
      ),
      ElectrostaticTable = Import[ElectrostaticTablefname];
    ];
  );

  LoadV1k::usage = "LoadV1k[] loads into session the matrix elements of V1k.";
  LoadV1k[] := (
    If[ValueQ[ReducedV1kTable], Return[]];
    PrintTemporary["Loading the association of matrix elements for V1k ..."];
    ReducedV1kTableFname = FileNameJoin[{moduleDir, "data", "ReducedV1kTable.m"}];
    If[!FileExistsQ[ReducedV1kTableFname],
      (PrintTemporary[">> ReducedV1kTable.m not found, generating ..."];
        ReducedV1kTable = GenerateReducedV1kTable[7];
      ),
      ReducedV1kTable = Import[ReducedV1kTableFname];
    ]
  );

  LoadSpinOrbit::usage = "LoadSpinOrbit[] loads into session the matrix elements of the spin-orbit interaction.";
  LoadSpinOrbit[] := (
    If[ValueQ[SpinOrbitTable], Return[]];
    PrintTemporary["Loading the association of matrix elements for spin-orbit ..."];
    SpinOrbitTableFname = FileNameJoin[{moduleDir, "data", "SpinOrbitTable.m"}];
    If[!FileExistsQ[SpinOrbitTableFname],
      (
        PrintTemporary[">> SpinOrbitTable.m not found, generating ..."];
        SpinOrbitTable = GenerateSpinOrbitTable[7, "Export" -> True];
      ),
      SpinOrbitTable = Import[SpinOrbitTableFname];
    ]
  );

  LoadSOOandECSOLS::usage = "LoadSOOandECSOLS[] loads into session the LS reduced matrix elements of the SOO-ECSO interaction.";
  LoadSOOandECSOLS[] := (
    If[ValueQ[SOOandECSOLSTable], Return[]];
    PrintTemporary["Loading the association of LS reduced matrix elements for SOO-ECSO ..."];
    SOOandECSOLSTableFname = FileNameJoin[{moduleDir, "data", "ReducedSOOandECSOLSTable.m"}];
    If[!FileExistsQ[SOOandECSOLSTableFname],
      (PrintTemporary[">> ReducedSOOandECSOLSTable.m not found, generating ..."];
        SOOandECSOLSTable = GenerateSOOandECSOLSTable[7];
      ),
      SOOandECSOLSTable = Import[SOOandECSOLSTableFname];
    ];
  );

  LoadSOOandECSO::usage = "LoadSOOandECSO[] loads into session the LSJ reduced matrix elements of spin-other-orbit and electrostatically-correlated-spin-orbit.";
  LoadSOOandECSO[] := (
    If[ValueQ[SOOandECSOTableFname], Return[]];
    PrintTemporary["Loading the association of matrix elements for spin-other-orbit and electrostatically-correlated-spin-orbit ..."];
    SOOandECSOTableFname = FileNameJoin[{moduleDir, "data", "SOOandECSOTable.m"}];
    If[!FileExistsQ[SOOandECSOTableFname],
      (PrintTemporary[">> SOOandECSOTable.m not found, generating ..."];
        SOOandECSOTable = GenerateSOOandECSOTable[7, "Export"->True];
      ),
      SOOandECSOTable = Import[SOOandECSOTableFname];
    ];
  );

  LoadT22::usage = "LoadT22[] loads into session the matrix elements of the double tensor operator T22.";
  LoadT22[] := (
    If[ValueQ[T22Table], Return[]];
    PrintTemporary["Loading the association of reduced T22 matrix elements ..."];
    T22TableFname = FileNameJoin[{moduleDir, "data", "ReducedT22Table.m"}];
    If[!FileExistsQ[T22TableFname],
      (PrintTemporary[">> ReducedT22Table.m not found, generating ..."];
        T22Table = GenerateT22Table[7];
      ),
      T22Table = Import[T22TableFname];
    ];
  );

  LoadSpinSpin::usage = "LoadSpinSpin[] loads into session the matrix elements of spin-spin.";
  LoadSpinSpin[] := (
    If[ValueQ[SpinSpinTable], Return[]];
    PrintTemporary["Loading the association of matrix elements for spin-spin ..."];
    SpinSpinTableFname = FileNameJoin[{moduleDir, "data", "SpinSpinTable.m"}];
    If[!FileExistsQ[SpinSpinTableFname],
      (PrintTemporary[">> SpinSpinTable.m not found, generating ..."];
        SpinSpinTable = GenerateSpinSpinTable[7, "Export" -> True];
      ),
      SpinSpinTable = Import[SpinSpinTableFname];
    ];
  );

  LoadThreeBody::usage = "LoadThreeBody[] loads into session the matrix elements of three-body configuration-interaction effects.";
  LoadThreeBody[] := (
    If[ValueQ[ThreeBodyTable], Return[]];
    PrintTemporary["Loading the association of matrix elements for three-body configuration-interaction effects ..."];
    ThreeBodyFname   = FileNameJoin[{moduleDir, "data", "ThreeBodyTable.m"}];
    ThreeBodiesFname = FileNameJoin[{moduleDir, "data", "ThreeBodyTables.m"}];
    If[!FileExistsQ[ThreeBodyFname],
      (PrintTemporary[">> ThreeBodyTable.m not found, generating ..."];
        {ThreeBodyTable, ThreeBodyTables} = GenerateThreeBodyTables["Export" -> True];
      ),
      ThreeBodyTable  = Import[ThreeBodyFname];
      ThreeBodyTables = Import[ThreeBodiesFname];
    ];
  );

  (* ##################### Load Functions ###################### *)
  (* ########################################################### *)

End[];

LoadTermLabels[];
LoadCFP[];

EndPackage[];
