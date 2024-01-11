#!/usr/bin/env wolframscript

(*This is simple command line interface to qlanth, it takes the name of a file with the input parameters to the model Hamiltonian. This file should be configured as an association using the adequate symbols for the model parameters. Every parameter that is not given will be coerced to be zero by ParamPad.

The number of f-electrons is given by the parameter "numE", and is a mandatory parameter.

The Association also needs to have a key called "saveRoot" which will be used to save the results of the calculation. The results will be saved in two files, one with the eigenenergies and one with the eigenvectors. The names of these files will be the value of "saveRoot" with the suffix "_energies.csv" and "_states.csv" respectively. Each of the lines in the energies file will be a single eigenenergy, and each of the lines in the states file will be a single eigenvector, with the elements of the vector separated by commas.

No specific ratios are assumed for M2/M0 or M4/M0, their values must be given explicitly, or they will be coerced to zero. Similarly for P4/P2, P6/P2.*)

Get["qlanth.m"];
Get["misc.m"];

argv = Rest @ $ScriptCommandLine;
argc = Length @ argv;
configFileFname = argv[[1]];

params = Import[configFileFname];

numE = params["numE"];
saveRoot = params["saveRoot"];

LoadSymbolicHamiltonians[{If[numE>7, 14-numE, numE]}];
params = ParamPad[params, "Print"->False];
If[numE > 7,
    params = HoleElectronConjugation[params]
];

symHam = symbolicHamiltonians[If[numE>7, 14-numE, numE]];
numHam = ReplaceInSparseArray[symHam, Normal[params]];

eigenSys = Eigensystem[numHam];
eigenEnergies = eigenSys[[1]];
eigenOrder = Ordering[eigenEnergies];
eigenEnergies = eigenEnergies[[eigenOrder]];
eigenVectors = eigenSys[[2]];
eigenVectors = eigenVectors[[eigenOrder]];

energiesFname = saveRoot<>"_energies.csv";
eigenvecsFname = saveRoot<>"_states.csv";
Export[energiesFname, eigenEnergies];
Export[eigenvecsFname, eigenVectors];
