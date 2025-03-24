(* ::Package:: *)

(* ::Section:: *)
(*QED model for FeynArts*)


Quit[]


(* ::Subsection:: *)
(*Load FeynRules*)


FR$Parallel=False;
$FeynRulesPath="/Users/vjhirsch/MG5/MG5_aMC_v3_5_7/ComputationalHEP/mathematica/FeynRules";
(*$FeynRulesPath=FileNameJoin[{$UserBaseDirectory,"Applications","FeynRules"}];*)
Import[FileNameJoin[$FeynRulesPath,"FeynRules.m"]];


(* ::Subsection:: *)
(*Load FeynRules model*)


If[$FrontEnd===Null,
nbDir=DirectoryName[$InputFileName],
nbDir=NotebookDirectory[]
];


frModelPath=FileNameJoin[{"/Users/vjhirsch/MG5/MG5_aMC_v3_5_7/ComputationalHEP/mathematica/FeynCalc/Examples/FeynRules/QED/Mathematica/","QED.fr"}];
LoadModel[frModelPath];


(* ::Subsection:: *)
(*Create FeynArts model*)


FR$Loop=True;
SetDirectory[FileNameJoin[{"/Users/vjhirsch/MG5/MG5_aMC_v3_5_7/ComputationalHEP/mathematica/FeynCalc","Models"}]];
WriteFeynArtsOutput[LQED,Output->"QED",CouplingRename->False];



