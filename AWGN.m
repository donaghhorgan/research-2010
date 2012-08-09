(* ::Package:: *)

(* ::Title:: *)
(*AWGN channel functions*)


(* ::Subsection::Closed:: *)
(*Copyright notice*)


(* ::Text:: *)
(*Mathematica function definitions for cooperative energy detection in AWGN channels.*)
(*Copyright (C) 2012 Donagh Horgan.*)
(*Email: donaghh@rennes.ucc.ie.*)
(**)
(*This program is free software : you can redistribute it and/or modify*)
(*it under the terms of the GNU General Public License as published by*)
(*the Free Software Foundation, either version 3 of the License, or*)
(*(at your option) any later version.*)
(**)
(*This program is distributed in the hope that it will be useful,*)
(*but WITHOUT ANY WARRANTY; without even the implied warranty of*)
(*MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See *)
(*COPYING for more details.*)
(**)
(*You should have received a copy of the GNU General Public License*)
(*along with this program. If not, see http://www.gnu.org/licenses.*)


(* ::Subsection::Closed:: *)
(*Version information*)


(* ::Text:: *)
(*27/06/2012*)
(*1.0*)


(* ::Subsection::Closed:: *)
(*Changelog*)


(* ::Text:: *)
(*Version 1.0: First working version, minor bug fixes to follow.*)


(* ::Section:: *)
(*Public*)


BeginPackage["AWGN`"];


(* ::Subsection::Closed:: *)
(*PDF of the signal to noise ratio*)


AWGNPDF;


(* ::Subsection:: *)
(*Decision probabilities*)


(* ::Subsubsection::Closed:: *)
(*Probability of false alarm*)


AWGNProbabilityOfFalseAlarm;


(* ::Subsubsection::Closed:: *)
(*Probability of detection*)


AWGNProbabilityOfDetection;


(* ::Subsubsection::Closed:: *)
(*Threshold*)


\[Lambda];


(* ::Subsection::Closed:: *)
(*Sample complexity*)


AWGNSampleComplexity;


(* ::Section:: *)
(*Private*)


Begin["`Private`"];


Needs["QFunction`"];


(* ::Subsection:: *)
(*PDF of the recieved energy*)


Options[AWGNPDF] = {Method->"Exact"};
AWGNPDF::usage="AWGNPDF[M, \!\(\*OverscriptBox[\(\[Gamma]\), \(_\)]\), x] evaluates the probability density function of the recieved energy at a single energy detector operating on an AWGN channel at x.
AWGNPDF[M, \!\(\*OverscriptBox[\(\[Gamma]\), \(_\)]\), x, n] evaluates the probability density function of the recieved energy at the fusion center of a cooperative network operating on an AWGN fading channel at x.

The following methods can be given:

Method\[Rule]\"Approximate\"
Method\[Rule]\"Exact\"

By default, Method\[Rule]\""<>ToString[Method/.Options[AWGNPDF]]<>"\".";
AWGNPDF[M_,\[Gamma]_,x_,OptionsPattern[]]:=Module[{n = 1},AWGNPDF[M,\[Gamma],x,n,Method->OptionValue[Method]]]
AWGNPDF[M_,\[Gamma]_,x_,n_,OptionsPattern[]]:=Switch[OptionValue[Method],
	"Exact",
	PDF[NoncentralChiSquareDistribution[M n, M n \[Gamma]], x],
	"Approximate",
	PDF[NormalDistribution[M n(1+\[Gamma]),Sqrt[2M n(1+2\[Gamma])]],x],
	_,
	AWGNPDF[M,\[Gamma],x,n,Method->"Exact"]
]


(* ::Subsection:: *)
(*Decision probabilities*)


(* ::Subsubsection::Closed:: *)
(*Probability of false alarm*)


Options[AWGNProbabilityOfFalseAlarm]={Method->"Approximate"};
AWGNProbabilityOfFalseAlarm::usage="AWGNProbabilityOfFalseAlarm[M, \[Lambda]] calculates the probability of false alarm for a single energy detector operating on an AWGN channel.
AWGNProbabilityOfFalseAlarm[M, \[Lambda], n] calculates the probability of false alarm for a cooperative network operating on an AWGN channel.

The following methods may be specified:

Method\[Rule]\"Approximate\"
Method\[Rule]\"Exact\"

By default, Method\[Rule]\""<>ToString[Method/.Options[AWGNProbabilityOfFalseAlarm]]<>"\".";
AWGNProbabilityOfFalseAlarm[M_,\[Lambda]_,OptionsPattern[]]:=Module[{n = 1},AWGNProbabilityOfFalseAlarm[M,\[Lambda],n,Method->OptionValue[Method]]]
AWGNProbabilityOfFalseAlarm[M_,\[Lambda]_,n_,OptionsPattern[]]:=Switch[OptionValue[Method],
	"Approximate",
	Q[(\[Lambda]-M n) / Sqrt[2M n]],
	"Exact",
	GammaRegularized[M n/2, \[Lambda] / 2],
	_,
	AWGNProbabilityOfFalseAlarm[M,\[Lambda],n]
]


(* ::Subsubsection:: *)
(*Probability of detection*)


Options[AWGNProbabilityOfDetection]={Method->"Approximate",LowSNR->True}
AWGNProbabilityOfDetection::usage="AWGNProbabilityOfDetection[M, \[Gamma], \[Lambda]] calculates the approximate probability of detection for a single energy detector operating on an AWGN channel.
AWGNProbabilityOfDetection[M, \[Gamma], \[Lambda], n] calculates the approximate probability of detection for a cooperative network operating on an AWGN channel.

The following methods may be specified:

Method\[Rule]\"Approximate\"
Method\[Rule]\"Exact\"

By default, Method\[Rule]\""<>ToString[Method/.Options[AWGNProbabilityOfDetection]]<>"\".

In addition, the LowSNR option can be used to specify whether to use a low signal to noise ratio approximation. By default, LowSNR\[Rule]"<>ToString[LowSNR/.Options[AWGNProbabilityOfDetection]]<>".";
AWGNProbabilityOfDetection[M_,\[Gamma]_,\[Lambda]_,OptionsPattern[]]:=Module[{n = 1},AWGNProbabilityOfDetection[M,\[Gamma],\[Lambda],n,Method->OptionValue[Method],LowSNR->OptionValue[LowSNR]]]
AWGNProbabilityOfDetection[M_,\[Gamma]_,\[Lambda]_,n_,OptionsPattern[]]:=Switch[OptionValue[Method],
	"Approximate",
	If[OptionValue[LowSNR],
		Q[(\[Lambda]-M n(1+\[Gamma])) / Sqrt[2M n]],
		Q[(\[Lambda]-M n(1+\[Gamma])) / Sqrt[2M n(1+2\[Gamma])]]
	],
	"Exact",
	MarcumQ[M n / 2, Sqrt[M n \[Gamma]], Sqrt[\[Lambda]]],
	_,
	AWGNProbabilityOfDetection[M,\[Gamma],\[Lambda],n]
]


(* ::Subsubsection::Closed:: *)
(*Threshold*)


Options[\[Lambda]]={Method->"Approximate"};
\[Lambda]::usage="\[Lambda][M, Pf] calculates a threshold suitable for use in the calculation of the decision probabilities for a single energy detector.
\[Lambda][M, Pf, n] calculates a threshold suitable for use in the calculation of the fusion center decision probabilities when Nb = \[Infinity].

The following methods may be specified:

Method\[Rule]\"Approximate\"
Method\[Rule]\"Exact\"

By default, Method\[Rule]\""<>ToString[Method/.Options[\[Lambda]]]<>"\".";
\[Lambda][M_,Pf_,OptionsPattern[]]:=Module[{n = 1},\[Lambda][M,Pf,n,Method->OptionValue[Method]]]
\[Lambda][M_,Pf_,n_,OptionsPattern[]]:=Switch[OptionValue[Method],
	"Approximate",
	Sqrt[2M n] InverseQ[Pf] + M n,
	"Exact",
	2 InverseGammaRegularized[M n / 2, Pf],
	_,
	\[Lambda][M,Pf,n]
]


(* ::Subsection::Closed:: *)
(*Sample complexity*)


AWGNSampleComplexity::usage="AWGNSampleComplexity[\[Gamma], Pf, Pd] calculates the approximate number of samples required for a single energy detector to operate with the specified decision probabilities at a given signal to noise ratio in an AWGN channel.
AWGNSampleComplexity[\[Gamma], Pf, Pd, n] calculates the approximate number of samples required for a cooperative network to operate with the specified decision probabilities at a given signal to noise ratio in an AWGN channel.";
AWGNSampleComplexity[\[Gamma]_,Pf_,Pd_,n_:1]:= (2 / n) * ((InverseQ[Pf] - InverseQ[Pd]) / \[Gamma])^2


End[]


EndPackage[];
