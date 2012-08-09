(* ::Package:: *)

(* ::Title:: *)
(*Rayleigh channel functions*)


(* ::Subsection::Closed:: *)
(*Copyright notice*)


(* ::Text:: *)
(*Mathematica function definitions for cooperative energy detection in Rayleigh channels.*)
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
(*02/08/2012*)
(*1.11*)


(* ::Subsection::Closed:: *)
(*Changelog*)


(* ::Text:: *)
(*Version 1.11: Moved database logging functions to the Network package.*)
(*Version 1.1: Introduced RelevantOptions function and changed function definitions, so that child options are inherited from parents.*)
(*Version 1.0: First working version, minor bug fixes to follow.*)


(* ::Section:: *)
(*Public*)


BeginPackage["Rayleigh`"];


(* ::Subsection::Closed:: *)
(*PDF of the signal to noise ratio*)


RayleighPDF;


(* ::Subsection::Closed:: *)
(*Detection probability*)


NRayleighProbabilityOfDetection;


(* ::Subsection::Closed:: *)
(*Sample complexity*)


NRayleighSampleComplexity;


(* ::Section:: *)
(*Private*)


Begin["`Private`"];


Needs["Nakagami`"];


(* ::Subsection::Closed:: *)
(*PDF of the signal to noise ratio*)


Options[RayleighPDF] = Options[NakagamiPDF];
RayleighPDF::usage="RayleighPDF[\!\(\*OverscriptBox[\(\[Gamma]\), \(_\)]\), x] evaluates the probability density function of the instantaneous signal to noise ratio at a single energy detector operating on a Rayleigh fading channel at x.\nRayleighPDF[\!\(\*OverscriptBox[\(\[Gamma]\), \(_\)]\), x, n] evaluates the probability density function of the average instantaneous signal to noise ratio at the fusion center of a cooperative network operating on a Rayleigh fading channel at x.\n\nThe following methods can be given:\n\nMethod\[Rule]\"Approximate\"\nMethod\[Rule]\"Exact\"\n\nBy default, Method\[Rule]\"Exact\".";
RayleighPDF[\[Gamma]_,x_,OptionsPattern[]]:=Module[{n = 1},RayleighPDF[\[Gamma],x,n,Method->OptionValue[Method]]]
RayleighPDF[\[Gamma]_,x_,n_,OptionsPattern[]]:=Module[{m = 1},NakagamiPDF[\[Gamma],m,x,n,Method->OptionValue[Method]]]


(* ::Subsection::Closed:: *)
(*Detection probability*)


Options[NRayleighProbabilityOfDetection]=Options[NNakagamiProbabilityOfDetection];
NRayleighProbabilityOfDetection::usage="NRayleighProbabilityOfDetection[M, \[Gamma], \[Lambda], m] calculates the probability of detection for a single energy detector operating on a Rayleigh fading channel.
NRayleighProbabilityOfDetection[M, \[Gamma], \[Lambda], m, n] calculates the probability of detection for the fusion center of a cooperative network operating on a Rayleigh fading channel.

The following methods can be given:

Method\[Rule]\"Approximate\"
Method\[Rule]{\"Approximate\", Algorithm\[Rule]...}
Method\[Rule]\"Exact\"
Method\[Rule]{\"Exact\", Algorithm\[Rule]...}

By default, Method\[Rule]\""<>ToString[Method/.Options[NRayleighProbabilityOfDetection]]<>"\".

For a given method, an algorithm may be specified. If Method\[Rule]\"Approximate\", then the following algorithms may be specified:

Algorithm\[Rule]\"Gaussian\"
Algorithm\[Rule]\"NGaussian\"

By default, Algorithm\[Rule]\""<>ToString[Algorithm/.Options[NRayleighProbabilityOfDetection]]<>"\". If Algorithm\[Rule]\"Gaussian\", then the switching point between the small and large mn approximations may be specified so that Method\[Rule]{\"Approximate\", Algorithm\[Rule]\"Gaussian\", LargeMN\[Rule]"<>ToString[LargeMN/.Options[NRayleighProbabilityOfDetection]]<>"}. By default, LargeMN\[Rule]"<>ToString[LargeMN/.Options[NRayleighProbabilityOfDetection]]<>". Similarly, if Algorithm\[Rule]\"NGaussian\", then a LowSNR boolean option may also be specified so that Method\[Rule]{\"Approximate\", Algorithm\[Rule]\"NGaussian\", LowSNR->"<>ToString[LowSNR/.Options[NRayleighProbabilityOfDetection]]<>"}. By default, LowSNR\[Rule]"<>ToString[LowSNR/.Options[NRayleighProbabilityOfDetection]]<>".

For the exact method, the following algorithms may be specified:

Algorithm\[Rule]\"Annamalai\"
Algorithm\[Rule]\"Digham\"
Algorithm\[Rule]\"Herath\"

By default, Algorithm\[Rule]\"Annamalai\".

In addition, timing options may be specified. The timing option is specified by:

Timed\[Rule]"<>ToString[Timed/.Options[ProbabilityOfDetection]]<>"
MaxIterations\[Rule]"<>ToString[MaxIterations/.Options[ProbabilityOfDetection]]<>"
MaxTime\[Rule]"<>ToString[MaxTime/.Options[ProbabilityOfDetection]]<>"

where the above options are the defaults, if not specified. If Timed\[Rule]True, then a {Pd, time} list of values will be returned.";
NRayleighProbabilityOfDetection[M_?NumericQ,\[Gamma]_?NumericQ,\[Lambda]_,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[NRayleighProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	NRayleighProbabilityOfDetection[M,\[Gamma],\[Lambda],n,RelevantOptions[NRayleighProbabilityOfDetection]]
]
NRayleighProbabilityOfDetection[M_?NumericQ,\[Gamma]_?NumericQ,\[Lambda]_,n_?IntegerQ,OptionsPattern[]]:=Module[{RelevantOptions, m = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[NRayleighProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	NNakagamiProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,RelevantOptions[NNakagamiProbabilityOfDetection]]
]


(* ::Subsection::Closed:: *)
(*Sample complexity*)


Options[NRayleighSampleComplexity]=Options[NNakagamiSampleComplexity];
NRayleighSampleComplexity::usage="NRayleighSampleComplexity[\[Gamma], Pf, Pd, m] calculates the sample complexity for a single energy detector operating on a Rayleigh fading channel.
NRayleighSampleComplexity[\[Gamma], Pf, Pd, m, n] calculates the sample complexity for a cooperative network operating on a Rayleigh fading channel.

The following methods can be given:

Method\[Rule]\"Approximate\"
Method\[Rule]\"Exact\"

By default, Method\[Rule]\""<>ToString[Method/.Options[NRayleighSampleComplexity]]<>"\".

If Method\[Rule]\"Approximate\", the LowSNR option may be specified. By default, LowSNR\[Rule]"<>ToString[LowSNR/.Options[NRayleighSampleComplexity]]<>".

Numerical tolerance can be specified using the Tolerance option. By default, Tolerance\[Rule]"<>ToString[Tolerance/.Options[NNakagamiSampleComplexity]//N//InputForm]<>".";
NRayleighSampleComplexity[\[Gamma]_?NumericQ,Pf_?NumericQ,Pd_?NumericQ,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[NRayleighSampleComplexity][[All,1]]],Options[target][[All,1]]];
	NRayleighSampleComplexity[\[Gamma],Pf,Pd,n,RelevantOptions[NRayleighSampleComplexity]]
]
NRayleighSampleComplexity[\[Gamma]_?NumericQ,Pf_?NumericQ,Pd_?NumericQ,n_?IntegerQ,OptionsPattern[]]:=Module[{RelevantOptions, m = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[NRayleighSampleComplexity][[All,1]]],Options[target][[All,1]]];
	NNakagamiSampleComplexity[\[Gamma],Pf,Pd,m,n,RelevantOptions[NNakagamiSampleComplexity]]
]


End[];


EndPackage[];
