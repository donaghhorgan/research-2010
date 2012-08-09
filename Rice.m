(* ::Package:: *)

(* ::Title:: *)
(*Rice channel functions*)


(* ::Subsection::Closed:: *)
(*Copyright notice*)


(* ::Text:: *)
(*Mathematica function definitions for cooperative energy detection in Rice channels.*)
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
(*1.02*)


(* ::Subsection::Closed:: *)
(*Changelog*)


(* ::Text:: *)
(*Version 1.02: Moved database logging functions to the Network package.*)
(*Version 1.01: Added sample complexity function.*)
(*Version 1.0: First working version, minor bug fixes to follow.*)


(* ::Section:: *)
(*Public*)


BeginPackage["Rice`"];


(* ::Subsection::Closed:: *)
(*PDF of the signal to noise ratio*)


RicePDF;


(* ::Subsection:: *)
(*Probabiity of detection*)


(* ::Subsubsection::Closed:: *)
(*Main function*)


NRiceProbabilityOfDetection;


(* ::Subsubsection::Closed:: *)
(*Annamalai's method*)


NAnnamalaiRiceProbabilityOfDetection;


(* ::Subsubsection::Closed:: *)
(*Horgan's method*)


NHorganRiceProbabilityOfDetection;


(* ::Subsubsection::Closed:: *)
(*Numerical Gaussian method*)


NGaussianRiceProbabilityOfDetection;


(* ::Subsection::Closed:: *)
(*Sample complexity*)


NRiceSampleComplexity


(* ::Section:: *)
(*Private*)


Begin["`Private`"];


<<Network`;
<<AWGN`;
<<ErfApprox`;


(* ::Subsection::Closed:: *)
(*PDF of the signal to noise ratio*)


Options[RicePDF] = {Method->"Exact"};
RicePDF::usage="RicePDF[\!\(\*OverscriptBox[\(\[Gamma]\), \(_\)]\), m, x] evaluates the probability density function of the instantaneous signal to noise ratio at a single energy detector operating on a Rice fading channel at x.
RicePDF[\!\(\*OverscriptBox[\(\[Gamma]\), \(_\)]\), m, x, n] evaluates the probability density function of the average instantaneous signal to noise ratio at the fusion center of a cooperative network operating on a Rice fading channel at x.

The following methods can be given:

Method\[Rule]\"Approximate\"
Method\[Rule]\"Exact\"

By default, Method\[Rule]\""<>ToString[Method/.Options[RicePDF]]<>"\".";
RicePDF[\[Gamma]_,m_,x_,OptionsPattern[]]:=Module[{n = 1},RicePDF[\[Gamma],m,x,n,Method->OptionValue[Method]]]
RicePDF[\[Gamma]_,m_,x_,n_,OptionsPattern[]]:=Switch[OptionValue[Method],
	"Exact",
	(m + 1) / \[Gamma] Exp[-m n - (m + 1) x / \[Gamma]]((m + 1) x / (m n \[Gamma]))^((n - 1) / 2) BesselI[n - 1, 2 Sqrt[m n (m + 1)x / \[Gamma]]],
	"Approximate",
	2 (m + 1) / \[Gamma] PDF[NormalDistribution[2 n + 2 m n, Sqrt[4 n + 8 m n]], 2 (m + 1) x / \[Gamma]],
	_,
	RicePDF[\[Gamma],m,x,n,Method->"Exact"]
]


(* ::Subsection:: *)
(*Probabiity of detection*)


(* ::Subsubsection::Closed:: *)
(*Main function*)


Options[NRiceProbabilityOfDetection]={Method->OptionValue[ProbabilityOfDetection,Method],Algorithm->OptionValue[ProbabilityOfDetection,Algorithm],LargeMN->OptionValue[ProbabilityOfDetection,LargeMN],LowSNR->OptionValue[ProbabilityOfDetection,LowSNR],Timed->OptionValue[ProbabilityOfDetection,Timed],MaxTime->OptionValue[ProbabilityOfDetection,MaxTime],MaxIterations->OptionValue[ProbabilityOfDetection,MaxIterations]};
NRiceProbabilityOfDetection::usage="NRiceProbabilityOfDetection[M, \[Gamma], \[Lambda], K] calculates the probability of detection for a single energy detector operating on a Rice-K fading channel.
NRiceProbabilityOfDetection[M, \[Gamma], \[Lambda], K, n] calculates the probability of detection for the fusion center of a cooperative network operating on a Rice-K fading channel.

The following methods can be given:

Method\[Rule]\"Approximate\"
Method\[Rule]{\"Approximate\", Algorithm\[Rule]...}
Method\[Rule]\"Exact\"
Method\[Rule]{\"Exact\", Algorithm\[Rule]...}

By default, Method\[Rule]\""<>ToString[Method/.Options[NRiceProbabilityOfDetection]]<>"\".

For a given method, an algorithm may be specified. If Method\[Rule]\"Approximate\", then the following algorithms may be specified:

Algorithm\[Rule]\"Horgan\"
Algorithm\[Rule]\"NGaussian\"

By default, Algorithm\[Rule]\""<>ToString[Algorithm/.Options[NRiceProbabilityOfDetection]]<>"\". If Algorithm\[Rule]\"Horgan\", then the switching point between the small and large mn approximations may be specified so that Method\[Rule]{\"Approximate\", Algorithm\[Rule]\"Horgan\", LargeMN\[Rule]"<>ToString[LargeMN/.Options[NRiceProbabilityOfDetection]]<>"}. By default, LargeMN\[Rule]"<>ToString[LargeMN/.Options[NRiceProbabilityOfDetection]]<>". Similarly, if Algorithm\[Rule]\"NGaussian\", then a LowSNR boolean option may also be specified so that Method\[Rule]{\"Approximate\", Algorithm\[Rule]\"NGaussian\", LowSNR->"<>ToString[LowSNR/.Options[NRiceProbabilityOfDetection]]<>"}. By default, LowSNR\[Rule]"<>ToString[LowSNR/.Options[NRiceProbabilityOfDetection]]<>".

For the exact method, the following algorithms may be specified:

Algorithm\[Rule]\"Annamalai\"
Algorithm\[Rule]\"Herath\"

By default, Algorithm\[Rule]\"Annamalai\".

In addition, timing options may be specified. The timing option is specified by:

Timed\[Rule]"<>ToString[Timed/.Options[ProbabilityOfDetection]]<>"
MaxIterations\[Rule]"<>ToString[MaxIterations/.Options[ProbabilityOfDetection]]<>"
MaxTime\[Rule]"<>ToString[MaxTime/.Options[ProbabilityOfDetection]]<>"

where the above options are the defaults, if not specified. If Timed\[Rule]True, then a {Pd, time} list of values will be returned.";
NRiceProbabilityOfDetection[M_,\[Gamma]_,\[Lambda]_,m_,OptionsPattern[]]:=Module[{n = 1, RelevantOptions},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[NRiceProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	NRiceProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,RelevantOptions[NRiceProbabilityOfDetection]]
]
NRiceProbabilityOfDetection[M_,\[Gamma]_,\[Lambda]_,m_,n_,OptionsPattern[]]:=Module[{f, result, time = 0, totaltime = 0, iterations = 0, channelType = "Rice", rationalPf, RelevantOptions},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[NRiceProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	Switch[OptionValue[Method],
		"Exact",
		Switch[OptionValue[Algorithm],
			"Annamalai",
			NAnnamalaiRiceProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,RelevantOptions[NAnnamalaiRiceProbabilityOfDetection]],
			"Herath",
			NHerathRiceProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,RelevantOptions[NHerathRiceProbabilityOfDetection]],
			_,
			NRiceProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,Algorithm->"Annamalai",RelevantOptions[NRiceProbabilityOfDetection]]
		],
		"Approximate",
		Switch[OptionValue[Algorithm],
			"Horgan",
			NHorganRiceProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,RelevantOptions[NHorganRiceProbabilityOfDetection]],
			"NGaussian",
			NGaussianRiceProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,RelevantOptions[NGaussianRiceProbabilityOfDetection]],
			_,
			NRiceProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,Algorithm->"NGaussian",RelevantOptions[NRiceProbabilityOfDetection]]
		],
		_,
		Undefined
	]
]


(* ::Subsubsection::Closed:: *)
(*Annamalai's method*)


Options[AnnamalaiLimit]={Tolerance->10^-8};
AnnamalaiLimit::usage="AnnamalaiLimit[M, \[Gamma], \[Lambda], m] calculates the truncation point for use in Herath's algorithm using the default tolerance for a single energy detector.
AnnamalaiLimit[M, \[Gamma], \[Lambda], m, n] calculates the truncation point for use in Herath's algorithm using the default tolerance for a cooperative network.

The calculation tolerance may be specified using the Tolerance option. By default, Tolerance\[Rule]"<>ToString[Tolerance/.Options[AnnamalaiLimit]//N//InputForm]<>".";
AnnamalaiLimit[M_?NumericQ,\[Gamma]_?NumericQ,\[Lambda]_,m_?NumericQ,OptionsPattern[]]:=Module[{n = 1},AnnamalaiLimit[M,\[Gamma],\[Lambda],m,n,Tolerance->OptionValue[Tolerance]]]
AnnamalaiLimit[M_?NumericQ,\[Gamma]_?NumericQ,\[Lambda]_,m_?NumericQ,n_?IntegerQ,OptionsPattern[]]:=Module[{tol=OptionValue[Tolerance]},
	1/2 (2(\[Lambda] / Sqrt[2M n] - InverseCDF[NormalDistribution[],tol])^2 - M n)//N//Ceiling
]


Options[NAnnamalaiRiceProbabilityOfDetection]={Timed->OptionValue[ProbabilityOfDetection,Timed],MaxTime->OptionValue[ProbabilityOfDetection,MaxTime],MaxIterations->OptionValue[ProbabilityOfDetection,MaxIterations]};
NAnnamalaiRiceProbabilityOfDetection::usage="NAnnamalaiRiceProbabilityOfDetection[M, \[Gamma], \[Lambda], m, lim] calculates the exact probability of detection for a single energy detector operating in a Rice-K fading channel using Annamalai's algorithm.
NAnnamalaiRiceProbabilityOfDetection[M, \[Gamma], \[Lambda], m, n, lim] calculates the exact probability of detection for a cooperative network operating in a Rice-m fading channel using Annamalai's algorithm.

Function timing may be specified using the following options:

Timed\[Rule]"<>ToString[Timed/.Options[ProbabilityOfDetection]]<>"
MaxIterations\[Rule]"<>ToString[MaxIterations/.Options[ProbabilityOfDetection]]<>"
MaxTime\[Rule]"<>ToString[MaxTime/.Options[ProbabilityOfDetection]]<>"

where the above settings are the defaults, if not specified. If Timed\[Rule]True, then a {Pd, time} list of values will be returned.";
NAnnamalaiRiceProbabilityOfDetection[M_?NumericQ,\[Gamma]_?NumericQ,\[Lambda]_,m_?NumericQ,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[NAnnamalaiRiceProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	NAnnamalaiRiceProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,RelevantOptions[NAnnamalaiRiceProbabilityOfDetection]]
]
NAnnamalaiRiceProbabilityOfDetection[M_?NumericQ,\[Gamma]_?NumericQ,\[Lambda]_,m_?NumericQ,n_?IntegerQ,OptionsPattern[]]:=Module[{lim, f, totaltime = 0, iterations = 0, time, result},
	lim = AnnamalaiLimit[M,\[Gamma],\[Lambda],m,n];
	f:=Block[{$MaxExtraPrecision=\[Infinity]},
		N[
			1 - ((2 (1 + m) / (2 + 2m + M \[Gamma]))^n) Exp[-m n M \[Gamma] / (2 + 2m + M \[Gamma])] (1 - GammaRegularized[M n/2, \[Lambda] / 2]) - Total[Table[(n + k - 1)! ((2 (1 + m) / (2 + 2m + M \[Gamma]))^n) (((M \[Gamma]) / (2 + 2m + M \[Gamma]))^k) Exp[-m n M \[Gamma] / (2 + 2m + M \[Gamma])] (1 / ((n - 1)! k!) + Total[Table[((2 m n (1 + m) / (2 + 2m + M \[Gamma]))^i) / (i! (n + i - 1)! (k - i)!),{i , 1, k}]]) (1 - GammaRegularized[M n / 2 + k, \[Lambda] / 2]),{k, 1, lim}]],
		20]
	];
	If[OptionValue[Timed],
		(* Evaluate result until MaxTime seconds of CPU time have been used or MaxIterations have been performed, whichever comes first *)
		While[totaltime < OptionValue[MaxTime] && iterations < OptionValue[MaxIterations],
			ClearSystemCache[];
			{time, result} = TimeConstrained[Timing[f],OptionValue[MaxTime],{OptionValue[MaxTime],Null}];
			totaltime += time;
			iterations++;
		];
		{result,totaltime/iterations},
		f
	]
]


(* ::Subsubsection::Closed:: *)
(*Herath' s method*)


Options[HerathLimit]={Tolerance->10^-6};
HerathLimit::usage="HerathLimit[M, \[Gamma], \[Lambda], m] calculates the truncation point for use in Herath's algorithm using the default tolerance for a single energy detector.
HerathLimit[M, \[Gamma], \[Lambda], m, n] calculates the truncation point for use in Herath's algorithm using the default tolerance for a cooperative network.

The calculation tolerance may be specified using the Tolerance option. By default, Tolerance\[Rule]"<>ToString[Tolerance/.Options[HerathLimit]//N//InputForm]<>".";
HerathLimit[M_?NumericQ,\[Gamma]_?NumericQ,\[Lambda]_,m_?NumericQ,OptionsPattern[]]:=Module[{n = 1},HerathLimit[M,\[Gamma],\[Lambda],m,n,Tolerance->OptionValue[Tolerance]]]
HerathLimit[M_?NumericQ,\[Gamma]_?NumericQ,\[Lambda]_,m_?NumericQ,n_?IntegerQ,OptionsPattern[]]:=Module[{x0,x1,f,fx0,tol=OptionValue[Tolerance]},
	Block[{$MaxExtraPrecision=\[Infinity]},
		f[j_?NumericQ]:=N[(m / ((M / 2) \[Gamma] + m))^(n m) Hypergeometric1F1[m n, j + 1, \[Lambda] (M / 2) \[Gamma] / (2 ((M / 2) \[Gamma] + m))] (1 - GammaRegularized[j + 1, \[Lambda] / 2]),20];
		(* Use Newton's method to find root *)
		Ceiling[j/.FindRoot[f[j]==tol,{j,Ceiling[M n/2]}]]
	]
]


Options[NHerathRiceProbabilityOfDetection]={Timed->OptionValue[ProbabilityOfDetection,Timed],MaxTime->OptionValue[ProbabilityOfDetection,MaxTime],MaxIterations->OptionValue[ProbabilityOfDetection,MaxIterations]};
NHerathRiceProbabilityOfDetection::usage="NHerathRiceProbabilityOfDetection[M, \[Gamma], \[Lambda], m] calculates the exact probability of detection for a single energy detector operating in a Rice-K fading channel using Herath's algorithm.
NHerathRiceProbabilityOfDetection[M, \[Gamma], \[Lambda], m, n] calculates the exact probability of detection for a cooperative network operating in a Rice-K fading channel using Herath's algorithm.

Function timing may be specified using the following options:

Timed\[Rule]"<>ToString[Timed/.Options[ProbabilityOfDetection]]<>"
MaxIterations\[Rule]"<>ToString[MaxIterations/.Options[ProbabilityOfDetection]]<>"
MaxTime\[Rule]"<>ToString[MaxTime/.Options[ProbabilityOfDetection]]<>"

where the above settings are the defaults, if not specified. If Timed\[Rule]True, then a {Pd, time} list of values will be returned.";
NHerathRiceProbabilityOfDetection[M_?NumericQ,\[Gamma]_?NumericQ,\[Lambda]_,m_?NumericQ,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[NHerathRiceProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	NHerathRiceProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,RelevantOptions[NHerathRiceProbabilityOfDetection]]
]
NHerathRiceProbabilityOfDetection[M_?NumericQ,\[Gamma]_?NumericQ,\[Lambda]_,m_?NumericQ,n_?IntegerQ,OptionsPattern[]]:=Module[{lim, f, totaltime = 0, iterations = 0, time, result},
	lim = HerathLimit[M,\[Gamma],\[Lambda],m,n];
	f:=Block[{$MaxExtraPrecision=\[Infinity]},
		N[
			Exp[-\[Lambda] / 2](((1 + m)Exp[-m \[Gamma] / (\[Gamma] + m + 1)])/(\[Gamma] + 1 + m))^n Total[Table[(m (\[Gamma] / (\[Gamma] + m + 1)) (1 - (\[Gamma] / (\[Gamma] + m + 1))))^(k-1)/(k - 1)! (1/(u - n - 1) (D[Exp[\[Lambda] z / 2]/((1 - z)(z - (\[Gamma] / (\[Gamma] + m + 1)))^(n + k - 1)),{z, u - L - 1}]/.z->0)+1/(n + k - 2)! (D[Exp[\[Lambda] z / 2]/(z^(u-n) (1 - z)),{z, n + k - 2}]/.z->(\[Gamma] / (\[Gamma] + m + 1)))),{k,1,lim}]]
		20]
	];
	If[OptionValue[Timed],
		(* Evaluate result until MaxTime seconds of CPU time have been used or MaxIterations have been performed, whichever comes first *)
		While[totaltime < OptionValue[MaxTime] && iterations < OptionValue[MaxIterations],
			ClearSystemCache[];
			{time, result} = TimeConstrained[Timing[f],OptionValue[MaxTime],{OptionValue[MaxTime],Null}];
			totaltime += time;
			iterations++;
		];
		{result,totaltime/iterations},
		f
	]
];


(* ::Subsubsection::Closed:: *)
(*Horgan' s method*)


Options[NHorganRiceProbabilityOfDetection]={LargeMN->OptionValue[ProbabilityOfDetection,LargeMN],Timed->OptionValue[ProbabilityOfDetection,Timed],MaxTime->OptionValue[ProbabilityOfDetection,MaxTime],MaxIterations->OptionValue[ProbabilityOfDetection,MaxIterations]};
NHorganRiceProbabilityOfDetection::usage="NHorganRiceProbabilityOfDetection[M, \[Gamma], \[Lambda], m] calculates the approximate probability of detection for a single energy detector operating in a Rice-m fading channel using Horgan's algorithm.
NHorganRiceProbabilityOfDetection[M, \[Gamma], \[Lambda], m, n] calculates the approximate probability of detection for a cooperative network operating in a Rice-m fading channel using Horgan's algorithm.

The switching point between the large and small mn approximations can be specified using the LargeMN option. By default, LargeMN\[Rule]"<>ToString[LargeMN/.Options[NHorganRiceProbabilityOfDetection]]<>".

Function timing may be specified using the following options:

Timed\[Rule]"<>ToString[Timed/.Options[ProbabilityOfDetection]]<>"
MaxIterations\[Rule]"<>ToString[MaxIterations/.Options[ProbabilityOfDetection]]<>"
MaxTime\[Rule]"<>ToString[MaxTime/.Options[ProbabilityOfDetection]]<>"

where the above settings are the defaults, if not specified. If Timed\[Rule]True, then a {Pd, time} list of values will be returned.";
NHorganRiceProbabilityOfDetection[M_?NumericQ,\[Gamma]_?NumericQ,\[Lambda]_,m_?NumericQ,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[NHorganRiceProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	NHorganRiceProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,RelevantOptions[NHorganRiceProbabilityOfDetection]]
]
NHorganRiceProbabilityOfDetection[M_?NumericQ,\[Gamma]_?NumericQ,\[Lambda]_,m_?NumericQ,n_?NumericQ,OptionsPattern[]]:=Module[{lim = Ceiling[Max[m E, n E]], f, totaltime = 0, iterations = 0, time, result},
	f:=Block[{$MaxExtraPrecision=\[Infinity]},
		N[
			Module[{x = Round[m n], tol = 10^-6},
				(* This method can only be used when m * n is an integer *)
				If[Abs[m n - x] <= tol,
					If[x < OptionValue[LargeMN],
						(1 / 2)E^(-m n) Total[Table[((m n)^k / (k!Gamma[n+k])) (-1)^(n+k-1) (D[(1+Erf[(M n-\[Lambda])/(2 Sqrt[M n])]-E^(((1+m) t (n (t+m t+M \[Gamma])-\[Gamma] \[Lambda]))/(M \[Gamma]^2)) (-1+Erf[(n (2 (1+m) t+M \[Gamma])-\[Gamma] \[Lambda])/(2 Sqrt[M n] \[Gamma])]))/t,{t,n+k-1}]/.t->1),{k,0,lim}]],
						(1/2 (1+Erf[((1+m) (M n (1+\[Gamma])-\[Lambda]))/(Sqrt[2] M Sqrt[n+2 m n] \[Gamma])])+(E^(-((b^2 (1+2 m) M \[Gamma]^2+4 c ((1+m)^2-a (1+2 m) M \[Gamma]^2)-(2 Sqrt[2] b (1+m)^2 (M n (1+\[Gamma])-\[Lambda]))/Sqrt[M n]+(2 a (1+m)^2 (-M n (1+\[Gamma])+\[Lambda])^2)/(M n))/(4 (-(1+m)^2+a (1+2 m) M \[Gamma]^2)))) (1+m) (-1+Erf[(-((b M)/(Sqrt[2] Sqrt[M n]))+(1+m)^2/(\[Gamma]+2 m \[Gamma])+a (M-\[Lambda]/n))/(Sqrt[2] Sqrt[-((a M)/n)+(1+m)^2/((n+2 m n) \[Gamma]^2)])]+Erfc[(M (-Sqrt[2] b (1+2 m) Sqrt[M n] \[Gamma]^2+2 (1+m)^2 n (1+\[Gamma]))-2 (1+m)^2 \[Lambda])/(2 Sqrt[2] M \[Gamma] Sqrt[(1+2 m) n ((1+m)^2-a (1+2 m) M \[Gamma]^2)])]))/(2 Sqrt[(1+m)^2-a (1+2 m) M \[Gamma]^2])+(E^(-((b^2 (1+2 m) M \[Gamma]^2+4 c ((1+m)^2-a (1+2 m) M \[Gamma]^2)+(2 Sqrt[2] b (1+m)^2 (M n (1+\[Gamma])-\[Lambda]))/Sqrt[M n]+(2 a (1+m)^2 (-M n (1+\[Gamma])+\[Lambda])^2)/(M n))/(4 (-(1+m)^2+a (1+2 m) M \[Gamma]^2)))) (1+m) (-2+Erfc[(M (Sqrt[2] b (1+2 m) Sqrt[M n] \[Gamma]^2+2 (1+m)^2 n (1+\[Gamma]))-2 (1+m)^2 \[Lambda])/(2 Sqrt[2] M \[Gamma] Sqrt[(1+2 m) n ((1+m)^2-a (1+2 m) M \[Gamma]^2)])]))/(2 Sqrt[(1+m)^2-a (1+2 m) M \[Gamma]^2]))/.LopezBenitezParameters[(-M (n+\[Gamma])+\[Lambda])/(2 Sqrt[M n])]
					],
					Undefined
				]
			],
		20]
	];
	If[OptionValue[Timed],
		(* Evaluate result until MaxTime seconds of CPU time have been used or MaxIterations have been performed, whichever comes first *)
		While[totaltime < OptionValue[MaxTime] && iterations < OptionValue[MaxIterations],
			ClearSystemCache[];
			{time, result} = TimeConstrained[Timing[f],OptionValue[MaxTime],{OptionValue[MaxTime],Null}];
			totaltime += time;
			iterations++;
		];
		{result,totaltime/iterations},
		f//N
	]
]


(* ::Subsubsection::Closed:: *)
(*Numerical Gaussian method*)


Options[NGaussianRiceProbabilityOfDetection] = {LowSNR->OptionValue[ProbabilityOfDetection,LowSNR],Timed->OptionValue[ProbabilityOfDetection,Timed],MaxTime->OptionValue[ProbabilityOfDetection,MaxTime],MaxIterations->OptionValue[ProbabilityOfDetection,MaxIterations]};
NGaussianRiceProbabilityOfDetection::usage="NGaussianRiceProbabilityOfDetection[M, \[Gamma], \[Lambda], m] calculates the approximate probability of detection for a single energy detector operating in a Rice-m fading channel using a numerical algorithm
NGaussianRiceProbabilityOfDetection[M, \[Gamma], \[Lambda], m, n] calculates the approximate probability of detection for a cooperative network operating in a Rice-m fading channel using a numerical algorithm.

A LowSNR option may be specified. By default, LowSNR\[Rule]"<>ToString[LowSNR/.Options[NGaussianRiceProbabilityOfDetection]]<>".

Function timing may be specified using the following options:

Timed\[Rule]"<>ToString[Timed/.Options[ProbabilityOfDetection]]<>"
MaxIterations\[Rule]"<>ToString[MaxIterations/.Options[ProbabilityOfDetection]]<>"
MaxTime\[Rule]"<>ToString[MaxTime/.Options[ProbabilityOfDetection]]<>"

where the above settings are the defaults, if not specified. If Timed\[Rule]True, then a {Pd, time} list of values will be returned.";
NGaussianRiceProbabilityOfDetection[M_?NumericQ,\[Gamma]_?NumericQ,\[Lambda]_,m_?NumericQ,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[NGaussianRiceProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	NGaussianRiceProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,RelevantOptions[NGaussianRiceProbabilityOfDetection]]
]
NGaussianRiceProbabilityOfDetection[M_?NumericQ,\[Gamma]_?NumericQ,\[Lambda]_,m_?NumericQ,n_?IntegerQ,OptionsPattern[]]:=Module[{lim, f, totaltime = 0, iterations = 0, time, result, RelevantOptions},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[NGaussianRiceProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	f:=Block[{$MaxExtraPrecision=\[Infinity]},
		N[
			NIntegrate[AWGNProbabilityOfDetection[M,x/n,\[Lambda],n,RelevantOptions[AWGNProbabilityOfDetection]]RicePDF[\[Gamma],m,x,n],{x,0,\[Infinity]}],
		20]
	];
	If[OptionValue[Timed],
		(* Evaluate result until MaxTime seconds of CPU time have been used or MaxIterations have been performed, whichever comes first *)
		While[totaltime < OptionValue[MaxTime] && iterations < OptionValue[MaxIterations],
			ClearSystemCache[];
			{time, result} = TimeConstrained[Timing[f],OptionValue[MaxTime],{OptionValue[MaxTime],Null}];
			totaltime += time;
			iterations++;
		];
		{result,totaltime/iterations},
		f
	]
]


(* ::Subsection::Closed:: *)
(*Sample complexity*)


Options[NRiceSampleComplexity]={Method->OptionValue[SampleComplexity,Method],LowSNR->OptionValue[SampleComplexity,LowSNR],Tolerance->OptionValue[SampleComplexity,Tolerance]};
NRiceSampleComplexity::usage="NRiceSampleComplexity[\[Gamma], Pf, Pd, m] calculates the sample complexity for a single energy detector operating on a Rice-m fading channel.
NRiceSampleComplexity[\[Gamma], Pf, Pd, m, n] calculates the sample complexity for a cooperative network operating on a Rice-m fading channel.

The following methods can be given:

Method\[Rule]\"Approximate\"
Method\[Rule]\"Exact\"

By default, Method\[Rule]\""<>ToString[Method/.Options[NRiceSampleComplexity]]<>"\".

If Method\[Rule]\"Approximate\", the LowSNR option may be specified. By default, LowSNR\[Rule]"<>ToString[LowSNR/.Options[NRiceSampleComplexity]]<>".

Numerical tolerance can be specified using the Tolerance option. By default, Tolerance\[Rule]"<>ToString[Tolerance/.Options[NRiceSampleComplexity]//N//InputForm]<>".";
NRiceSampleComplexity::tol="The difference between the result `1` and the constraint `2` was greater than the specified tolerance `3`.";
NRiceSampleComplexity[\[Gamma]_?NumericQ,Pf_?NumericQ,Pd_?NumericQ,m_?NumericQ,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[NRiceSampleComplexity][[All,1]]],Options[target][[All,1]]];
	NRiceSampleComplexity[\[Gamma],Pf,Pd,m,n,RelevantOptions[NRiceSampleComplexity]]
]
NRiceSampleComplexity[\[Gamma]_?NumericQ,Pf_?NumericQ,Pd_?NumericQ,m_?NumericQ,n_?IntegerQ,OptionsPattern[]]:=Module[{RelevantOptions, tol = OptionValue[Tolerance], intialGuess = Max[(20 / (n m^2)), 1] * SampleComplexity[\[Gamma],Pf,Pd,n], courseGuess, fineGuess, result},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[NRiceSampleComplexity][[All,1]]],Options[target][[All,1]]];
	(* Temporarily disable error checking - we'll do our own *)
	Off[FindRoot::reged,FindRoot::lstol];
	Switch[OptionValue[Method],
		"Approximate",
		(* Only use Gaussian method if it is valid *)
		If[intialGuess <= 250,
			result = NRiceSampleComplexity[\[Gamma],Pf,Pd,m,n,tol,Method->"Exact",Tolerance->OptionValue[Tolerance]];,
			fineGuess = M/.FindRoot[NRiceProbabilityOfDetection[M,\[Gamma],\[Lambda][M,Pf,n],m,n,RelevantOptions[NRiceProbabilityOfDetection]] == Pd, {M, intialGuess, 1, \[Infinity]}];
			result = NRiceProbabilityOfDetection[fineGuess,\[Gamma],\[Lambda][fineGuess,Pf,n],m,n,RelevantOptions[NRiceProbabilityOfDetection]];
		];,
		"Exact",
		(* If Gaussian approximation is valid, then use it to speed up the calculation *)
		If[intialGuess <= 250,
			courseGuess = intialGuess;,
			courseGuess = M/.FindRoot[NRiceProbabilityOfDetection[M,\[Gamma],\[Lambda][M,Pf,n],m,n,RelevantOptions[NRiceProbabilityOfDetection]] == Pd, {M, intialGuess, 1, \[Infinity]}];
		];
		fineGuess = M/.FindRoot[NRiceProbabilityOfDetection[courseGuess,\[Gamma],\[Lambda][courseGuess,Pf,n,RelevantOptions[\[Lambda]]],m,n] == Pd, {M, courseGuess, 1, \[Infinity]}];
		result = NRiceProbabilityOfDetection[fineGuess,\[Gamma],\[Lambda][fineGuess,Pf,n,RelevantOptions[\[Lambda]]],m,n];
	];
	On[FindRoot::reged,FindRoot::lstol];
	If[Abs[result - Pd] <= tol//TrueQ,
		fineGuess,
		Message[NRiceSampleComplexity::tol, result//N, Pd//N, tol//N]
	]
]


End[];


EndPackage[];
