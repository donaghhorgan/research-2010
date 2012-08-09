(* ::Package:: *)

(* ::Title:: *)
(*Error function approximations*)


(* ::Subsection::Closed:: *)
(*Copyright notice*)


(* ::Text:: *)
(*Mathematica approximations for the error function and complementary error function.*)
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


(* ::Subsection::Closed:: *)
(*Note*)


(* ::Text:: *)
(*This package implements some approximations for the error function and its complement as defined in "Versatile, Accurate, and Analytically Tractable Approximation for the Gaussian Q-Function" by L\[OAcute]pez-Ben\[IAcute]tez and Casadevall (IEEE TRANSACTIONS ON COMMUNICATIONS, VOL. 59, NO. 4, APRIL 2011).*)


(* ::Section:: *)
(*Public*)


BeginPackage["ErfApprox`"];


(* ::Subsection::Closed:: *)
(*Erf*)


ErfApprox;


NErfApprox;


(* ::Subsection::Closed:: *)
(*Erfc*)


ErfcApprox;


NErfcApprox;


(* ::Subsection::Closed:: *)
(*Diagnostics*)


NErfApproxError;


NErfcApproxError;


(* ::Subsection::Closed:: *)
(*Misc*)


LopezBenitezParameters;


a;


b;


c;


erfApproxAssumptions;


(* ::Section:: *)
(*Private*)


Begin["`Private`"];


(* ::Subsection::Closed:: *)
(*Erf*)


Options[ErfApprox] = {Method->"Lopez-Benitez"};
ErfApprox::usage="ErfApprox[x] calculates the approximate value of the error function at x.

The following methods can be given:

Method->\"Chiani\"
Method->\"Lopez-Benitez\"
Method->\"Loskot\"
Method->\"Karagiannidis\"

By default, Method->"<>ToString[Method/.Options[ErfApprox]]<>".";
ErfApprox[x_,OptionsPattern[]]:=1 - ErfcApprox[x,Method->OptionValue[Method]]


Options[NErfApprox] = Options[ErfApprox];
NErfApprox::usage="NErfcApprox[x] numerically evaluates the approximate value of the error function at x.

The following methods can be given:

Method->\"Chiani\"
Method->\"Lopez-Benitez\"
Method->\"Loskot\"
Method->\"Karagiannidis\"

By default, Method->"<>ToString[Method/.Options[NErfApprox]]<>".";
NErfApprox[x_,OptionsPattern[]]:=1 - NErfcApprox[x,Method->OptionValue[Method]]


(* ::Subsection::Closed:: *)
(*Erfc*)


Options[ErfcApprox] = Options[ErfApprox];
ErfcApprox::usage="ErfcApprox[x] calculates the approximate value of the complementary error function at x.

The following methods can be given:

Method->\"Chiani\"
Method->\"Lopez-Benitez\"
Method->\"Loskot\"
Method->\"Karagiannidis\"

By default, Method->"<>ToString[Method/.Options[ErfcApprox]]<>".";
ErfcApprox[x_,OptionsPattern[]]:=Switch[OptionValue[Method],
	"Chiani",
	erfApproxAssumptions = {a > 0, b > 0, c > 0};
	a Exp[-x^2] + b Exp[-c x^2],
	"Lopez-Benitez",
	erfApproxAssumptions = {a < 0, b < 0, c < 0};
	2Exp[2a x^2+Sqrt[2]b x+c],
	"Loskot",
	erfApproxAssumptions = {a > 0, b > 0, c > 0, d > 0, e > 0, f > 0};
	2(a*Exp[-2b x^2]+c*Exp[-2d x^2]+e Exp[-2f x^2]),
	"Karagiannidis",
	erfApproxAssumptions = {a > 0, b > 0};
	(1-Exp[-a x])Exp[-x^2]/(b Sqrt[Pi]x),
	_,
	ErfcApprox[x,Method->"Lopez-Benitez"]
]


Options[NErfcApprox] = Options[ErfApprox];
NErfcApprox::usage="NErfcApprox[x] numerically evaluates the approximate value of the complementary error function at x.

The following methods can be given:

Method->\"Chiani\"
Method->\"Lopez-Benitez\"
Method->\"Loskot\"
Method->\"Karagiannidis\"

By default, Method->"<>ToString[Method/.Options[NErfcApprox]]<>".";
NErfcApprox[x_,OptionsPattern[]]:=If[x < 0,
	2 - ErfcApprox[-x,Method->OptionValue[Method]],
	Switch[OptionValue[Method],
		"Chiani",
		ErfcApprox[x,Method->OptionValue[Method]]/.{a -> 1 / 6, b -> 1 / 2, c -> 4 / 3},
		"Lopez-Benitez",
		ErfcApprox[x,Method->OptionValue[Method]]/.LopezBenitezParameters[x],
		"Loskot",
		ErfcApprox[x,Method->OptionValue[Method]]/.{a -> 0.168, b -> 0.876, c -> 0.144, d -> 0.525, e -> 0.002, f -> 0.603},
		"Karagiannidis",
		ErfcApprox[x,Method->OptionValue[Method]]/.{a -> 1.98, b -> 1.135},
		_,
		NErfcApprox[x,Method->"Lopez-Benitez"]
	]
]


(* ::Subsection::Closed:: *)
(*Diagnostics*)


Options[NErfApproxError]={Method->"Lopez-Benitez"};
NErfApproxError::usage="NErfApproxError[x] calculates the error between the error function and its approximation.

The following methods can be given:

Method->\"Chiani\"
Method->\"Lopez-Benitez\"
Method->\"Loskot\"
Method->\"Karagiannidis\"

By default, Method->"<>ToString[Method/.Options[NErfApproxError]]<>".";
NErfcApproxError[x_,OptionsPattern[]]:=(Erfc[x] - NErfcApprox[x,Method->OptionValue[Method]]) / Erfc[x];


Options[NErfcApproxError]=Options[NErfApproxError];
NErfcApproxError::usage="NErfcApproxError[x] calculates the error between the complementary error function and its approximation.

The following methods can be given:

Method->\"Chiani\"
Method->\"Lopez-Benitez\"
Method->\"Loskot\"
Method->\"Karagiannidis\"

By default, Method->"<>ToString[Method/.Options[NErfcApproxError]]<>".";
NErfApproxError[x_,OptionsPattern[]]:=(Erf[x] - NErfApprox[x,Method->OptionValue[Method]]) / Erf[x];


(* ::Subsection::Closed:: *)
(*Misc*)


(* ::Text:: *)
(*Here, the min-MARE parameters are used:*)


LopezBenitezParameters::usage="LopezBenitezParameters[x] returns Lopez-Benitez's min-MARE parameters for x.";
LopezBenitezParameters[x_]:=If[x < 0,
	LopezBenitezParameters[-x],
	Which[0 <= x < 2,
		{a->-0.3976, b->-0.7418, c->-0.7019},
		2 <= x < 4,
		{a->-0.4369, b->-0.6511, c->-0.7358},
		4 <= x < 6,
		{a->-0.4577, b->-0.5695, c->-0.7864},
		6 <= x < 8,
		{a->-0.4698, b->-0.5026, c->-0.8444},
		8 <= x < 10,
		{a->-0.4774, b->-0.4484, c->-0.9049},
		10 <= x,
		{a->-0.4920, b->-0.2887, c->-1.1893}
	]
]


End[];


EndPackage[];
