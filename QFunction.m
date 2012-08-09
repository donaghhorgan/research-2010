(* ::Package:: *)

(* ::Title:: *)
(*Q-function*)


(* ::Subsection::Closed:: *)
(*Copyright notice*)


(* ::Text:: *)
(*Mathematica implementation of the Gaussian Q-function.*)
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
(*Version 1.0: Basic implemenation of the Gaussian Q-function.*)


(* ::Section:: *)
(*Public*)


BeginPackage["QFunction`"];


(* ::Subsection::Closed:: *)
(*Q-function*)


Q;


(* ::Subsection::Closed:: *)
(*Inverse Q-function*)


InverseQ;


(* ::Section:: *)
(*Private*)


Begin["`Private`"];


(* ::Subsection::Closed:: *)
(*Q-function*)


Q::usage="Q[x] calculates the value of Gaussian Q-function at x.";
Q[x_]:=1-CDF[NormalDistribution[], x]


(* ::Subsection::Closed:: *)
(*Inverse Q-function*)


InverseQ::usage="InverseQ[P] calculates the value of the inverse of the Gaussian Q-function at P."
InverseQ[P_]:=InverseCDF[NormalDistribution[], 1 - P]


End[];


EndPackage[];
