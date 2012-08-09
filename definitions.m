(* ::Package:: *)

(* ::Title:: *)
(*Definitions*)


(* ::Subsection::Closed:: *)
(*Copyright notice*)


(* ::Text:: *)
(*Mathematica cooperative energy detection and fading channel function definitions.*)
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


(* ::Subsection:: *)
(*Version information*)


(* ::Text:: *)
(*27/06/2012*)
(*1.01*)


(* ::Subsection::Closed:: *)
(*Changelog*)


(* ::Text:: *)
(*Version 1.01: Grouped similar functionality into individual packages.*)
(*Version 1.0: First working version, bug fixes to follow.*)


(* ::Subsection::Closed:: *)
(*Notation guide*)


(* ::Text:: *)
(*Pf: Probability of false alarm*)
(*Pd: Probability of detection*)
(*M: Number of samples*)
(*\[Gamma]: Instantaneous signal to noise ratio*)
(*\bar{\[Gamma]}: Average instantaneous signal to noise ratio*)
(*n: Number of nodes*)
(*m: Nakagami fading parameter*)
(*\[Lambda]: Decision threshold at the local node or fusion center (the context will make clear which)*)
(*k: Fusion center voting rule*)
(*\[Rho]: Average correlation coefficient between nodes*)


(* ::Subsection:: *)
(*Load packages*)


<<QFunction`
<<AWGN`
<<Nakagami`
<<Rayleigh`
<<Rice`
<<Network`
