(* ::Package:: *)

(* ::Title:: *)
(*Database logging functions*)


(* ::Subsection::Closed:: *)
(*Copyright notice*)


(* ::Text:: *)
(*Mathematica database logging functions.*)
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
(*Version 1.0: Basic database logging functionality, minor bug fixes to follow.*)


(* ::Section:: *)
(*Public*)


BeginPackage["DBLogging`"];


(* ::Subsection::Closed:: *)
(*Result fetching*)


GetProbabilityOfDetection;


GetSampleComplexity;


(* ::Subsection::Closed:: *)
(*Result caching*)


CacheProbabilityOfDetection;


CacheSampleComplexity;


(* ::Section:: *)
(*Private*)


Begin["`Private`"];


Needs["SQLite`"];


(* ::Subsection::Closed:: *)
(*Database parameters*)


databaseName="data.sqlite";


(* ::Subsection::Closed:: *)
(*Result fetching*)


(* ::Subsubsection::Closed:: *)
(*Probability of detection*)


GetProbabilityOfDetection::usage="GetProbabilityOfDetection[algorithm, channelType, M, \[Gamma], Pf, n, m] fetches the specified record from the database.";
GetProbabilityOfDetection[algorithm_?StringQ,channelType_?StringQ,M_?NumericQ,\[Gamma]_?NumericQ,Pf_?NumericQ,n_?IntegerQ,m_?NumericQ]:=Module[{db,result,tableName="data",columnNames={"algorithm","channelType","sampleComplexity","n","m","snrdb","pf","pd"},columnTypes={"TEXT","TEXT","NUMERIC","INTEGER","NUMERIC","INTEGER","NUMERIC","NUMERIC"}},
	If[FileExistsQ[databaseName],
		db = SQLiteOpenDatabase[databaseName];
		result = SQLiteLookupRecord[db,tableName,columnNames,{algorithm,channelType,M,n,m,10Log[10,\[Gamma]]//Round,Pf//N,Null}];
		SQLiteCloseDatabase[db];,
		result = Undefined;
	];
	If[result == {}//TrueQ,
		Null,
		result[[1]][[1]]
	]
]


(* ::Subsubsection::Closed:: *)
(*Sample complexity*)


GetSampleComplexity::usage="GetSampleComplexity[channel, precision, n, \[Gamma], Pf, Pm] fetches the specified record from the database.";
GetSampleComplexity[channel_,precision_,n_?IntegerQ,\[Gamma]_?NumericQ,Pf_?NumericQ,Pm_?NumericQ]:=Module[{db,result,p,channelType,m,tableName="sample_complexity",columnNames={"channel_type","precision","number_of_nodes","fading_parameter","snr_db","probability_of_false_alarm","probability_of_missed_detection","sample_complexity"},columnTypes={"TEXT","TEXT","INTEGER","REAL","REAL","REAL","REAL","REAL"}},
	If[FileExistsQ[databaseName],
		db = SQLiteOpenDatabase[databaseName];
		If[ListQ[channel],
			{channelType,m}=channel,
			{channelType,m}={channel,1}
		];
		p = If[precision==\[Infinity],"Infinity",precision];
		result = SQLiteLookupRecord[db,tableName,columnNames,{channelType,p,n,m,Round[10Log[10,\[Gamma]]],Pf//N,Pm//N,Null}];
		SQLiteCloseDatabase[db];,
		result = Undefined;
	];
	If[result == {}//TrueQ,
		Null,
		result[[1]][[1]]
	]
]


(* ::Subsection::Closed:: *)
(*Result caching*)


(* ::Subsubsection::Closed:: *)
(*Probability of detection*)


CacheProbabilityOfDetection::usage="CacheProbabilityOfDetection[algorithm, channelType, M, \[Gamma], Pf, n, m, result] caches the specified record in the database.";
CacheProbabilityOfDetection[algorithm_?StringQ,channelType_?StringQ,M_?NumericQ,\[Gamma]_?NumericQ,Pf_?NumericQ,n_?IntegerQ,m_?NumericQ,result_?NumericQ]:=Module[{db,tableName="data",columnNames={"algorithm","channelType","sampleComplexity","n","m","snrdb","pf","pd"},columnTypes={"TEXT","TEXT","NUMERIC","INTEGER","NUMERIC","INTEGER","NUMERIC","NUMERIC"}},
	If[!FileExistsQ[databaseName],
		db = SQLiteOpenDatabase[databaseName];
		SQLiteCreateTable[db, tableName, columnNames, columnTypes],
		db = SQLiteOpenDatabase[databaseName];
	];
	SQLiteInsertRecord[db,tableName,columnNames,{algorithm,channelType,M,n,m,10Log[10,\[Gamma]]//Round,Pf//N,result,algorithm}];
	SQLiteCloseDatabase[db];
]


(* ::Subsubsection::Closed:: *)
(*Sample complexity*)


CacheSampleComplexity::usage="CacheSampleComplexity[channel, precision, n, \[Gamma], Pf, Pm, result] caches the specified record in the database.";
CacheSampleComplexity[channel_,precision_,n_?IntegerQ,\[Gamma]_?NumericQ,Pf_?NumericQ,Pm_?NumericQ,result_?NumericQ]:=Module[{db,p,channelType,m,tableName="sample_complexity",columnNames={"channel_type","precision","number_of_nodes","fading_parameter","snr_db","probability_of_false_alarm","probability_of_missed_detection","sample_complexity"},columnTypes={"TEXT","TEXT","INTEGER","REAL","REAL","REAL","REAL","REAL"}},
	If[!FileExistsQ[databaseName],
		db = SQLiteOpenDatabase[databaseName];
		SQLiteCreateTable[db, tableName, columnNames, columnTypes],
		db = SQLiteOpenDatabase[databaseName];
	];
	If[ListQ[channel],
		{channelType,m}=channel,
		{channelType,m}={channel,1}
	];
	p = If[precision==\[Infinity],"Infinity",precision];
	SQLiteInsertRecord[db,tableName,columnNames,{channelType,p,n,m,Round[10Log[10,\[Gamma]]],Pf//N,Pm//N,result}];
	SQLiteCloseDatabase[db];
]


End[];


EndPackage[];
