(* ::Package:: *)

(* ::Title:: *)
(*SQLite database connectivity package*)


(* ::Subsection::Closed:: *)
(*Copyright notice*)


(* ::Text:: *)
(*Mathematica SQLite connectivity package,*)
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
(*Version 1.0: Basic SQLite functionality, minor bug fixes and functionality improvements to follow.*)


(* ::Section:: *)
(*Public*)


BeginPackage["SQLite`"];


SQLiteOpenDatabase;


SQLiteCloseDatabase;


SQLiteCreateTable;


SQLiteInsertRecord;


SQLiteLookupRecord;


(* ::Section:: *)
(*Private*)


Begin["`Private`"];


SQLiteOpenDatabase::usage="SQLiteOpenDatabase[databaseName] opens the database specified by the string databaseName and returns a file handle. If no database exists, then one is created.";
SQLiteOpenDatabase[db_String]:=Database`OpenDatabase[db];


SQLiteCloseDatabase::usage="SQLiteCloseDatabase[db] closes the database specified by the file handle db.";
SQLiteCloseDatabase[db_Database`Database]:=Database`CloseDatabase[db];


SQLiteCreateTable::usage="SQLiteCreateTable[db, tableName, columnNames, columnTypes] creates a table named tableName with columns specified by columnNames and columnTypes in the database specified by the handle db.";
SQLiteCreateTable[db_Database`Database, tableName_String, columnNames_List, columnTypes_List]:=Module[{query},
	query = "CREATE TABLE "<>tableName<>" (id INTEGER PRIMARY KEY, "<>StringJoin[Riffle[Riffle[ToString/@columnNames,ToString/@columnTypes],{" ",", "}]]<>");";
	Database`QueryDatabase[db, query]
];


SQLiteInsertRecord::usage="SQLiteInsertRecord[db, tableName, columnNames, values] inserts a new record into the specified table.";
SQLiteInsertRecord[db_Database`Database,tableName_String,columnNames_List,values_List]:=Module[{query},
	query = "INSERT INTO "<>tableName<>" ("<>StringJoin[Riffle[ToString/@columnNames,", "]]<>") values ("<>StringJoin[Riffle[ToString/@Table["?",{Length[columnNames]}],", "]]<>");";
    Database`QueryDatabase[db, query, values]
];


SQLiteLookupRecord::usage="SQLiteLookupRecord[db, tableName, columnNames, values] returns the value of the entries specified by tableName, columnNames and values.";
SQLiteLookupRecord[db_Database`Database,tableName_String,columnNames_List,values_List]:=Module[{query,validColumns,validValues,returnedColumns},
	validColumns = Delete[columnNames,Position[values,Null]];
	returnedColumns = Delete[columnNames,Complement[Table[{i},{i,Length[values]}],Position[values,Null]]];
	If[Length[returnedColumns]==0, returnedColumns = {"*"}];
	validValues = Delete[values,Position[values,Null]];
	validValues = Table[If[TrueQ[Head[validValues[[i]]]==String],"\""<>validValues[[i]]<>"\"",validValues[[i]]],{i,Length[validValues]}];
	query = "SELECT "<>StringJoin[Riffle[ToString/@returnedColumns,", "]]<>" FROM "<>tableName<>" WHERE "<>StringJoin[Riffle[Riffle[ToString/@validColumns,ToString/@validValues],If[Length[validColumns]==1,{" = "},{" = "," and "}]]]<>";";
    Database`QueryDatabase[db, query]
];


End[];


EndPackage[];
