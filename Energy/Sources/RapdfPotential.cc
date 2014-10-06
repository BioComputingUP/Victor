/*  This file is part of Victor.

    Victor is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Victor is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Victor.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
 *@Class               RapdfPotential 
 *@Project   Victor
 *@Description 
*    This class implements the all-atom residue specific probability 
*    discriminatory function from Samudrala & Moult (JMB 1998).
*/
// Includes:
#include <RapdfPotential.h>
#include <AminoAcidCode.h>
#include <Spacer.h>
#include <cstring>

using namespace Biopool;

// Global constants, typedefs, etc. (to avoid):
 

// CONSTRUCTORS/DESTRUCTOR:

/**
 *@Description    basic constructor, allocate the information from the ram.par file 
 */
RapdfPotential::RapdfPotential(){
     
     char *victor=getenv("VICTOR_ROOT");
	if (victor  == NULL)
		ERROR("Environment variable VICTOR_ROOT was not found.\n Use the command:\n export VICTOR_ROOT=......", exception);
     path = "data/ram.par";
    string inputFile = getenv("VICTOR_ROOT");
    if (inputFile.length() < 3)
       ERROR("Environment variable VICTOR_ROOT was not found.", exception);
    inputFile += path;
  ifstream in(inputFile.c_str());
  if (!in)
    ERROR("Could not read data file.", exception);
  for (unsigned int i = 0; i < MAX_BINS; i++)
    for (unsigned int j = 0; j < MAX_TYPES; j++)      {
	string dummy;
	in >> dummy;
	for (unsigned int k = 0; k < MAX_TYPES; k++)
	  in >> prob[i][j][k];
      }
  in.close();
}

// PREDICATES:

/**
 *@Description Calculates the energy for the amino acids in the Spacer
*@param   spacer reference(Spacer&),
 *@return energy value (long double)
*/
long double RapdfPotential::calculateEnergy(Spacer& sp){
  long double en = 0.0;
  unsigned int size = sp.sizeAmino();
  for (unsigned int i = 0; i < size-1; i++){
	AminoAcid& aa = sp.getAmino(i);
	string aaType = aa.getType();
	for (unsigned int ii = i+1; ii < size; ii++) {
	    AminoAcid& aa2 = sp.getAmino(ii);
	    string aaType2 = aa2.getType();
	    for (unsigned int j = 0; j < aa.size(); j++)
	      for (unsigned int k = 0; k < aa2.size(); k++)
		en += calculateEnergy(aa[j], aa2[k], aaType, aaType2);
	  }
      }
  
  return en;
}
/**
 *@Description Calculates the energy for a group of amino acids in the Spacer limited by the index values as start and end
*@param   spacer reference(Spacer&), start and end positions (unsigned int,unsigned int)
 *@return energy value (long double)
*/
long double RapdfPotential::calculateEnergy(Spacer& sp, unsigned int index1, unsigned int index2){
  long double en = 0.0;
  for (unsigned int i = index1; i < index2; i++) {
      AminoAcid& aa = sp.getAmino(i);
      string aaType = aa.getType();
      for (unsigned int ii = i+1; ii < index2; ii++){
	  AminoAcid& aa2 = sp.getAmino(ii);
	  string aaType2 = aa2.getType();
	  for (unsigned int j = 0; j < aa.size(); j++)
	    for (unsigned int k = 0; k < aa2.size(); k++)
	      en += calculateEnergy(aa[j], aa2[k], aaType, aaType2);
	}
    }
  return en;
}

/**
 *@Description Calculates the energy for the amino acids in the Spacer that are different to the given aa
*@param  amino acid reference(AminoAcid&),spacer reference(Spacer&),
 *@return energy value (long double)
*/
long double RapdfPotential::calculateEnergy(AminoAcid& aa, Spacer& sp){
  long double en = 0.0;
  string aaType = aa.getType();
  for (unsigned int i = 0; i < sp.sizeAmino(); i++){
      AminoAcid& aa2 = sp.getAmino(i);
      if (aa == aa2) // exclude self-energy
	continue;
      string aaType2 = aa2.getType();
      for (unsigned int j = 0; j < aa.size(); j++)
	for (unsigned int k = 0; k < aa2.size(); k++)
	  en += calculateEnergy(aa[j], aa2[k], aaType, aaType2);
    }
  return en;
}

/**
 *@Description Calculates the energy between 2 amino acids
*@param  amino acids references( AminoAcid&,AminoAcid&)
 *@return energy value (long double)
*/
long double RapdfPotential::calculateEnergy(AminoAcid& aa, AminoAcid& aa2){
  if (aa == aa2) // exclude self-energy
    return 0.0;
  long double en = 0.0;
  string aaType = aa.getType();
  string aaType2 = aa2.getType();
  for (unsigned int j = 0; j < aa.size(); j++)
    for (unsigned int k = 0; k < aa2.size(); k++)
      en += calculateEnergy(aa[j], aa2[k], aaType, aaType2);

  return en;
}

// MODIFIERS:

/******************************************************************/
// HELPERS:
/**
 *@Description Retrieves the corresponding index for the given distance
*@param  distance (double)
 *@return corresponding index (unsigned int)
*/
unsigned int RapdfPotential::pGetDistanceBinOne(double distance)
/* No zeros for the distance bin */{
  if (distance < 10.0){
      if (distance < 3.0) return 0;
      if (distance < 4.0) return 1;
      if (distance < 5.0) return 2;
      if (distance < 6.0) return 3;
      if (distance < 7.0) return 4;
      if (distance < 8.0) return 5;
      if (distance < 9.0) return 6;
      return 7;
    }
  else    {
      if (distance < 11.0) return 8;
      if (distance < 12.0) return 9;
      if (distance < 13.0) return 10;
      if (distance < 14.0) return 11;
      if (distance < 15.0) return 12;
      if (distance < 16.0) return 13;
      if (distance < 17.0) return 14;
      if (distance < 18.0) return 15;
      if (distance < 19.0) return 16;
      if (distance < 20.0) return 17;
    }
  return 999;
}

/******************************************************************/
/**
 *@Description Retrieves the corresponding index for the given group(AAtype_AtomType)
*@param  group name (char*)
 *@return corresponding index (unsigned int)
*/
unsigned int RapdfPotential::pGetGroupBin(const char* group_name){
  switch (group_name[0]){
    case 'A':
      if (strcmp(group_name, "AN") == 0) return 1;
      if (strcmp(group_name, "ACA") == 0) return 2;
      if (strcmp(group_name, "AC") == 0) return 3;
      if (strcmp(group_name, "AO") == 0) return 4;
      if (strcmp(group_name, "ACB") == 0) return 5;
      break;
    case 'C':
      if (strcmp(group_name, "CN") == 0) return 6;
      if (strcmp(group_name, "CCA") == 0) return 7;
      if (strcmp(group_name, "CC") == 0) return 8;
      if (strcmp(group_name, "CO") == 0) return 9;
      if (strcmp(group_name, "CCB") == 0) return 10;
      if (strcmp(group_name, "CSG") == 0) return 11;
      break;
    case 'D':
      if (strcmp(group_name, "DN") == 0) return 12;
      if (strcmp(group_name, "DCA") == 0) return 13;
      if (strcmp(group_name, "DC") == 0) return 14;
      if (strcmp(group_name, "DO") == 0) return 15;
      if (strcmp(group_name, "DCB") == 0) return 16;
      if (strcmp(group_name, "DCG") == 0) return 17;
      if (strcmp(group_name, "DOD1") == 0) return 18;
      if (strcmp(group_name, "DOD2") == 0) return 19;
      break;
    case 'E':
      if (strcmp(group_name, "EN") == 0) return 20;
      if (strcmp(group_name, "ECA") == 0) return 21;
      if (strcmp(group_name, "EC") == 0) return 22;
      if (strcmp(group_name, "EO") == 0) return 23;
      if (strcmp(group_name, "ECB") == 0) return 24;
      if (strcmp(group_name, "ECG") == 0) return 25;
      if (strcmp(group_name, "ECD") == 0) return 26;
      if (strcmp(group_name, "EOE1") == 0) return 27;
      if (strcmp(group_name, "EOE2") == 0) return 28;
      break;
    case 'F':
      if (strcmp(group_name, "FN") == 0) return 29;
      if (strcmp(group_name, "FCA") == 0) return 30;
      if (strcmp(group_name, "FC") == 0) return 31;
      if (strcmp(group_name, "FO") == 0) return 32;
      if (strcmp(group_name, "FCB") == 0) return 33;
      if (strcmp(group_name, "FCG") == 0) return 34;
      if (strcmp(group_name, "FCD1") == 0) return 35;
      if (strcmp(group_name, "FCD2") == 0) return 36;
      if (strcmp(group_name, "FCE1") == 0) return 37;
      if (strcmp(group_name, "FCE2") == 0) return 38;
      if (strcmp(group_name, "FCZ") == 0) return 39;
      break;
    case 'G':
      if (strcmp(group_name, "GN") == 0) return 40;
      if (strcmp(group_name, "GCA") == 0) return 41;
      if (strcmp(group_name, "GC") == 0) return 42;
      if (strcmp(group_name, "GO") == 0) return 43;
      break;
    case 'H':
      if (strcmp(group_name, "HN") == 0) return 44;
      if (strcmp(group_name, "HCA") == 0) return 45;
      if (strcmp(group_name, "HC") == 0) return 46;
      if (strcmp(group_name, "HO") == 0) return 47;
      if (strcmp(group_name, "HCB") == 0) return 48;
      if (strcmp(group_name, "HCG") == 0) return 49;
      if (strcmp(group_name, "HND1") == 0) return 50;
      if (strcmp(group_name, "HCD2") == 0) return 51;
      if (strcmp(group_name, "HCE1") == 0) return 52;
      if (strcmp(group_name, "HNE2") == 0) return 53;
      break;
    case 'I':
      if (strcmp(group_name, "IN") == 0) return 54;
      if (strcmp(group_name, "ICA") == 0) return 55;
      if (strcmp(group_name, "IC") == 0) return 56;
      if (strcmp(group_name, "IO") == 0) return 57;
      if (strcmp(group_name, "ICB") == 0) return 58;
      if (strcmp(group_name, "ICG1") == 0) return 59;
      if (strcmp(group_name, "ICG2") == 0) return 60;
      if (strcmp(group_name, "ICD1") == 0) return 61;
      break;
    case 'K':
      if (strcmp(group_name, "KN") == 0) return 62;
      if (strcmp(group_name, "KCA") == 0) return 63;
      if (strcmp(group_name, "KC") == 0) return 64;
      if (strcmp(group_name, "KO") == 0) return 65;
      if (strcmp(group_name, "KCB") == 0) return 66;
      if (strcmp(group_name, "KCG") == 0) return 67;
      if (strcmp(group_name, "KCD") == 0) return 68;
      if (strcmp(group_name, "KCE") == 0) return 69;
      if (strcmp(group_name, "KNZ") == 0) return 70;
      break;
    case 'L':
      if (strcmp(group_name, "LN") == 0) return 71;
      if (strcmp(group_name, "LCA") == 0) return 72;
      if (strcmp(group_name, "LC") == 0) return 73;
      if (strcmp(group_name, "LO") == 0) return 74;
      if (strcmp(group_name, "LCB") == 0) return 75;
      if (strcmp(group_name, "LCG") == 0) return 76;
      if (strcmp(group_name, "LCD1") == 0) return 77;
      if (strcmp(group_name, "LCD2") == 0) return 78;
      break;
    case 'M':
      if (strcmp(group_name, "MN") == 0) return 79;
      if (strcmp(group_name, "MCA") == 0) return 80;
      if (strcmp(group_name, "MC") == 0) return 81;
      if (strcmp(group_name, "MO") == 0) return 82;
      if (strcmp(group_name, "MCB") == 0) return 83;
      if (strcmp(group_name, "MCG") == 0) return 84;
      if (strcmp(group_name, "MSD") == 0) return 85;
      if (strcmp(group_name, "MCE") == 0) return 86;
      break;
    case 'N':
      if (strcmp(group_name, "NN") == 0) return 87;
      if (strcmp(group_name, "NCA") == 0) return 88;
      if (strcmp(group_name, "NC") == 0) return 89;
      if (strcmp(group_name, "NO") == 0) return 90;
      if (strcmp(group_name, "NCB") == 0) return 91;
      if (strcmp(group_name, "NCG") == 0) return 92;
      if (strcmp(group_name, "NOD1") == 0) return 93;
      if (strcmp(group_name, "NND2") == 0) return 94;
      break;
    case 'P':
      if (strcmp(group_name, "PN") == 0) return 95;
      if (strcmp(group_name, "PCA") == 0) return 96;
      if (strcmp(group_name, "PC") == 0) return 97;
      if (strcmp(group_name, "PO") == 0) return 98;
      if (strcmp(group_name, "PCB") == 0) return 99;
      if (strcmp(group_name, "PCG") == 0) return 100;
      if (strcmp(group_name, "PCD") == 0) return 101;
      break;
    case 'Q':
      if (strcmp(group_name, "QN") == 0) return 102;
      if (strcmp(group_name, "QCA") == 0) return 103;
      if (strcmp(group_name, "QC") == 0) return 104;
      if (strcmp(group_name, "QO") == 0) return 105;
      if (strcmp(group_name, "QCB") == 0) return 106;
      if (strcmp(group_name, "QCG") == 0) return 107;
      if (strcmp(group_name, "QCD") == 0) return 108;
      if (strcmp(group_name, "QOE1") == 0) return 109;
      if (strcmp(group_name, "QNE2") == 0) return 110;
      break;
    case 'R':
      if (strcmp(group_name, "RN") == 0) return 111;
      if (strcmp(group_name, "RCA") == 0) return 112;
      if (strcmp(group_name, "RC") == 0) return 113;
      if (strcmp(group_name, "RO") == 0) return 114;
      if (strcmp(group_name, "RCB") == 0) return 115;
      if (strcmp(group_name, "RCG") == 0) return 116;
      if (strcmp(group_name, "RCD") == 0) return 117;
      if (strcmp(group_name, "RNE") == 0) return 118;
      if (strcmp(group_name, "RCZ") == 0) return 119;
      if (strcmp(group_name, "RNH1") == 0) return 120;
      if (strcmp(group_name, "RNH2") == 0) return 121;
      break;
    case 'S':
      if (strcmp(group_name, "SN") == 0) return 122;
      if (strcmp(group_name, "SCA") == 0) return 123;
      if (strcmp(group_name, "SC") == 0) return 124;
      if (strcmp(group_name, "SO") == 0) return 125;
      if (strcmp(group_name, "SCB") == 0) return 126;
      if (strcmp(group_name, "SOG") == 0) return 127;
      break;
    case 'T':
      if (strcmp(group_name, "TN") == 0) return 128;
      if (strcmp(group_name, "TCA") == 0) return 129;
      if (strcmp(group_name, "TC") == 0) return 130;
      if (strcmp(group_name, "TO") == 0) return 131;
      if (strcmp(group_name, "TCB") == 0) return 132;
      if (strcmp(group_name, "TOG1") == 0) return 133;
      if (strcmp(group_name, "TCG2") == 0) return 134;
      break;
    case 'V':
      if (strcmp(group_name, "VN") == 0) return 135;
      if (strcmp(group_name, "VCA") == 0) return 136;
      if (strcmp(group_name, "VC") == 0) return 137;
      if (strcmp(group_name, "VO") == 0) return 138;
      if (strcmp(group_name, "VCB") == 0) return 139;
      if (strcmp(group_name, "VCG1") == 0) return 140;
      if (strcmp(group_name, "VCG2") == 0) return 141;
      break;
    case 'W':
      if (strcmp(group_name, "WN") == 0) return 142;
      if (strcmp(group_name, "WCA") == 0) return 143;
      if (strcmp(group_name, "WC") == 0) return 144;
      if (strcmp(group_name, "WO") == 0) return 145;
      if (strcmp(group_name, "WCB") == 0) return 146;
      if (strcmp(group_name, "WCG") == 0) return 147;
      if (strcmp(group_name, "WCD1") == 0) return 148;
      if (strcmp(group_name, "WCD2") == 0) return 149;
      if (strcmp(group_name, "WNE1") == 0) return 150;
      if (strcmp(group_name, "WCE2") == 0) return 151;
      if (strcmp(group_name, "WCE3") == 0) return 152;
      if (strcmp(group_name, "WCZ2") == 0) return 153;
      if (strcmp(group_name, "WCZ3") == 0) return 154;
      if (strcmp(group_name, "WCH2") == 0) return 155;
      break;
    case 'Y':
      if (strcmp(group_name, "YN") == 0) return 156;
      if (strcmp(group_name, "YCA") == 0) return 157;
      if (strcmp(group_name, "YC") == 0) return 158;
      if (strcmp(group_name, "YO") == 0) return 159;
      if (strcmp(group_name, "YCB") == 0) return 160;
      if (strcmp(group_name, "YCG") == 0) return 161;
      if (strcmp(group_name, "YCD1") == 0) return 162;
      if (strcmp(group_name, "YCD2") == 0) return 163;
      if (strcmp(group_name, "YCE1") == 0) return 164;
      if (strcmp(group_name, "YCE2") == 0) return 165;
      if (strcmp(group_name, "YCZ") == 0) return 166;
      if (strcmp(group_name, "YOH") == 0) return 167;
      break;
    }

  return 999;
}

/******************************************************************/
