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
/*
 
*@Description:     Translator: PDB names to internal one-word-code and 
*                  vice versa. Provides some simple predicates dealing with 
*                  one-word-code. 
*
*@NB: AtomTranslator converts some "illegal" entries (eg. 1HD, 2HE, ...)
* to their "legal" equivalent (eg. HD, HE, ...)
*@Copyright:       This file contains information from the Bioinformatics 
*                  Template Library (BTL).
*
*                  Copyright (C) 1997,1998 Birkbeck College, Malet Street, 
*                  London WC1E 7HX, U.K. (classlib@mail.cryst.bbk.ac.uk)
*
*/

#ifndef __ATOM_CODE_H__
#define __ATOM_CODE_H__

// Includes

#include <string>
#include <Debug.h>

/** Internal one-word-code for PDB names.

    This table is taken from the Bioinformatics Template Library (BTL).

    The order of this table is the standard order in which atoms of an
    AminoAcid are stored. The idea is that if you point at atom i then
    all atoms with indices > i "follow" atom[i] in a natural sense. Of
    course there are some special cases, but in general I think this
    is a good idea. ;-)

    The default AtomCode X is moved to the very back of this list.

    */
enum AtomCode 
{ 
  // Backbone atoms
  N,   
  CA,
  C,   
  O,   
  // Side chain atoms
  CB,
  SG,  
  OG,  
  CG,  
  OG1, 
  CG1, 
  CG2,
  CD,  
  OD,  
  SD,  
  CD1, 
  OD1, 
  ND1, 
  CD2, 
  OD2, 
  ND2, 
  CE,  
  NE,  
  CE1, 
  OE1, 
  NE1, 
  CE2, 
  OE2, 
  NE2, 
  CE3, 
  CZ,  
  NZ,  
  CZ2, 
  CZ3, 
  OH,  
  NH1, 
  NH2, 
  CH2, 
  OXT, 
  P,   
  O1P, 
  O2P, 
  O5S, 
  C5S, 
  C4S, 
  O4S, 
  C3S, 
  O3S, 
  C2S, 
  O2S, 
  C1S, 
  N9,  
  C8,  
  N7,  
  C5,  
  C6,  
  O6,  
  N6,  
  N1,  
  C2,  
  O2,  
  N2,  
  N3,  
  C4,  
  O4,  
  N4,  
  C5M, 

  // Hydrogens
  H,
  H1,
  H2,
  H3,
  HA,
  HA2,
  HA3,
  HB,
  HB1,
  HB2,
  HB3,
  HD1,
  HD2,
  HD3,
  HD11,
  HD12,
  HD13,
  HD21,
  HD22,
  HD23,
  HE,
  HE1,
  HE2,
  HE3,
  HE21,
  HE22,
  HG,
  HG1,
  HG2,
  HG11,
  HG12,
  HG13,
  HG21,
  HG22,
  HG23,
  HG3,
  HH,
  HH2,
  HH11,
  HH12,
  HH21,
  HH22,
  HZ,
  HZ1,
  HZ2,
  HZ3,
  
  // side chain pseudo atoms:
  XA,
  XR,
  XN,
  XD,
  XC,
  XQ,
  XE,
  XG,
  XH,
  XI,
  XL,
  XK,
  XM,
  XF,
  XP,
  XS,
  XT,
  XW,
  XY,
  XV,

  X,     //  Please leave as last element! Corresponds to unknown atom type 
  ATOM_CODE_SIZE // number of atom types
};

// ---------------------------------------------------------------------------
//                                AtomTranslator
// -----------------x-------------------x-------------------x-----------------

/** Is second behind first in amino acid? This predicate is used to
    determine if atom second is to be moved if atom first changed its
    position. 

    NOT COMPLETE YET! SPECIAL HANDLING FOR SIDE CHAIN REQUIRED! */
inline
bool
follows(AtomCode first, AtomCode second)
{
  if (first == ATOM_CODE_SIZE || second == ATOM_CODE_SIZE)
    {
      return false;
    }
  else
    {
      return first < second;
    }
}

///check if atom is side chain beta atom
inline
bool 
isBetaAtom(AtomCode code) 
{
    return ( code == CB  );
}

///check if atom is side chain gamma atom
inline
bool 
isGammaAtom(AtomCode code) 
{
    return (( code == SG  ) || ( code == OG ) || ( code == CG )
	    || ( code == OG1 ) || ( code == CG1 ) || ( code == CG2 ));
}

///check if atom is side chain delta atom
inline
bool 
isDeltaAtom(AtomCode code) 
{
    return (( code == CD  ) || ( code == OD ) || ( code == SD )
	    || ( code == CD1 ) || ( code == OD1 ) || ( code == ND1 )
	    || ( code == CD2 ) || ( code == OD2 ) || ( code == ND2 ));
}

///check if atom is side chain epsilon atom
inline
bool 
isEpsilonAtom(AtomCode code) 
{
    return (( code == CE  ) || ( code == NE ) || ( code == CE1 )
	    || ( code == OE1 ) || ( code == NE1 ) || ( code == CE2 )
	    || ( code == OE2 ) || ( code == NE2 ) || ( code = CE3 ));
}

///check if atom is side chain zeta atom
inline
bool 
isZetaAtom(AtomCode code) 
{
    return (( code == CZ  ) || ( code == NZ ) || ( code == CZ2 )
	    || ( code == CZ3 ));
}

///check if atom is side chain eta atom
inline
bool 
isEtaAtom(AtomCode code) 
{
    return (( code == OH  ) || ( code == NH1 ) || ( code == NH2 )
	    || ( code == CH2 ));
}

/** is atom a C atom ?  */
// To be verified! 
inline
bool 
isCAtom(AtomCode code) 
{
  return (( code == C   ) || ( code == CB  ) || ( code == CG  ) 
	  || ( code == CG1 ) || ( code == CG2 ) || ( code == CD  ) 
	  || ( code == CD1 ) || ( code == CD2 ) || ( code == CE  ) 
	  || ( code == CE1 ) || ( code == CE2 ) || ( code == CE3 ) 
	  || ( code == CZ  ) || ( code == CZ2 ) || ( code == CZ3 ) 
	  || ( code == C5S ) || ( code == C4S ) || ( code == C3S ) 
	  || ( code == C2S ) || ( code == C1S ) || ( code == C8  ) 
	  || ( code == C5  ) || ( code == C6  ) || ( code == C2  ) 
	  || ( code == C4  ) || ( code == C5M ) );
}

/** is atom a C atom ?  */
// To be verified! 
inline
bool 
isNAtom(AtomCode code) 
{
  return ( (code == N9) || (code == N7) || (code == N6)
	   || (code == N1) || (code == N2) || (code == N3)
	   || (code == N4) || (code == ND1) || (code == ND2)
	   || (code == NE) || (code == NE1) || (code == NE2)
	   || (code == NZ) || (code == NH1) || (code == NH2)
	   || (code == N) );
}

/** is atom a C atom ?  */
// To be verified! 
inline
bool 
isOAtom(AtomCode code) 
{
  return ( (code == OXT) || (code == O1P) || (code == O2P)
	   || (code == O5S) || (code == O4S) || (code == O3S)
	   || (code == O2S) || (code == O6) || (code == O2)
	   || (code == O4) || (code == OD1) || (code == OD2)
	   || (code == OE1) || (code == OE2) || (code == OH)
	   || (code == O) || (code == OG) || (code == OG1)
	   || (code == OD) );  
}

/** is code specifying an H atom ? */
// Updated 2014 by Damiano Piovesan
inline
bool 
isHAtom(AtomCode code) 
{
  return (  (code == H) || (code == H1) || (code == H2) ||  (code == H3) ||  (code == HA) ||
            (code == HA2) || (code == HA3) || (code == HB) || (code == HB1) ||
            (code == HB2) || (code == HB3) || (code == HD1) || (code == HD2) ||
            (code == HD3) || (code == HD11) || (code == HD12) || (code == HD13) ||
            (code == HD21) || (code == HD22) || (code == HD23) || (code == HE) ||
            (code == HE1) || (code == HE2) || (code == HE3) || (code == HE21) ||
            (code == HE22) || (code == HG) || (code == HG1) ||  (code == HG2) ||
            (code == HG11) || (code == HG12) || (code == HG13) || 
            (code == HG21) || (code == HG22) ||
            (code == HG23) || (code == HG3) || (code == HH) || (code == HH2) ||
            (code == HH11) || (code == HH12) || (code == HH21) || (code == HH22) ||
            (code == HZ) || (code == HZ1) || (code == HZ2) || (code == HZ3)
	 );
}

/**
   true, if atom type name is known 
   magic code X corresponds to an unknown atom code
*/
inline
bool
isKnownAtom(AtomCode code)
{
  return !(code == X);
}


/** Is code a non-H atom? */
inline
bool
isHeavyAtom(AtomCode code)
{
  return (isKnownAtom(code))&&(!isHAtom(code));
}

// To be verified! 
inline
bool 
isBackboneAtom(AtomCode code) 
{
  return code == C || code == CA || code == N || code == O || code == OXT;
}

/** Translate string into atom code enum. */
inline
AtomCode
AtomTranslator(const string& name)
{
  DUMP(name);
  if (name == "")
    {
      return X;
    }
  if (name == "X")
    {
      return X;
    }
  if (name == "N")
    {
      return N;
    }
  else if (name == "CA")
    {
      return CA;
    }
  else if (name == "C")
    {
      return C;
    }
  else if (name == "O")
    {
      return O;
    }
  else if (name == "CB")
    {
      return CB;
    }
  else if (name == "SG")
    {
      return SG;
    }
  else if (name == "OG")
    {
      return OG;
    }
  else if (name == "CG")
    {
      return CG;
    }
  else if (name == "OG1")
    {
      return OG1;
    }
  else if (name == "CG1")
    {
      return CG1;
    }
  else if (name == "CG2")
    {
      return CG2;
    }
  else if (name == "CD")
    {
      return CD;
    }
  else if (name == "OD")
    {
      return OD;
    }
  else if (name == "SD")
    {
      return SD;
    }
  else if (name == "CD1")
    {
      return CD1;
    }
  else if (name == "OD1")
    {
      return OD1;
    }
  else if (name == "ND1")
    {
      return ND1;
    }
  else if (name == "CD2")
    {
      return CD2;
    }
  else if (name == "OD2")
    {
      return OD2;
    }
  else if (name == "ND2")
    {
      return ND2;
    }
  else if (name == "CE")
    {
      return CE;
    }
  else if (name == "NE")
    {
      return NE;
    }
  else if (name == "CE1")
    {
      return CE1;
    }
  else if (name == "OE1")
    {
      return OE1;
    }
  else if (name == "NE1")
    {
      return NE1;
    }
  else if (name == "CE2")
    {
      return CE2;
    }
  else if (name == "OE2")
    {
      return OE2;
    }
  else if (name == "NE2")
    {
      return NE2;
    }
  else if (name == "CE3")
    {
      return CE3;
    }
  else if (name == "CZ")
    {
      return CZ;
    }
  else if (name == "NZ")
    {
      return NZ;
    }
  else if (name == "CZ2")
    {
      return CZ2;
    }
  else if (name == "CZ3")
    {
      return CZ3;
    }
  else if (name == "OH")
    {
      return OH;
    }
  else if (name == "NH1")
    {
      return NH1;
    }
  else if (name == "NH2")
    {
      return NH2;
    }
  else if (name == "CH2")
    {
      return CH2;
    }
  else if (name == "OXT")
    {
      return OXT;
    }
  else if (name == "P")
    {
      return P;
    }
  else if (name == "O1P")
    {
      return O1P;
    }
  else if (name == "O2P")
    {
      return O2P;
    }
  else if (name == "O5S")
    {
      return O5S;
    }
  else if (name == "C5S")
    {
      return C5S;
    }
  else if (name == "C4S")
    {
      return C4S;
    }
  else if (name == "O4S")
    {
      return O4S;
    }
  else if (name == "C3S")
    {
      return C3S;
    }
  else if (name == "O3S")
    {
      return O3S;
    }
  else if (name == "C2S")
    {
      return C2S;
    }
  else if (name == "O2S")
    {
      return O2S;
    }
  else if (name == "C1S")
    {
      return C1S;
    }
  else if (name == "N9")
    {
      return N9;
    }
  else if (name == "C8")
    {
      return C8;
    }
  else if (name == "N7")
    {
      return N7;
    }
  else if (name == "C5")
    {
      return C5;
    }
  else if (name == "C6")
    {
      return C6;
    }
  else if (name == "O6")
    {
      return O6;
    }
  else if (name == "N6")
    {
      return N6;
    }
  else if (name == "N1")
    {
      return N1;
    }
  else if (name == "C2")
    {
      return C2;
    }
  else if (name == "O2")
    {
      return O2;
    }
  else if (name == "N2")
    {
      return N2;
    }
  else if (name == "N3")
    {
      return N3;
    }
  else if (name == "C4")
    {
      return C4;
    }
  else if (name == "O4")
    {
      return O4;
    }
  else if (name == "N4")
    {
      return N4;
    }
  else if (name == "C5M")
    {
      return C5M;
    }
  
  
  // Hydrogens
  else if ( (name == "H" ) || (name == "HN") )
    {
      return H;
    }
  else if ( (name == "H1" ) || (name == "HN1") ) // 1H ???
    {
      return H1;
    }
  else if ( (name == "H2" ) || (name == "HN2") ) // 2H ???
    {
      return H2;
    }
  else if ( (name == "H3") || (name == "HN3") ) // 3H ???
    {
      return H3;
    }
  else if (name == "HA")
    {
      return HA;
    }
   else if ( (name == "HA2") || (name == "1HA") ) // HA1 conflict old GLY
    {
      return HA2;
    }
   else if ( (name == "HA3") || (name == "2HA") ) // HA2 conflict old GLY
    {
      return HA3;
    }
  else if (name == "HB") 
    {
      return HB;
    }
  else if (name == "HB1" ) // 1HB conflict
    {
      return HB1;
    }
  else if (name == "HB2" ) // 2HB conflict
    {
      return HB2;
    }
  else if ( (name == "HB3") || (name == "3HB") )
    {
      return HB3;
    }
  else if (name == "HD1")
    {
      return HD1;
    }
  else if ( (name == "HD2") || (name == "1HD") )
    {
      return HD2;
    }
  else if ( (name == "HD3") || (name == "2HD") )
    {
      return HD3;
    }
  else if ( (name == "HD11") || (name == "1HD1") )
    {
      return HD11;
    }
  else if ( (name == "HD12") || (name == "2HD1") )
    {
      return HD12;
    }
  else if ( (name == "HD13") || (name == "3HD1") )
    {
      return HD13;
    }
  else if ( (name == "HD21") || (name == "1HD2") )
    {
      return HD21;
    }
  else if ( (name == "HD22") || (name == "2HD2") )
    {
      return HD22;
    }
  else if ( (name == "HD23") || (name == "3HD2") )
    {
      return HD23;
    }
  else if (name == "HE") // 1HE conflict
    {
      return HE;
    }
  else if (name == "HE1") // 1HE conflict
    {
      return HE1;
    }
  else if (name == "HE2") // 1HE, 2HE conflict
    {
      return HE2;
    }
  else if ( (name == "HE3") || (name == "3HE") ) // 2HE conflict
    {
      return HE3;
    }
  else if ((name == "HE21") || (name == "1HE2") )  
    {
      return HE21;
    }
  else if ((name == "HE22") || (name == "2HE2") )  
    {
      return HE22;
    }
  else if ( name == "HG" )
    {
      return HG;
    }
  else if (name == "HG1")
    {
      return HG1;
    }
  else if ( (name == "1HG") || (name == "HG2") )
    {
      return HG2;
    }
  else if (name == "HG11") // 1HG1 conflict
    {
      return HG11;
    }
  else if (name == "HG12")  // 2HG1 conflict
    {
      return HG12;
    }
  else if ( (name == "3HG1") || (name == "HG13") ) // 2HG1 conflict
    {
      return HG13;
    }
  else if ( (name == "HG21") || (name == "1HG2") )
    {
      return HG21;
    }
  else if ( (name == "HG22") || (name == "2HG2") )
    {
      return HG22;
    }
  else if ((name == "HG23") || (name == "3HG2") )
    {
      return HG23;
    }
  else if ((name == "HG3") || (name == "2HG"))
    {
      return HG3;
    }
  else if (name == "HH")
    {
      return HH;
    }
  else if (name == "HH2")
    {
      return HH2;
    }
  else if ( (name == "HH11") || (name == "1HH1") )
    {
      return HH11;
    }
  else if ( (name == "HH12") || (name == "2HH1") )
    {
      return HH12;
    }
  else if ( (name == "HH21") || (name == "1HH2") )
    {
      return HH21;
    }
  else if ( (name == "HH22") || (name == "2HH2") )
    {
      return HH22;
    }
  else if (name == "HZ")
    {
      return HZ;
    }
  else if ( (name == "HZ1") || (name == "1HZ"))
    {
      return HZ1;
    }
  else if ( (name == "HZ2") || (name == "2HZ"))
    {
      return HZ2;
    }
  else if ( (name == "HZ3") || (name == "3HZ"))
    {
      return HZ3;
    }

  
  
  // side chain pseudo atoms:
  else if (name == "XA")
    {
      return XA;
    }
  else if (name == "XD")
    {
      return XD;
    }
  else if (name == "XR")
    {
      return XR;
    }
  else if (name == "XN")
    {
      return XN;
    }
  else if (name == "XC")
    {
      return XC;
    }
  else if (name == "XQ")
    {
      return XQ;
    }
  else if (name == "XE")
    {
      return XE;
    }
  else if (name == "XG")
    {
      return XG;
    }
  else if (name == "XG")
    {
      return XG;
    }
  else if (name == "XH")
    {
      return XH;
    }
  else if (name == "XI")
    {
      return XI;
    }
  else if (name == "XL")
    {
      return XL;
    }
  else if (name == "XK")
    {
      return XK;
    }
  else if (name == "XM")
    {
      return XM;
    }
  else if (name == "XF")
    {
      return XF;
    }
  else if (name == "XP")
    {
      return XP;
    }
  else if (name == "XS")
    {
      return XS;
    }
  else if (name == "XT")
    {
      return XT;
    }
  else if (name == "XW")
    {
      return XW;
    }
  else if (name == "XY")
    {
      return XY;
    }
  else if (name == "XV")
    {
      return XV;
    }

  return X; // atom type name is unknown
}

inline string AtomTranslator(AtomCode code) 
{
  switch (code)
    {
    case X:
      return "X";
    case N:
      return "N";
    case CA:
      return "CA";
    case C:
      return "C";
    case O:
      return "O";
    case CB:
      return "CB";
    case SG:
      return "SG";
    case OG:
      return "OG";
    case CG:
      return "CG";
    case OG1:
      return "OG1";
    case CG1:
      return "CG1";
    case CG2:
      return "CG2";
    case CD:
      return "CD";
    case OD:
      return "OD";
    case SD:
      return "SD";
    case CD1:
      return "CD1";
    case OD1:
      return "OD1";
    case ND1:
      return "ND1";
    case CD2:
      return "CD2";
    case OD2:
      return "OD2";
    case ND2:
      return "ND2";
    case CE:
      return "CE";
    case NE:
      return "NE";
    case CE1:
      return "CE1";
    case OE1:
      return "OE1";
    case NE1:
      return "NE1";
    case CE2:
      return "CE2";
    case OE2:
      return "OE2";
    case NE2:
      return "NE2";
    case CE3:
      return "CE3";
    case CZ:
      return "CZ";
    case NZ:
      return "NZ";
    case CZ2:
      return "CZ2";
    case CZ3:
      return "CZ3";
    case OH:
      return "OH";
    case NH1:
      return "NH1";
    case NH2:
      return "NH2";
    case CH2:
      return "CH2";
    case OXT:
      return "OXT";
    case P:
      return "P";
    case O1P:
      return "O1P";
    case O2P:
      return "O2P";
    case O5S:
      return "O5S";
    case C5S:
      return "C5S";
    case C4S:
      return "C4S";
    case O4S:
      return "O4S";
    case C3S:
      return "C3S";
    case O3S:
      return "O3S";
    case C2S:
      return "C2S";
    case O2S:
      return "O2S";
    case C1S:
      return "C1S";
    case N9:
      return "N9";
    case C8:
      return "C8";
    case N7:
      return "N7";
    case C5:
      return "C5";
    case C6:
      return "C6";
    case O6:
      return "O6";
    case N6:
      return "N6";
    case N1:
      return "N1";
    case C2:
      return "C2";
    case O2:
      return "O2";
    case N2:
      return "N2";
    case N3:
      return "N3";
    case C4:
      return "C4";
    case O4:
      return "O4";
    case N4:
      return "N4";
    case C5M:
      return "C5M";
    
    // Hydrogens by Damiano Piovesan 2014
    case H:
      return "H";
    case H1:
      return "H1";
    case H2:
      return "H2";
    case H3:
      return "H3";
    case HA:
      return "HA";
    case HA2:
      return "HA2";
    case HA3:
      return "HA3";
    case HB:
      return "HB";
    case HB1:
      return "HB1";
    case HB2:
      return "HB2";
    case HB3:
      return "HB3";
    case HD1:
      return "HD1";
    case HD2:
      return "HD2";
    case HD3:
      return "HD3";
    case HD11:
      return "HD11";
    case HD12:
      return "HD12";
    case HD13:
      return "HD13";
    case HD21:
      return "HD21";
    case HD22:
      return "HD22";
    case HD23:
      return "HD23";
    case HE:
      return "HE";
    case HE1:
      return "HE1";
    case HE2:
      return "HE2";
    case HE3:
      return "HE3";
    case HE21:
      return "HE21";
    case HE22:
      return "HE22";
    case HG:
      return "HG";
    case HG1:
      return "HG1";
    case HG2:
      return "HG2";
    case HG11:
      return "HG11";
    case HG12:
      return "HG12";
    case HG13:
      return "HG13";
    case HG21:
      return "HG21";
    case HG22:
      return "HG22";
    case HG23:
      return "HG23";
    case HG3:
      return "HG3";
    case HH:
      return "HH";
    case HH2:
      return "HH2";
    case HH11:
      return "HH11";
    case HH12:
      return "HH12";
    case HH21:
      return "HH21";
    case HH22:
      return "HH22";
    case HZ:
      return "HZ";
    case HZ1:
      return "HZ1";
    case HZ2:
      return "HZ2";
    case HZ3:
      return "HZ3";
      
    // side chain pseudo atoms:
    case XA:
      return "XA";
    case XR:
      return "XR";
    case XN:
      return "XN";
    case XD:
      return "XD";
    case XC:
      return "XC";
    case XQ:
      return "XQ";
    case XE:
      return "XE";
    case XG:
      return "XG";
    case XH:
      return "XH";
    case XI:
      return "XI";
    case XL:
      return "XL";
    case XK:
      return "XK";
    case XM:
      return "XM";
    case XF:
      return "XF";
    case XP:
      return "XP";
    case XS:
      return "XS";
    case XT:
      return "XT";
    case XW:
      return "XW";
    case XY:
      return "XY";
    case XV:
      return "XV";


    case ATOM_CODE_SIZE: 
      ERROR("AtomTranslator(AtomCode code): unknown code", exception);
    }

  ERROR("AtomTranslator(AtomCode code): unknown code", exception);

  return "X";
}

//********

inline ostream& operator<<(ostream& os, const AtomCode& rval)
{
  os << AtomTranslator(rval);
  return os;
}

inline istream& operator>>(istream& is, AtomCode& rval)
{
  string name;
  is >> name;
  rval = AtomTranslator(name);
  return is;
}

//*******
// imported from Qmean (aa_map.cpp)
//*******
inline int get_all_atom_bin(string group_name) {
  switch (group_name[0])
    {
    case 'A':
      if (group_name == "AN") return 1;
      if (group_name == "ACA") return 2;
      if (group_name == "AC") return 3;
      if (group_name == "AO") return 4;
      if (group_name == "ACB") return 5;
      break;
    case 'C':
      if (group_name == "CN") return 6;
      if (group_name == "CCA") return 7;
      if (group_name == "CC") return 8;
      if (group_name == "CO") return 9;
      if (group_name == "CCB") return 10;
      if (group_name == "CSG") return 11;
      break;
    case 'D':
      if (group_name == "DN") return 12;
      if (group_name == "DCA") return 13;
      if (group_name == "DC") return 14;
      if (group_name == "DO") return 15;
      if (group_name == "DCB") return 16;
      if (group_name == "DCG") return 17;
      if (group_name == "DOD1") return 18;
      if (group_name == "DOD2") return 19;
      break;
    case 'E':
      if (group_name == "EN") return 20;
      if (group_name == "ECA") return 21;
      if (group_name == "EC") return 22;
      if (group_name == "EO") return 23;
      if (group_name == "ECB") return 24;
      if (group_name == "ECG") return 25;
      if (group_name == "ECD") return 26;
      if (group_name == "EOE1") return 27;
      if (group_name == "EOE2") return 28;
      break;
    case 'F':
      if (group_name == "FN") return 29;
      if (group_name == "FCA") return 30;
      if (group_name == "FC") return 31;
      if (group_name == "FO") return 32;
      if (group_name == "FCB") return 33;
      if (group_name == "FCG") return 34;
      if (group_name == "FCD1") return 35;
      if (group_name == "FCD2") return 36;
      if (group_name == "FCE1") return 37;
      if (group_name == "FCE2") return 38;
      if (group_name == "FCZ") return 39;
      break;
    case 'G':
      if (group_name == "GN") return 40;
      if (group_name == "GCA") return 41;
      if (group_name == "GC") return 42;
      if (group_name == "GO") return 43;
      break;
    case 'H':
      if (group_name == "HN") return 44;
      if (group_name == "HCA") return 45;
      if (group_name == "HC") return 46;
      if (group_name == "HO") return 47;
      if (group_name == "HCB") return 48;
      if (group_name == "HCG") return 49;
      if (group_name == "HND1") return 50;
      if (group_name == "HCD2") return 51;
      if (group_name == "HCE1") return 52;
      if (group_name == "HNE2") return 53;
      break;
    case 'I':
      if (group_name == "IN") return 54;
      if (group_name == "ICA") return 55;
      if (group_name == "IC") return 56;
      if (group_name == "IO") return 57;
      if (group_name == "ICB") return 58;
      if (group_name == "ICG1") return 59;
      if (group_name == "ICG2") return 60;
      if (group_name == "ICD1") return 61;
      break;
    case 'K':
      if (group_name == "KN") return 62;
      if (group_name == "KCA") return 63;
      if (group_name == "KC") return 64;
      if (group_name == "KO") return 65;
      if (group_name == "KCB") return 66;
      if (group_name == "KCG") return 67;
      if (group_name == "KCD") return 68;
      if (group_name == "KCE") return 69;
      if (group_name == "KNZ") return 70;
      break;
    case 'L':
      if (group_name == "LN") return 71;
      if (group_name == "LCA") return 72;
      if (group_name == "LC") return 73;
      if (group_name == "LO") return 74;
      if (group_name == "LCB") return 75;
      if (group_name == "LCG") return 76;
      if (group_name == "LCD1") return 77;
      if (group_name == "LCD2") return 78;
      break;
    case 'M':
      if (group_name == "MN") return 79;
      if (group_name == "MCA") return 80;
      if (group_name == "MC") return 81;
      if (group_name == "MO") return 82;
      if (group_name == "MCB") return 83;
      if (group_name == "MCG") return 84;
      if (group_name == "MSD") return 85;
      if (group_name == "MCE") return 86;
      break;
    case 'N':
      if (group_name == "NN") return 87;
      if (group_name == "NCA") return 88;
      if (group_name == "NC") return 89;
      if (group_name == "NO") return 90;
      if (group_name == "NCB") return 91;
      if (group_name == "NCG") return 92;
      if (group_name == "NOD1") return 93;
      if (group_name == "NND2") return 94;
      break;
    case 'P':
      if (group_name == "PN") return 95;
      if (group_name == "PCA") return 96;
      if (group_name == "PC") return 97;
      if (group_name == "PO") return 98;
      if (group_name == "PCB") return 99;
      if (group_name == "PCG") return 100;
      if (group_name == "PCD") return 101;
      break;
    case 'Q':
      if (group_name == "QN") return 102;
      if (group_name == "QCA") return 103;
      if (group_name == "QC") return 104;
      if (group_name == "QO") return 105;
      if (group_name == "QCB") return 106;
      if (group_name == "QCG") return 107;
      if (group_name == "QCD") return 108;
      if (group_name == "QOE1") return 109;
      if (group_name == "QNE2") return 110;
      break;
    case 'R':
      if (group_name == "RN") return 111;
      if (group_name == "RCA") return 112;
      if (group_name == "RC") return 113;
      if (group_name == "RO") return 114;
      if (group_name == "RCB") return 115;
      if (group_name == "RCG") return 116;
      if (group_name == "RCD") return 117;
      if (group_name == "RNE") return 118;
      if (group_name == "RCZ") return 119;
      if (group_name == "RNH1") return 120;
      if (group_name == "RNH2") return 121;
      break;
    case 'S':
      if (group_name == "SN") return 122;
      if (group_name == "SCA") return 123;
      if (group_name == "SC") return 124;
      if (group_name == "SO") return 125;
      if (group_name == "SCB") return 126;
      if (group_name == "SOG") return 127;
      break;
    case 'T':
      if (group_name == "TN") return 128;
      if (group_name == "TCA") return 129;
      if (group_name == "TC") return 130;
      if (group_name == "TO") return 131;
      if (group_name == "TCB") return 132;
      if (group_name == "TOG1") return 133;
      if (group_name == "TCG2") return 134;
      break;
    case 'V':
      if (group_name == "VN") return 135;
      if (group_name == "VCA") return 136;
      if (group_name == "VC") return 137;
      if (group_name == "VO") return 138;
      if (group_name == "VCB") return 139;
      if (group_name == "VCG1") return 140;
      if (group_name == "VCG2") return 141;
      break;
    case 'W':
      if (group_name == "WN") return 142;
      if (group_name == "WCA") return 143;
      if (group_name == "WC") return 144;
      if (group_name == "WO") return 145;
      if (group_name == "WCB") return 146;
      if (group_name == "WCG") return 147;
      if (group_name == "WCD1") return 148;
      if (group_name == "WCD2") return 149;
      if (group_name == "WNE1") return 150;
      if (group_name == "WCE2") return 151;
      if (group_name == "WCE3") return 152;
      if (group_name == "WCZ2") return 153;
      if (group_name == "WCZ3") return 154;
      if (group_name == "WCH2") return 155;
      break;
    case 'Y':
      if (group_name == "YN") return 156;
      if (group_name == "YCA") return 157;
      if (group_name == "YC") return 158;
      if (group_name == "YO") return 159;
      if (group_name == "YCB") return 160;
      if (group_name == "YCG") return 161;
      if (group_name == "YCD1") return 162;
      if (group_name == "YCD2") return 163;
      if (group_name == "YCE1") return 164;
      if (group_name == "YCE2") return 165;
      if (group_name == "YCZ") return 166;
      if (group_name == "YOH") return 167;
      break;
    }

	return 0;	// e.g. if hydrogen
}
#endif /* __ATOM_CODE_H__ */

