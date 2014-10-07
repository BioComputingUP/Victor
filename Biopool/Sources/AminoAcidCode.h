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
 
 *@Description:     Translator: PDB names to internal one-word-code and 
 *                  vice versa. Provides some simple predicates dealing with 
 *                  one-word-code. 
 *@ Copyright:       This file contains information from the Bioinformatics 
 *                  Template Library (BTL).
 *
 *                  Copyright (C) 1997,1998 Birkbeck College, Malet Street, 
 *                  London WC1E 7HX, U.K. (classlib@mail.cryst.bbk.ac.uk)
 */

#ifndef __AminoAcid_CODE_H__
#define __AminoAcid_CODE_H__

// Includes

#include <string>
#include <Debug.h>

enum AminoAcidCode {
    ALA = 0,
    CYS,
    ASP,
    GLU,
    PHE,
    GLY,
    HIS,
    ILE,
    LYS,
    LEU,
    MET,
    ASN,
    PRO,
    GLN,
    ARG,
    SER,
    THR,
    VAL,
    TRP,
    TYR,
    XXX, //  Corresponds to unknown AminoAcid type 
    AminoAcid_CODE_SIZE // number of AminoAcid types
};

/**
   true, if AminoAcid type name is known 
   magic code XXX corresponds to an unknown AminoAcid code
 */
inline
bool
isKnownAminoAcid(AminoAcidCode code) {
    return !(code == XXX);
}

/** Translate string into AminoAcid code enum. */
inline
AminoAcidCode
aminoAcidOneLetterTranslator(char name) {
    DUMP(name);
    if (name == ' ') {
        return XXX;
    }
    if (name == 'X') {
        return XXX;
    }
    if (name == 'A') {
        return ALA;
    } else if (name == 'C') {
        return CYS;
    } else if (name == 'D') {
        return ASP;
    } else if (name == 'E') {
        return GLU;
    } else if (name == 'F') {
        return PHE;
    } else if (name == 'G') {
        return GLY;
    } else if (name == 'H') {
        return HIS;
    } else if (name == 'I') {
        return ILE;
    } else if (name == 'K') {
        return LYS;
    } else if (name == 'L') {
        return LEU;
    } else if (name == 'M') {
        return MET;
    } else if (name == 'N') {
        return ASN;
    } else if (name == 'P') {
        return PRO;
    } else if (name == 'Q') {
        return GLN;
    } else if (name == 'R') {
        return ARG;
    } else if (name == 'S') {
        return SER;
    } else if (name == 'T') {
        return THR;
    } else if (name == 'V') {
        return VAL;
    } else if (name == 'W') {
        return TRP;
    } else if (name == 'Y') {
        return TYR;
    }

    return XXX; // AminoAcid type name is unknown
}

/** Translate string into AminoAcid code enum. */
inline
AminoAcidCode
aminoAcidThreeLetterTranslator(const string& name) {
    DUMP(name);
    if (name == "") {
        return XXX;
    }
    if (name == "XXX") {
        return XXX;
    }
    if (name == "ALA") {
        return ALA;
    } else if (name == "CYS") {
        return CYS;
    } else if (name == "ASP") {
        return ASP;
    } else if (name == "GLU") {
        return GLU;
    } else if (name == "PHE") {
        return PHE;
    } else if (name == "GLY") {
        return GLY;
    } else if (name == "HIS") {
        return HIS;
    } else if (name == "ILE") {
        return ILE;
    } else if (name == "LYS") {
        return LYS;
    } else if (name == "LEU") {
        return LEU;
    } else if (name == "MET") {
        return MET;
    } else if (name == "ASN") {
        return ASN;
    } else if (name == "PRO") {
        return PRO;
    } else if (name == "GLN") {
        return GLN;
    } else if (name == "ARG") {
        return ARG;
    } else if (name == "SER") {
        return SER;
    } else if (name == "THR") {
        return THR;
    } else if (name == "VAL") {
        return VAL;
    } else if (name == "TRP") {
        return TRP;
    } else if (name == "TYR") {
        return TYR;
    }

    return XXX; // AminoAcid type name is unknown
}

/**
 *@Description Returns the corresponding three letter code as a string.
 *@param aminoacidCode 
 *@return string
 */
inline string aminoAcidThreeLetterTranslator(AminoAcidCode code) {
    switch (code) {
        case XXX:
            return "XXX";
        case ALA:
            return "ALA";
        case CYS:
            return "CYS";
        case ASP:
            return "ASP";
        case GLU:
            return "GLU";
        case PHE:
            return "PHE";
        case GLY:
            return "GLY";
        case HIS:
            return "HIS";
        case ILE:
            return "ILE";
        case LYS:
            return "LYS";
        case LEU:
            return "LEU";
        case MET:
            return "MET";
        case ASN:
            return "ASN";
        case PRO:
            return "PRO";
        case GLN:
            return "GLN";
        case ARG:
            return "ARG";
        case SER:
            return "SER";
        case THR:
            return "THR";
        case VAL:
            return "VAL";
        case TRP:
            return "TRP";
        case TYR:
            return "TYR";

        case AminoAcid_CODE_SIZE:
            ERROR("AminoAcidTranslator(AminoAcidCode code): unknown code", exception);
    }

    ERROR("AminoAcidTranslator(AminoAcidCode code): unknown code", eXXXception);

    return "X";
}

/**
 *@Description Returns the corresponding one letter code as a string.
 *@param aminoacidCode 
 *@return char
 */
inline char aminoAcidOneLetterTranslator(AminoAcidCode code) {
    switch (code) {
        case XXX:
            return 'X';
        case ALA:
            return 'A';
        case CYS:
            return 'C';
        case ASP:
            return 'D';
        case GLU:
            return 'E';
        case PHE:
            return 'F';
        case GLY:
            return 'G';
        case HIS:
            return 'H';
        case ILE:
            return 'I';
        case LYS:
            return 'K';
        case LEU:
            return 'L';
        case MET:
            return 'M';
        case ASN:
            return 'N';
        case PRO:
            return 'P';
        case GLN:
            return 'Q';
        case ARG:
            return 'R';
        case SER:
            return 'S';
        case THR:
            return 'T';
        case VAL:
            return 'V';
        case TRP:
            return 'W';
        case TYR:
            return 'Y';

        case AminoAcid_CODE_SIZE:
            ERROR("AminoAcidTranslator(AminoAcidCode code): unknown code", exception);
    }


    ERROR("AminoAcidTranslator(AminoAcidCode code): unknown code", exception);

    return 'X';
}

/**
 *@Description Returns the corresponding aminoacid name as a string.
 *@param aminoacidCode 
 *@return string
 */

inline string aminoAcidFullTranslator(AminoAcidCode code) {
    switch (code) {
        case XXX:
            return "Unknown";
        case ALA:
            return "Alanine";
        case CYS:
            return "Cysteine";
        case ASP:
            return "Aspartate";
        case GLU:
            return "Glutamate";
        case PHE:
            return "Phenylalanine";
        case GLY:
            return "Glycine";
        case HIS:
            return "Histidine";
        case ILE:
            return "Isoleucine";
        case LYS:
            return "Lysine";
        case LEU:
            return "Leucine";
        case MET:
            return "Methionine";
        case ASN:
            return "Asparagine";
        case PRO:
            return "Proline";
        case GLN:
            return "Glutamine";
        case ARG:
            return "Arginine";
        case SER:
            return "Serine";
        case THR:
            return "Threonine";
        case VAL:
            return "Valine";
        case TRP:
            return "Tryptophane";
        case TYR:
            return "Tyrosine";

        case AminoAcid_CODE_SIZE:
            ERROR("AminoAcidFullTranslator(AminoAcidCode code): unknown code", exception);
    }

    ERROR("AminoAcidFullTranslator(AminoAcidCode code): unknown code", exception);

    return "Unknown";
}

/**
 *@Description Returns the corresponding three letter code as a string.
 *@param one letter code (char)
 *@return string
 */
inline string oneLetter2ThreeLetter(char oneLetter) {
    AminoAcidCode code;
    code = aminoAcidOneLetterTranslator(oneLetter);
    return aminoAcidThreeLetterTranslator(code);
}

/**
 *@Description Returns the corresponding one letter code.
 *@param three letter code (string)
 *@return char
 */
inline char threeLetter2OneLetter(const string& threeLetter) {
    AminoAcidCode code;
    code = aminoAcidThreeLetterTranslator(threeLetter);
    return aminoAcidOneLetterTranslator(code);
}

/**
 *@Description Returns one aminoacid code.
 *@param aminoacid code 
 *@return aminoacid code
 */
inline AminoAcidCode& operator++(AminoAcidCode& ac, int) {
    return ac = ((ac == XXX) ? ALA : AminoAcidCode(ac + 1));
}

/**
 *@Description Returns one aminoacid code.
 *@param aminoacid code 
 *@return aminoacid code
 */
inline AminoAcidCode& operator--(AminoAcidCode& ac, int) {
    return ac = ((ac == ALA) ? XXX : AminoAcidCode(ac - 1));
}

/**
 *@Description verifyes if the aminoacid is a polar one.
 *@param aminoacid code 
 *@return boolean
 */
inline bool isPolar(AminoAcidCode code) {
    if ((code == CYS) ||
            (code == ASP) ||
            (code == GLU) ||
            (code == HIS) ||
            (code == LYS) ||
            (code == ASN) ||
            (code == GLN) ||
            (code == ARG) ||
            (code == SER) ||
            (code == THR) ||
            (code == TRP) ||
            (code == TYR))
        return true;

    return false;
}

/**
 *@Description verifyes if the aminoacid is a  Hydrophobic one.
 *@param aminoacid code 
 *@return boolean
 */
inline bool isHydrophobic(AminoAcidCode code) {
    if ((code == ALA) ||
            (code == PHE) ||
            (code == GLY) ||
            (code == ILE) ||
            (code == LEU) ||
            (code == MET) ||
            (code == PRO) ||
            (code == VAL))
        return true;

    return false;
}

enum StateCode {
    HELIX,
    STRAND,
    TURN,
    COIL,
    YYY,
    State_CODE_SIZE // number of State types
};

/**
 *@Description states  name for a translator.
 *@param translator  one letter code
 *@return statecode(char)
 */
inline StateCode stateCodeOneLetterTranslator(char name) {
    DUMP(name);
    if (name == ' ') {
        return YYY;
    }
    if (name == 'Y') {
        return YYY;
    }
    if (name == 'H') {
        return HELIX;
    } else if (name == 'E') {
        return STRAND;
    } else if (name == 'T') {
        return TURN;
    } else if (name == 'C') {
        return COIL;
    }

    return YYY; // StateCode type name is unknown
}

/**
 *@Description states the one letter name for a translator.
 *@param translator name
 *@return statecode
 */
inline char stateCodeOneLetterTranslator(StateCode code) {
    switch (code) {
        case YYY:
            return 'Y';
        case HELIX:
            return 'H';
        case STRAND:
            return 'E';
        case TURN:
            return 'T';
        case COIL:
            return 'C';

        case State_CODE_SIZE:
            ERROR("StateCodeTranslator(StateCode code): unknown code", exception);
    }


    ERROR("StateCodeTranslator(StateCode code): unknown code", exception);

    return 'X';
}

/**
 *@Description Returns one state code.
 *@param statecode code 
 *@return statecode code
 */
inline StateCode& operator++(StateCode& sc, int) {
    return sc = ((sc == YYY) ? HELIX : StateCode(sc + 1));
}

#endif /* __AminoAcid_CODE_H__ */

