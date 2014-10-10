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


#ifndef __Nucleotide_CODE_H__
#define __Nucleotide_CODE_H__

// Includes

#include <string>
#include <Debug.h>


/**
 *      Translator: PDB names to internal one-word-code and 
 *                  vice versa. Provides some simple predicates dealing with 
 *                  one-word-code. 
 * Copyright:       This file contains information from the Bioinformatics Template Library (BTL).
 *
 *                  Copyright (C) 1997,1998 Birkbeck College, Malet Street, 
 *                  London WC1E 7HX, U.K. (classlib@mail.cryst.bbk.ac.uk)
 */

enum NucleotideCode {
    ADENINE = 0,
    THYMINE,
    CYTOSINE,
    GUANINE,
    URACIL,
    XX, //  Corresponds to unknown Nucleotide type 
    Nucleotide_CODE_SIZE // number of Nucleotide types
};

/**
 *   Return true, if Nucleotide type name is known.
 *  The magic code XXX corresponds to an unknown Nucleotide code
 */
inline
bool
isKnownNucleotide(NucleotideCode code) {
    return !(code == XX);
}

/** 
 *   Translate string into Nucleotide code enum. 
 */
inline
NucleotideCode
nucleotideOneLetterTranslator(char name) {
    DUMP(name);
    if (name == ' ') {
        return XX;
    }
    if (name == 'X') {
        return XX;
    }
    if (name == 'A') {
        return ADENINE;
    } else if (name == 'T') {
        return THYMINE;
    } else if (name == 'C') {
        return CYTOSINE;
    } else if (name == 'G') {
        return GUANINE;
    } else if (name == 'U') {
        return URACIL;
    }

    return XX; // Nucleotide type name is unknown
}

/** 
 *   Translate string into Nucleotide code enum. 
 */
inline
NucleotideCode
nucleotideThreeLetterTranslator(const string& name) {
    DUMP(name);
    if (name == "") {
        return XX;
    }
    if (name == "X") {
        return XX;
    }
    if (name == "A") {
        return ADENINE;
    } else if (name == "T") {
        return THYMINE;
    } else if (name == "C") {
        return CYTOSINE;
    } else if (name == "G") {
        return GUANINE;
    } else if (name == "U") {
        return URACIL;
    } else if (name == "DA") {
        return ADENINE;
    } else if (name == "DT") {
        return THYMINE;
    } else if (name == "DC") {
        return CYTOSINE;
    } else if (name == "DG") {
        return GUANINE;
    } else if (name == "DU") {
        return URACIL;
    }

    return XX; // Nucleotide type name is unknown
}

/**
 *  Returns the corresponding three letter code as a string.
 *@param aminoacidCode 
 *@return string
 */
inline string nucleotideThreeLetterTranslator(NucleotideCode code) {
    switch (code) {
        case XX:
            return "X";
        case ADENINE:
            return "A";
        case THYMINE:
            return "T";
        case CYTOSINE:
            return "C";
        case GUANINE:
            return "G";
        case URACIL:
            return "U";
        case Nucleotide_CODE_SIZE:
            ERROR("NucleotideTranslator(NucleotideCode code): unknown code", exception);
    }

    ERROR("NucleotideTranslator(NucleotideCode code): unknown code", eXXXception);

    return "X";
}

/**
 *  Returns the corresponding three letter code as a string.
 *@param one letter code (char)
 *@return string
 */
inline string nucleotideOneLetter2ThreeLetter(char oneLetter) {
    NucleotideCode code;
    code = nucleotideOneLetterTranslator(oneLetter);
    return nucleotideThreeLetterTranslator(code);
}

/**
 *  Returns the corresponding one letter code.
 *@param three letter code (string)
 *@return char
 */
inline char nucleotideThreeLetter2OneLetter(const string& threeLetter) {
    NucleotideCode code;
    code = nucleotideThreeLetterTranslator(threeLetter);
    return nucleotideOneLetterTranslator(code);
}

/**
 *  Returns one aminoacid code.
 *@param aminoacid code 
 *@return aminoacid code
 */
inline NucleotideCode& operator++(NucleotideCode& ac, int) {
    return ac = ((ac == XX) ? ADENINE : NucleotideCode(ac + 1));
}

/**
 *  Returns one aminoacid code.
 *@param aminoacid code 
 *@return aminoacid code
 */
inline NucleotideCode& operator--(NucleotideCode& ac, int) {
    return ac = ((ac == ADENINE) ? XX : NucleotideCode(ac - 1));
}

#endif /* __Nucleotide_CODE_H__ */

