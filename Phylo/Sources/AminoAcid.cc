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

#include <AminoAcid.h>

namespace Victor {
namespace Phylo {
    typedef AminoAcid::Code Code;

    const char AminoAcid::letter1s[] = {
        'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M',
        'F', 'P', 'S', 'T', 'W', 'Y', 'V',
        'U', 'O',
        'B', 'Z', 'J', 'X'
    };
    const string AminoAcid::letter3s[] = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
        "TYR", "VAL",
        "SEC", "PYL",
        "ASX", "GLX", "XLE", "XAA"
    };
    const string AminoAcid::names[] = {
        "Alanine", "Arginine", "Asparagine", "Aspartic acid", "Cysteine",
        "Glutamic acid", "Glutamine", "Glycine", "Histeine", "Isoleucine",
        "Leucine", "Lysine", "Methionine", "Phenylalanine", "Proline",
        "Serine", "Threonine", "Tryptophan", "Tyrosine", "Valine",
        "Selenocysteine", "Pyrrolisine",
        "Asparagine / Aspartic acid", "Glutamine / Glutamic acid",
        "Leucine / Isoleucine", "Unspecified / Unknown"
    };



    Code
    AminoAcid::letter1ToCode(const char code) {
        switch (code) {
            case 'A': return A;
            case 'R': return R;
            case 'N': return N;
            case 'D': return D;
            case 'C': return C;
            case 'E': return E;
            case 'Q': return Q;
            case 'G': return G;
            case 'H': return H;
            case 'I': return I;
            case 'L': return L;
            case 'K': return K;
            case 'M': return M;
            case 'F': return F;
            case 'P': return P;
            case 'S': return S;
            case 'T': return T;
            case 'W': return W;
            case 'Y': return Y;
            case 'V': return V;

            case 'U': return U;
            case 'O': return O;

            case 'B': return B;
            case 'Z': return Z;
            case 'J': return J;
            default:  return X;
        }
    }



    Code
    AminoAcid::letter3ToCode(const string &code) {
        size_t index = 0;
        bool found   = false;

        while (!found && index < X) {
            found = code == letter3s[index];
            if (!found) {
                index++;
            }
        }

        return (Code) index;
    }



    AminoAcid::AminoAcid(const Code code)
    : code(code) {
    }



    AminoAcid::AminoAcid(const char code)
    : code(letter1ToCode(code)) {
    }



    AminoAcid::AminoAcid(const string &code)
    : code(letter3ToCode(code)) {
    }



    AminoAcid::AminoAcid(const AminoAcid &origin)
    : AminoAcid(origin.code) {
    }



    AminoAcid::~AminoAcid() {
    }



    Code
    AminoAcid::getCode() const {
        return code;
    }



    char
    AminoAcid::getSingleLetter() const {
        return letter1s[code];
    }



    string
    AminoAcid::getThreeLetter() const {
        return letter3s[code];
    }



    string
    AminoAcid::getFullName() const {
        return names[code];
    }

}  // namespace Phylo
}  // namespace Victor