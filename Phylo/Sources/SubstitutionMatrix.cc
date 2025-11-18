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

#include <SubstitutionMatrix.h>
#include <iostream>
#include <sstream>
#include <string.h>
#include <AminoAcid.h>
#include <precalculated_matrices.h>

namespace Victor {
namespace Phylo {
    typedef SubstitutionMatrix::Identifier Identifier;

    /** @brief Returns a precalculate matrix as a map.
     * @param[in]           identifier Identifier of the matrix
     * @return Precalculated matrix
     */
    static map<char, map<char, int>> &precalculated_matrices(const Identifier identifier);
    
    const string SubstitutionMatrix::residues = "ARNDCEQGHILKMFPSTWYVUOBZJX";


    SubstitutionMatrix::SubstitutionMatrix(const map<char, map<char, int>> &scores) {
        const size_t size = 26;
        int new_scores[26 * 26],
            max_score = 0,
            min_score = 0;
        double avg_mismatch_score = 0.0;
        map<char, map<char, int>> scores_map = scores;

        for (size_t i = 0; i < size; i++) {
            const char r1 = residues[i];
            for (size_t j = 0; j < size; j++) {
                const char r2 = residues[j];
                const int score = scores_map[r1][r2];

                new_scores[hash(r1) * size + hash(r2)]  = score;
                min_score           = std::min(score, min_score);
                max_score           = std::max(score, max_score);
                avg_mismatch_score += score;
            }
        }
        avg_mismatch_score /= (size * size);

        memcpy(&this->scores, &new_scores, size * size * sizeof(int));
        this->min_score          = min_score;
        this->max_score          = max_score;
        this->avg_mismatch_score = avg_mismatch_score;
    }



    SubstitutionMatrix::SubstitutionMatrix(const Identifier identifier)
    : SubstitutionMatrix(precalculated_matrices(identifier)) {
    }



    SubstitutionMatrix::SubstitutionMatrix(const SubMatrix &matrix)
    : SubstitutionMatrix(SubMatrixToMap(matrix)) {
    }



    SubstitutionMatrix::SubstitutionMatrix(const SubstitutionMatrix &matrix)
    : min_score(matrix.min_score), max_score(matrix.max_score),
    avg_mismatch_score(matrix.avg_mismatch_score) {
        memcpy(&scores, &matrix.scores, matrix.getScoresNumber() * sizeof(int));
    }



    SubstitutionMatrix::~SubstitutionMatrix() {
    }



    bool
    SubstitutionMatrix::operator==(const SubstitutionMatrix &other) const {
        return isEqual(other);
    }



    int
    SubstitutionMatrix::operator()(const char a, const char b) const {
        return getScore(a, b);
    }



    SubstitutionMatrix &
    SubstitutionMatrix::operator+(const int value) const {
        return add(value);
    }



    SubstitutionMatrix &
    SubstitutionMatrix::operator-(const int value) const {
        return add(-value);
    }



    SubstitutionMatrix &
    SubstitutionMatrix::operator-() const {
        return multiply(-1);
    }



    SubstitutionMatrix &
    SubstitutionMatrix::operator*(const int value) const {
        return multiply(value);
    }



    bool
    SubstitutionMatrix::isEqual(const SubstitutionMatrix &other) const {
        const size_t size = getScoresNumber();
        bool is_equal = true;
        size_t i = 0;

        // If matrices differ in sizes, they are not equal
        if (size != other.getScoresNumber()) {
            return false;
        }

        // Checks whether every score is equal
        while (is_equal && i < size) {
            is_equal &= scores[i] == other.scores[i];
            i++;
        }

        return is_equal;
    }



    size_t
    SubstitutionMatrix::getSize() const {
        return getResidues().size();
    }



    size_t
    SubstitutionMatrix::getScoresNumber() const {
        return getSize() * getSize();
    }



    string
    SubstitutionMatrix::getResidues() const {
        return residues;
    }



    int
    SubstitutionMatrix::getScore(const char a, const char b) const {
        return scores[hash(a) * getSize() + hash(b)];
    }



    int
    SubstitutionMatrix::getMinScore() const {
        return min_score;
    }



    int
    SubstitutionMatrix::getMaxScore() const {
        return max_score;
    }



    double
    SubstitutionMatrix::getAvgMismatchScore() const {
        return avg_mismatch_score;
    }



    map<char, map<char, int>> &
    SubstitutionMatrix::getScoresAsMap() const {
        map<char, map<char, int>> *scores_map = new map<char, map<char, int>>();

        for (char r1 : residues) {
            for (char r2: residues) {
                (*scores_map)[r1][r2] = getScore(r1, r2);
            }
        }

        return *scores_map;
    }



    SubMatrix &
    SubstitutionMatrix::asSubMatrix() const {
        stringstream stream;

        // Saves matrix into a stream
        stream << getResidues() << "\n" << getSize() << "\n";
        for (char r1 : getResidues()) {
            stream << getSize();
            for (char r2 : getResidues()) {
                if (getScore(r1, r2) >= 0 && getScore(r1, r2) < 10) {
                    stream << " ";
                }
                stream << " " << getScore(r1, r2) << " ";
            }
            stream << "\n";
        }
        stream << "#" << "\n";

        // Build SubMatrix from the stream
        return *(new SubMatrix(stream));
    }



    SubstitutionMatrix &
    SubstitutionMatrix::add(const int value) const {
        return fma(value, 1);
    }



    SubstitutionMatrix &
    SubstitutionMatrix::multiply(const int value) const {
        return fma(0, value);
    }



    SubstitutionMatrix &
    SubstitutionMatrix::fma(const int addend, const int factor) const {
        map<char, map<char, int>> scores_map;

        for (char r1 : residues) {
            for (char r2: residues) {
                scores_map[r1][r2] = addend + factor * getScore(r1, r2);
            }
        }

        return *(new SubstitutionMatrix(scores_map));
    }



    Identifier
    SubstitutionMatrix::stringToIdentifier(const string &identifier) {
        if (identifier == "BLOSUM30") return BLOSUM30;
        if (identifier == "BLOSUM35") return BLOSUM35;
        if (identifier == "BLOSUM40") return BLOSUM40;
        if (identifier == "BLOSUM45") return BLOSUM45;
        if (identifier == "BLOSUM50") return BLOSUM50;
        if (identifier == "BLOSUM55") return BLOSUM55;
        if (identifier == "BLOSUM62") return BLOSUM62;
        if (identifier == "BLOSUM65") return BLOSUM65;
        if (identifier == "BLOSUM70") return BLOSUM70;
        if (identifier == "BLOSUM75") return BLOSUM75;
        if (identifier == "BLOSUM80") return BLOSUM80;
        if (identifier == "BLOSUM90") return BLOSUM90;

        if (identifier == "PAM20")  return PAM20;
        if (identifier == "PAM60")  return PAM60;
        if (identifier == "PAM120") return PAM120;
        if (identifier == "PAM160") return PAM160;
        if (identifier == "PAM250") return PAM250;
        if (identifier == "PAM350") return PAM350;

        if (identifier == "MD40")  return MD40;
        if (identifier == "MD120") return MD120;
        if (identifier == "MD250") return MD250;
        if (identifier == "MD350") return MD350;

        if (identifier == "IDENTITY") return IDENTITY;

        if (identifier == "GON40")  return GON40;
        if (identifier == "GON80")  return GON80;
        if (identifier == "GON120") return GON120;
        if (identifier == "GON160") return GON160;
        if (identifier == "GON250") return GON250;
        if (identifier == "GON300") return GON300;
        if (identifier == "GON350") return GON350;

        return IDENTITY;
    }



    size_t
    SubstitutionMatrix::hash(const char amino) {
        return (size_t) AminoAcid::letter1ToCode(amino);
    }



    map<char, map<char, int>> &
    SubstitutionMatrix::SubMatrixToMap(const SubMatrix &submatrix) {
        map<char, map<char, int>> *scores_map = new map<char, map<char, int>>();

        for (char r1 : residues) {
            for (char r2 : residues) {
                (*scores_map)[r1][r2] = submatrix.score[r1][r2];
            }
        }

        return *scores_map;
    }



////////////////////////////////////////////////////////////////////////
// Precalculated matrices
static map<char, map<char, int>> &precalculated_matrices(const Identifier identifier) {
    const string residues = _residues;;
    const size_t size     = residues.size();
    map<char, map<char, int>> *matrix = new map<char, map<char, int>>();
    int *ptr;

    // Finds correct matrix
    switch (identifier) {
        // BLOSUM family
        case SubstitutionMatrix::BLOSUM30: ptr = _blosum30; break;
        case SubstitutionMatrix::BLOSUM35: ptr = _blosum35; break;
        case SubstitutionMatrix::BLOSUM40: ptr = _blosum40; break;
        case SubstitutionMatrix::BLOSUM45: ptr = _blosum45; break;
        case SubstitutionMatrix::BLOSUM50: ptr = _blosum50; break;
        case SubstitutionMatrix::BLOSUM55: ptr = _blosum55; break;
        case SubstitutionMatrix::BLOSUM62: ptr = _blosum62; break;
        case SubstitutionMatrix::BLOSUM65: ptr = _blosum65; break;
        case SubstitutionMatrix::BLOSUM70: ptr = _blosum70; break;
        case SubstitutionMatrix::BLOSUM75: ptr = _blosum75; break;
        case SubstitutionMatrix::BLOSUM80: ptr = _blosum80; break;
        case SubstitutionMatrix::BLOSUM90: ptr = _blosum90; break;

        // PAM family
        case SubstitutionMatrix::PAM20:  ptr = _pam20;  break;
        case SubstitutionMatrix::PAM60:  ptr = _pam60;  break;
        case SubstitutionMatrix::PAM120: ptr = _pam120; break;
        case SubstitutionMatrix::PAM160: ptr = _pam160; break;
        case SubstitutionMatrix::PAM250: ptr = _pam250; break;
        case SubstitutionMatrix::PAM350: ptr = _pam350; break;

        // MD family
        case SubstitutionMatrix::MD40:  ptr = _md40;  break;
        case SubstitutionMatrix::MD120: ptr = _md120; break;
        case SubstitutionMatrix::MD250: ptr = _md250; break;
        case SubstitutionMatrix::MD350: ptr = _md350; break;

        // Identity
        case SubstitutionMatrix::IDENTITY: ptr = _idmat; break;

        // GON family
        case SubstitutionMatrix::GON40:  ptr = _gon40;  break;
        case SubstitutionMatrix::GON80:  ptr = _gon80;  break;
        case SubstitutionMatrix::GON120: ptr = _gon120; break;
        case SubstitutionMatrix::GON160: ptr = _gon160; break;
        case SubstitutionMatrix::GON250: ptr = _gon250; break;
        case SubstitutionMatrix::GON300: ptr = _gon300; break;
        case SubstitutionMatrix::GON350: ptr = _gon350; break;

        // Default: identity
        default: ptr = _idmat;
    }

    // Converts raw matrix to map
    for (size_t i = 0; i < size; i++) {
        const char r_i = residues[i];
        for (size_t j = 0; j <= i; j++) {
            const char r_j = residues[j];
            (*matrix)[r_i][r_j] = ptr[(i * i + i) / 2 + j];
            (*matrix)[r_j][r_i] = ptr[(i * i + i) / 2 + j];
        }
    }

    return *matrix;
}

}  // namespace Phylo
}  // namespace Victor