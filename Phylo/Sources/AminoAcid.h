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

#ifndef _VICTOR_PHYLO_AMINOACID_H_
#define _VICTOR_PHYLO_AMINOACID_H_

/** Number of amino acids. */
#define _AMINOACID_SIZE 26
#include <string>

using std::string;

namespace Victor{
namespace Phylo {

/** @brief Implements an amino acid. 
 * 
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 */
class AminoAcid {
 public:
    /** Hashed values for known amino acids. */
    typedef enum {
        A,    ///< Alanine
        R,    ///< Arginine
        N,    ///< Asparagine
        D,    ///< Aspartic acid
        C,    ///< Cysteine
        E,    ///< Glutamic acid
        Q,    ///< Glutamine
        G,    ///< Glycine
        H,    ///< Histidine
        I,    ///< Isoleucine
        L,    ///< Leucine
        K,    ///< Lysine
        M,    ///< Methionine
        F,    ///< Phenylalanine
        P,    ///< Proline
        S,    ///< Serine
        T,    ///< Threonine
        W,    ///< Tryptophan
        Y,    ///< Tyrosine
        V,    ///< Valine

        U,    ///< Selenocysteine
        O,    ///< Pyrrolisine

        B,    ///< Asparagine or aspartic acid
        Z,    ///< Glutamine or glutamic acid
        J,    ///< Leucine or isoleucine
        X     ///< Unspecified or unknown
    } Code;



    /** @brief Returns code of given amino acid.
     * @param[in]       code      Single letter code
     * @return Hash value
     */
    static Code letter1ToCode(const char code);



    /** @brief Returns code of given amino acid.
     * @param[in]       code      Three letter code
     * @return Hash value
     */
    static Code letter3ToCode(const string &code);



    /** @brief Builds an amino acid of given code.
     * @param[in]       code      Code of the amino acid
     */
    explicit AminoAcid(const Code code);

    /** @brief Builds an amino acid of given code.
     * @param[in]       code      Single letter code
     */
    explicit AminoAcid(const char code);

    /** @brief Builds an amino acid of given code.
     * @param[in]       code      Three letter code
     */
    explicit AminoAcid(const string &code);

    /** @brief Copy constructor.
     * @param[in]       origin    Original object
     */
    explicit AminoAcid(const AminoAcid &origin);

    /** @brief Destructor. */
    virtual ~AminoAcid();



    /** @brief Returns code of this amino acid.
     * @return Code of this amino acid
     */
    Code getCode() const;

    /** @brief Returns single letter code of this amino acid.
     * @return Single letter code
     */
    char getSingleLetter() const;

    /** @brief Returns thre letters code of this amino acid.
     * @return Three letters code
     */
    string getThreeLetter() const;

    /** @brief Returns fullname of this amino acid.
     * @return Full name of this amino acid
     */
    string getFullName() const;



 private:
    static const char   letter1s[_AMINOACID_SIZE]; ///< 1-letter codes of amino acids
    static const string letter3s[_AMINOACID_SIZE]; ///< 3-letters code of amino acids
    static const string names[_AMINOACID_SIZE];    ///< Full names of amino acids

    const Code code;    ///< Code of this amino acid
};

}
}

#endif  // _VICTOR_PHYLO_AMINOACID_H_