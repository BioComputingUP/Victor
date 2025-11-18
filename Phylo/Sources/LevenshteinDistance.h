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

#ifndef _VICTOR_PHYLO_LEVENSHTEINDISTANCE_H_
#define _VICTOR_PHYLO_LEVENSHTEINDISTANCE_H_

#include <DistanceMatrixBuilder.h>

namespace Victor {
namespace Phylo {

/** @brief Builds a distance matrix by considering Levenshtein distance.
 *
 * Compute pairwise alignment for every pair of sequences, then creates
 * the distance matrix using Levenshtein distance.
 *
 * This class follows the Strategy Design Pattern.
 * 
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 */
class LevenshteinDistance: public DistanceMatrixBuilder {
 public:
    /** @brief Destructor. */
    virtual ~LevenshteinDistance();



 private:
    /** @brief Returns distance between two sequences.
     * @param[in]       seq1      First sequence
     * @param[in]       seq2      Second sequence
     * @return Distance between seq1 and seq2
     */
    virtual double computeDistance(const string &seq1, const string &seq2) const;
};

}  // namespace Phylo
}  // namespace Victor

#endif  // _VICTOR_PHYLO_LEVENSHTEINDISTANCE_H_