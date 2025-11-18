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

#ifndef _VICTOR_PHYLO_MULTIPLEALIGNMENTALGORITHM_H_
#define _VICTOR_PHYLO_MULTIPLEALIGNMENTALGORITHM_H_

#include <MultipleAlignment.h>
#include <AlignmentBase.h>

using namespace Victor::Align2;

namespace Victor {
namespace Phylo {

/** @brief Implements a multiple sequence alignment algorithm.
 *
 * Multiple sequence alignment (MSA) algorithms produce MSAs by considering
 * a set of sequences.
 *
 * MSA algorithms are implemented as functors from an alignment to a MSA.
 *
 * This class follows the Strategy and the Facade Design Patterns.
 * 
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 */
class MultipleAlignmentAlgorithm {
 public:
    /** @brief  Builds an MSA from an alignment.
     * @param[in]       alignment The alignment
     * @return Multiple sequence alignment
     */
    inline MultipleAlignment &operator()(const AlignmentBase &alignment) {
        return apply(alignment);
    }


    
    /** @brief  Builds an MSA from an alignment.
     * @param[in]       alignment The alignment
     * @return Multiple sequence alignment
     */
    virtual MultipleAlignment &apply(const AlignmentBase &alignment) const = 0;
};

}  // namespace Phylo
}  // namespace Victor

#endif  // _VICTOR_PHYLO_MULTIPLEALIGNMENTALGORITHM_H_