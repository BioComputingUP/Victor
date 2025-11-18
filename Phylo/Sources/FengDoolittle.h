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

#ifndef _VICTOR_PHYLO_FENGDOOLITTLE_H_
#define _VICTOR_PHYLO_FENGDOOLITTLE_H_

#define _GLIBCXX_IOMANIP
#include <FengDoolittleDistance.h>
#include <DistanceMatrix.h>
#include <FitchMargoliash.h>
#include <RootedTree.h>
#include <Sequence.h>
#include <MultipleAlignmentAlgorithm.h>

#include <AGPFunction.h>
#include <AlignmentBase.h>
#include <Alignment.h>
#include <AlignmentData.h>
#include <SequenceData.h>
#include <Structure.h>
#include <NWAlign.h>
#include <ScoringS2S.h>

using namespace Victor::Align2;

namespace Victor {
namespace Phylo {

/** @brief Implements the Feng and Doolittle multiple sequence
 * alignment algorithm.
 * 
 * Implements Feng-Doolittle as a functor.
 * Feng-Doolittle is a progressive multiple sequence alignment program.
 * For more information see:
 * - Progressive Sequence Alignment as a Prerequisite to Correct
 *   Phylogenetic Trees
 *
 * This class follows the Strategy and the Facade Design Patterns.
 * 
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 */
class FengDoolittle: public MultipleAlignmentAlgorithm {
 public:
    /** @brief Builds a Feng-Doolittle aligner with given configuration.
     * @param[in]       gap_open  Gap opening penalty (default: 10.0)
     * @param[in]       gap_extension  Gap extension penalty (default: 0.1)
     */
    explicit FengDoolittle(
        const double gap_open = 10.0,
        const double gap_extension = 0.1);

    /** @brief Copy constructor.
     * @param[in]       origin    Original Feng-Doolittle aligner
     */
    explicit FengDoolittle(const FengDoolittle &origin);

    /** @brief Destructor. */
    ~FengDoolittle();



    /** @brief Builds an MSA from a set of sequences.
     * @param[in]       alignment Set of sequences
     * @return Multiple sequence alignment
     */
    virtual MultipleAlignment &apply(const AlignmentBase &alignment) const;



 private:
    typedef map<string, string> Sequences;  ///< Set of sequences

    static SubMatrix matrix;                ///< Substitution matix
    static const FitchMargoliash build_tree; ///< Clustering algorithm

    AGPFunction &gap_function;              ///< Gap function


    /** @brief Aligns next nodes in the guide tree.
     * @param[in]       sequences Sequences to align
     * @param[in]       node      Next node of the guide tree
     * @return Multiple sequence alignment
     */
    MultipleAlignment *align(
        Sequences &sequences,
        const RootedTree &node) const;
};

}  // namespace Phylo
}  // namespace Victor

#endif  // _VICTOR_PHYLO_FENGDOOLITTLE_H_