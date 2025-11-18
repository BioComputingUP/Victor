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

#ifndef _VICTOR_PHYLO_CLUSTALW_H_
#define _VICTOR_PHYLO_CLUSTALW_H_

#define _GLIBCXX_IOMANIP
#include <DistanceMatrix.h>
#include <ClusteringAlgorithm.h>
#include <PhylogeneticTree.h>
#include <RootedTree.h>
#include <DistanceMatrixBuilder.h>
#include <Sequence.h>
#include <SubstitutionMatrix.h>
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

/** @brief Implements the ClustalW algorithm.
 *
 * Implements ClustalW as a functor.
 * ClustalW is a progressive multiple sequence alignment program. It is
 * an improved version of the original Clustal which uses sequence
 * weighting, position-specific gap penalties and weight matrix choice.
 *
 * For more information see:
 * - CLUSTAL W: improving the sensitivity of progressive multiple sequence
 *   alignment through sequence weighting, position-specific gap penalties
 *   and weight matrix choice
 *
 * This class follows the Strategy and Facade Design Patterns.
 *
 * Please note that neither the original paper nor other available
 * resources provide a good/serious description of the (many) implementation
 * details. Moreover, the freely available source code of Clustal W is
 * far too complex and too poorly documented. For these reason, this
 * work (had to) follow a reverse-engeneering approach.
 * 
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 */
class ClustalW: public MultipleAlignmentAlgorithm {
 public:
    /** Available weight matrix families. */
    typedef enum {
        PAM,            ///< PAM series
        BLOSUM          ///< BLOSUM series
    } WeightMatrix;



    /** @brief Builds a ClustalW aligner configured with given paramters.
     * @param[in]       distance_matrix_builder  Distance matrix builder
     * @param[in]       clustering_algorithm     Clustering algorithm
     * @param[in]       weight_matrix            Weight matrix family (progressive alignment)
     * @param[in]       gap_open                 Gap open penalty (progressive alignment)
     * @param[in]       gap_extension            Gap extension penalty (progressive alignment)
     */
    explicit ClustalW(
        const DistanceMatrixBuilder &distance_matrix_builder,
        const ClusteringAlgorithm &clustering_algorithm,
        const WeightMatrix weight_matrix,
        const double gap_open = 10.0,
        const double gap_extension = 0.1);

    /** @brief Copy constructor.
     * @param[in]       origin    Original object
     */
    explicit ClustalW(const ClustalW &origin);

    /** @brief Destructor. */
    virtual ~ClustalW();



    /** @brief  Builds an MSA from an alignment.
     * @param[in]       alignment The alignment
     * @return Multiple sequence alignment
     */
    virtual MultipleAlignment &apply(const AlignmentBase &alignment) const;



 private:
    typedef map<string, double> Weights;    ///< Weight map
    typedef map<string, string> Sequences;  ///< Set of sequences

    /** Residue-specific gap modification factors (Pascarella and Argos). */
    static map<char, double> gap_modifiers;

    const DistanceMatrixBuilder &build_matrix;   ///< Distance matrix builder
    const ClusteringAlgorithm &build_tree;  ///< Phylogenetic tree builder
    const WeightMatrix weight_matrix;       ///< Weight matrix family
    const double gap_open;                  ///< Gap open penalty
    const double gap_extension;             ///< Gap extension penalty


    /** @brief Returns weights.
     *
     * Weights are normalized.
     * @param[in]       guide_tree Guide tree of the alignment
     * @return Weights of the sequences
     */
    Weights getWeights(const RootedTree &guide_tree) const;



    /** @brief Performs multiple sequence alignment.
     *
     * Order of alignments is determined by topology of guide tree.
     * @param[in]       sequences Sequences to align
     * @param[in]       node      Node on the guide tree
     * @param[in]       weights   Weights
     * @return Multiple sequence alignment
     */
    MultipleAlignment *align(
        Sequences &sequences,
        const RootedTree &node,
        Weights &weights) const;



    /** @brief Returns initial gap opening penalty.
     * 
     * Initial gap opening penalty is computed considering dependence on
     * the weight matrix, dependence on the similarity of the sequences
     * and dependences on the length of the sequences:
     * GOP = GOP + log(min(N, M)) * (average mismatch score) * (percentage
     * of identity)
     * where M, N are the length of the sequences.
     * @param[in]       A         First (group of) sequence
     * @param[in]       B         Second (group of) sequence
     * @param[in]       matrix    Substitution matrix
     * @return Initial gap opening penalty
     */
    double getInitialGOP(
        const MultipleAlignment &A,
        const MultipleAlignment &B,
        const SubstitutionMatrix &matrix) const;

    /** @brief Returns initial gap extension penalty.
     *
     * Initial gap extension penalty is computed considering dependence
     * on the difference in the lengths of the sequences:
     * GEP = GEP * (1.0 + |log(N/M)|)
     * where M, N are the length of the sequences.
     * @param[in]       A         First (group of) sequence
     * @param[in]       B         Second (group of) sequence
     * @return Initial gap extension penalty
     */
    double getInitialGEP(
        const MultipleAlignment &A,
        const MultipleAlignment &B) const;

    /** @brief Returns position-specific gap open penalty.
     * 
     * Position-specific gap open penalty is computed considering
     * presence of gaps, nearby gaps and hydrophilic stretches.
     * @param[in]       A         First (group of) sequence
     * @param[in]       B         Second (group of) sequence
     * @param[in]       matrix    Substitution matrix
     * @param[in]       position  Current position
     * @return Position-specific gap open penalty
     */
    double getPositionSpecificGOP(
        const MultipleAlignment &A,
        const MultipleAlignment &B,
        const SubstitutionMatrix &matrix,
        const size_t position) const;

    /** @brief Returns position-specific gap extension penalty.
     *
     * Position specific gap extension penalty is computed considering
     * presence of gaps: if there is a gap, GEP is halved.
     * @param[in]       A         First (group of) sequence
     * @param[in]       B         Second (group of) sequence
     * @param[in]       position  Current positon
     * @return Position-specific gap extension penalty
     */
    double getPositionSpecificGEP(
        const MultipleAlignment &A,
        const MultipleAlignment &B,
        const size_t position) const;

    /** @brief Returns a substitution matrix.
     *
     * Returns substitution matrix depending on the distance of the 
     * (groups of) sequences to align. Substitution matrix is shifted so
     * that its minimum value is 0.
     * @param[in]       guide_tree Guide tree
     * @return Substitution matrix
     */
    SubstitutionMatrix &getSubstitutionMatrix(const RootedTree &guide_tree) const;



    /** @brief Returns match score.
     *
     * Computes match score between given MSAs, at given positions, using
     * given substitutioin matrix and given weights for the sequences.
     * @param[in]           A         First MSA
     * @param[in]           B         Second MSA
     * @param[in]           i         Row
     * @param[in]           j         Column
     * @param[in]           matrix    Substitution matrix
     * @param[in]           weights   Weights of the sequences
     * @return Score of the alignment.
     */
    static double S(
        const MultipleAlignment &A,
        const MultipleAlignment &B,
        const size_t i,
        const size_t j,
        const SubstitutionMatrix &matrix,
        Weights &weights);
};

}  // namespace Phylo
}  // namespace Victor

#endif  // _VICTOR_PHYLO_CLUSTALW_H_