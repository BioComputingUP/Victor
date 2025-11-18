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

/** Maximum distance to search gaps (position-specific GOP). */
#define GAP_SEARCH_DISTANCE 8
#include <ClustalW.h>

using namespace std;
using namespace Victor::Phylo;

/** Direction for a traceback. */
typedef enum {
    NONE,     ///< No direction
    UP,       ///< Up
    LEFT,     ///< Left
    DIAG      ///< Diagonal
} Direction;

/** @brief Pascarella and Argos residue specific gap modification factors.
 * @return Map from residues to factors.
 */
static map<char, double> PascarellaArgos();

/** @brief Selects best direction.
 * @param[in]           diagonal  Score from the diagonal
 * @param[in]           up        Score from top
 * @param[in]           left      Score from bottom
 * @return Best direction
 */
static inline Direction select_direction(
    const double diagonal,
    const double up,
    const double left);

/** @brief Performs the traceback operation.
 * @param[in]           A         First MSA
 * @param[in]           B         Second MSA
 * @param[in]           direction Traceback matrix
 * @return Vector of residues
 */
static vector<string> traceback(
    const MultipleAlignment &A,
    const MultipleAlignment &B,
    vector<vector<Direction>> &direction);

/** @brief Builds a vector of sequences.
 *
 * Joins sequences from two MSAs and a vector of residues.
 * @param[in]           A         First MSA
 * @param[in]           B         Second MSA
 * @param[in]           residues  Vector of residues
 * @return Vector of sequences
 */
static vector<Sequence> join(
    const MultipleAlignment &A,
    const MultipleAlignment &B,
    vector<string> &residues);





namespace Victor {
namespace Phylo {
    typedef map<string, double> Weights;
    typedef vector<Sequence> Sequences;

    map<char, double> ClustalW::gap_modifiers = PascarellaArgos();

    ClustalW::ClustalW(
        const DistanceMatrixBuilder &distance_matrix_builder,
        const ClusteringAlgorithm &clustering_algorithm,
        const WeightMatrix weight_matrix,
        const double gap_open,
        const double gap_extension)
    : build_matrix(distance_matrix_builder),
      build_tree(clustering_algorithm),
      weight_matrix(weight_matrix),
      gap_open(gap_open),
      gap_extension(gap_extension) {
    }



    ClustalW::ClustalW(const ClustalW &origin)
    : ClustalW(
        origin.build_matrix,
        origin.build_tree,
        origin.weight_matrix,
        origin.gap_open,
        origin.gap_extension
    ) {
    }



    ClustalW::~ClustalW() {
    }



    MultipleAlignment &
    ClustalW::apply(const AlignmentBase &alignment) const {
        // Reads sequences
        Sequences sequences = Sequence::toMap(alignment);
        
        // Computes distance matrix (pairwise alignments)
        DistanceMatrix d = build_matrix(alignment);

        // Computes guide tree and weights
        RootedTree guide_tree = build_tree(d).asRootedTree();
        Weights weights = getWeights(guide_tree);

        // Progressive alignments
        return *align(sequences, guide_tree, weights);
    }



    Weights
    ClustalW::getWeights(const RootedTree &guide_tree) const {
        map<string, double> weights;

        // Computes weighs for every sequence
        double max_weight = 0.0;
        for (auto &leaf : guide_tree.getLeaves()) {
            RootedTree node = *leaf;
            double weight   = node.getDistance();
            const string id = node.getLabel();

            while (!node.isRoot()) {
                weight += node.getDistance() / node.getLeavesNumber();
                node    = node.getParent();
            }

            weights[id] = weight;
            max_weight  = max(weight, max_weight);
        }

        // Normalizes weights
        for (auto &entry : weights) {
            weights[entry.first] = entry.second / max_weight;
        }

        return weights;
    }



    MultipleAlignment *
    ClustalW::align(
        Sequences &sequences,
        const RootedTree &node,
        Weights &weights) const {

        ////////////////////////////////////////////////////////////////
        // Single sequence: just returns it
        if (node.isLeaf()) {
            const string identifier = node.getLabel();
            const string residues   = sequences[identifier];
            return new MultipleAlignment(Sequence(identifier, residues));
        }


        ////////////////////////////////////////////////////////////////
        // Single inheritance: straightforward
        if (node.getChildrenNumber() == 1) {
            return align(sequences, node[0], weights);
        }


        ////////////////////////////////////////////////////////////////
        // Two sequences: performs alignment
        const MultipleAlignment *A   = align(sequences, node[0], weights),
                                *B   = align(sequences, node[1], weights);
        const SubstitutionMatrix ss(getSubstitutionMatrix(node));
        const size_t M = A->getLength(),
                     N = B->getLength();


#if VERBOSE==1
        cout << "Now aligning:" << endl;
        for (size_t i = 0; i < A->getSize(); i++) {
            cout << A->getSequence(i).getIdentifier() << endl;
        }
        cout << "---------- against ----------" << endl;
        for (size_t j = 0; j < B->getSize(); j++) {
            cout << B->getSequence(j).getIdentifier() << endl;
        }
        cout << endl << endl;
#endif


        // Allocates score and traceback matrices
        vector<vector<double>> score(N + 1);
        vector<vector<Direction>> direction(N + 1);
        vector<vector<size_t>> hgaps(N + 1);
        vector<vector<size_t>> vgaps(N + 1);
        for (size_t i = 0; i <= N; i++) {
            score[i]     = vector<double>(M + 1);
            direction[i] = vector<Direction>(M + 1);

            hgaps[i] = vector<size_t>(M + 1);
            vgaps[i] = vector<size_t>(M + 1);
        }


        // Initializes first row
        score[0][0]     = 0;
        direction[0][0] = NONE;
        for (size_t j = 1; j <= M; j++) {
            const double GOP = getPositionSpecificGOP(*A, *B, ss, 0);
            const double GEP = getPositionSpecificGEP(*A, *B, j - 1);
            score[0][j]      = -(GOP + GEP * (j - 1));
            direction[0][j]  = LEFT;
            hgaps[0][j] = 0;
            vgaps[0][j] = 1;
        }


        // Initializes first column
        for (size_t i = 1; i <= N; i++) {
            const double GOP = getPositionSpecificGOP(*B, *A, ss, 0),
                         GEP = getPositionSpecificGEP(*B, *A, i - 1);
            score[i][0]      = -(GOP + GEP * (i - 1));
            direction[i][0]  = UP;
            hgaps[i][0] = 1;
            vgaps[i][0] = 0;
        }


        // Calculates score matrix
        for (size_t i = 1; i <= N; i++) {
            for (size_t j = 1; j <= M; j++) {
                double vgap = 0.0;
                double hgap = 0.0;

                hgap = (direction[i][j - 1] != LEFT)
                     ? getPositionSpecificGOP(*A, *B, ss, j - 1)
                     : getPositionSpecificGEP(*A, *B, j - 1) * hgaps[i][j - 1];
                vgap = (direction[i - 1][j] != UP)
                     ? getPositionSpecificGOP(*B, *A, ss, i - 1)
                     : getPositionSpecificGEP(*B, *A, i - 1) * vgaps[i - 1][j];

                const double d = score[i-1][j-1] + S(*A, *B, i, j, ss, weights),
                             u = score[i - 1][j] - vgap,
                             l = score[i][j - 1] - hgap;

                Direction dir   = select_direction(d, u, l);
                score[i][j]     = max(d, max(u, l));
                direction[i][j] = dir;

                // Written in term of Linear programming for performance
                hgaps[i][j] = 1 + (dir == LEFT) * hgaps[i][j - 1];
                vgaps[i][j] = 1 + (dir == UP)   * vgaps[i - 1][j];
            }
        }


        // Traceback
        vector<string> residues = traceback(*A, *B, direction);


        // Builds sequences
        vector<Sequence> joint = join(*A, *B, residues);


        // Returns multiple alignment
        delete A;
        delete B;
        return new MultipleAlignment(joint);
    }

    

    double
    ClustalW::getInitialGOP(
        const MultipleAlignment &A,
        const MultipleAlignment &B,
        const SubstitutionMatrix &matrix) const {
        const double M = A.getLength(),
                     N = B.getLength(),
                     avg = matrix.getAvgMismatchScore(),
                     identity = A.getPercentIdentity();
        return gap_open + log10(min(M, N)) * avg * identity;
    }



    double
    ClustalW::getInitialGEP(
        const MultipleAlignment &A,
        const MultipleAlignment &B) const {
        const double N = A.getLength(),
                     M = B.getLength();
        return gap_extension * (1.0 + abs(log10(N / M)));
    }



    double
    ClustalW::getPositionSpecificGOP(
        const MultipleAlignment &A,
        const MultipleAlignment &B,
        const SubstitutionMatrix &matrix,
        const size_t position) const {
        double GOP = getInitialGOP(A, B, matrix);

        // If there are sequences with gaps, ignore other rules
        if (A.getGapNumber(position) > 0) {
            return GOP * 0.3 * (1.0 - A.getGapNumber(position) / A.getSize());
        }


        // Checks for gaps within 8 spaces
        size_t distance = 0;
        bool gaps_nearby = false;
        while (!gaps_nearby && distance <= GAP_SEARCH_DISTANCE) {
            distance++;
            const size_t gaps_before = A.getGapNumber(position - distance),
                         gaps_after  = A.getGapNumber(position + distance);
            gaps_nearby = (gaps_before + gaps_after) > 0;
        }
        if (gaps_nearby) {
            GOP = GOP * (4.0 - distance / 4.0);
        }


        // Checks for hydrophilic stretches
        if (A.hasHydrophilicStretch(position)) {
            GOP = GOP * 2.0 / 3.0;
        }


        // Uses Pascarella-Argos table, if no hydropilic stretches
        else {
            const string &residues = A.getResidues(position);
            double average = 0.0;

            for (char residue : residues) {
                average += gap_modifiers[residue];
            }
            average /= A.getSize();

            GOP = GOP * average;
        }

        return GOP;
    }



    double
    ClustalW::getPositionSpecificGEP(
        const MultipleAlignment &A,
        const MultipleAlignment &B,
        const size_t position) const {
        const double GEP = getInitialGEP(A, B);

        return (A.getGapNumber(position) > 0) ? GEP * 0.5 : GEP;
    }



    SubstitutionMatrix &
    ClustalW::getSubstitutionMatrix(const RootedTree &guide_tree) const {
        const double max_distance = guide_tree.getRoot().getMaxDistance();
        const double raw_distance = guide_tree[0].getDistance(guide_tree[1]);
        const double distance     = raw_distance / max_distance;
        SubstitutionMatrix::Identifier identifier;

        // PAM family
        if (weight_matrix == PAM && distance >= 0.8)
            identifier = SubstitutionMatrix::PAM20;
        else if (weight_matrix == PAM && distance >= 0.6)
            identifier = SubstitutionMatrix::PAM60;
        else if (weight_matrix == PAM && distance >= 0.4)
            identifier = SubstitutionMatrix::PAM120;
        else if (weight_matrix == PAM)
            identifier = SubstitutionMatrix::PAM350;

        // BLOSUM family
        else if (weight_matrix == BLOSUM && distance >= 0.8)
            identifier = SubstitutionMatrix::BLOSUM80;
        else if (weight_matrix == BLOSUM && distance >= 0.6)
            identifier = SubstitutionMatrix::BLOSUM62;
        else if (weight_matrix == BLOSUM && distance >= 0.3)
            identifier = SubstitutionMatrix::BLOSUM45;
        else if (weight_matrix == BLOSUM)
            identifier = SubstitutionMatrix::BLOSUM30;

        // Default
        else
            identifier = SubstitutionMatrix::BLOSUM62;


        SubstitutionMatrix matrix(identifier);

        // Shifts matrix so that minimum value is 0
        return matrix - matrix.getMinScore();
    }



    double
    ClustalW::S(
        const MultipleAlignment &A,
        const MultipleAlignment &B,
        const size_t i,
        const size_t j,
        const SubstitutionMatrix &matrix,
        Weights &weights) {
        const size_t M = A.getSize(),
                     N = B.getSize();
        double score = 0.0;
        
        // Full version
        for (size_t a = 0; a < M; a++) {
            const Sequence seqA = A.getSequence(a);
            const char rA       = seqA.getResidues()[j - 1];
            const double wA     = weights[seqA.getIdentifier()];
            for (size_t b = 0; b < N; b++) {
                const Sequence seqB = B.getSequence(b);
                const char rB       = seqB.getResidues()[i - 1];
                const double wB     = weights[seqB.getIdentifier()];

                score += (1 - (rA == '-' || rB == '-')) * matrix(rA, rB) * wA * wB;
            }
        }
        score /= (M * N);

        // Consensus only version
        /*
        const char a = A.getConsensus(i),
                   b = B.getConsensus(j);
        score = (1 - (a == '-' || b == '-')) * matrix(a, b);
        */

        return score;
    }
}  // namespace Phylo
}  // namespace Victor





////////////////////////////////////////////////////////////////////////
// Support (local) functions.
static map<char, double> PascarellaArgos() {
    map<char, double> modifiers;

    modifiers['A'] = 1.13;
    modifiers['C'] = 1.13;
    modifiers['D'] = 0.96;
    modifiers['E'] = 1.31;
    modifiers['F'] = 1.20;
    modifiers['G'] = 0.61;
    modifiers['H'] = 1.00;
    modifiers['I'] = 1.32;
    modifiers['K'] = 0.96;
    modifiers['L'] = 1.21;
    modifiers['M'] = 1.29;
    modifiers['N'] = 0.63;
    modifiers['P'] = 0.74;
    modifiers['Q'] = 1.07;
    modifiers['R'] = 0.72;
    modifiers['S'] = 0.76;
    modifiers['T'] = 0.89;
    modifiers['V'] = 1.25;
    modifiers['Y'] = 1.00;
    modifiers['W'] = 1.23;

    return modifiers;
}



static inline Direction select_direction(
    const double diagonal,
    const double up,
    const double left) {
    return (diagonal > up && diagonal > left)
         ? DIAG
         : ((up > left) ? UP : LEFT);
}



static vector<string> traceback(
    const MultipleAlignment &A,
    const MultipleAlignment &B,
    vector<vector<Direction>> &direction) {
    const size_t M = A.getSize(),
                 N = B.getSize();
    size_t i = B.getLength(),
           j = A.getLength();
    vector<string> residues(M + N);
    Direction next = direction[i][j];

    while (next != NONE) {
        if (next == DIAG) {
            i--; j--;
            for (size_t k = 0; k < M; k++)
                residues[k] = A.getResidues(j)[k] + residues[k];

            for (size_t k = 0; k < N; k++)
                residues[M + k] = B.getResidues(i)[k] + residues[M + k];
        } else if (next == UP) {
            i--;
            for (size_t k = 0; k < M; k++)
                residues[k] = "-" + residues[k];

            for (size_t k = 0; k < N; k++)
                residues[M + k] = B.getResidues(i)[k] + residues[M + k];
        } else {
            j--;
            for (size_t k = 0; k < M; k++)
                residues[k] = A.getResidues(j)[k] + residues[k];

            for (size_t k = 0; k < N; k++)
                residues[M + k] = "-" + residues[M + k];
        }

        next = direction[i][j];
    }

    return residues;
}



static vector<Sequence> join(
    const MultipleAlignment &A,
    const MultipleAlignment &B,
    vector<string> &residues) {
    const size_t M = A.getSize(),
                 N = B.getSize();
    vector<Sequence> sequences;

    for (size_t i = 0; i < M; i++) {
        sequences.push_back(Sequence(
            A.getSequence(i).getIdentifier(),
            residues[i]
        ));
    }

    for (size_t j = 0; j < N; j++) {
        sequences.push_back(Sequence(
            B.getSequence(j).getIdentifier(),
            residues[M + j]
        ));
    }

    return sequences;
}
