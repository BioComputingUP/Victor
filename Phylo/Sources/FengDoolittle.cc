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

#include <FengDoolittle.h>
#include <SubstitutionMatrix.h>
#include <phylo.h>

using namespace std;
using namespace Victor::Phylo;

/** @brief Replaces gaps in a sequence.
 * @param[in]           sequence  Sequence
 * @return Sequence where gaps have been replaced
 */
static string &replace_gaps(const string &sequence);

/** @brief Returns the PAM 250 substitution matrix.
 * @return PAM 250 matrix
 */
static SubMatrix &init_matrix();



namespace Victor {
namespace Phylo {
    typedef map<string, string> Sequences;

    // For some unknown reasons, matrix cannot be initialized statically (althought
    // it is a static member). A good place is the constructor.
    SubMatrix FengDoolittle::matrix;
    const FitchMargoliash FengDoolittle::build_tree;


    FengDoolittle::FengDoolittle(const double gap_open, const double gap_extension)
    : gap_function(*new AGPFunction(gap_open, gap_extension)) {
        matrix = init_matrix();
    }



    FengDoolittle::FengDoolittle(const FengDoolittle &origin)
    : gap_function(origin.gap_function) {
    }



    FengDoolittle::~FengDoolittle() {
    }



    MultipleAlignment &
    FengDoolittle::apply(const AlignmentBase &alignment) const {
        // Reads sequences
        Sequences sequences = Sequence::toMap(alignment);

        // Computes distance matrix (pairwise alignment)
        FengDoolittleDistance build_matrix(matrix, gap_function);
        DistanceMatrix d = build_matrix(alignment);

        // Computes guide tree
        RootedTree guide_tree = build_tree(d).asRootedTree();

        // Progressive alignments
        return *align(sequences, guide_tree);
    }



    MultipleAlignment *
    FengDoolittle::align(
        Sequences &sequences,
        const RootedTree &node) const {
        ////////////////////////////////////////////////////////////////
        // Single sequence: trivially aligned
        if (node.isLeaf()) {
            const string identifier = node.getLabel(),
                         residues   = sequences[identifier];
            return new MultipleAlignment(Sequence(identifier, residues));
        }


        ////////////////////////////////////////////////////////////////
        // Single inheritance: straightforward
        if (node.getChildrenNumber() == 1) {
            return align(sequences, node[0]);
        }


        ////////////////////////////////////////////////////////////////
        // Two MSAs: performs alignment
        const MultipleAlignment *A = align(sequences, node[0]),
                                *B = align(sequences, node[1]);
        double max_score = -999999.0;
        string horizontal, vertical;


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


        // Finds best alignment
        for (size_t i = 0; i < A->getSize(); i++) {
            const string seqA = replace_gaps(A->getSequence(i).getResidues());
            for (size_t j = 0; j < B->getSize(); j++) {
                const string seqB = replace_gaps(B->getSequence(j).getResidues());

                AlignmentData *data = new SequenceData(2, seqA, seqB, "", "");
                ScoringS2S scoring_scheme(&matrix, data, nullptr, 1.0);
                NWAlign nw_align(data, &gap_function, &scoring_scheme);
                const double score = nw_align.getScore();

                if (score > max_score) {
                    max_score  = score;
                    horizontal = nw_align.getMatch()[0];
                    vertical   = nw_align.getMatch()[1];
                }
            }
        }


        // Builds a new MSA using best alignment
        vector<Sequence> new_sequences;
        for (size_t i = 0; i < A->getSize(); i++) {
            const string identifier = A->getSequence(i).getIdentifier(),
                         sequence   = A->getSequence(i).getResidues();
            size_t cursor   = 0;
            string residues = "";
            for (char residue : horizontal) {
                residues += (residue != '-') ? sequence[cursor++] : '-';
            }
            new_sequences.push_back(Sequence(identifier, residues));
        }
        
        for (size_t i = 0; i < B->getSize(); i++) {
            const string identifier = B->getSequence(i).getIdentifier(),
                         sequence   = B->getSequence(i).getResidues();
            size_t cursor   = 0;
            string residues = "";
            for (char residue : vertical) {
                residues += (residue != '-') ? sequence[cursor++] : '-';
            }
            new_sequences.push_back(Sequence(identifier, residues));
        }


        // Returns multiple alignmnet
        return new MultipleAlignment(new_sequences);
    }

}  // namespace Phylo
}  // namespace Victor



////////////////////////////////////////////////////////////////////////
// Support functions
static string &replace_gaps(const string &sequence) {
    string *new_sequence = new string();

    for (char residue : sequence) {
        (*new_sequence) += (residue != '-') ? residue : 'X';
    }

    return *new_sequence;
}



static SubMatrix &init_matrix() {
    SubstitutionMatrix m(SubstitutionMatrix::PAM250);
    return (m - m.getMinScore()).asSubMatrix();
}