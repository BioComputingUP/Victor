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

#include <IdentityPercentage.h>
#include <Sequence.h>
#include <SequenceData.h>
#include <ScoringS2S.h>
#include <NWAlign.h>

namespace Victor {
namespace Phylo {
    IdentityPercentage::IdentityPercentage(
        SubMatrix &matrix,
        GapFunction &gap)
    : substitution_matrix(matrix), gap_function(gap) {
    }



    IdentityPercentage::IdentityPercentage(const IdentityPercentage &origin)
    : IdentityPercentage(origin.substitution_matrix, origin.gap_function) {
    }



    IdentityPercentage::~IdentityPercentage() {
    }



    double
    IdentityPercentage::computeDistance(const string &seq1, const string &seq2) const {
        return 1.0 - identityPercentage(seq1, seq2);
    }



    double
    IdentityPercentage::identityPercentage(const string &seq1, const string &seq2) const {
        AlignmentData *alignment_data = new SequenceData(2, seq1, seq2, "", "");
        ScoringS2S    scoring_scheme(&substitution_matrix, alignment_data, nullptr, 1.0);
        AlignmentBase alignment;

        NWAlign align(alignment_data, &gap_function, &scoring_scheme);

        vector<string> result = align.getMatch();

        return alignment.calculatePairwiseIdentity(result[0], result[1]);
    }
}  // namespace Phylo
}  // namespace Victor