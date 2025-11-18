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

#include <FengDoolittleDistance.h>
#include <algorithm>
#include <SequenceData.h>
#include <ScoringS2S.h>
#include <NWAlign.h>
#include <Sequence.h>

namespace Victor {
namespace Phylo {
    FengDoolittleDistance::FengDoolittleDistance(
        SubMatrix &matrix,
        GapFunction &gap)
    : substitution_matrix(matrix), gap_function(gap) {
    }



    FengDoolittleDistance::FengDoolittleDistance(const FengDoolittleDistance &origin)
    : FengDoolittleDistance(origin.substitution_matrix, origin.gap_function) {
    }



    FengDoolittleDistance::~FengDoolittleDistance() {
    }



    double
    FengDoolittleDistance::computeDistance(const string &seq1, const string &seq2) const {
        string r1 = seq1, r2 = seq2;

        std::random_shuffle(r1.begin(), r1.end());
        std::random_shuffle(r2.begin(), r2.end());

        const double S_x  = computeScore(seq1, seq1),
                     S_y  = computeScore(seq2, seq2),
                     S_xy = computeScore(seq1, seq2),
                     S_r  = computeScore(r1, r2),
                     d    = 2.0 * (S_xy - S_r) / (S_x + S_y + 2.0 * S_r);

        return -log((d > 0.0) ? d : 1.0);
    }



    double
    FengDoolittleDistance::computeScore(const string &seq1, const string &seq2) const {
        AlignmentData *alignment_data = new SequenceData(2, seq1, seq2, "", "");
        Structure     *structure      = nullptr;
        const double  cSeq            = 1.0;
        ScoringS2S    scoring_scheme(&substitution_matrix, alignment_data, structure, cSeq);
        NWAlign align(alignment_data, &gap_function, &scoring_scheme);

        return align.getScore();
    }
}  // namespace Phylo
}  // namespace Victor