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

#include <DistanceMatrixBuilder.h>
#include <Sequence.h>

namespace Victor {
namespace Phylo {
    DistanceMatrixBuilder::~DistanceMatrixBuilder() {
    }



    DistanceMatrix &
    DistanceMatrixBuilder::operator()(const AlignmentBase &alignment) const {
        return apply(alignment);
    }



    DistanceMatrix &
    DistanceMatrixBuilder::apply(const AlignmentBase &alignment) const {
        const size_t n                   = alignment.size() + 1;
        const vector<Sequence> sequences = Sequence::toSequences(alignment);
        DistanceMatrix * const d         = new DistanceMatrix();
        
        // Adds OTUs to the matrix
        for (auto &sequence : sequences) {
            d->addOTU(sequence.getIdentifier());
        }

        // Computes distance matrix
        for (size_t i = 0; i < n; i++) {
            for (size_t j = i + 1; j < n; j++) {
                (*d)(
                    sequences[i].getIdentifier(),
                    sequences[j].getIdentifier(),
                    computeDistance(sequences[i].getResidues(), sequences[j].getResidues())
                );
            }
        }

        return *d;
    }
}  // namespace Phylo
}  // namespace Victor