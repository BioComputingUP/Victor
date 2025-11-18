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

#include <LevenshteinDistance.h>

namespace Victor {
namespace Phylo {
    LevenshteinDistance::~LevenshteinDistance() {
    }



    double
    LevenshteinDistance::computeDistance(const string &seq1, const string &seq2) const {
        const size_t m = seq1.size();
        const size_t n = seq2.size();

        if (m == 0) return n;
        if (n == 0) return m;

        size_t *costs = new size_t[n + 1];

        // Initializes cost vector
        for (size_t k = 0; k <= n; k++)
            costs[k] = k;

        // Computes distances
        for (size_t i = 0; i < m; i++) {
            size_t corner = i;
            costs[0] = i + 1;

            for (size_t j = 0; j < n; j++) {
                size_t upper = costs[j + 1];

                costs[j + 1] = (seq1[i] == seq2[j])
                             ? corner
                             : min(costs[j], min(upper, corner)) + 1;
                
                corner = upper;
            }
        }

        double distance = (double) costs[n];
        delete [] costs;

        return distance;
    }
}  // namespace Phylo
}  // namespace Victor