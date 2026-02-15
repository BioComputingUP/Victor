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

#include "CWScoringScheme.h"

namespace Victor { namespace Phylogenesis{

	CWScoringScheme::CWScoringScheme(SubMatrix *sub, std::vector<Sequence*> vertical, std::vector<Sequence*> horizontal) :
			vertical(vertical), horizontal(horizontal), sub(sub){
	}

	double CWScoringScheme::scoring(int i, int j){

		//notice that in the alignment algorithm sequences indices start at 1 and not 0
		double score = 0;
		for(unsigned int ii=0; ii<vertical.size(); ii++){
			char v_i = vertical[ii]->sequence[i-1];
			for(unsigned int jj=0; jj<horizontal.size(); jj++){
				//v_i with h_j
				char h_j = horizontal[jj]->sequence[j-1];
				//score with gap equals 0
				if(v_i != '-' && h_j != '-'){
					score += (sub->score[v_i][h_j] * vertical[ii]->normalizedWeight * horizontal[jj]->normalizedWeight);
				}
			}
		}
		return score/(double)(horizontal.size()*vertical.size());
	}

}}
