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

#include <SequenceData.h>
#include "MsaNWG.h"

namespace Victor { namespace Phylogenesis {

static int DIAG = 0;
static int LEFT = 1;
static int UP = 2;

MsaNWG::MsaNWG(MsaScoringScheme* ss, MsaGapFunction* gf) :
		ss(ss), gf(gf){
}

std::vector<Sequence*> MsaNWG::getMatch(const vector<Sequence*>& vertical, const vector<Sequence*>& horizontal){

	n = vertical[0]->sequence.length();
	m = horizontal[0]->sequence.length();

	F  = vector<vector<double> >(n + 1);
	F1 = vector<vector<double> >(n + 1);
	F2 = vector<vector<double> >(n + 1);
	D  = vector<vector<int> >(n + 1);
	D1 = vector<vector<int> >(n + 1);
	D2 = vector<vector<int> >(n + 1);

	vector<double> frow(m + 1, -std::numeric_limits<double>::max());
	vector<int> drow(m + 1, -1);

	for (unsigned int i = 0; i < F.size(); ++i) {
		F[i] =  frow;
		F1[i] = frow;
		F2[i] = frow;
		D[i]  = drow;
		D1[i] = drow;
		D2[i] = drow;
	}

	calculateMatrix(vertical, horizontal);

	//Traceback

	string res1;
	string res2;

	int i = n;
	int j = m;
	int nextI = -1;
	int nextJ = -1;

	double val = max(max(F[n][m], F1[n][m]), F2[n][m]);
	int K = -1;

	if(val == F[n][m]){
		K = DIAG;
	}
	else if(val == F2[n][m]){
		K = LEFT;
	}
	else if(val == F1[n][m]){
		K = UP;
	}
	else{
		ERROR("Error in NWGAlign", exception);
	}

	while (i > 0 || j > 0)  {
		if(K == DIAG){
			nextI = i-1;
			nextJ = j-1;
			K = D[i][j];
		}
		else if(K == UP){
			nextI = i-1;
			nextJ = j;
			K = D1[i][j];
		}
		else if(K == LEFT){
			nextI = i;
			nextJ = j-1;
			K = D2[i][j];
		}
		else{
			ERROR("Error in NWGAlign", exception);
		}

		// Here reports only where there are gaps in the profile alignment
		if (i == nextI) {
			res1 = res1 + "-";
		} else {
			res1 = res1 + "*";
		}
		if (j == nextJ) {
			res2 = res2 + "-";
		} else {
			res2 = res2 + "*";
		}

		i = nextI;
		j = nextJ;
	}

	reverse(res1.begin(), res1.end());
	reverse(res2.begin(), res2.end());

	//Insert gaps

	for (unsigned int i = 0; i < res1.length(); ++i){
		if(res1[i] == '-'){
			for(unsigned int j = 0; j < vertical.size(); ++j){
				vertical[j]->sequence.insert(i, 1, '-');
			}
		}
	}
	for (unsigned int i = 0; i < res2.length(); ++i){
		if(res2[i] == '-'){
			for(unsigned int j = 0; j < horizontal.size(); ++j){
				horizontal[j]->sequence.insert(i, 1, '-');
			}
		}
	}

	std::vector<Sequence*> msa(vertical);
	msa.insert(msa.end(),horizontal.begin(),horizontal.end());
	return msa;
}

void MsaNWG::calculateMatrix(const vector<Sequence*>& vertical, const vector<Sequence*>& horizontal)
{
	F[0][0] = F1[0][0] = F2[0][0] = 0;

	//Init first column
	F1[1][0] = 0;
	D1[1][0] = UP;
	F[1][0] = F2[1][0] = -std::numeric_limits<double>::max();
	for (int i = 2; i <= static_cast<int> (n); i++) {
		F1[i][0] = F1[i-1][0] -gf->getVerticalExtensionPenalty(i);
		D1[i][0] = UP;
		F[i][0] = F2[i][0] = -std::numeric_limits<double>::max();
	}

	//Init first row
	F2[0][1] = 0;
	D2[0][1] = LEFT;
	F[0][1] = F1[0][1] = -std::numeric_limits<double>::max();
	for (int j = 2; j <= static_cast<int> (m); j++) {
		F2[0][j] = F2[0][j-1] - gf->getHorizontalExtensionPenalty(j);
		D2[0][j] = LEFT;
		F[0][j] = F1[0][j] = -std::numeric_limits<double>::max();
	}

    for (int i = 1; i <= static_cast<int> (n); i++){
        for (int j = 1; j <= static_cast<int> (m); j++) {

            //F - diag (i-1,j-1)
        	double s  = ss->scoring(i, j);
        	double m  = F[i-1][j-1] + s;
			double m1 = F1[i-1][j-1] + s;
			double m2 = F2[i-1][j-1] + s;
			double mv = max(max(m, m1), m2);

			if(mv == m){
				F[i][j] = mv;
				D[i][j] = DIAG;
			}
			else if(mv == m2){
				F[i][j] = m2;
				D[i][j] = LEFT;
			}
			else if(mv == m1){
				F[i][j] = m1;
				D[i][j] = UP;
			}
			else{
				ERROR("Error in NWGAlign", exception);
			}

            //F1 - up (i-1,j)
			m  = F[i-1][j] - gf->getVerticalOpenPenalty(i);
			m1 = F1[i-1][j] - gf->getVerticalExtensionPenalty(i);
			mv = max(m, m1);

			if(mv == m){
				F1[i][j] = mv;
				D1[i][j] = DIAG;
			}
			else if(mv == m1){
				F1[i][j] = m1;
				D1[i][j] = UP;
			}
			else{
				ERROR("Error in NWGAlign", exception);
			}

            //F2 - left (i,j-1)
			m  = F[i][j-1] - gf->getHorizontalOpenPenalty(j);
			m2 = F2[i][j-1] - gf->getHorizontalExtensionPenalty(j);
			mv = max(m, m2);

			if(mv == m){
				F2[i][j] = mv;
				D2[i][j] = DIAG;
			}
			else if(mv == m2){
				F2[i][j] = m2;
				D2[i][j] = LEFT;
			}
			else{
				ERROR("Error in NWGAlign", exception);
			}
        }
   }
}

double MsaNWG::getScore(){
	return max(max(F[n][m], F1[n][m]), F2[n][m]);
}

}}
