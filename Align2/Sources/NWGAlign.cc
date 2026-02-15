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
// --*- C++ -*------x-----------------------------------------------------------
//
//
// Description:     Implement Needleman-Wunsch-Gotoh global alignment.
//
// -----------------x-----------------------------------------------------------

#include <NWGAlign.h>

namespace Victor { namespace Align2{

	static int DIAG = 0;
	static int LEFT = 1;
	static int UP = 2;

    // CONSTRUCTORS:
    /**
     * 
     * @param ad
     * @param gf
     * @param ss
     */
	NWGAlign::NWGAlign(AlignmentData *ad, GapFunction *gf, ScoringScheme *ss)
    : Align(ad, gf, ss), F1((ad->getSequence(1)).size() + 1),
	  	  	  	  	     F2((ad->getSequence(1)).size() + 1),
						 D((ad->getSequence(1)).size() + 1),
						 D1((ad->getSequence(1)).size() + 1),
						 D2((ad->getSequence(1)).size() + 1){

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

        pCalculateMatrix(true);
    }
    /**
     * 
     * @param ad
     * @param gf
     * @param ss
     * @param v1
     * @param v2
     */
	NWGAlign::NWGAlign(AlignmentData *ad, GapFunction *gf, ScoringScheme *ss,
            const vector<unsigned int> &v1, const vector<unsigned int> &v2)
    : Align(ad, gf, ss) {
        pCalculateMatrix(v1, v2, true);
    }

	NWGAlign::NWGAlign(const NWGAlign &orig) : Align(orig) {
    }

	NWGAlign::~NWGAlign() {
    }


    // OPERATORS:
    /**
     * 
     * @param orig
     * @return 
     */
	NWGAlign&
	NWGAlign::operator =(const NWGAlign &orig) {
        if (&orig != this)
            copy(orig);
        POSTCOND((orig == *this), exception);
        return *this;
    }


    // PREDICATES:
    /**
     * 
     */
    void
	NWGAlign::getMultiMatch() {
    	//TODO
    }

    vector<string>
    NWGAlign::getMatch() const {

        string res1, res2;
        res1Pos.clear();
        res2Pos.clear();

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

        	if (i == nextI) {
                res1 = res1 + "-";
                res1Pos.push_back(INVALID_POS);
            } else {
                res1 = res1 + (ad->getSequence(1))[i - 1];
                res1Pos.push_back(i - 1);
            }
            if (j == nextJ) {
                res2 = res2 + "-";
                res2Pos.push_back(INVALID_POS);
            } else {
                res2 = res2 + (ad->getSequence(2))[j - 1];
                res2Pos.push_back(j - 1);
            }
            ad->calculateMatch(i, nextI, j, nextJ);

            i = nextI;
            j = nextJ;
        }

        reverse(res1.begin(), res1.end());
        reverse(res2.begin(), res2.end());
        reverse(res1Pos.begin(), res1Pos.end());
        reverse(res2Pos.begin(), res2Pos.end());
        vector<string> res(2);
        res[0] = res1;
        res[1] = res2;

        return res;
    }

    double NWGAlign::getScore() const {
    	return max(max(F[n][m], F1[n][m]), F2[n][m]);
    }

	/**
	*
	* @param num
	* @return
	*/
	vector<Alignment>
	NWGAlign::generateMultiMatch(unsigned int num) {

	  //TODO draft version for one optimal alignment

	  vector<Alignment> va;
	  this->getMatch();
	  double tmpScore = this->getScore();
	  ad->getMatch();
	  Alignment *tmp = new Alignment;
	  *tmp = ad->generateMatch(tmpScore);
	  va.push_back(*tmp);

	  return va;
	}

    // MODIFIERS:

    void
	NWGAlign::copy(const NWGAlign &orig) {
        Align::copy(orig);
    }
    /**
     * 
     * @return 
     */
    NWGAlign*
	NWGAlign::newCopy() {
    	NWGAlign *tmp = new NWGAlign(*this);
        return tmp;
    }

    // HELPERS:
    /**
     * 
     * @param update
     */
    void
	NWGAlign::pCalculateMatrix(bool update) {

        if (update){
        	F[0][0] = F1[0][0] = F2[0][0] = 0;

        	//Init first column
        	for (int i = 1; i <= static_cast<int> (n); i++) {
				F1[i][0] = -gf->getOpenPenalty(1) - gf->getExtensionPenalty(i) * (i - 1);
				D1[i][0] = UP;
				F[i][0] = F2[i][0] = -std::numeric_limits<double>::max();
			}

        	//Init first row
			for (int j = 1; j <= static_cast<int> (m); j++) {
				F2[0][j] = -gf->getOpenPenalty(1) - gf->getExtensionPenalty(j) * (j - 1);
				D2[0][j] = LEFT;
				F[0][j] = F1[0][j] = -std::numeric_limits<double>::max();
			}
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
				m  = F[i-1][j] - gf->getOpenPenalty(i);
				m1 = F1[i-1][j] - gf->getExtensionPenalty(i);
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
				m  = F[i][j-1] - gf->getOpenPenalty(j);
				m2 = F2[i][j-1] - gf->getExtensionPenalty(j);
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

    // SSEA variant
    /**
     * 
     * @param v1
     * @param v2
     * @param update
     */
    void
	NWGAlign::pCalculateMatrix(const vector<unsigned int> &v1,
            const vector<unsigned int> &v2, bool update) {

    	//TODO
    }

}} // namespace
