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

#ifndef _MSANWG_H_
#define _MSANWG_H_

#include <MsaGapFunction.h>
#include <MsaScoringScheme.h>
#include <Traceback.h>
#include <Sequence.h>
#include <vector>

namespace Victor { namespace Phylogenesis {

using namespace Victor::Align2;

	/** @brief  Implement Needleman-Wunsch-Gotoh global alignment for MSA.
	 * 	@author Andrea Nasari (andrea.nasari\@studenti.unipd.it)
	 *
	 **/
	class MsaNWG{

	public:

		/// Default constructor.
		MsaNWG(MsaScoringScheme* ss, MsaGapFunction* gf);

		/// Returns a sequence vector containing a multiple alignment with maximal score.
		std::vector<Sequence*> getMatch(const vector<Sequence*>& vertical, const vector<Sequence*>& horizontal);

		 /// Returns alignment score.
		double getScore();

	private:

		// ATTRIBUTES

		MsaScoringScheme* ss;
		MsaGapFunction* gf;

		vector< vector<double> > F;  ///< Gotoh diagonal score matrix.
		vector< vector<double> > F1; ///< Gotoh left score matrix.
		vector< vector<double> > F2; ///< Gotoh up score matrix.
		vector< vector<int> > D;  ///< Gotoh diag traceback matrix.
		vector< vector<int> > D1; ///< Gotoh left traceback matrix.
		vector< vector<int> > D2; ///< Gotoh up traceback matrix.

		unsigned int n; ///< Length of vertical sequences.
		unsigned int m; ///< Length of horizontal sequences.

		void calculateMatrix(const vector<Sequence*>& vertical, const vector<Sequence*>& horizontal);
	};

}}

#endif
