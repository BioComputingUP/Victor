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

#ifndef __ClustalWScoringScheme_h__
#define __ClustalWScoringScheme_h__

#include <MsaScoringScheme.h>
#include <Sequence.h>
#include <SubMatrix.h>

using namespace Victor::Align2;

namespace Victor { namespace Phylogenesis{

	/** @brief    Implements ClustalW Multiple Sequence Alignment (MSA) scoring function.
	 *  @author   Andrea Nasari (andrea.nasari\@studenti.unipd.it)
	 *
	 **/
	class CWScoringScheme : public MsaScoringScheme{

		public:

			/// Default constructor.
			CWScoringScheme(SubMatrix *sub, std::vector<Sequence*> vertical, std::vector<Sequence*> horizontal);

			/// Scoring function.
			virtual double scoring(int i, int j);

		private:

			// ATTRIBUTES

			SubMatrix* sub;

			std::vector<Sequence*> horizontal, vertical;
	};

}}

#endif
