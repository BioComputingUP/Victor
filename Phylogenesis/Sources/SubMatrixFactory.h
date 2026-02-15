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

#ifndef __SubMatrixFactory_h__
#define __SubMatrixFactory_h__

#include <SubMatrix.h>

namespace Victor { namespace Phylogenesis {

	using namespace Victor::Align2;

	enum SubMatrixFamily{
		BLOSUM,
		PAM
	};

	enum SubMatrixName {
		BLOSUM30,
		BLOSUM45,
		BLOSUM62,
		BLOSUM80,
		PAM20,
		PAM60,
		PAM120,
		PAM350
	};

	/** @brief  Implement a factory for SumMatrix objects
	 *  @author Andrea Nasari (andrea.nasari\@studenti.unipd.it)
	 *
	 **/
	class SubMatrixFactory {
		private:
			SubMatrix* blosum30;
			SubMatrix* blosum45;
			SubMatrix* blosum62;
			SubMatrix* blosum80;
			SubMatrix* pam20;
			SubMatrix* pam60;
			SubMatrix* pam120;
			SubMatrix* pam350;
		public:

			/// Default constructor.
			SubMatrixFactory(bool positive);

			/// Deconstructor.
			virtual ~SubMatrixFactory();

			/// Returns the SubMatrix by name.
			SubMatrix* getSubMatrix(SubMatrixName name);

			/// Rescore the SubMatrix with positive values and multiplied by the scalar.
			void rescore(SubMatrix* matrix, double scalar = 1.0);
	};

}}

#endif
