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

#ifndef __CWGapFunction_h__
#define __CWGapFunction_h__

#include <MsaGapFunction.h>
#include <Sequence.h>
#include <string>
#include <vector>
#include <map>

namespace Victor { namespace Phylogenesis {

	/** @brief Implements ClustalW Multiple Sequence Alignment (MSA) gap function
	 *  @author Andrea Nasari (andrea.nasari\@studenti.unipd.it)
	 *
	 **/
	class CWGapFunction: public MsaGapFunction {
	public:

		/// Default constructor.
		CWGapFunction(double initialGop, double initiapGep, const std::vector<Sequence*>& vertical, const std::vector<Sequence*>& horizontal);

		/// Destructor.
		virtual ~CWGapFunction();

		/// Return open gap penalty for horizontal seq at position p.
		double getHorizontalOpenPenalty(int pos);

		/// Return extension gap penalty for  horizontal seq at position p.
		double getHorizontalExtensionPenalty(int pos);

		/// Return open gap penalty for vertical seq at position p.
		double getVerticalOpenPenalty(int pos);

		/// Return extension gap penalty for vertical seq at position p.
		double getVerticalExtensionPenalty(int pos);

	private:

		// ATTRIBUTES:
		double initialGop, initialGep;
		std::vector<double> horizontalGops;
		std::vector<double> verticalGops;
		std::vector<double> horizontalGeps;
		std::vector<double> verticalGeps;

		/// Init initial gap values
		void initGapValues(const std::vector<Sequence*>& seqs, std::vector<double>& gops, std::vector<double>& geps);
	};
}}

#endif
