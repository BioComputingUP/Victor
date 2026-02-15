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

#ifndef PHYLOGENESIS_SOURCES_MULTIPLEALIGNMENT_H_
#define PHYLOGENESIS_SOURCES_MULTIPLEALIGNMENT_H_

#include <Sequence.h>
#include <SubMatrix.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;
using namespace Victor::Align2;

namespace Victor {
namespace Phylogenesis {

	/** @brief  Multiple Sequence Alignment result
	 * 	@author Andrea Nasari (andrea.nasari\@studenti.unipd.it)
	 *
	 **/
	class MultipleAlignment {
	public:

		/// Default constructor.
		MultipleAlignment(vector<Sequence*> sequences);

		/// Deconstructor.
		virtual ~MultipleAlignment();

		/// Save single sequence in FASTA format.
		static void saveFasta(string t, string tName, ostream &output);

		/// Save as FASTA like output.
		virtual void saveFasta(ostream &output) const;

		/// Save single line in CLUSTAL format.
		static void saveClustal(string t, string tName, ostream &output, unsigned int from, unsigned int maxNameLength);

		/// Save as CLUSTAL like output.
		virtual void saveClustal(ostream &output) const;

		/// Sequences
		vector<Sequence*> sequences;

	private:

		unsigned int  maxNameLength;
	};

}}

#endif
