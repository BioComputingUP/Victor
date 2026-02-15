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

#include "MultipleAlignment.h"

namespace Victor {
namespace Phylogenesis {

	MultipleAlignment::MultipleAlignment(vector<Sequence*> sequences) : sequences(sequences) {
		maxNameLength = 0;
		for (unsigned int i = 0; i < sequences.size(); ++i) {
			string tName = sequences[i]->name;
			tName = tName.substr(0,tName.find(' '));
			if(tName.length() > maxNameLength)
				maxNameLength = tName.length();
		}
	}

	MultipleAlignment::~MultipleAlignment() {
	}

	void MultipleAlignment::saveFasta(string t, string tName, ostream &output) {
		output << ">" << tName << "\n";
		for (unsigned int i = 0; i < t.length(); i++) {
			if ((i != 0) && ((i % 60) == 0))
				output << "\n";
			output << t[i];
		}
		output << "\n";
	}

	void MultipleAlignment::saveClustal(string t, string tName, ostream &output, unsigned int from, unsigned int maxNameLength) {
		tName = tName.substr(0,tName.find(' '));
		output << tName;
		for (int i = 0; i < (static_cast<int>(maxNameLength) + 10 - static_cast<int> (tName.length())); ++i)
			output << " ";
		unsigned int max = ((from + 60) < t.length()) ? from + 60 : t.length();
		for (unsigned int i = from; i < max; i++)
			output << t[i];
		output << "\n";
	}

	void MultipleAlignment::saveFasta(ostream &output) const {
	   for (unsigned int j = 0; j < sequences.size(); j++)
		   saveFasta(sequences[j]->sequence, sequences[j]->name, output);
	}

	void MultipleAlignment::saveClustal(ostream &output) const {
	   for (unsigned int from = 0; from < sequences[0]->sequence.length(); from += 60) {
		   for (unsigned int j = 0; j < sequences.size(); j++)
			   saveClustal(sequences[j]->sequence, sequences[j]->name, output, from, maxNameLength);
		   output << "\n";
	   }
	}

}}
