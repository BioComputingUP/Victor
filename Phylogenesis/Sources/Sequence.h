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

#ifndef __Sequence_H__
#define __Sequence_H__

#include <string>
#include <vector>
#include <fstream>

namespace Victor { namespace Phylogenesis{

	/** @brief    It represents a generic sequence with name.
	 *  @author   Andrea Nasari (andrea.nasari\@studenti.unipd.it)
	 *
	 **/
	class Sequence{
	public:
		/// Default constructor.
		Sequence();

		/// Constructor.
		Sequence(std::string name, std::string sequence, int id);

		/// Check if there is a gap at position
		bool isGap(int pos);

		/// Check if amino acid at pos is hydrophilic
		bool isHydrophilic(int pos);

		/// Check if there is an hydrophilic stretch at position
		bool isHydrophilicStretch(int pos);

		/// Read FASTA format input.
		static std::vector<Sequence*> LoadFasta(std::ifstream& inputFile);

		/// Get sequence without gaps
		static std::string GetPureSequence(const std::string &s);

		/// Get average length of the sequence
		static double GetAvgLength(std::vector<Sequence*> seq);

		// ATTRIBUTES:

		/// Sequence's name
		std::string name;
		/// Generic sequence
		std::string sequence;
		/// Id
		int id;
		/// Sequence weight
		double weight;
		/// Normalized sequence weight
		double normalizedWeight;
		/// Is divergent
		bool divergent;


	};

}}

#endif
