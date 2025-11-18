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

#define HYDROPHILIC_RANGE 5
#define HYDROPHILIC_RESIDUES "DEGKNQPRS"
#include <MultipleAlignment.h>
#include <algorithm>

/** @brief Tells whether a residue is a hydrophiic residue.
 * @param[in]           residue   Residue to check
 * @retval true         Residue is hydrophilic
 * @retval false        Residue is not hydrophilic
 */
static bool is_hydrophiic_residue(const char residue) {
    for (char r : HYDROPHILIC_RESIDUES) {
        if (residue == r) {
            return true;
        }
    }

    return false;
}



/** @brief Tells whether a sequence has a hydrophilic stretch.
 * @param[in]           sequence  Sequence to check
 * @param[in]           position  Position to check
 * @retval true         The sequence has a hydrophilic stretch
 * @retval false        The sequence has not hydrophilic stretches
 */
static bool has_hydrophilic_stretch(const string &sequence, const size_t position) {
    size_t start = position - HYDROPHILIC_RANGE,
           end   = position + HYDROPHILIC_RANGE;
    bool found = false;

    if (start > sequence.size()) start = 0;
    if (end   > sequence.size()) end   = 0;

    bool previous = false;
    size_t p = start, consecutives = 0;
    while (!found && p < end) {
        if (is_hydrophiic_residue(sequence[p])) {
            if (previous) {
                consecutives++;
            }
            previous = true;
        } else {
            consecutives = 0;
            previous = false;
        }

        found = consecutives == HYDROPHILIC_RANGE;
        p++;
    }

    return found;
}



namespace Victor {
namespace Phylo {
    const string MultipleAlignment::residues = "ARNDCQEGHILKMFPSTWYVUBZX-";


    MultipleAlignment::MultipleAlignment(const vector<Sequence> &sequences)
    : sequences(sequences) {
    }



    MultipleAlignment::MultipleAlignment(const Sequence &sequence) {
        sequences.push_back(sequence);
    }



    MultipleAlignment::MultipleAlignment() {
        MultipleAlignment(vector<Sequence>());
    }



    MultipleAlignment::~MultipleAlignment() {
    }



    void
    MultipleAlignment::operator=(const MultipleAlignment &other) {
        sequences = other.sequences;
    }



    bool
    MultipleAlignment::isEmpty() const {
        return getSize() == 0;
    }



    bool
    MultipleAlignment::hasSequence(const string &identifier) const {
        return getSequence(identifier) < getSize();
    }



    bool
    MultipleAlignment::hasHydrophilicStretch(const size_t position) const {
        bool found = false;

        size_t i = 0;
        while (!found && i < getSize()) {
            const string &sequence = getSequence(i).getResidues();
            found = has_hydrophilic_stretch(sequence, position);
            i++;
        }

        return found;
    }



    size_t
    MultipleAlignment::getSize() const {
        return sequences.size();
    }



    size_t
    MultipleAlignment::getLength() const {
        return (!isEmpty()) ? getSequence(0).getResidues().size() : 0;
    }



    Sequence
    MultipleAlignment::getSequence(const size_t index) const {
        return (index < getSize()) ? sequences[index] : Sequence("", "");
    }



    size_t
    MultipleAlignment::getSequence(const string &identifier) const {
        size_t index = 0;

        for (index = 0; index < getSize(); index++) {
            if (getSequence(index).getIdentifier() == identifier) {
                return index;
            }
        }

        return index;
    }



    string
    MultipleAlignment::getResidues(const size_t position) const {
        char *s = new char[getSize()];
        size_t idx = 0;

        if (position < getLength()) {
            for (auto &sequence : sequences) {
                s[idx++] = sequence.getResidue(position);
            }
        }

        string result(s);
        delete [] s;
        return result;
    }



    size_t
    MultipleAlignment::getResiduesNumber(const size_t position) const {
        return getResidues(position).size();
    }



    size_t
    MultipleAlignment::getNonGapNumber(const size_t position) const {
        return getResiduesNumber(position) - getGapNumber(position);
    }



    size_t
    MultipleAlignment::getGapNumber(const size_t position) const {
        return getCount('-', position);
    }



    size_t
    MultipleAlignment::getCount(
        const char residue,
        const size_t position) const {
        size_t count = 0;
        
        for (auto &sequence : sequences) {
            count += sequence.getResidue(position) == residue;
        }
        
        return count;
    }



    double
    MultipleAlignment::getFrequency(
        const char residue,
        const size_t position) const {
        return (getResiduesNumber(position) != 0)
             ? (double) getCount(residue, position) / (double) getResiduesNumber(position)
             : 0.0;
    }



    map<char, double> &
    MultipleAlignment::getFrequencies(const size_t position) const {
        map<char, double> *frequencies = new map<char, double>();

        for (char residue : residues) {
            (*frequencies)[residue] = getFrequency(residue, position);
        }

        return *frequencies;
    }



    double
    MultipleAlignment::getPercentIdentity() const {
        size_t conserved = 0;

        for (size_t j = 0; j < getSize(); j++) {
            const string column = getResidues(j);
            size_t i            = 1;
            bool is_conserved   = true;
            while (is_conserved && i < getSize()) {
                is_conserved &= column[i] == column[0];
                i++;
            }
            conserved += is_conserved ? 1 : 0;
        }

        return (double) conserved / getLength();
    }



    char
    MultipleAlignment::getConsensus(const size_t position) const {
        char consensus   = '-';
        size_t max_count = 0;

        for (char residue : residues) {
            if (residue == '-') continue;
            const size_t count = getCount(residue, position);

            if (count > max_count) {
                max_count = count;
                consensus = residue;
            }
        }

        return consensus;
    }

    

    string &
    MultipleAlignment::getConsensus() const {
        string *consensus = new string();

        for (size_t i = 0; i < getLength(); i++) {
            (*consensus) += getConsensus(i);
        }

        return *consensus;
    }



    string &
    MultipleAlignment::asString() const {
        size_t longest_name = 0;
        string *result = new string();
        vector<string> names(getSize());
        vector<size_t> positions(getSize());

        // Finds longest name of sequence
        for (size_t i = 0; i < getSize(); i++) {
            const size_t length = getSequence(i).getIdentifier().size();
            longest_name = std::max(longest_name, length);
        }

        // Pads names
        for (size_t i = 0; i < getSize(); i++) {
            const string identifier = getSequence(i).getIdentifier();
            string name             = identifier;
            for (size_t j = identifier.size(); j < longest_name; j++) {
                name += " ";
            }
            name += "    ";
            names[i] = name;
        }

        // Prints data
        for (size_t j = 0; j < getLength(); j += 50) {
            for (size_t i = 0; i < getSize(); i++) {
                bool print_positions = false;
                (*result) += names[i];
                for (size_t k = 0; k < 50 && j + k < getLength(); k++) {
                    const char symbol = getSequence(i).getResidues()[j + k];
                    (*result)       += symbol;
                    positions[i]    += (symbol != '-');
                    print_positions |= (symbol != '-');
                }
                if (print_positions) {
                    (*result) += " " + to_string(positions[i]);
                }
                (*result) += "\n";
            }
            (*result) += "\n\n";
        }

        return *result;
    }



    void
    MultipleAlignment::saveClustalW(ostream &output) const {
        output << asString() << endl;
    }
}  // namespace Phylo
}  // namespace Victor