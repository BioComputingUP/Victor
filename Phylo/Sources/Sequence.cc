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

#include <Sequence.h>

namespace Victor {
namespace Phylo {
    Sequence::Sequence(const string &identifier, const string &residues)
    : residues(residues) {
        this->identifier = Sequence::formatIdentifier(identifier);
    }


    Sequence::~Sequence() {
    }



    vector<Sequence> &
    Sequence::toSequences(const AlignmentBase &alignment) {
        string identifier, residues;
        vector<Sequence> *sequences = new vector<Sequence>();

        // Inserts "target" sequence
        identifier = alignment.getTargetName();
        residues  = AlignmentBase::getPureSequence(alignment.getTarget());
        sequences->push_back(Sequence(identifier, residues));

        // Inserts every template sequence
        for (size_t i = 0; i < alignment.size(); i++) {
            identifier = alignment.getTemplateName(i);
            residues  = AlignmentBase::getPureSequence(alignment.getTemplate(i));
            sequences->push_back(Sequence(identifier, residues));            
        }

        return *sequences;
    }



    map<string, string> &
    Sequence::toMap(const AlignmentBase &alignment) {
        vector<Sequence> list = toSequences(alignment);
        map<string, string> *sequences = new map<string, string>;

        for (auto &sequence : list) {
            (*sequences)[sequence.getIdentifier()] = sequence.getResidues();
        }

        return *sequences;
    }



    bool
    Sequence::isEmpty() const {
        return residues.empty();
    }



    string
    Sequence::getIdentifier() const {
        return identifier;
    }



    string
    Sequence::getID() const {
        return getIdentifier();
    }



    char
    Sequence::getResidue(const size_t position) const {
        return (position < residues.size()) ? residues[position] : '\0';
    }



    string
    Sequence::getResidues() const {
        return residues;
    }



    string
    Sequence::formatIdentifier(const string &identifier) {
        const string::size_type pos = identifier.find(' ', 0);
        return (pos < identifier.size())
             ? identifier.substr(0, pos)
             : identifier;
    }

}  // namespace Phylo
}  // namespace Victor
