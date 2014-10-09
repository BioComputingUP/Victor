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

// Description:     Base class for scoring schemes.
//
// -----------------x-----------------------------------------------------------

#include <ScoringScheme.h>

namespace Biopool {

    // CONSTRUCTORS:
    /**
     * @description
     * @param sub
     * @param ad
     * @param str
     */
    ScoringScheme::ScoringScheme(SubMatrix *sub, AlignmentData *ad, Structure *str)
    : sub(sub), ad(ad), str(str) {
        if ((!checkSequence(ad->getSequence(1))) ||
                (!checkSequence(ad->getSequence(2)))) {
            cout << "Illegal sequence:\n"
                    << ad->getSequence(1) << "\n"
                    << ad->getSequence(2) << "\n"
                    << sub->getResidues() << endl;
            ERROR("Error checking sequence.", exception);
        }
    }

    ScoringScheme::ScoringScheme(const ScoringScheme &orig) {
        copy(orig);
    }

    ScoringScheme::~ScoringScheme() {
    }


    // OPERATORS:
    /**
     * @description
     * @param orig
     * @return 
     */
    ScoringScheme&
            ScoringScheme::operator =(const ScoringScheme &orig) {
        if (&orig != this)
            copy(orig);
        POSTCOND((orig == *this), exception);
        return *this;
    }


    // PREDICATES:
    /**
     * @description
     * @param s
     * @return 
     */
    bool
    ScoringScheme::checkSequence(const string &s) const {
        string residues = sub->getResidues();
        bool found = false;

        for (unsigned int i = 0; i < s.size(); i++) {
            found = false;

            for (unsigned int j = 0; j < residues.size(); j++)
                if (residues[j] == s[i]) {
                    found = true;
                    break;
                }

            if (!found) {
                cout << "checkSequence unsuccessful!" << endl;
                return false;
            }
        }

        return true;
    }


    // MODIFIERS:
    /**
     * @description
     * @param orig
     */
    void    
    ScoringScheme::copy(const ScoringScheme &orig) {
        sub = orig.sub->newCopy();
        ad = orig.ad->newCopy();
        str = orig.str->newCopy();
    }

    void
    ScoringScheme::reverse() {
        if (str != 0)
            str->reverse();
    }

} // namespace
