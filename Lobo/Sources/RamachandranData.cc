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


// Includes:
#include <RamachandranData.h>
#include <IoTools.h>

// Global constants, typedefs, etc. (to avoid):

namespace Victor { namespace Lobo {
   
    float RamachandranData::PHI_ANGLE_TOL = 0.0;
    float RamachandranData::PSI_ANGLE_TOL = 0.0;

    // CONSTRUCTORS/DESTRUCTOR:

    /**
     *@Description Basic constructor
     */
    RamachandranData::RamachandranData() : nextRama(0), ramaPhi(), ramaPsi() {
    }

    /**
     *@Description Constructor based on an original object copy
     */
    RamachandranData::RamachandranData(const RamachandranData& orig) {
        this->copy(orig);
    }

    /**
     *@Description Basic destructor
     */
    RamachandranData::~RamachandranData() {
        PRINT_NAME;
    }


    // PREDICATES:

    /**
     *@Description obtains a random phi value,  
     *@param  flag to consider an advance(bool)
     *@return  corresponding random value(double)
     */
    double RamachandranData::getRandomPhi(bool noAdvance) {
        if (!noAdvance) {
            nextRama++;
            if (nextRama >= ramaPhi.size())
                nextRama = 0;
        }

        INVARIANT(nextRama < ramaPhi.size(), exception);
        return ramaPhi[nextRama] + PHI_ANGLE_TOL * pGetRand();
    }

    /**
     *@Description obtains a random psi value,  
     *@param  flag to consider an advance(bool)
     *@return  corresponding random value(double)
     */
    double RamachandranData::getRandomPsi(bool noAdvance) {
        if (!noAdvance) {
            nextRama++;
            if (nextRama >= ramaPsi.size())
                nextRama = 0;
        }

        INVARIANT(nextRama < ramaPsi.size(), exception);
        return ramaPsi[nextRama] + PSI_ANGLE_TOL * pGetRand();
    }


    // MODIFIERS:

    /**
     *@Description copies the original object
     *@param  reference to the original object(const RamachandranData&)
     *@return  changes are made internally(void)
     */
    void RamachandranData::copy(const RamachandranData& orig) {
        PRINT_NAME;
        nextRama = orig.nextRama;
        ramaPhi = orig.ramaPhi;
        ramaPsi = orig.ramaPsi;
    }

    /**
     *@Description loads data of the object from an input file
     *@param  reference to the input file(istream&)
     *@return  changes are made internally(void)
     */
    void RamachandranData::load(istream& input) {
        unsigned long nEntry;
        PRECOND(input, exception);

        input >> nEntry;

        ramaPhi.reserve(nEntry);
        ramaPsi.reserve(nEntry);

        if (nEntry <= 0)
            ERROR("Invalid number of entries for Rama data.", exception);

        double tempPhi, tempPsi;
        for (unsigned long i = 0; i < nEntry; i++) {
            INVARIANT(input, exception);
            input >> tempPhi >> tempPsi;
            skipToNewLine(input);
            ramaPhi.push_back(tempPhi);
            ramaPsi.push_back(tempPsi);
        }

        nextRama = rand() % nEntry;
    }

    /**
     *@Description Saves the data from the object in an output file
     *@param  reference to the output file(ostream)
     */
    void RamachandranData::save(ostream& output) {
        PRECOND(output, exception);

        output << ramaPhi.size() << "\n";

        for (unsigned long i = 0; i < ramaPhi.size(); i++) {
            INVARIANT(output, exception);
            output << setw(7) << setprecision(4) << ramaPhi[i] << "   "
                    << setw(7) << setprecision(4) << ramaPsi[i] << "   GLY\n";
        }
    }

    /**
     *@Description defines a cluster
     *@param  value for the cutoff (double)
     *@return  changes are made internally(void)
     */
    void RamachandranData::cluster(double cutoff) {
        if (ramaPhi.size() == 0)
            return;

        vector<double> tmpPhi, tmpPsi;

        tmpPhi.push_back(ramaPhi[0]);
        tmpPsi.push_back(ramaPsi[0]);

        cutoff = sqr(cutoff);

        for (unsigned long i = 0; i < ramaPhi.size(); i++) {
            bool contained = false;

            for (unsigned long j = 0; j < tmpPhi.size(); j++) if (sqr(ramaPhi[i] - tmpPhi[j]) + sqr(ramaPsi[i] - tmpPsi[j])
                        < cutoff) {
                    contained = true;
                    break;
                }

            if (!contained) {
                tmpPhi.push_back(ramaPhi[i]);
                tmpPsi.push_back(ramaPsi[i]);
            }
        }

        ramaPhi = tmpPhi;
        ramaPsi = tmpPsi;
    }


    // OPERATORS:

    /**
     *@Description Allows to assign and original object into another
     *@param  the original object reference(const RamachandranData& )
     *@return  reference to the new object( RamachandranData& )
     */
    RamachandranData& RamachandranData::operator=(const RamachandranData& orig) {
        if (&orig != this)
            copy(orig);
        return *this;
    }

}}
