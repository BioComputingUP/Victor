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

#ifndef _POLARSOLVATIONPOTENTIAL_H_
#define _POLARSOLVATIONPOTENTIAL_H_


// Includes:
#include <vector>
#include <Spacer.h>
#include <Potential.h>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {

    const double SOLVATION_CUTOFF_DISTANCE_POLAR = 7.0;

    /**
     * @Brief   This class implements a knowledge-based solvation with polar/hydrophobic information potential.
     * 
     * @Description    The likelihood for each of the 20 amino acids to adopt a given relative solvent accessibility is used to derive the pseudo-energy. The relative solvent accessibility is estimated as the number of other
     *   Cβ atoms within a sphere of radius 7 Å centered on the residue’s Cβ atom. See also SolvationPotential.
     *
    */
    
    class PolarSolvationPotential : public Potential {
    public:

        // CONSTRUCTORS/DESTRUCTOR:
        PolarSolvationPotential();

        virtual ~PolarSolvationPotential() {
            PRINT_NAME;
        }

        // PREDICATES:

        virtual long double calculateEnergy(Spacer& sp) {
            return calculateSolvation(sp);
        }
        virtual long double calculateEnergy(Spacer& sp, unsigned int index1,
                unsigned int index2);

        virtual long double calculateEnergy(AminoAcid& aa, Spacer& sp) {
            return calculateSolvation(aa, sp);
        }
        long double calculateSolvation(Spacer& sp);
        long double calculateSolvation(AminoAcid& aa, Spacer& sp,
                unsigned int start = 0,
                unsigned int end = 9999);

        // MODIFIERS:

        // OPERATORS:

    protected:

        // HELPERS:

    private:

        // ATTRIBUTES:
        vector<vector<vector<int> > > sumPolar;

        static unsigned int MAX_BINS;
        static unsigned int BIN_POLAR;
        static string SOLV_PARAM_FILE;
    };

} // namespace
#endif //_POLARSOLVATIONPOTENTIAL_H_
