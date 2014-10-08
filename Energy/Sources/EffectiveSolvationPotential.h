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

#ifndef _EFFECTIVESOLVATIONPOTENTIAL_H_
#define _EFFECTIVESOLVATIONPOTENTIAL_H_


// Includes:
#include <vector>
#include <Spacer.h>
#include <Potential.h>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {

    const double SOLVATION_CUTOFF_DISTANCE_EFFECTIVE = 10.0;

    /**
     * @brief Implements a knowledge-based solvation with polar/hydrophobic information potential. 
     * A coefficient is used to normalize the propensity.
     * 
     *@Description  This class implements a knowledge-based solvation with polar/hydrophobic 
     *    information potential. The final values are normalized by an hardcoded coefficient (see source)
     * */
    class EffectiveSolvationPotential : public Potential {
    public:

        // CONSTRUCTORS/DESTRUCTOR:
        EffectiveSolvationPotential();

        virtual ~EffectiveSolvationPotential() {
            PRINT_NAME;
        }

        // PREDICATES:
        virtual long double calculateEnergy(Spacer& sp);
        virtual long double calculateEnergy(Spacer& sp, unsigned int index1,
                unsigned int index2);
        virtual long double calculateEnergy(AminoAcid& aa, Spacer& sp);

        // MODIFIERS:

        // OPERATORS:

    protected:

        // HELPERS:
        bool isPolar(AminoAcid& aa);
        double pCalcFracBuried(unsigned int index, Spacer& sp);

    private:

        // ATTRIBUTES:

        vector<double> solvCoeff;

    };

    // ---------------------------------------------------------------------------
    //                          EffectiveSolvationPotential
    // -----------------x-------------------x-------------------x-----------------

    /**
     *@Description Verifies if the amino acid is a polar one
     *@param reference of the aa(aminoAcid&)
     *@return  result of the validation(bool)
     */

    inline bool EffectiveSolvationPotential::isPolar(AminoAcid& aa) {
        AminoAcidCode c = static_cast<AminoAcidCode> (aa.getCode());
        if ((c == ARG) ||
                (c == ASN) ||
                (c == ASP) ||
                (c == GLN) ||
                (c == GLU) ||
                (c == LYS) ||
                (c == PRO))
            return true;

        return false;
    }


} // namespace
#endif //_EFFECTIVESOLVATIONPOTENTIAL_H_
