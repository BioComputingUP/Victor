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
/**
 *@Class               EnergyFeatures 
 *@Project        Victor
 *@Description 
 *    Interface/wrapper for energy feature calculation, e.g. in FRST2.
 */
#ifndef _ENERGYFEATURE_H_
#define _ENERGYFEATURE_H_


// Includes:
#include <vector>
#include <Spacer.h>
#include <RapdfPotential.h>
#include <SolvationPotential.h>
#include <PhiPsi.h>
#include <PhiPsiOmegaChi1Chi2.h>


// Global constants, typedefs, etc. (to avoid):

namespace Biopool {

    /**@brief  Interface/wrapper for energy feature calculation, e.g. in FRST2.
     * 
     *@Description 
     *@This 
     * */
    class EnergyFeatures {
    public:

        // CONSTRUCTORS/DESTRUCTOR:
        EnergyFeatures();

        virtual ~EnergyFeatures() {
            PRINT_NAME;
        }

        // PREDICATES:
        static void showFeatures();
        static double calculateBackboneHydrogenBonds(Spacer& sp);
        static vector<double> calculateAaComposition(Spacer& sp);
        static vector<double> calculateSecondaryComposition(Spacer& sp);
        static vector<double> calculateMesoStateComposition(Spacer& sp);
        static vector<double> calculateChainBreaks(Spacer& sp);
        static double calculateClashes(Spacer& sp);

        /// Main function wrapping up feature calculation for FRST2
        vector<double> calculateFeatures(Spacer& sp);

        // MODIFIERS:

        // OPERATORS:

    protected:

    private:

        // HELPERS:

        // ATTRIBUTES:

        SolvationPotential solv;
        RapdfPotential rapdf;
        TorsionPotential* tors;
        TorsionPotential* tors2;

    };

    // ---------------------------------------------------------------------------
    //                                EnergyFeatures
    // -----------------x-------------------x-------------------x-----------------


} // namespace
#endif //_ENERGYFEATURE_H_
