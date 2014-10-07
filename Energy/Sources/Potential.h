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
 *@Class             Potential 
 *@Description 
 *    Base class for energy calculation
 */
#ifndef _POTENTIAL_H_
#define _POTENTIAL_H_


// Includes:
#include <vector>
#include <Spacer.h>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {

    /**@brief Abstract class for the energy potential.
     * 
     *@Description Includes methods that allow to obtain propensity. 
     * It also defines a "struct" variable with angles and energy associated to an aminoacid.
     * 
     * */
    
    class Potential {
    public:

        struct ANGLES {
            double phi;
            double psi;
            double omega;
            double energy;
        };

        // CONSTRUCTORS/DESTRUCTOR:

        Potential() {
        }

        Potential(string inputFile) {
        }

        virtual ~Potential() {
            PRINT_NAME;
        }

        // PREDICATES:
        virtual long double calculateEnergy(Spacer& sp) = 0;
        virtual long double calculateEnergy(Spacer& sp, unsigned int index1,
                unsigned int index2) = 0;
        virtual long double calculateEnergy(AminoAcid& aa, Spacer& sp) = 0;

        virtual long double calculateEnergy(AminoAcid& resid, AminoAcidCode type, Spacer& sp) {
            ERROR("ERROR. NOT IMPLEMENTED FOR THIS CLASS.", exception)
        }

        virtual long double pReturnMaxPropensity(const AminoAcidCode type) const {
            ERROR("ERROR. NOT IMPLEMENTED FOR THIS CLASS.", exception);
        }

        virtual long double pReturnMinPropensity(const AminoAcidCode type) const {
            ERROR("ERROR. NOT IMPLEMENTED FOR THIS CLASS.", exception);
        }

        virtual vector< vector<ANGLES> >* orderedEnergy() {
            ERROR("ERROR. NOT IMPLEMENTED FOR THIS CLASS.", exception);
        }
        // MODIFIERS:

        // OPERATORS:

    protected:

        // HELPERS:

        // ATTRIBUTES:

    private:

    };

    // ---------------------------------------------------------------------------
    //                                  Potential
    // -----------------x-------------------x-------------------x-----------------

    /**@Description Energy operator. It compares enrgies of two Potential::ANGLES structures
     * @param    Two ANGLES structures associated to two aminoacids
     * @return    True, if the energy of the object 1 is greater than the energy of the object 2 (bool)  
     * 
     * */
            
    struct EnergyGreater : public binary_function<const Potential::ANGLES &, const Potential::ANGLES &, bool> {

        bool operator()(const Potential::ANGLES &ref1, const Potential::ANGLES &ref2) {
            return ref1.energy > ref2.energy;
        }
    };

} // namespace
#endif //_POTENTIAL_H_
