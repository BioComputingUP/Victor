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
 *@Class              TorsionPotential
 *@Project        Victor
 */

#ifndef _TORSIONPOTENTIAL_H_
#define _TORSIONPOTENTIAL_H_
// Includes:
#include <vector>
#include <Spacer.h>
#include <AminoAcidCode.h>
#include <Potential.h>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {

    /** @brief This class implements a simple torsion potential based on the 
     *    statistical preference of aminoacid types for certain angles.
     * 
     * @Description Is a generic class for administrate different angles. Includes methods that allow to calculate max and min energy, also contains virtual methods that will be defined in the differents angles classes.
     * @This 
     * */
    class TorsionPotential : public Potential {
    public:
        // CONSTRUCTORS/DESTRUCTOR:

        TorsionPotential() {
        }

        virtual ~TorsionPotential() {
            PRINT_NAME;
        }

        // PREDICATES:
        virtual long double calculateEnergy(Spacer& sp) = 0;
        virtual long double calculateEnergy(Spacer& sp, unsigned int index1,
                unsigned int index2) = 0;
        virtual long double calculateEnergy(AminoAcid& aa) = 0;
        virtual long double calculateEnergy(AminoAcid& aa, Spacer& sp) = 0;
        virtual long double calculateEnergy(AminoAcid& diheds, AminoAcidCode code) = 0;

        //Max propensities environement independent
        virtual double calculateMaxEnergy(Spacer& sp);
        virtual double calculateMinEnergy(Spacer& sp);
        virtual double calculateMaxEnergy(unsigned int amino);
        virtual double calculateMinEnergy(unsigned int amino);
        virtual double pReturnMaxPropensities(int amino) = 0;

        virtual double pReturnMinPropensities(int amino) {
            ERROR("ERROR. NOT IMPLEMENTED FOR THIS CLASS.", exception)
        }

        //Max propensities pre-amino angle phi/psi depend (not implemented in all version)

        virtual double pReturnMaxPropensitiesPreAngle(int amino, int prephi, int prepsi) {
            ERROR("ERROR. NOT IMPLEMENTED FOR THIS CLASS.", exception)
        }

        virtual double pReturnMinPropensitiesPreAngle(int amino, int prephi, int prepsi) {
            ERROR("ERROR. NOT IMPLEMENTED FOR THIS CLASS.", exception)
        }
        //Usefull here only for calculate relative energy

        virtual int sGetPropBin2(double p) {
            ERROR("ERROR. NOT IMPLEMENTED FOR THIS CLASS.", exception)
        }

        virtual vector< vector<ANGLES> >* orderedEnergy() {
            ERROR("ERROR. NOT IMPLEMENTED FOR THIS CLASS.", exception);
        }

        virtual string getLabel() {
            return "torsion potential";
        }

        // MODIFIERS:

        // OPERATORS:

    protected:

        // HELPERS:
        virtual void pConstructData() = 0;
        virtual void pResetData() = 0;

        virtual double getOmegaAngle(int prop, long RANGE_OMEGA) {
            ERROR("ERROR. NOT IMPLEMENTED FOR THIS CLASS.", exception)
        }

        virtual double getPhiPsiAngle(int prop, long SIZE_TABLE, int ARC_STEP) {
            ERROR("ERROR. NOT IMPLEMENTED FOR THIS CLASS.", exception)
        }

        // ATTRIBUTES:

    private:

    };

    // ---------------------------------------------------------------------------
    //                            TorsionPotential
    // ---------------------------------------------------------------------------

    /**
     * @Description Calculates the maximum energy for an amino acid
     * @param  index of the amino acid, consider the enum list(unsigned int)
     * @return  the corresponding value( double)
     */
    inline double TorsionPotential::calculateMaxEnergy(unsigned int amino) {
        return -log(pReturnMaxPropensities(static_cast<int> (amino)));
    }

    /**
     * @Description Calculates the minimum energy for an amino acid
     * @param  index of the amino acid, consider the enum list(unsigned int)
     * @return  the corresponding value( double)
     */
    inline double TorsionPotential::calculateMinEnergy(unsigned int amino) {
        return -log(pReturnMinPropensities(static_cast<int> (amino)));
    }

    /**
     * @Description Calculates the maximum energy for the amino acids in the spacer
     * @param  spacer reference(Spacer&)
     * @return  the corresponding value( double)
     */
    inline double TorsionPotential::calculateMaxEnergy(Spacer& sp) {
        long double tmp = 0.0;
        for (unsigned int i = 1; i < sp.sizeAmino() - 1; i++)
            tmp += calculateMaxEnergy(sp.getAmino(i).getCode());
        return tmp;
    }

    /**
     * @Description Calculates the minimum energy for the amino acids in the spacer
     * @param  spacer reference(Spacer&)
     * @return  the corresponding value( double)
     */
    inline double TorsionPotential::calculateMinEnergy(Spacer& sp) {
        long double tmp = 0.0;
        for (unsigned int i = 1; i < sp.sizeAmino() - 1; i++)
            tmp += calculateMinEnergy(sp.getAmino(i).getCode());
        return tmp;
    }


} // namespace
#endif// _TORSIONPOTENTIAL_H_






