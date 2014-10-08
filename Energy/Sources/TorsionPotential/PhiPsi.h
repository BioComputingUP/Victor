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


#ifndef _PHIPSI_H_
#define _PHIPSI_H_


// Includes:
#include <vector>
#include <Spacer.h>
#include <AminoAcidCode.h>
#include <TorsionPotential.h>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {

    /** @brief class manages the angle qualities and the energy 
     * 
     * @Description Includes methods that allow to obtain information about the angle and the energy. This class implements a simple torsion potential based on the 
     *    statistical preference of aminoacid types for phi and psi angles.

     * */
    class PhiPsi : public TorsionPotential {
    public:

        // CONSTRUCTORS/DESTRUCTOR:
        PhiPsi(int SET_ARC1 = 10,
                string knownledge = "data/tor.par"); //default knownledge TOP500

        virtual ~PhiPsi() {
            PRINT_NAME;
        }

        // PREDICATES:
        virtual long double calculateEnergy(Spacer& sp);
        virtual long double calculateEnergy(Spacer& sp, unsigned int index1,
                unsigned int index2);
        virtual long double calculateEnergySmooth(AminoAcid& aa);
        virtual long double calculateEnergy(AminoAcid& aa);
        virtual long double calculateEnergy(AminoAcid& aa, double phi, double psi);

        virtual long double calculateEnergy(AminoAcid& aa, Spacer& sp) {
            return calculateEnergy(aa);
        }
        virtual long double calculateEnergy(AminoAcid& diheds, AminoAcidCode code);
        virtual double pReturnMaxPropensities(int amino);
        virtual double pReturnMinPropensities(int amino);
        int sGetPropBinDiff(double p, double& diff);

        string getLabel() {
            return "phi-psi";
        }

        // MODIFIERS:
        virtual void setArcStep(int n);

        int getArcStep() {
            return ARC_STEP;
        }
        // OPERATORS:
        int getPhiPsiIndexTable(double angle);
        int getPropensity(AminoAcidCode aa, int phi, int psi);
        double getEnergyFromPhiPsi(AminoAcidCode code, double phi, double psi);

    protected:

        // HELPERS:
        virtual void pConstructData();
        virtual void pResetData();
        virtual double pGetMaxPropensities(int amino);
        virtual double pGetMinPropensities(int amino);
        int sGetPropBin(double p);
        void sAddProp(int code, int x, int y);
        virtual void pConstructMaxPropensities();
        virtual void pConstructMinPropensities();


        // ATTRIBUTES:

        string TOR_PARAM_FILE; // File with prop torsion angles
        int ARC_STEP; // important: must be a divisior of 360 !!!! 
        int SIZE_OF_TABLE; // "granularity" props.
        int amino_count[AminoAcid_CODE_SIZE]; //number of amino.
        vector<vector<vector<int>* >* > propensities; // the propensities table.
        vector<vector<int>* > all_propensities; // the sum of propropensities table.
        double total; //total numer of ammino considered.
        vector<double> amino_max_propensities; //vector with max amino propensities
        // according to knowledge.
        vector<double> amino_min_propensities; //vector with min amino propensities
        // according to knowledge.

    private:

    };

    // ---------------------------------------------------------------------------
    //                            PhiPsi

} // namespace
#endif// _PHIPSI_H_






