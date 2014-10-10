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


#ifndef _RAPDFPOTENTIAL_H_
#define _RAPDFPOTENTIAL_H_
// Includes:
#include <vector>
#include <Potential.h>

// Global constants, typedefs, etc. (to avoid):
const unsigned int MAX_BINS = 18;
const unsigned int MAX_TYPES = 168;

 namespace Victor { namespace Energy {

    /**
     * @brief Distance-dependent residue-specific all-atom probability discriminatory function.
     * 
     *@Description
     * The RAPDF potential (Samudrala and Moult, 1998) is a distance-dependent residue-specific all-atom probability discriminatory function. It discriminates between 167 protein heavy atom types, meaning that
     * different types are assigned to each common atom of the 20 amino acids. E.g. the Cα of an isoleucine is a different type from the Cα of a glycine. Interactions are divided based on distance (d) into 18
     * bins (b), covering distances of up to 20 Å (d → b := [0,..., 3]→ 1, (3,..., 4]→ 2,(4,..., 5]→ 3,...,(19,..., 20]→ 18). The originally published parameters (Samudrala et al., 1998), as downloaded from the ProStar website (URL: http://prostar.carb.nist.gov/) were used. 
     * 
     * */
    class RapdfPotential : public Potential {
    public:

        // CONSTRUCTORS/DESTRUCTOR:
        RapdfPotential();

        virtual ~RapdfPotential() {
            PRINT_NAME;
        }

        // PREDICATES:
        virtual long double calculateEnergy(Spacer& sp);
        virtual long double calculateEnergy(Spacer& sp, unsigned int index1,
                unsigned int index2);
        virtual long double calculateEnergy(AminoAcid& aa, Spacer& sp);
        virtual long double calculateEnergy(AminoAcid& aa, AminoAcid& aa2);
        virtual long double calculateEnergy(Atom& at1, Atom& at2, string aaType,
                string aaType2);

        // MODIFIERS: 
        // OPERATORS:

    protected:

        // HELPERS:
        unsigned int pGetDistanceBinOne(double distance);
        unsigned int pGetGroupBin(const char* group_name);

        // ATTRIBUTES:
        double prob[MAX_BINS][MAX_TYPES][MAX_TYPES];

    public:

        string path;
    private:

    };

    // ---------------------------------------------------------------------------
    //                               RapdfPotential
    // -----------------x-------------------x-------------------x-----------------

        /**
     *@Description calculates the energy between two atoms 
     *@param   the references to the atoms (Atom&,Atom&)the amino acid types(string,string)
     *@return    the value of maximum propensity(long double)
     */
    inline long double RapdfPotential::calculateEnergy(Atom& at1, Atom& at2, string aaType,
            string aaType2) {
        double d = at1.distance(at2);
        if ((d >= 20.0) || (at1.getType() == "OXT") || (at2.getType() == "OXT")
                || ((at1.getType() == "CB") && (aaType == "GLY"))
                || ((at2.getType() == "CB") && (aaType2 == "GLY")))
            return 0.0;
        string tmp1 = threeLetter2OneLetter(aaType) + at1.getType();
        string tmp2 = threeLetter2OneLetter(aaType2) + at2.getType();
        unsigned int dist = pGetDistanceBinOne(d);
        unsigned int grp1 = pGetGroupBin(tmp1.c_str());
        unsigned int grp2 = pGetGroupBin(tmp2.c_str());
        if (dist + grp1 + grp2 < 999)
            return prob[dist][grp1][grp2];
        else // ignore errors
            return 0;
    }

}} // namespace
#endif //_RAPDFPOTENTIAL_H_
