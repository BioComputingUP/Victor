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


#ifndef _RAMACHANDRANDATA_H_
#define _RAMACHANDRANDATA_H_

// Includes:
#include <vector>
#include <Debug.h>

namespace Victor { namespace Lobo {

    // Global constants, typedefs, etc. (to avoid):

    /** 
     *@Description Implements the container for ramachandran plot-like 
    *    phi/psi angle combinations for the LoopTable class(es).
    *      
    */
    inline double sqr(double x) {
        return x*x;
    }

    /**
     *@brief This struct implements the container for ramachandran plot-like 
     *    phi/psi angle combinations for the LoopTable class(es).
     * */
    struct RamachandranData {
    public:
        // CONSTRUCTORS/DESTRUCTOR:
        RamachandranData();
        RamachandranData(const RamachandranData& orig);
        virtual ~RamachandranData();

        // PREDICATES:
        double getRandomPhi(bool noAdvance = false);
        double getRandomPsi(bool noAdvance = false);

        static double getAngleTol() {
            return (PHI_ANGLE_TOL + PSI_ANGLE_TOL) / 2.0;
        }

        // MODIFIERS:
        void copy(const RamachandranData& orig);
        void load(istream& input);
        void save(ostream& output);
        void cluster(double cutoff);

        static void setAngleTol(float t) {
            PHI_ANGLE_TOL = t;
            PSI_ANGLE_TOL = t;
        }

        static void setPhiAngleTol(float t) {
            PHI_ANGLE_TOL = t;
        }

        static void setPsiAngleTol(float t) {
            PSI_ANGLE_TOL = t;
        }

        // OPERATORS:
        RamachandranData& operator=(const RamachandranData& orig);

    private:

        // ATTRIBUTES:
        static float PHI_ANGLE_TOL;
        static float PSI_ANGLE_TOL;

        unsigned long nextRama;
        vector<double> ramaPhi;
        vector<double> ramaPsi;

        // HELPERS:

        double pGetRand();

    };


    // ---------------------------------------------------------------------------
    //                               RamachandranData
    // -----------------x-------------------x-------------------x-----------------

    /**
     *@Description returns a random value
     *@param  none
     *@return  the corresponding value( double)
     */
    inline double RamachandranData::pGetRand() {
        double tmp = 0.0;
        for (unsigned int i = 0; i < 12; i++)
            tmp += static_cast<double> (rand()) / RAND_MAX;

        return tmp - 6;
    }

}} // namespace

#endif //_RAMACHANDRANDATA_H_
