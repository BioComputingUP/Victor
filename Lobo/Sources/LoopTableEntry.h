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

#ifndef _LOOPTABLEENTRY_H_
#define _LOOPTABLEENTRY_H_

// Includes:
#include <VectorTransformation.h>
#include <vector3.h>
#include <matrix3.h>
#include <Debug.h>

namespace Victor { namespace Lobo {

    // Global constants, typedefs, etc. (to avoid):

    /**@brief Implements an entry for the loop table (equivalent to the
     *    old "ProteinTableEntry" from the previous version. 
     * 
     *@Description  The two structs below serve to store a single entry in the protein table.
     *	To this end they contain endPoint, endDirection and normal, which are 
     *      required for successful concatenation. The main difference between these 
     *      two structs is that CompressedProteinTableEntry contains a table entry 
     *      which had to be packed due to memory restrictions, hence the name.
     *	ProteinTableEntry contains two helper functions to concatenate a single
     *	amino acid combination and to rotate its endPoint back into the XY plane.
     * */
    class LoopTableEntry {
    public:

        // CONSTRUCTORS/DESTRUCTOR:
        LoopTableEntry();
        LoopTableEntry(const LoopTableEntry& orig);
        virtual ~LoopTableEntry();

        // PREDICATES:
        // 
        double calculateDeviation(const LoopTableEntry&, unsigned int nAmino = 0) const;

        float getBOND_LENGTH_CALPHA_TO_CPRIME();
        float getBOND_LENGTH_CPRIME_TO_N();
        float getBOND_ANGLE_AT_CALPHA_TO_CPRIME();
        float getBOND_ANGLE_AT_CPRIME_TO_N();
        float getBOND_LENGTH_N_TO_CALPHA();
        float getBOND_ANGLE_AT_N_TO_CALPHA();

        static double getBondLengthTol() {
            return BOND_LENGTH_TOL;
        }

        static double getBondAngleTol() {
            return BOND_ANGLE_TOL;
        }

        // MODIFIERS:
        void copy(const LoopTableEntry& orig);

        void setToSingleAminoAcid(); // initializes the entry for single aminoacid
        void rotate(const vgVector3<float>& axis, const double angle);
        void rotate(const vgMatrix3<float>& rotationMatrix);
        // 

        LoopTableEntry concatenate(const LoopTableEntry&, unsigned int);
        // 
        LoopTableEntry setToOrigin(LoopTableEntry&, unsigned int,
                VectorTransformation& vt);
        // 
        float rotateIntoXYPlane(VectorTransformation& vt);
        // 

        static void setBondLengthTol(float t) {
            BOND_LENGTH_TOL = t;
        }

        static void setBondAngleTol(float t) {
            BOND_ANGLE_TOL = t;
        }

        // OPERATORS:
        LoopTableEntry& operator=(const LoopTableEntry& orig);
        vgVector3<float>& operator[](unsigned int n);
        const vgVector3<float>& operator[](unsigned int n) const;

        // ATTRIBUTES:
        vgVector3<float> endPoint;
        vgVector3<float> endDirection;
        vgVector3<float> endNormal;
        vgVector3<float> midPoint;
        vgVector3<float> midDirection;
        vgVector3<float> midNormal;

        // NB: contrary to usual conventions the attributes above are _public_,
        // for simplicity's sake, because this whole class is essentially a
        // glorified container element (eg. like vgVector3)

        // TESTERS:
        void printTable(unsigned int k = 0) const;

        // ATTRIBUTES:
        static float BOND_LENGTH_TOL;
        static float BOND_ANGLE_TOL;
        static float BOND_LENGTH_CALPHA_TO_CPRIME;
        static float BOND_LENGTH_CPRIME_TO_N;
        static float BOND_ANGLE_AT_CALPHA_TO_CPRIME;
        static float BOND_ANGLE_AT_CPRIME_TO_N;
        static float BOND_LENGTH_N_TO_CALPHA;
        static float BOND_ANGLE_AT_N_TO_CALPHA;
        static float BOND_LENGTH_CALPHA_TO_CPRIME_SD;
        static float BOND_LENGTH_CPRIME_TO_N_SD;
        static float BOND_ANGLE_AT_CALPHA_TO_CPRIME_SD;
        static float BOND_ANGLE_AT_CPRIME_TO_N_SD;
        static float BOND_LENGTH_N_TO_CALPHA_SD;
        static float BOND_ANGLE_AT_N_TO_CALPHA_SD;
        static float LAMBDA_EP;
        static float LAMBDA_ED;
        static float LAMBDA_EN;

    private:

        double pGetRand();

    };

    // for internal use (in LoopTable::read()):

/**@brief This structure allows to manage a compressed loop table information in vectors.
     */
    struct CompressedLoopTableEntry {
    public:

        CompressedLoopTableEntry() : endPoint(0, 0, 0), endDirection(0, 0, 0),
        endNormal(0, 0, 0), midPoint(0, 0, 0), midDirection(0, 0, 0),
        midNormal(0, 0, 0) {
        }

        vgVector3<unsigned short>& operator[](unsigned int n) {
            switch (n) {
                case 0: return endPoint;
                case 1: return endDirection;
                case 2: return endNormal;
                case 3: return midPoint;
                case 4: return midDirection;
                case 5: return midNormal;
            }
            ERROR("Invalid entry.", exception);
            return endPoint;
        }

        vgVector3<unsigned short> endPoint;
        vgVector3<unsigned short> endDirection;
        vgVector3<unsigned short> endNormal;
        vgVector3<unsigned short> midPoint;
        vgVector3<unsigned short> midDirection;
        vgVector3<unsigned short> midNormal;
    };

    // ---------------------------------------------------------------------------
    //                               LoopTableEntry
    // -----------------x-------------------x-------------------x-----------------

    /**
     *@Description returns a random number
     *@param  none
     *@return  corresponding value ( double)
     */
    inline double LoopTableEntry::pGetRand() {
        double tmp = 0.0;
        for (unsigned int i = 0; i < 12; i++)
            tmp += static_cast<double> (rand()) / RAND_MAX;

        return tmp - 6;
    }

    /**
     *@Description returns the BOND_LENGTH_CALPHA_TO_CPRIME
     *@param  none
     *@return  corresponding value ( float)
     */
    inline float LoopTableEntry::getBOND_LENGTH_CALPHA_TO_CPRIME() {
        return BOND_LENGTH_CALPHA_TO_CPRIME
                + BOND_LENGTH_TOL * BOND_LENGTH_CALPHA_TO_CPRIME_SD * pGetRand();
    }

    /**
     *@Description returns the value of BOND_LENGTH_CPRIME_TO_N
     *@param  none
     *@return  corresponding value ( float)
     */
    inline float LoopTableEntry::getBOND_LENGTH_CPRIME_TO_N() {
        return BOND_LENGTH_CPRIME_TO_N
                + BOND_LENGTH_TOL * BOND_LENGTH_CPRIME_TO_N_SD * pGetRand();
    }

    /**
     *@Description returns the value of  BOND_ANGLE_AT_CALPHA_TO_CPRIME
     *@param  none
     *@return  corresponding value ( float)
     */

    inline float LoopTableEntry::getBOND_ANGLE_AT_CALPHA_TO_CPRIME() {
        return BOND_ANGLE_AT_CALPHA_TO_CPRIME
                + BOND_ANGLE_TOL * BOND_ANGLE_AT_CALPHA_TO_CPRIME_SD * pGetRand();
    }

    /**
     *@Description returns the value of  BOND_ANGLE_AT_CPRIME_TO_N
     *@param  none
     *@return  corresponding value ( float)
     */
    inline float LoopTableEntry::getBOND_ANGLE_AT_CPRIME_TO_N() {
        return BOND_ANGLE_AT_CPRIME_TO_N
                + BOND_ANGLE_TOL * BOND_ANGLE_AT_CPRIME_TO_N_SD * pGetRand();
    }

    /**
     *@Description returns the value of BOND_LENGTH_N_TO_CALPHA
     *@param  none
     *@return  corresponding value ( float)
     */
    inline float LoopTableEntry::getBOND_LENGTH_N_TO_CALPHA() {
        return BOND_LENGTH_N_TO_CALPHA
                + BOND_LENGTH_TOL * BOND_LENGTH_N_TO_CALPHA_SD * pGetRand();
    }

    /**
     *@Description returns the value of  BOND_ANGLE_AT_N_TO_CALPHA
     *@param  none
     *@return  corresponding value ( float)
     */
    inline float LoopTableEntry::getBOND_ANGLE_AT_N_TO_CALPHA() {
        return BOND_ANGLE_AT_N_TO_CALPHA
                + BOND_ANGLE_TOL * BOND_ANGLE_AT_N_TO_CALPHA_SD * pGetRand();
    }


}} // namespace

#endif //_LOOPTABLEENTRY_H_


