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

#ifndef _INTCOORDCONVERTER_H_
#define _INTCOORDCONVERTER_H_

// Includes
#include <iostream>
#include <Atom.h>
#include <Debug.h>
#include <IoTools.h>
#include <vglMath.h>

namespace Victor { namespace Biopool { 
    
    class AminoAcid;



    const static double EPSILON = 0.0001;

    const double RAD2DEG = 180 / M_PI;
    const double DEG2RAD = M_PI / 180;

    const double UNDEFINED_ANGLE = DEG2RAD * 999;

    /**
     * @brief Implements methods to manage translation vector and angles between vectors
     * 
     *   Includes methods that allow to get, calculate angels and traslation vector, rotation matrix, etc .
     * */
    
    class IntCoordConverter {
    public:
        // CONSTRUCTORS/DESTRUCTOR:
        IntCoordConverter();
        virtual ~IntCoordConverter();

        // PREDICATES:
        static double getBondAngle(Atom& a1, Atom& a2, Atom& a3);
        static double getTorsionAngle(Atom& a1, Atom& a2, Atom& a3, Atom& a4);
        static double getBondLength(Atom& a1, Atom& a2);
        ///returns angle between vectors detected by couple of points
        double getAngle(vgVector3<double>& p11, vgVector3<double>& p12,
                vgVector3<double>& p21, vgVector3<double>& p22);
        ///calcluates angle between two vectors
        double getAngle(vgVector3<double>& v1, vgVector3<double>& v2);
        ///returns translation vector that translates a point into another
        vgVector3<double> calculateTrans(vgVector3<double>& p1, vgVector3<double>& p2);
        ///calculates rotation matrix that rotate a vector into another 
        vgMatrix3<double> calculateRot(vgVector3<double>& p11, vgVector3<double>& p12,
                vgVector3<double>& p21, vgVector3<double>& p22);

        // MODIFIERS:
        void setBondAngle(Atom& a1, Atom& a2, Atom& a3, double angle);
        void setBondAngle(Atom& a1, Atom& a2, Atom& a3, double angle,
                double old_angle);
        void setTorsionAngle(Atom& a1, Atom& a2, Atom& a3, Atom& a4, double angle);
        void setTorsionAngle(Atom& a1, Atom& a2, Atom& a3, Atom& a4,
                double angle, double old_angle);
        void setBondLength(Atom& a1, Atom& a2, double length);

        void connectStructure(AminoAcid& aa, AminoAcid& prev);
        void connectReverseStructure(AminoAcid& aa, AminoAcid& prev);
        void zAtomToCartesian(Atom& atbLP, const double bondLength, Atom& atbAP,
                const double bondAngle, Atom& attAP,
                const double torsionAngle, const int chiral, Atom& at);

        // OPERATORS:

        //protected:
        // HELPERS:
        double normalize(double angle) const;
        /** Get new angle in [0,2pi). */
        double getAngleDifference(double a, double b) const;
        /** Get difference between two angles. */

        // ATTRIBUTES:
    };


    // ---------------------------------------------------------------------------
    //                              IntCoordConverter
    // -----------------x-------------------x-------------------x-----------------




    // PREDICATES:

    // MODIFIERS:

    // OPERATORS:

    // HELPERS:

    /**
     *   Returns the angle between two vectors of coordinates
     * @param a and b contain the 3D coordinates (vgVector3<double>)
     * @return The angle in RADIANTS (double)
     */
    inline
    double
    angle(const vgVector3<double>& a, const vgVector3<double>& b) {
        PRECOND(a.length() * b.length() != 0.0,
                exception);

        long double d = (a * b) / (a.length() * b.length());

        // Cope with numeric inaccurracy.
        if (d > 1.0)
            d = 1.0;
        else if (d < -1.0)
            d = -1.0;

        INVARIANT(fabs(d) <= 1.0, logic_error);

        return acos(d);
    }
    /**
     *   Returns the angle between two vectors of coordinates
     * @param a and b contain the 3D coordinates (vgVector3<float>)
     * @return The angle in RADIANTS (double)
     */
    inline
    double
    angle(const vgVector3<float>& a, const vgVector3<float>& b) {
        PRECOND(a.length() * b.length() != 0.0,
                exception);

        long double d = (a * b) / (a.length() * b.length());

        // Cope with numeric inaccurracy.
        if (d > 1.0)
            d = 1.0;
        else if (d < -1.0)
            d = -1.0;

        INVARIANT(fabs(d) <= 1.0, logic_error);

        return acos(d);
    }

    /**
     *   Align (ie. move) v2 with (to) v1.
     * @param v1 and v2 contain the 3D coordinates (vgVector3<double>)
     */
    

    inline void alignVectors(vgVector3<double> v1, vgVector3<double> v2,
            vgMatrix3<double>& res) {
        PRECOND(v1.length() * v2.length() != 0.0, exception);

        // special case if v1 = -v2 would not rotate at all
        if ((v1 + v2).length() < EPSILON)
            v2.z += EPSILON;
        vgVector3<double> axis = (v1.normalize()).cross(v2.normalize());
        res = vgMatrix3<double>::createRotationMatrix(axis, -angle(v1.normalize(),
                v2.normalize()));

#ifndef NDEBUG
        v2 = res * v2;
        POSTCOND((v2.normalize() - v1.normalize()).length() < EPSILON, exception);
#endif
    }

    /**
     *   Align (ie. move) v2 with (to) v1.
     * @param v1 and v2 contain the 3D coordinates (vgVector3<float>)
     */
    inline void
    alignVectors(vgVector3<float> v1, vgVector3<float> v2,
            vgMatrix3<float>& res) {
        PRECOND(v1.length() * v2.length() != 0.0, exception);

        // special case if v1 = -v2 would not rotate at all
        if ((v1 + v2).length() < EPSILON)
            v2.z += EPSILON;
        vgVector3<float> axis = (v1.normalize()).cross(v2.normalize());
        res = vgMatrix3<float>::createRotationMatrix(axis, -angle(v1.normalize(),
                v2.normalize()));

#ifndef NDEBUG
        v2 = res * v2;
        POSTCOND((v2.normalize() - v1.normalize()).length() < EPSILON, exception);
#endif
    }

}} //namespace
#endif //_INTCOORDCONVERTER_H_

