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

// Includes:
#include <IntCoordConverter.h>
#include <vector3.h>
#include <matrix3.h>
#include <AminoAcid.h>

// Global constants, typedefs, etc. (to avoid):
using namespace Biopool;

// CONSTRUCTORS/DESTRUCTOR:

IntCoordConverter::IntCoordConverter() {
}

IntCoordConverter::~IntCoordConverter() {
    PRINT_NAME;
}

// PREDICATES:

double
IntCoordConverter::getAngle(vgVector3<double>& p11, vgVector3<double>& p12,
        vgVector3<double>& p21, vgVector3<double>& p22) {
    vgVector3 <double> v1 = p12 - p11;
    vgVector3 <double> v2 = p21 - p22;

    return angle(v1, v2);
}

double
IntCoordConverter::getAngle(vgVector3<double>& v1, vgVector3<double>& v2) {
    return angle(v1, v2);
}

vgVector3<double>
IntCoordConverter::calculateTrans(vgVector3<double>& p1, vgVector3<double>& p2) {
    vgVector3<double> t = p2 - p1;
    return t;
}

vgMatrix3<double>
IntCoordConverter::calculateRot(vgVector3<double>& p11, vgVector3<double>& p12,
        vgVector3<double>& p21, vgVector3<double>& p22) {
    IntCoordConverter icc;

    vgVector3 <double> v1 = p12 - p11;
    vgVector3 <double> v2 = p21 - p22;

    double ang = angle(v1, v2);
    vgVector3<double> axis = v1.cross(v2);

    return vgMatrix3<double>::createRotationMatrix(axis, ang);
}

double
IntCoordConverter::getBondAngle(Atom& a1, Atom& a2, Atom& a3) {
    PRECOND((&a1 != &a2) && (&a2 != &a3) && (a1.isBond(a2))
            && (a2.isBond(a3)) && (!a1.isBond(a3)), exception);

    vgVector3<double> v1 = a1.getCoords() - a2.getCoords();
    vgVector3<double> v2 = a3.getCoords() - a2.getCoords();

    INVARIANT(fabs(angle(v1, v2)) <= M_PI, exception);

    return angle(v1, v2);
}

double
IntCoordConverter::getBondLength(Atom& a1, Atom& a2) {
    PRECOND((&a1 != &a2) && (a1.isBond(a2)), exception);
    return (a1.getCoords() - a2.getCoords()).length();
}

/* I took the definition of a torsion angle from the Babel package. 
   Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl */
double
IntCoordConverter::getTorsionAngle(Atom& a1, Atom& a2, Atom& a3, Atom& a4) {
    PRECOND((&a1 != &a2) && (&a2 != &a3) && (&a3 != &a4) && (a1.isBond(a2))
            && (a2.isBond(a3)) && (a3.isBond(a4)) && (!a1.isBond(a3))
            && (!a2.isBond(a4)), exception);

    vgVector3<double> b1 = a1.getCoords() - a2.getCoords();
    vgVector3<double> b2 = a2.getCoords() - a3.getCoords();
    vgVector3<double> b3 = a3.getCoords() - a4.getCoords();

    vgVector3<double> c1 = b1.cross(b2);
    vgVector3<double> c2 = b2.cross(b3);
    vgVector3<double> c3 = c1.cross(c2);

    double ang = angle(c1, c2);

    INVARIANT(fabs(ang) <= M_PI, exception);

    if (b2 * c3 > 0.0)
        return -ang;
    else
        return ang;
}

double
IntCoordConverter::normalize(double angle) const {
    while (angle < -M_PI)
        angle += 2 * M_PI;

    while (angle >= M_PI)
        angle -= 2 * M_PI;

    POSTCOND(-M_PI <= angle && angle < M_PI, exception);
    return angle;
}

double
IntCoordConverter::getAngleDifference(double a, double b) const {
    normalize(a);
    normalize(b);
    double diff = fabs(a - b);
    return min(diff, fabs(diff - 2 * M_PI));
}

// MODIFIERS:

void
IntCoordConverter::setBondAngle(Atom& a1, Atom& a2, Atom& a3, double angle) {
    vgVector3<double> v1 = a1.getCoords() - a2.getCoords();
    vgVector3<double> v2 = a3.getCoords() - a2.getCoords();
    setBondAngle(a1, a2, a3, angle, -(::angle(v1, v2)));
}

void
IntCoordConverter::setBondAngle(Atom& a1, Atom& a2, Atom& a3, double angle,
        double old_angle) {
    PRECOND((&a1 != &a2) && (&a2 != &a3) && (a1.isBond(a2))
            && (a2.isBond(a3)) && (!a1.isBond(a3)), exception);

    if (fabs(angle - old_angle) < EPSILON)
        return;

    vgVector3<double> v1 = a1.getCoords() - a2.getCoords();
    vgVector3<double> v2 = a3.getCoords() - a2.getCoords();

    if (a3.isInBond(a2)) // case: 1 --> 2 --> 3
        a3.addRot(vgMatrix3<double>::createRotationMatrix(v1.cross(v2),
            getAngleDifference(angle, old_angle)));
    else // case: 3 --> 2 --> 1
        a1.addRot(vgMatrix3<double>::createRotationMatrix(v2.cross(v1),
            getAngleDifference(angle, old_angle)));

    POSTCOND(fabs(getBondAngle(a1, a2, a3) - angle) < EPSILON, exception);
}

void
IntCoordConverter::setTorsionAngle(Atom& a1, Atom& a2, Atom& a3, Atom& a4,
        double angle) {
    a3.sync();
    setTorsionAngle(a1, a2, a3, a4, angle, getTorsionAngle(a1, a2, a3, a4));
}

void
IntCoordConverter::setTorsionAngle(Atom& a1, Atom& a2, Atom& a3, Atom& a4,
        double angle, double old_angle) {
    PRECOND((&a1 != &a2) && (&a2 != &a3) && (&a3 != &a4) && (a1.isBond(a2))
            && (a2.isBond(a3)) && (a3.isBond(a4)) && (!a1.isBond(a3))
            && (!a2.isBond(a4)), exception);

    if (fabs(angle - old_angle) < EPSILON)
        return;

    if (a3.isInBond(a2)) // case: 1 --> 2 --> 3 --> 4
    {
        for (unsigned int i = 0; i < a3.sizeOutBonds(); i++)
            a3.getOutBond(i).addRot(vgMatrix3<double>::createRotationMatrix(
                a3.getCoords() - a2.getCoords(), angle - old_angle));
    } else // case: 4 --> 3 --> 2 --> 1
    {
        for (unsigned int i = 0; i < a2.sizeOutBonds(); i++)
            a2.getOutBond(i).addRot(vgMatrix3<double>::createRotationMatrix(
                a2.getCoords() - a3.getCoords(), angle - old_angle));
    }

    POSTCOND((fabs(getTorsionAngle(a1, a2, a3, a4) - angle) < EPSILON)
            || (fabs(getTorsionAngle(a1, a2, a3, a4) - angle) - DEG2RAD * 360.0
            < EPSILON), exception);
}

void
IntCoordConverter::setBondLength(Atom& a1, Atom& a2, double _length) {
    PRECOND((&a1 != &a2) && (a1.isBond(a2)), exception);

    if (a1.isInBond(a2)) {
        if (a1.getTrans().length() != _length)
            a1.setTrans(a1.getTrans().normalize() * _length);
    } else {
        if (a2.getTrans().length() != _length)
            a2.setTrans(a2.getTrans().normalize() * _length);
    }

    POSTCOND(fabs(getBondLength(a1, a2) - _length) < EPSILON, exception);
}

void
IntCoordConverter::connectStructure(AminoAcid& aa, AminoAcid& prev) {
    PRECOND((prev.isMember(OXT)), exception);
    // concatenate the new amino acid with its predecessor
    // by removing the previous OXT and replacing it 
    // with this' N-atom
    aa.bindIn(aa[N], prev, prev[C]);
    vgMatrix3<double> res(1);
    alignVectors(prev[OXT].getTrans(), aa[N].getTrans(), res);
    aa[N].setRot(res);
    aa[N].setTrans(prev[OXT].getTrans().length()
            * aa[N].getTrans().normalize());
    //    prev[OXT].setTrans(aa[N].getTrans().length() 
    //  		       * prev[OXT].getTrans().normalize());
    // align length of C to OXT with C to N
    aa.sync();
    prev.removeAtom(prev[OXT]);
    aa.removeHAtomsfromLeadingNH3();
    prev.setPsi(prev.getPsi()); // apply prev's torsion angles as they ought
    prev.setOmega(prev.getOmega()); // to be (psi & omega are only truly defined
    aa.setPhi(aa.getPhi());
    aa.sync(); // on bound structures)
}

void
IntCoordConverter::connectReverseStructure(AminoAcid& aa, AminoAcid& prev) {
    PRECOND((prev.isMember(OXT)), exception);
    aa.sync();
    prev.sync();

    // concatenate the new amino acid with its predecessor
    // by removing the previous OXT and replacing it 
    // with this' N-atom
    // Align the two structures, such that prev and not aa changes coords:
    vgMatrix3<double> res(1);
    alignVectors(aa[N].getTrans(), prev[OXT].getTrans(), res);
    prev.setRot(res);
    prev[OXT].setTrans(aa[N].getTrans().length()
            * prev[OXT].getTrans().normalize());
    // align length of C to OXT with C to N

    prev.addTrans(aa[N].getCoords() - prev[OXT].getCoords());
    prev.sync();

    prev.removeAtom(prev[OXT]);
    aa.removeHAtomsfromLeadingNH3();

    aa.bindIn(aa[N], prev, prev[C]);
    vgVector3<double> tmp(0.0, 0.0, 0.0);
    aa.setTrans(tmp);


    // pass omega angle on to prev and update its coords:
    if (prev.getOmega() > 990)
        res = vgMatrix3<double>::createRotationMatrix(aa[N].getTrans(),
            (getTorsionAngle(prev[CA], prev[C], aa[N], aa[CA])
            - DEG2RAD * 180));
    else
        res = vgMatrix3<double>::createRotationMatrix(aa[N].getTrans(),
            (getTorsionAngle(prev[CA], prev[C], aa[N], aa[CA])
            - DEG2RAD * prev.getOmega()));
    for (unsigned int i = prev.size(); i > 0; i--)
        prev[i - 1].setCoords(res *
            (prev[i - 1].getCoords() - prev[C].getCoords()) + prev[C].getCoords());

    prev.sync();
    // pass phi angle on to prev and update its coords:
    INVARIANT(aa.getPhi() < 990, exception);
    res = vgMatrix3<double>::createRotationMatrix(-aa[CA].getTrans(),
            (getTorsionAngle(prev[N], prev[CA], prev[C], aa[N])
            - DEG2RAD * prev.getPhi()));
    for (unsigned int i = prev.size(); i > 0; i--)
        prev[i - 1].setCoords(res *
            (prev[i - 1].getCoords() - aa[N].getCoords()) + aa[N].getCoords());

}

/**
 *@Description Converts from Tinker format in cartesian. WARNING: All angles for this function are in *DEGREES* not radians!
 *  @param atblP  : atombondLengthPartner
 *                 atbAP  : atombondAnglePartner
 *                 attAP  : atomtorsionAnglePartner
 *                 at     : result atom
 *@return void
 */

void
IntCoordConverter::zAtomToCartesian(Atom& atbLP, const double bondLength,
        Atom& atbAP, const double bondAngle, Atom& attAP, const double torsionAngle,
        const int chiral, Atom& at) {
    double sine, sine2, a, b, c;
    double cosine;
    vgVector3<double> u1, u2, u3, u4;

    // mathematical and geometrical constants
    double eps = 1e-8;
    double sin_1 = sin(bondAngle / 57.29577951308232088);
    double cos_1 = cos(bondAngle / 57.29577951308232088);
    double sin_2 = sin(torsionAngle / 57.29577951308232088);
    double cos_2 = cos(torsionAngle / 57.29577951308232088);

    vgVector3<double> aCatbAP = atbAP.getCoords();
    vgVector3<double> aCattAP = attAP.getCoords();
    vgVector3<double> aCatbLP = atbLP.getCoords();

    // first, the case where second angle is a dihedral angle
    if (chiral == 0) {

        u1 = (aCatbAP - aCattAP).normalize();
        u2 = (aCatbLP - aCatbAP).normalize();

        u3[0] = u1[1] * u2[2] - u1[2] * u2[1];
        u3[1] = u1[2] * u2[0] - u1[0] * u2[2];
        u3[2] = u1[0] * u2[1] - u1[1] * u2[0];
        cosine = u1[0] * u2[0] + u1[1] * u2[1] + u1[2] * u2[2];

        if (fabs(cosine) < 1.0)
            sine = sqrt(1.0 - (cosine * cosine));
        else {
            ERROR("IntCoordConverter::zAtomToCartesian : Undefined Dihedral Angle",
                    exception);
            sine = sqrt((cosine * cosine) - 1.0);
        }

        u3 /= sine;
        u4[0] = u3[1] * u2[2] - u3[2] * u2[1];
        u4[1] = u3[2] * u2[0] - u3[0] * u2[2];
        u4[2] = u3[0] * u2[1] - u3[1] * u2[0];

        at.setCoords(aCatbLP[0] + bondLength
                * (-u2[0] * cos_1 + u4[0] * sin_1 * cos_2 + u3[0] * sin_1 * sin_2),
                aCatbLP[1] + bondLength
                * (-u2[1] * cos_1 + u4[1] * sin_1 * cos_2 + u3[1] * sin_1 * sin_2),
                aCatbLP[2] + bondLength
                * (-u2[2] * cos_1 + u4[2] * sin_1 * cos_2 + u3[2] * sin_1 * sin_2));


    }// now, the case where second angle is a bondLength angle
    else if (abs(chiral) == 1) {
        u1 = (aCatbLP - aCattAP).normalize();
        u2 = (aCatbAP - aCatbLP).normalize();

        u3[0] = u1[1] * u2[2] - u1[2] * u2[1];
        u3[1] = u1[2] * u2[0] - u1[0] * u2[2];
        u3[2] = u1[0] * u2[1] - u1[1] * u2[0];
        cosine = u1[0] * u2[0] + u1[1] * u2[1] + u1[2] * u2[2];

        if (fabs(cosine) < 1.0)
            sine2 = 1.0 - (cosine * cosine);
        else {
            cout << cosine << endl;
            ERROR("IntLoader::zAtomToCartesian : Defining Atoms Colinear at Atom.", exception);
        }
        a = (-cos_2 - cosine * cos_1) / sine2;
        b = (cos_1 + cosine * cos_2) / sine2;
        c = (a * cos_2 + 1. - b * cos_1) / sine2;

        if (c > eps)
            c = chiral * sqrt(c);
        else if (c < -eps) {
            double temp1 = a * u1[0] + b * u2[0];
            double temp2 = a * u1[1] + b * u2[1];
            double temp3 = a * u1[2] + b * u2[2];
            c = sqrt(temp1 * temp1 + temp2 * temp2 + temp3 * temp3);
            a /= c;
            b /= c;
            c = 0;
#ifdef DEBUG
            ERROR("IntLoader::zAtomToCartesian : Sum of bondLength Angles Too Large.",
                    exception);
#endif
        }
        else
            c = 0;



        at.setCoords(aCatbLP[0] + bondLength
                * (a * u1[0] + b * u2[0] + c * u3[0]),
                aCatbLP[1] + bondLength
                * (a * u1[1] + b * u2[1] + c * u3[1]),
                aCatbLP[2] + bondLength
                * (a * u1[2] + b * u2[2] + c * u3[2]));
    }
}
