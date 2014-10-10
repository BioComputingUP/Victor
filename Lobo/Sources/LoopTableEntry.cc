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
#include <LoopTableEntry.h>
#include <IntCoordConverter.h>

// Global constants, typedefs, etc. (to avoid):

//using namespace Victor;
namespace Victor { namespace Lobo {
    float LoopTableEntry::BOND_LENGTH_TOL = 2.0;
    float LoopTableEntry::BOND_ANGLE_TOL = 2.0;
    float LoopTableEntry::BOND_LENGTH_CALPHA_TO_CPRIME = 1.52;
    float LoopTableEntry::BOND_LENGTH_CPRIME_TO_N = 1.33;
    float LoopTableEntry::BOND_ANGLE_AT_CALPHA_TO_CPRIME = 111.6;
    float LoopTableEntry::BOND_ANGLE_AT_CPRIME_TO_N = 116.4;
    float LoopTableEntry::BOND_LENGTH_N_TO_CALPHA = 1.458;
    float LoopTableEntry::BOND_ANGLE_AT_N_TO_CALPHA = 121.7;
    float LoopTableEntry::BOND_LENGTH_CALPHA_TO_CPRIME_SD = 0.02;
    float LoopTableEntry::BOND_LENGTH_CPRIME_TO_N_SD = 0.015;
    float LoopTableEntry::BOND_ANGLE_AT_CALPHA_TO_CPRIME_SD = 2.5;
    float LoopTableEntry::BOND_ANGLE_AT_CPRIME_TO_N_SD = 2.0;
    float LoopTableEntry::BOND_LENGTH_N_TO_CALPHA_SD = 0.016;
    float LoopTableEntry::BOND_ANGLE_AT_N_TO_CALPHA_SD = 1.7;
    float LoopTableEntry::LAMBDA_EP = 0.5;
    float LoopTableEntry::LAMBDA_ED = 1.0;
    float LoopTableEntry::LAMBDA_EN = 1.0;

    /**
     *@Description Calculates the square of a number
     *@param  value to multiply to by itself(float)
     *@return  corresponding value(float)
     */
    static float sqr(float x) {
        return x*x;
    }


    // CONSTRUCTORS/DESTRUCTOR:

    /**
     *@Description basic constructor
     */
    LoopTableEntry::LoopTableEntry() : endPoint(0, 0, 0),
    endDirection(1, 0, 0), endNormal(0, 0, 1), midPoint(0, 0, 0),
    midDirection(0, 0, 0), midNormal(0, 0, 0) {
    }

    /**
     *@Description Constructor based on the copy of another original object
     *@param reference to the original object(const LoopTableEntry&)  
     */
    LoopTableEntry::LoopTableEntry(const LoopTableEntry& orig) {
        this->copy(orig);
    }

    /**
     *@Description basic destructor
     */
    LoopTableEntry::~LoopTableEntry() {
        PRINT_NAME;
    }


    // PREDICATES:

    /**
     *@Description Calculates the deviation between src1 and src2 (the object data and the one given 
     * parameter )
     *@param  reference to the Loop Table entry(const LoopTableEntry&) , number of amino acids(unsigned int)
     *@return  corresponding value(double)
     */
    double LoopTableEntry::calculateDeviation(const LoopTableEntry& other,
            unsigned int nAmino) const {
        double tmpLambda = LAMBDA_EP;

        return tmpLambda * (sqr(endPoint.x - other.endPoint.x)
                + sqr(endPoint.y - other.endPoint.y)
                + sqr(endPoint.z - other.endPoint.z))
                + LAMBDA_ED * (sqr(endDirection.x - other.endDirection.x)
                + sqr(endDirection.y - other.endDirection.y)
                + sqr(endDirection.z - other.endDirection.z))
                + LAMBDA_EN * (sqr(endNormal.x - other.endNormal.x)
                + sqr(endNormal.y - other.endNormal.y)
                + sqr(endNormal.z - other.endNormal.z))
                ;
    }

    // MODIFIERS:

    /**
     *@Description copies an original object into the structure
     *@param  reference to the original object(const LoopTableEntry&)
     *@return  changes are made internally(void)
     */
    void LoopTableEntry::copy(const LoopTableEntry& orig) {
        PRINT_NAME;
        endPoint = orig.endPoint;
        endDirection = orig.endDirection;
        endNormal = orig.endNormal;
        midPoint = orig.midPoint;
        midDirection = orig.midDirection;
        midNormal = orig.midNormal;
    }

    /**
     *@Description Defaults the protein table entry to the simple case where it is a single amino acid.
     *@param  none
     *@return  changes are made internally(void)
     */
    void LoopTableEntry::setToSingleAminoAcid() {
        vgVector3<float> tmpV(0.0, getBOND_LENGTH_N_TO_CALPHA(), 0.0);
        vgVector3<float> t1(getBOND_LENGTH_CALPHA_TO_CPRIME()
                * sin(getBOND_ANGLE_AT_CALPHA_TO_CPRIME()
                / 57.29577951308232088),
                tmpV[1] - (getBOND_LENGTH_CALPHA_TO_CPRIME()
                * cos(getBOND_ANGLE_AT_CALPHA_TO_CPRIME()
                / 57.29577951308232088)), 0.0);
        vgVector3<float> t3(0, 0, 1);
        vgVector3<float> t4(0, 0, 0);


        endPoint = t1;

        vgVector3<float> t0(0, getBOND_LENGTH_N_TO_CALPHA(), 0);

        endDirection = vgMatrix3<float>::createRotationMatrix(
                t3, DEG2RAD * getBOND_ANGLE_AT_CPRIME_TO_N())
                * (t0 - t1);
        endDirection = endDirection.normalize() * getBOND_LENGTH_CPRIME_TO_N();

        endNormal = t3;
        midPoint = t4;
        midDirection = t4;
        midNormal = t4;
    }

    /**
     *@Description Rotates this along axis by angle. Axis can be any valid
     *    vector, angle is indicated in radians.
     *@param  reference to the axis(const vgVector3<float>&), angle to rotate(const double)
     *@return  changes are made internally(void)
     */
    void LoopTableEntry::rotate(const vgVector3<float>& axis, const double _angle) {
        vgMatrix3<float> rotationMatrix =
                vgMatrix3<float>::createRotationMatrix(
                const_cast<vgVector3<float>&> (axis).normalize(), _angle);

        endPoint = rotationMatrix * endPoint;
        endDirection = rotationMatrix * endDirection;
        endNormal = rotationMatrix * endNormal;
        midPoint = rotationMatrix * midPoint;
        midDirection = rotationMatrix * midDirection;
        midNormal = rotationMatrix * midNormal;
    }

    /**
     *@Description rotates this by angle along axis
     *@param  reference to the rotation matrix( const vgMatrix3<float>&)
     *@return  changes are made internally(void)
     */
    void LoopTableEntry::rotate(const vgMatrix3<float>& rotationMatrix) {
        endPoint = rotationMatrix * endPoint;
        endDirection = rotationMatrix * endDirection;
        endNormal = rotationMatrix * endNormal;
        midPoint = rotationMatrix * midPoint;
        midDirection = rotationMatrix * midDirection;
        midNormal = rotationMatrix * midNormal;
    }

    /**
     *@Description Concatenates the table entry this with dest, the result
     *    is returned as function value. dest is aligned with this before 
     *    concatenation. The Omega torsion angle is to be twisted considered 
     *    for even nChain starting points of src.
     *@param  reference to the loop table entry(const LoopTableEntry&),number of chains(unsigned int)
     *@return  the resulting loop table entry(LoopTableEntry)
     */
    LoopTableEntry LoopTableEntry::concatenate(const LoopTableEntry& dest, unsigned int nChain) {

        PRECOND(nChain >= 2, exception);

        const LoopTableEntry tSource = *this;
        LoopTableEntry tDest = dest;

        vgVector3<float> FPdest(0, getBOND_LENGTH_N_TO_CALPHA(), 0);
        // the first point of dest

        // 1.A. for the two tables to be concatenated, tDest must be aligned 
        // (ie. same endNormal) to tSource. If tSource.endNormal is unchanged it 
        // points to (0,0,1)

        vgVector3<float> refNormalDest(0, 0, 1); // the reference normal

        if (nChain % 2 != 0)
            refNormalDest = -refNormalDest;
        // zNormal has to be reverted to compensate Omega angle, see 2. below

        // calculate the angle between both normals:
        vgMatrix3<float> rotationMatrix;
        alignVectors(tSource.endNormal, refNormalDest, rotationMatrix);

        refNormalDest = rotationMatrix * refNormalDest;

        tDest.rotate(rotationMatrix); // update tDest 
        FPdest = rotationMatrix * FPdest; // ...and its first point 

        // 1.B. create the reference direction ("should-be") for tDest's first point:

        rotationMatrix = vgMatrix3<float>::createRotationMatrix(
                tSource.endNormal, DEG2RAD * ((nChain % 2 != 0)
                ? - getBOND_ANGLE_AT_N_TO_CALPHA()
                : getBOND_ANGLE_AT_N_TO_CALPHA()));

        vgVector3<float> tempEndDir = -tSource.endDirection;
        tempEndDir = rotationMatrix * tempEndDir;
        tempEndDir.normalize();

        // 1.C. calculate the angle between tempEndDir ("should-be") and 
        // yNormal ("is") and rotate tDest to match the former:

        alignVectors(tempEndDir, FPdest, rotationMatrix);
        tDest.rotate(rotationMatrix);

        //  2. apply omega == 180 degrees, where needed:
        // omega has to be applied on every second entry, because the repeated 
        // rotations compensate each other.

        if (nChain % 2 == 0)
            tDest.rotate(
                -const_cast<vgVector3<float>&> (tSource.endDirection).normalize(),
                DEG2RAD * 180.0);
        // nChain % 2 == 0 is correct because it refers to *not* flipping over uneven
        // aminoacids, ie. 1st, 3rd, 5th, etc.

        // 3. translate tDest by tSource's length (ep + ed) and set 'mid'-data

        tDest.endPoint += (tSource.endPoint + tSource.endDirection);

        tDest.midPoint = tSource.endPoint;
        tDest.midDirection = tSource.endDirection;
        tDest.midNormal = tSource.endNormal;

        return tDest;
    }

    /**
     *@Description setToOrigin is the dual method to concatenate. It reverses the 
     *    concatenation held there, ie. the source is 'subtracted' from 
     *    destination.
     *@param  reference to the loop table entry(const LoopTableEntry&),number of chains(unsigned int)
     *@return   the resulting loop table entry(LoopTableEntry)
     */
    LoopTableEntry LoopTableEntry::setToOrigin(LoopTableEntry& dest, unsigned int nChain,
            VectorTransformation& vt) {

        PRECOND(nChain >= 1, exception);

        LoopTableEntry source = *this;
        LoopTableEntry result = dest;
        vgMatrix3<float> rotationMatrix;

        //  1. translate destination back to 'origin' by  source's length (ep + ed):

        result.endPoint -= (source.endPoint + source.endDirection);
        vt.addTrans(source.endPoint + source.endDirection);

        //  2. apply omega == 180 degrees, where needed:
        // omega has to be applied on every second entry, because the repeated 
        // rotations compensate each other.

        if ((nChain % 2) == 0) {
            if (source.endDirection.length() == 0)
                ERROR("LoopTableEntry::setToOrigin() : source.eD has length 0",
                    exception);

            result.rotate(-source.endDirection.normalize(), -(DEG2RAD * 180.0));
            vt.addRot(vgMatrix3<float>::createRotationMatrix(
                    -source.endDirection.normalize(),
                    (DEG2RAD * 180.0)));
        }

        // 3.A. for the second table to be set to origin, destination must be aligned
        // (ie. same normal) to back to (0,0,1). If source.normal is unchanged it 
        // points to (0,0,1) destination's normal indicating the 'Alignment' is 
        // identical to source.normal

        vgVector3<float> zNormal(0, 0, 1); // the reference normal
        if (nChain % 2 != 0)
            zNormal = -zNormal; // .....

        alignVectors(zNormal, source.endNormal, rotationMatrix);

        vt.addAlignVectors(source.endNormal, zNormal);

        result.rotate(rotationMatrix);
        source.endDirection = rotationMatrix * source.endDirection;
        source.endNormal = rotationMatrix * source.endNormal;

        // 3.B. create the reference direction ("should-be") for destination's first
        // point:

        rotationMatrix = vgMatrix3<float>::createRotationMatrix(source.endNormal,
                DEG2RAD * ((nChain % 2 != 0)
                ? - getBOND_ANGLE_AT_N_TO_CALPHA()
                : getBOND_ANGLE_AT_N_TO_CALPHA()));

        vgVector3<float> tempEndDir = -source.endDirection;
        tempEndDir = rotationMatrix * tempEndDir;
        tempEndDir.normalize();

        // 3.C. calculate the angle between tempEndDir ("should-be") and 
        // yNormal ("is") and rotate destination to match the former:

        vgVector3<float> yNormal(0, 1, 0);

        alignVectors(yNormal, tempEndDir, rotationMatrix);

        vt.addAlignVectors(tempEndDir, yNormal);

        result.rotate(rotationMatrix);

        return result;
        // NB: 'mid'-data for result is invalid.
    }

    /**
     *@Description  Rotates the endPoint of this back into the X,Y plane and
     *    updates endDirection, normal accordingly.
     *    rotate EP into the x,y plane (and update ED & N accordingly) by using 
     *    the arccos of the x-element of EP as reference for the required angle
     *@param  reference to the transformation vector(VectorTransformation&)
     *@return  corresponding value
     */
    float LoopTableEntry::rotateIntoXYPlane(VectorTransformation& vt) {
        double angle = -acos(endPoint.x
                / sqrt(endPoint.x * endPoint.x + endPoint.z * endPoint.z));

        // NB: this formula is correct (!) because we are using the y-axis as normal
        // the problem is hence reduced to the 2-D case of rotating (x,z) to (x,0)

        // check for axis symmetry
        if (endPoint.z > 0)
            angle = -angle;

        vgVector3<float> axis(0, 1, 0);

        rotate(axis, angle);
        vgMatrix3<float> rotMat =
                vgMatrix3<float>::createRotationMatrix(axis, -angle);

        POSTCOND(fabs(endPoint.z) < EPSILON, exception);
        POSTCOND(endPoint.x >= -EPSILON, exception);

        vt.addRot(rotMat);
        endPoint.z = 0.0;
        return angle;
    }


    // OPERATORS:

    /**
     *@Description Assigns one original Entry table into another
     *@param  the original loop entry table reference(const LoopTableEntry& )
     *@return   reference to the resulting loop table entry(LoopTableEntry)
     */
    LoopTableEntry& LoopTableEntry::operator=(const LoopTableEntry& orig) {
        if (&orig != this)
            copy(orig);
        return *this;
    }
    /**
     *@Description Obtains the corresponding value depending of the index
     *@param  index (unsigned int)
     *@return  reference to the resulting vector(vgVector3<float>&)
     */
    vgVector3<float>& LoopTableEntry::operator[](unsigned int i) {
        switch (i) {
            case 0: return endPoint;
            case 1: return endDirection;
            case 2: return endNormal;
            case 3: return midPoint;
            case 4: return midDirection;
            case 5: return midNormal;
            default: ERROR("Invalid argument for operator [].", exception);
        }
        return endPoint;
    }
    /**
     *@Description Obtains the corresponding value depending of the index
     *@param  index (unsigned int)
     *@return  reference to the resulting vector(vgVector3<float>&)
     */
    const vgVector3<float>& LoopTableEntry::operator[](unsigned int i) const {
        switch (i) {
            case 0: return endPoint;
            case 1: return endDirection;
            case 2: return endNormal;
            case 3: return midPoint;
            case 4: return midDirection;
            case 5: return midNormal;
            default: ERROR("Invalid argument for operator [].", exception);
        }
        return endPoint;
    }


    // TESTERS:

    /**
     *@Description Prints an entry data
     *@param  index of the entry
     *@return  the output is printed
     */
    void LoopTableEntry::printTable(unsigned int k) const {
        cout << "Entry " << setw(4) << k << "\t EP: " << setw(5) << setprecision(4)
                << endPoint.x << "\t ED: " << setw(5) << setprecision(4)
                << endDirection.x << "\t N: " << setw(5) << setprecision(4)
                << endNormal.x << "\t MP: " << setw(5) << setprecision(4)
                << midPoint.x << "\t MD: " << setw(5) << setprecision(4)
                << midDirection.x << "\t MN: " << setw(5) << setprecision(4)
                << midNormal.x
                << "\n"
                << "\t\t EP: " << setw(5) << setprecision(4) << endPoint.y
                << "\t ED: " << setw(5) << setprecision(4) << endDirection.y
                << "\t N: " << setw(5) << setprecision(4) << endNormal.y
                << "\t MP: " << setw(5) << setprecision(4) << midPoint.y
                << "\t MD: " << setw(5) << setprecision(4) << midDirection.y
                << "\t MN: " << setw(5) << setprecision(4) << midNormal.y
                << "\n"
                << "\t\t EP: " << setw(5) << setprecision(4) << endPoint.z
                << "\t ED: " << setw(5) << setprecision(4) << endDirection.z
                << "\t N: " << setw(5) << setprecision(4) << endNormal.z
                << "\t MP: " << setw(5) << setprecision(4) << midPoint.z
                << "\t MD: " << setw(5) << setprecision(4) << midDirection.z
                << "\t MN: " << setw(5) << setprecision(4) << midNormal.z
                << "\n";
    }

}}
