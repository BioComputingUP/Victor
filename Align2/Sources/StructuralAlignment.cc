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
// -*- C++ -*------------------------------------------------------------------

//
//  Description:
//    Class for structural alignments.
//
// ---------------------------------------------------------------------------

// Includes:
#include <StructuralAlignment.h>
#include <String2Number.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;

StructuralAlignment::StructuralAlignment() : AlignmentBase(), spTarget(),
spTemplate(), equivData(), fragData() {
    PRINT_NAME;
}

StructuralAlignment::StructuralAlignment(const StructuralAlignment& orig) {
    PRINT_NAME;
    this->copy(orig);
}

StructuralAlignment::~StructuralAlignment() {
    PRINT_NAME;
}

// MODIFIERS:
    /**
     * @description
     * @param input
     * @param spNew
     */
void
StructuralAlignment::loadCE(istream& input, Spacer& spNew) {
    if (spNew.sizeAmino() == 0)
        ERROR("PDB structure is invalid.", exception);

    spTemplate = spNew;

    vgMatrix3<double> rot(1);
    vgVector3<double> trans(0, 0, 0);

    while (input) {
        string tmp;
        input >> tmp;
        if (tmp == "X2") {
            // get transformation data
            string tmp2;

            tmp2 = readLine(input);
            rot[0] = stod(tmp2.substr(3, 9));
            rot[3] = stod(tmp2.substr(20, 9));
            rot[6] = stod(tmp2.substr(37, 9));
            trans[0] = stod(tmp2.substr(56, 10));

            input >> tmp;
            if (tmp != "Y2")
                ERROR("Unexpected end of CE file.", exception);

            tmp2 = readLine(input);
            rot[1] = stod(tmp2.substr(3, 9));
            rot[4] = stod(tmp2.substr(20, 9));
            rot[7] = stod(tmp2.substr(37, 9));
            trans[1] = stod(tmp2.substr(56, 10));

            input >> tmp;
            if (tmp != "Z2")
                ERROR("Unexpected end of CE file.", exception);

            tmp2 = readLine(input);
            rot[2] = stod(tmp2.substr(3, 9));
            rot[5] = stod(tmp2.substr(20, 9));
            rot[8] = stod(tmp2.substr(37, 9));
            trans[2] = stod(tmp2.substr(56, 10));

            cout << "\t" << rot[0] << "\t" << rot[3] << "\t" << rot[6]
                    << "\t" << trans[0] << "\n";
            cout << "\t" << rot[1] << "\t" << rot[4] << "\t" << rot[7]
                    << "\t" << trans[1] << "\n";
            cout << "\t" << rot[2] << "\t" << rot[5] << "\t" << rot[8]
                    << "\t" << trans[2] << "\n";

        } else
            skipToNewLine(input);

    }

    pExecStructAli(rot, trans);
}
/**
 * @description
 */
void
StructuralAlignment::buildEquivalenceNetwork() {
    unsigned int minIndex = 0;

    for (unsigned int i = 0; i < spTarget.sizeAmino(); i++) {
        double minD = 9999.9;
        double minD2 = 9999.9;
        unsigned int oldIndex = minIndex;
        minIndex = 0;

        for (unsigned int j = oldIndex; j < spTemplate.sizeAmino(); j++) {
            double tmpD = spTarget.getAmino(i)[CA].distance(
                    spTemplate.getAmino(j)[CA]) + 0.1 * (j - oldIndex);
            if (tmpD < minD) {
                minD = tmpD;
                minD2 = spTarget.getAmino(i)[CA].distance(
                        spTemplate.getAmino(j)[CA]);
                minIndex = j;
            }
        }
        EData tmpE(minIndex, minD2);
        equivData.push_back(tmpE);
    }

    cout << endl;

}
/**
 * @description
 * @param maxDist
 */
void
StructuralAlignment::buildFragmentNetwork(double maxDist) {
    unsigned int minE = 9999;
    //      unsigned int maxE = 9999;
    //      unsigned int currE = 9999;

    for (unsigned int i = 0; i < equivData.size(); i++) {
        if (equivData[i].dist <= maxDist) { // it's ok
            if (minE < 9999) // min open
                continue;
            else // min closed
            {
                minE = i;
            }
        } else { // it's not ok
            if (minE < 9999) // min open
            {
                if (minE != i - 1) // insert fragment only if longer than 1 AA
                {
                    double tmpD = 0;
                    for (unsigned int j = minE; j < i; j++)
                        tmpD += equivData[j].dist;
                    tmpD /= i - minE;

                    FData fd(minE, i - 1, tmpD);
                    fragData.push_back(fd);
                }

                minE = 9999;
            } else // min closed
                continue;

        }
    }

    if (minE < 9999) { // close last fragment
        double tmpD = 0;
        for (unsigned int j = minE; j < equivData.size() - 1; j++)
            tmpD += equivData[j].dist;
        tmpD /= equivData.size() - 1 - minE;

        FData fd(minE, equivData.size() - 1, tmpD);
        fragData.push_back(fd);
    }
}
/**
 * @description
 * @param maxDist
 */
void
StructuralAlignment::writeData(double maxDist) {
    // first, write equivalence data
    cout << "   0\t";

    for (unsigned int i = 0; i < equivData.size(); i++) {
        if (equivData[i].dist <= maxDist)
            cout << setw(5) << setprecision(2) << equivData[i].dist << " ("
                << setw(3) << equivData[i].other << ")\t";
        else
            cout << "  xxx    \t";

        if ((i + 1) % 10 == 0)
            cout << "\n" << setw(4) << i + 1 << "\t";
    }

    cout << "\n\n";

    // now write fragment data
    for (unsigned int i = 0; i < fragData.size(); i++) {
        cout << setw(3) << i << "\t" << setw(4) << fragData[i].min
                << "\t" << setw(4) << fragData[i].max << "\t" << setw(5)
                << setprecision(3) << fragData[i].dist << "\n";
    }
}
/**
 * @description
 * @param orig
 */
void
StructuralAlignment::copy(const StructuralAlignment& orig) {
    PRINT_NAME;

    AlignmentBase::copy(orig);

    spTarget = orig.spTarget;
    spTemplate = orig.spTemplate;

    equivData.clear();
    for (unsigned int i = 0; i < orig.equivData.size(); i++)
        equivData.push_back(orig.equivData[i]);

    fragData.clear();
    for (unsigned int i = 0; i < orig.fragData.size(); i++)
        fragData.push_back(orig.fragData[i]);

}


// HELPERS:
/**
 * @description
 * @param rot
 * @param trans
 */
void
StructuralAlignment::pExecStructAli(vgMatrix3<double> rot, vgVector3<double> trans) {
    for (unsigned int i = 0; i < spTemplate.sizeAmino(); i++)
        for (unsigned int j = 0; j < spTemplate.getAmino(i).size(); j++) {
            vgVector3<double> tmpV = spTemplate.getAmino(i)[j].getCoords();
            vgVector3<double> tmpV2;

            tmpV2[0] = rot[0] * tmpV[0] + rot[3] * tmpV[1]
                    + rot[6] * tmpV[2] + trans[0];
            tmpV2[1] = rot[1] * tmpV[0] + rot[4] * tmpV[1]
                    + rot[7] * tmpV[2] + trans[1];
            tmpV2[2] = rot[2] * tmpV[0] + rot[5] * tmpV[1]
                    + rot[8] * tmpV[2] + trans[2];

            spTemplate.getAmino(i)[j].setCoords(tmpV2);

        }

}

