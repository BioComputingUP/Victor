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
#include <EnergyFeatures.h>
#include <AminoAcidCode.h>
#include <Spacer.h>

using namespace Victor;

// Global constants, typedefs, etc. (to avoid):

// CONSTRUCTORS/DESTRUCTOR:

/**
 *  Basic constructor 
 */
EnergyFeatures::EnergyFeatures() : solv(), rapdf() {
    tors = new PhiPsi(10);
    tors2 = new PhiPsiOmegaChi1Chi2(10);
}

// PREDICATES:

/**
 *  print the possible output features
 */
void EnergyFeatures::showFeatures() {
    cout << "The ouput features (15 per protein model) are:\n"
            << "\t RAPDF energy\n"
            << "\t SOLV energy\n"
            << "\t HYDB energy\n"
            << "\t TORS energy\n"
            << "\t TAP score - normalized TORS score\n"
            << "\t length\n"
            << "\t secondary structure composition (3 values)\n"
            << "\t chain breaks, ie. Ca-Ca > 4.5 A at < 7.5, < 10, < 15, < 20, "
            << ">20 A (5 values)\n"
            << "\t hard sphere backbone Ca-Ca clashes, ie. distance < 2.75 A\n"
            << endl;
}

/**
 *    calculates the quantity of backbone hidrogen Bonds from the spacer 
 *@param reference of a Spacer(Spacer&)
 *@return  quantity of hydrogen bonds (double)
 */
double EnergyFeatures::calculateBackboneHydrogenBonds(Spacer& sp) {
    double count = 0;
    for (int i = 0; i < (int) sp.sizeAmino(); i++)
        for (int j = 0; j < (int) sp.sizeAmino(); j++)
            if (abs(i - j) > 1) {
                double dist = sp.getAmino(i)[N].distance(sp.getAmino(j)[O]);
                double dist2 = sp.getAmino(i)[N].distance(sp.getAmino(j)[C]);
                double dist3 = sp.getAmino(i)[CA].distance(sp.getAmino(j)[O]);

                if ((dist <= 4.0) && (dist >= 2.0) && (dist < dist2)
                        && (dist < dist3)) {
                    count++;
                    break;
                }
            }
    return count;
}

/**
 *  Returns the calculated  Aa Composition for all the aa in the spacer
 *@param reference of a Spacer(spacer&)
 *@return  a vector containing all the compositions for each of the aas (vector double)
 */
vector<double> EnergyFeatures::calculateAaComposition(Spacer& sp) {
    vector<double> aa;
    for (unsigned int i = 0; i < 20; i++)
        aa.push_back(0.0);

    for (unsigned int i = 0; i < sp.sizeAmino(); i++) {
        aa[static_cast<unsigned int> (sp.getAmino(i).getCode())]++;
    }

    for (unsigned int i = 0; i < 20; i++)
        aa[i] /= sp.sizeAmino();

    return aa;
}

/**
 *  Returns the calculated  Secondary Composition from the amino acids in the spacer
 *@param reference of a Spacer(Spacer&)
 *@return a vector of 3 containing the corresponding values (vector double)
 */

vector<double> EnergyFeatures::calculateSecondaryComposition(Spacer& sp) {
    vector<double> sec;
    for (unsigned int i = 0; i < 3; i++)
        sec.push_back(0.0);

    for (unsigned int i = 0; i < sp.sizeAmino(); i++) {
        if (sp.getAmino(i).getState() == HELIX)
            sec[0]++;
        else if (sp.getAmino(i).getState() == STRAND)
            sec[1]++;
        else
            sec[2]++;
    }
    for (unsigned int i = 0; i < 3; i++)
        sec[i] /= sp.sizeAmino();

    return sec;
}

/**
 *  Returns the calculated  Meso State Composition for the amino acids in the spacer
 *@param reference of a Spacer(Spacer&)
 *@return  vector of 36 containing the corresponding values (vector double)
 */

vector<double> calculateMesoStateComposition(Spacer& sp) {
    vector<double> meso;
    for (unsigned int i = 0; i < 36; i++)
        meso.push_back(0.0);
    unsigned int x, y;
    for (unsigned int i = 0; i < sp.sizeAmino(); i++) {
        if (sp.getAmino(i).getPhi() < -120.0)
            x = 0;
        else if (sp.getAmino(i).getPhi() < -60.0)
            x = 1;
        else if (sp.getAmino(i).getPhi() < 0.0)
            x = 2;
        else if (sp.getAmino(i).getPhi() < 60.0)
            x = 3;
        else if (sp.getAmino(i).getPhi() < 120.0)
            x = 4;
        else
            x = 5;

        if (sp.getAmino(i).getPsi() < -120.0)
            y = 0;
        else if (sp.getAmino(i).getPsi() < -60.0)
            y = 1;
        else if (sp.getAmino(i).getPsi() < 0.0)
            y = 2;
        else if (sp.getAmino(i).getPsi() < 60.0)
            y = 3;
        else if (sp.getAmino(i).getPsi() < 120.0)
            y = 4;
        else
            y = 5;

        meso[x * 6 + y]++;
    }

    for (unsigned int i = 0; i < 36; i++)
        meso[i] /= sp.sizeAmino();

    return meso;
}

/**
 *  Returns the calculated  chain breaks for the amino acids in the spacer
 *@param reference of a Spacer(Spacer&)
 *@return  vector of 5 containing the corresponding values (vector double)
 */
vector<double> EnergyFeatures::calculateChainBreaks(Spacer& sp) {
    vector<double> chain;
    for (unsigned int i = 0; i < 5; i++)
        chain.push_back(0.0);

    for (unsigned int i = 0; i < sp.sizeAmino() - 1; i++)
        if (sp.getAmino(i).getAtom(CA).distance(sp.getAmino(i + 1).getAtom(CA)) > 4.5) {
            if (sp.getAmino(i).getAtom(CA).distance(sp.getAmino(i + 1).getAtom(CA)) > 20.0)
                chain[4]++;
            else if (sp.getAmino(i).getAtom(CA).distance(sp.getAmino(i + 1).getAtom(CA)) > 15.0)
                chain[3]++;
            else if (sp.getAmino(i).getAtom(CA).distance(sp.getAmino(i + 1).getAtom(CA)) > 10.0)
                chain[2]++;
            else if (sp.getAmino(i).getAtom(CA).distance(sp.getAmino(i + 1).getAtom(CA)) > 7.5)
                chain[1]++;
            else
                chain[0]++;
        }
    return chain;
}

/**
 *  Returns the calculated  the clashes for the amino acids in the spacer
 *@param reference of a Spacer(Spacer&)
 *@return  double containing the corresponding value ( double)
 */
double EnergyFeatures::calculateClashes(Spacer& sp) {
    double clash = 0.0;
    for (unsigned int i = 0; i < sp.sizeAmino() - 2; i++)
        for (unsigned int j = i + 2; j < sp.sizeAmino(); j++) {
            if (sp.getAmino(i).getAtom(CA).distance(sp.getAmino(j).getAtom(CA))
                    < 2.75)
                clash++;
        }
    return clash;
}

/**
 *  Returns the calculated  all the features for the amino acids in the spacer
 *@param reference of a Spacer(Spacer&)
 *@return  vector  containing the corresponding values (vector double)
 */
vector<double> EnergyFeatures::calculateFeatures(Spacer& sp) {
    //Calls to the other functions to calculate all the other values  
    vector<double> tmp;

    tmp.push_back(rapdf.calculateEnergy(sp));
    tmp.push_back(solv.calculateSolvation(sp));
    tmp.push_back(-1.0
            * EnergyFeatures::calculateBackboneHydrogenBonds(sp));
    tmp.push_back(tors->calculateEnergy(sp));
    tmp.push_back(tors2->calculateEnergy(sp)
            / (tors2->calculateMaxEnergy(sp)));

    tmp.push_back(sp.sizeAmino());

    vector<double> secComposition =
            EnergyFeatures::calculateSecondaryComposition(sp);
    for (unsigned int i = 0; i < secComposition.size(); i++)
        tmp.push_back(secComposition[i]);

    vector<double> chainBreak = EnergyFeatures::calculateChainBreaks(sp);
    for (unsigned int i = 0; i < chainBreak.size(); i++)
        tmp.push_back(chainBreak[i]);

    tmp.push_back(EnergyFeatures::calculateClashes(sp));

    return tmp;
}

// MODIFIERS:


// HELPERS:

