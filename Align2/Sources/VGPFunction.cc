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
// --*- C++ -*------x-----------------------------------------------------------
//

//
// Description:     Implement VGP (Variable Gap Penalty) function. Some
//                  explanations can be found in:
//
//                  Madhusudhan MS., Marti-Renom MA., Sanchez R., Sali A.
//                  Variable gap penalty for protein sequence-structure alignment.
//                  Department of Biopharmaceutical Sciences and Pharmaceutical
//                  Chemistry, University of California at San Francisco, 94143, USA.
//                  PMID: 16423846 [PubMed - indexed for MEDLINE]
//
// -----------------x-----------------------------------------------------------

#include <VGPFunction.h>

namespace Biopool {

    // CONSTRUCTORS:

    VGPFunction::VGPFunction(string pdbFileName, string chainID) : o(14.00), e(1.00), extType(0),
    extCounter(0), wH(1.00), wS(1.00), wB(1.00), wC(1.00), wD(1.00) {
        pExtractPdbInfo(pdbFileName, chainID);
    }

    VGPFunction::VGPFunction(string pdbFileName, string chainID, double o, double e,
            unsigned int extType, double wH, double wS, double wB, double wC, double wD)
    : o(o), e(e), extType(extType), extCounter(0), wH(wH), wS(wS), wB(wB),
    wC(wC), wD(wD) {
        pExtractPdbInfo(pdbFileName, chainID);
    }

    VGPFunction::VGPFunction(const VGPFunction &orig) : GapFunction(orig) {
        copy(orig);
    }

    VGPFunction::~VGPFunction() {
    }


    // OPERATORS:

    VGPFunction&
            VGPFunction::operator =(const VGPFunction &orig) {
        if (&orig != this)
            copy(orig);
        POSTCOND((orig == *this), exception);
        return *this;
    }


    // PREDICATES:

    double
    VGPFunction::getOpenPenalty(int p) {
        double penH = hContent[p - 1];
        double penS = sContent[p - 1];
        double penB = 1 - solvAccess[p - 1];

        double penC = 0.00;
        if ((penH == 1.00) || (penS == 1.00))
            penC = 1.00;
        else
            penC = 1 - (min(180.00, max(0.00, bbStraight[p - 1])) / 180.00);

        const double GAMMA = 2.00;
        const double D0 = 3.00;
        double penD = pow(max(0.00, (spaceProx[p - 1] - D0)), GAMMA);

        extCounter = 0;
        return o + wH * penH + wS * penS + wB * penB + wC * penC + wD * penD;
    }

    double
    VGPFunction::getExtensionPenalty(int p) {
        const double STEP1 = 10.00;
        const double STEP2 = 1.00;

        extCounter++;

        if (extType == 1)
            return e - min(e, (e / STEP1 * (extCounter - 1)));
        if (extType == 2)
            return e * pow(1 / e, (extCounter * STEP2));

        return e;
    }


    // MODIFIERS:

    void
    VGPFunction::copy(const VGPFunction &orig) {
        GapFunction::copy(orig);
        o = orig.o;
        e = orig.e;
        extType = orig.extType;

        hContent.clear();
        for (unsigned int i = 0; i < orig.hContent.size(); i++)
            hContent.push_back(orig.hContent[i]);

        sContent.clear();
        for (unsigned int i = 0; i < orig.sContent.size(); i++)
            sContent.push_back(orig.sContent[i]);

        solvAccess.clear();
        for (unsigned int i = 0; i < orig.solvAccess.size(); i++)
            solvAccess.push_back(orig.solvAccess[i]);

        bbStraight.clear();
        for (unsigned int i = 0; i < orig.bbStraight.size(); i++)
            bbStraight.push_back(orig.bbStraight[i]);

        spaceProx.clear();
        for (unsigned int i = 0; i < orig.spaceProx.size(); i++)
            spaceProx.push_back(orig.spaceProx[i]);

        wH = orig.wH;
        wS = orig.wS;
        wB = orig.wB;
        wC = orig.wC;
        wD = orig.wD;
    }

    VGPFunction*
    VGPFunction::newCopy() {
        VGPFunction *tmp = new VGPFunction(*this);
        return tmp;
    }


    // HELPERS:

    void
    VGPFunction::pExtractPdbInfo(string pdbFileName, string chainID) {
        // --------------------------------------------------
        // 1. Load PDB file
        // --------------------------------------------------

        ifstream pdbFile(pdbFileName.c_str());
        if (!pdbFile)
            ERROR("Error opening template PDB file.", exception);
        PdbLoader pdb(pdbFile);
        pdb.setNoHAtoms();
        if (chainID != " ") {
            if (chainID.size() > 1)
                ERROR("You can choose only 1 chain", error);
            pdb.setChain(chainID[0]);
        }
        Spacer *sp;



        pdb.setChain(chainID[0]);


        Protein prot;
        prot.load(pdb);
        sp = prot.getSpacer(chainID[0]);

        // --------------------------------------------------
        // 2. Extract structural infos from PDB file
        // --------------------------------------------------

        string sec;
        vector< vgVector3<double> > CaCoords;

        for (unsigned int k = 0; k < sp->sizeAmino(); k++) {
            switch (sp->getAmino(k).getState()) {
                case HELIX:
                    sec += 'H';
                    break;
                case STRAND:
                    sec += 'E';
                    break;
                default:
                    sec += 'C';
                    break;
            };

            CaCoords.push_back(sp->getAmino(k)[CA].getCoords());
        }


        // --------------------------------------------------
        // 3. Calculate structural context from infos
        // --------------------------------------------------

        //TODO improve default assignment for unseen residues
        //define SA that will be used in case of gaps.
        //double previousSolvAcc = 0.5;

        for (unsigned int k = 0; k < sp->sizeAmino(); k++) {
            // Helical/strand content

            if (k == 0) {
                hContent.push_back(0.00);
                sContent.push_back(0.00);
            } else
                switch (sec[k]) {
                    case 'H':
                        sContent.push_back(0.00);
                        if (sec[k - 1] == 'H')
                            hContent.push_back(1.00);
                        else
                            hContent.push_back(0.00);
                        break;
                    case 'E':
                        hContent.push_back(0.00);
                        if (sec[k - 1] == 'E')
                            sContent.push_back(1.00);
                        else
                            sContent.push_back(0.00);
                        break;
                    default:
                        hContent.push_back(0.00);
                        sContent.push_back(0.00);
                        break;
                };

            // Solvent accessibility
            try {
                double solvAcc = getSolvAccess(*sp, k, 0, sp->sizeAmino());
                solvAccess.push_back(solvAcc);
                //previousSolvAcc = solvAcc;
            } catch (const char* exc) {
                ERROR("In the template there are missing atoms, fix it or use try another gap penalty function", exception);

            }

            // Backbone straightness

            unsigned int start = (k < 3) ? 0 : k - 3;
            unsigned int end = (k >= (sp->sizeAmino() - 3)) ? sp->sizeAmino() - 1 : k + 3;
            vgVector3<double> vec1 = CaCoords[k] - CaCoords[start];
            vgVector3<double> vec2 = CaCoords[end] - CaCoords[k];
            double mod1 = sqrt(vec1.square());
            double mod2 = sqrt(vec2.square());
            double dotP = vec1 * vec2;
            double arc = (dotP != 0.00) ? acos(dotP / (mod1 * mod2)) : acos(0.00);

            bbStraight.push_back(arc * 180.00 / 3.14);


            // Space proximity

            if (k == 0)
                spaceProx.push_back(0.00);
            else {
                vgVector3<double> vec = CaCoords[k] - CaCoords[k - 1];
                double mod = sqrt(vec.square());

                spaceProx.push_back(mod);
            }
        }

    }

} // namespace
