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

#include <VGPFunction2.h>

namespace Biopool {

    // CONSTRUCTORS:
    /**
     * @description
     * @param secFileName
     */
    VGPFunction2::VGPFunction2(string secFileName) : o(14.00), e(1.00), extType(0),
    extCounter(0), wH(1.00), wS(1.00) {
        pExtractSecInfo(secFileName);
    }
    /**
     * @description
     * @param secFileName
     * @param o
     * @param e
     * @param extType
     * @param wH
     * @param wS
     */
    VGPFunction2::VGPFunction2(string secFileName, double o, double e,
            unsigned int extType, double wH, double wS) : o(o), e(e), extType(extType),
    extCounter(0), wH(wH), wS(wS) {
        pExtractSecInfo(secFileName);
    }

    VGPFunction2::VGPFunction2(const VGPFunction2 &orig) : GapFunction(orig) {
        copy(orig);
    }

    VGPFunction2::~VGPFunction2() {
    }


    // OPERATORS:
    /**
     * @description
     * @param orig
     * @return 
     */
    VGPFunction2&
            VGPFunction2::operator =(const VGPFunction2 &orig) {
        if (&orig != this)
            copy(orig);
        POSTCOND((orig == *this), exception);
        return *this;
    }


    // PREDICATES:
    /**
     * @description
     * @param p
     * @return 
     */
    double
    VGPFunction2::getOpenPenalty(int p) {
        double penH = hContent[p - 1];
        double penS = sContent[p - 1];

        extCounter = 0;
        return o + wH * penH + wS * penS;
    }
    /**
     * @description
     * @param p
     * @return 
     */
    double
    VGPFunction2::getExtensionPenalty(int p) {
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
    /**
     * @description
     * @param orig
     */
    void
    VGPFunction2::copy(const VGPFunction2 &orig) {
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

        wH = orig.wH;
        wS = orig.wS;
    }
    /**
     * @description
     * @return 
     */
    VGPFunction2*
    VGPFunction2::newCopy() {
        VGPFunction2 *tmp = new VGPFunction2(*this);
        return tmp;
    }


    // HELPERS:
    /**
     * @description
     * @param secFileName
     */
    void
    VGPFunction2::pExtractSecInfo(string secFileName) {
        // --------------------------------------------------
        // 1. Load secondary structure FASTA file
        // --------------------------------------------------

        ifstream secFile(secFileName.c_str());
        if (!secFile)
            ERROR("Error opening secondary structure FASTA file.", exception);
        Alignment aliSec;
        aliSec.loadFasta(secFile);
        if (aliSec.size() < 1)
            ERROR("Secondary structure FASTA file must contain two sequences.", exception);
        string sec = Alignment::getPureSequence(aliSec.getTemplate());


        // --------------------------------------------------
        // 2. Calculate structural context from infos
        // --------------------------------------------------

        for (unsigned int k = 0; k < sec.length(); k++)
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
                };
    }

} // namespace
