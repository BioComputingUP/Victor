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

// -*- C++ -*-------x-----------------------------------------------------------
//
//
// Description:     Calculate a frequency profile or PSSM using Henikoff
//                  weighting scheme. Some explanations can be found in:
//
//                  Guoli Wang, Roland L. Dunbrack jr.
//                  Scoring profile-to-profile sequence alignments.
//                  Institute for Cancer Research, Fox Chase Cancer Center,
//                  Philadelphia, Pennsylvania 19111, USA. March 16, 2004.
//
// -----------------x-----------------------------------------------------------

#include <HenikoffProfile.h>
#include <ctime>
namespace Victor { namespace Align2{

    // CONSTRUCTORS:

    HenikoffProfile::HenikoffProfile() : Profile() {
    }

    HenikoffProfile::HenikoffProfile(const HenikoffProfile &orig) : Profile(orig) {
        copy(orig);
    }

    HenikoffProfile::~HenikoffProfile() {
    }


    // OPERATORS:
    /**
     *  
     * @param orig
     * @return 
     */
    HenikoffProfile&
            HenikoffProfile::operator =(const HenikoffProfile &orig) {
        if (&orig != this)
            copy(orig);
        POSTCOND((orig == *this), exception);
        return *this;
    }


    // MODIFIERS:
    /**
     *  
     * @param orig
     */
    void
    HenikoffProfile::copy(const HenikoffProfile &orig) {
        Profile::copy(orig);

        aliWeight.clear();
        for (unsigned int i = 0; i < orig.aliWeight.size(); i++) {
            vector<double> tmp;
            for (unsigned int j = 0; j < orig.aliWeight[0].size(); j++)
                tmp.push_back(orig.aliWeight[i][j]);
            aliWeight.push_back(tmp);
        }
    }
    /**
     *  
     * @return 
     */
    HenikoffProfile*
    HenikoffProfile::newCopy() {
        HenikoffProfile *tmp = new HenikoffProfile(*this);
        return tmp;
    }


    // HELPERS:
    /**
     * 
     * @param ali
     * @param cLen
     */
    void //suggested: use cLen=25 to save computational time
    HenikoffProfile::pCalculateWeight(Alignment &ali, unsigned int cLen) {
        const string residue_indices = "ARNDCQEGHILKMFPSTWYV";


        // --------------------------------------------------
        // 1. Calculate the positions of the first aminoacids
        //    for every sequence
        // --------------------------------------------------

        vector<unsigned int> firstAminoPos;

        firstAminoPos.push_back(0); // master sequence

        cout << "target \n" << ali.getTarget();
        cout << "\ntemplate \n" << ali.getTemplate(0) << "\n";
        for (unsigned int j = 0; j < (numSeq - 1); j++) // other sequences
            for (unsigned int i = 0; i < seqLen; i++) {

                if (ali.getTemplatePos(i, j) != '-') {
                    firstAminoPos.push_back(i);
                    break;
                }
            }

        // --------------------------------------------------
        // 2. Calculate the positions of the last aminoacids
        //    for every sequence
        // --------------------------------------------------

        vector<unsigned int> lastAminoPos;

        lastAminoPos.push_back(seqLen - 1); // master sequence

        for (unsigned int j = 0; j < (numSeq - 1); j++) // other sequences
            for (unsigned int i = seqLen; i > 0; i--)
                if (ali.getTemplatePos(i - 1, j) != '-') {
                    lastAminoPos.push_back(i - 1);
                    break;
                }


        // --------------------------------------------------
        // 3. Calculate weights for master sequence
        // --------------------------------------------------
        vector<double> wSeq;
        vector<unsigned int> seqSubset;

        for (unsigned int i = 0; i < seqLen; i++)
            if (ali.getTargetPos(i) == '-')
                wSeq.push_back(0.00);
            else {
                // Calculate the subset of sequences

                seqSubset.clear();

                seqSubset.push_back(0);
                for (unsigned int k = 0; k < (numSeq - 1); k++)
                    if (ali.getTemplatePos(i, k) != '-')
                        seqSubset.push_back(k + 1);

                unsigned int subsetSize = seqSubset.size();


                // Calculate Cleft and Cright

                unsigned int Cleft = firstAminoPos[seqSubset[0]];
                unsigned int Cright = lastAminoPos[seqSubset[0]];

                for (unsigned int k = 1; k < subsetSize; k++) {
                    if (firstAminoPos[seqSubset[k]] > Cleft)
                        Cleft = firstAminoPos[seqSubset[k]];
                    if (lastAminoPos[seqSubset[k]] < Cright)
                        Cright = lastAminoPos[seqSubset[k]];
                }

                double sum = 0.00;

                for (unsigned int p = Cleft; p <= Cright; p++) {
                    // Calculate the number of different aminoacids

                    unsigned int Ndiff = 0;

                    for (unsigned int index = 0; index < 20; index++) {
                        if ((seqSubset[0] == 0) && (ali.getTargetPos(p) == residue_indices[index])) {
                            Ndiff++;
                            continue;
                        }
                        for (unsigned int k = 1; k < subsetSize; k++)
                            if (ali.getTemplatePos(p, seqSubset[k] - 1) == residue_indices[index]) {
                                Ndiff++;
                                break;
                            }
                    }

                    // Calculate the number of the same aminoacid

                    unsigned int n = 0;

                    if (seqSubset[0] == 0)
                        n++;
                    for (unsigned int k = 1; k < subsetSize; k++)
                        if (ali.getTemplatePos(p, seqSubset[k] - 1) == ali.getTargetPos(p))
                            n++;


                    sum += (1 / (double) (Ndiff * n));
                }

                wSeq.push_back((1 / (double) (Cright - Cleft + 1)) * sum);
            }
        aliWeight.push_back(wSeq);

        // --------------------------------------------------
        // 4. Calculate weights for other sequences
        // --------------------------------------------------

        //4.1 calculate seqSubset, Cleft, Crigth ,Ndiff for each position, outside main loop
        //this is a time optimization of section 4. Francesco Lovo 2012
        vector <pair<int, int> > Cvalues; //Cleft and Cright for each subSeq
        vector <vector<unsigned int> > allseqSubseq;
        vector <vector<unsigned int> > allNdiff; //for each position i, we store Ndiff
        // for each column j between Cleft and Cright
        bool henikoffReduction = false;
        for (unsigned int i = 0; i < seqLen; i++) {
            // Calculate the subset of sequences
            seqSubset.clear();

            if (ali.getTargetPos(i) != '-')
                seqSubset.push_back(0);

            for (unsigned int k = 0; k < (numSeq - 1); k++)
                if (ali.getTemplatePos(i, k) != '-')
                    seqSubset.push_back(k + 1);

            unsigned int subsetSize = seqSubset.size();

            // Calculate Cleft and Cright
            unsigned int Cleft = firstAminoPos[seqSubset[0]];
            unsigned int Cright = lastAminoPos[seqSubset[0]];

            for (unsigned int k = 1; k < subsetSize; k++) {
                if (firstAminoPos[seqSubset[k]] > Cleft)
                    Cleft = firstAminoPos[seqSubset[k]];
                if (lastAminoPos[seqSubset[k]] < Cright)
                    Cright = lastAminoPos[seqSubset[k]];
            }
            //this is violation to Henikoff formula, but used to save some computational time
            if (i - Cleft > cLen) {
                (Cleft = i - cLen);
                henikoffReduction = true;
            }
            if (Cright - i > cLen) {
                (Cright = i + cLen);
                henikoffReduction = true;
            }

            //copy to global variables
            Cvalues.push_back(make_pair(Cleft, Cright));
            allseqSubseq.push_back(seqSubset);

            // Calculate the number of different aminoacids
            vector<unsigned int> pDiff;

            for (unsigned int p = Cleft; p <= Cright; p++) {
                unsigned int Ndiff = 0;

                for (unsigned int index = 0; index < 20; index++) {
                    if ((seqSubset[0] == 0) && (ali.getTargetPos(p) == residue_indices[index])) {
                        Ndiff++;
                        continue;
                    }
                    for (unsigned int k = 1; k < subsetSize; k++)
                        if (ali.getTemplatePos(p, seqSubset[k] - 1) == residue_indices[index]) {
                            Ndiff++;
                            break;
                        }
                }
                pDiff.push_back(Ndiff);
            }
            allNdiff.push_back(pDiff);
        }
        if (henikoffReduction) cout << "henikoff analisys reduced!\n"; //just a worning

        for (unsigned int j = 0; j < (numSeq - 1); j++) {
            wSeq.clear();
            for (unsigned int i = 0; i < seqLen; i++)
                if (ali.getTemplatePos(i, j) == '-')
                    wSeq.push_back(0.00);
                else {
                    unsigned int Cleft = Cvalues[i].first;
                    unsigned int Cright = Cvalues[i].second;
                    seqSubset = allseqSubseq[i];
                    unsigned int subsetSize = seqSubset.size();
                    double sum = 0.00;

                    for (unsigned int p = Cleft; p <= Cright; p++) {
                        // Calculate the number of different aminoacids
                        unsigned int Ndiff = allNdiff[i][p - Cleft];

                        // Calculate the number of the same aminoacid
                        unsigned int n = 0;
                        if ((seqSubset[0] == 0) && (ali.getTargetPos(p) == ali.getTemplatePos(p, j)))
                            n++;
                        for (unsigned int k = 1; k < subsetSize; k++)
                            if (ali.getTemplatePos(p, seqSubset[k] - 1) == ali.getTemplatePos(p, j))
                                n++;

                        sum += (1 / (double) (Ndiff * n));
                    }
                    wSeq.push_back((1 / (double) (Cright - Cleft + 1)) * sum);
                }
            aliWeight.push_back(wSeq);
        }

    }
    /**
     *  
     * @param freq
     * @param freqGap
     * @param ali
     * @param i
     */
    void
    HenikoffProfile::pCalculateRawFrequency(vector<double> &freq, double &freqGap,
            Alignment &ali, unsigned int i) {
        if (aminoAcidOneLetterTranslator(ali.getTargetPos(i)) != XXX)
            freq[aminoAcidOneLetterTranslator(ali.getTargetPos(i))] += aliWeight[0][i];
        else
            freqGap++;

        for (unsigned int j = 0; j < (numSeq - 1); j++)
            if (aminoAcidOneLetterTranslator(ali.getTemplatePos(i, j)) != XXX)
                freq[aminoAcidOneLetterTranslator(ali.getTemplatePos(i, j))] += aliWeight[j + 1][i];
            else
                freqGap++;
    }
    /**
     * 
     * @param ali
     */
    void
    HenikoffProfile::pConstructData(Alignment &ali) {
        if (!gap) {
            ali.purgeTargetInsertions();
            seqLen = ali.getTarget().size();
        }

        pCalculateWeight(ali);

        gapFreq.reserve(seqLen);
        for (unsigned int i = 0; i < seqLen; i++) {
            vector<double> tmp;
            for (AminoAcidCode i = ALA; i <= TYR; i++)
                tmp.push_back(0.00);
            profAliFrequency.push_back(tmp);
            gapFreq.push_back(0.00);
        }

        profAliFrequency.reserve(seqLen);
        for (unsigned int i = 0; i < seqLen; i++) {
            pCalculateRawFrequency(profAliFrequency[i], gapFreq[i], ali, i);

            double frequencySum = 0.00;
            for (AminoAcidCode j = ALA; j <= TYR; j++)
                frequencySum += profAliFrequency[i][j];

            for (AminoAcidCode j = ALA; j <= TYR; j++)
                profAliFrequency[i][j] /= frequencySum;
        }



        setSeq(ali.getTarget());
    }

}} // namespace
