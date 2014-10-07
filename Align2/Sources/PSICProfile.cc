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
// Description:     Calculate a frequency profile or PSSM using PSIC weighting
//                  scheme. Some explanations can be found in:
//
//                  Guoli Wang, Roland L. Dunbrack jr.
//                  Scoring profile-to-profile sequence alignments.
//                  Institute for Cancer Research, Fox Chase Cancer Center,
//                  Philadelphia, Pennsylvania 19111, USA. March 16, 2004.
//
// -----------------x-----------------------------------------------------------

#include <PSICProfile.h>

namespace Biopool {

    // CONSTRUCTORS:

    PSICProfile::PSICProfile() : Profile() {
    }

    PSICProfile::PSICProfile(const PSICProfile &orig) : Profile(orig) {
        copy(orig);
    }

    PSICProfile::~PSICProfile() {
    }


    // OPERATORS:

    PSICProfile&
            PSICProfile::operator =(const PSICProfile &orig) {
        if (&orig != this)
            copy(orig);
        POSTCOND((orig == *this), exception);
        return *this;
    }


    // MODIFIERS:

    void
    PSICProfile::copy(const PSICProfile &orig) {
        Profile::copy(orig);

        aliWeight.clear();
        for (unsigned int i = 0; i < orig.aliWeight.size(); i++) {
            vector<double> tmp;
            for (unsigned int j = 0; j < orig.aliWeight[0].size(); j++)
                tmp.push_back(orig.aliWeight[i][j]);
            aliWeight.push_back(tmp);
        }
    }

    PSICProfile*
    PSICProfile::newCopy() {
        PSICProfile *tmp = new PSICProfile(*this);
        return tmp;
    }


    // HELPERS:

    void
    PSICProfile::pCalculateWeight(Alignment &ali) {
        // --------------------------------------------------
        // 1. Calculate the positions of the first aminoacids
        //    for every sequence
        // --------------------------------------------------

        vector<unsigned int> firstAminoPos;

        firstAminoPos.push_back(0); // master sequence

        for (unsigned int j = 0; j < (numSeq - 1); j++) // other sequences
            for (unsigned int i = 0; i < seqLen; i++)
                if (ali.getTemplatePos(i, j) != '-') {
                    firstAminoPos.push_back(i);
                    break;
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

        vector<double> wMaster;

        for (unsigned int i = 0; i < seqLen; i++)
            if (ali.getTargetPos(i) == '-')
                wMaster.push_back(0.00);
            else {
                // Calculate the subset of sequences

                vector<unsigned int> seqSubset;

                seqSubset.push_back(0);
                for (unsigned int k = 0; k < (numSeq - 1); k++)
                    if (ali.getTemplatePos(i, k) == ali.getTargetPos(i))
                        seqSubset.push_back(k + 1);


                // Calculate Cleft and Cright

                unsigned int Cleft = firstAminoPos[seqSubset[0]];
                unsigned int Cright = lastAminoPos[seqSubset[0]];

                for (unsigned int k = 1; k < seqSubset.size(); k++) {
                    if (firstAminoPos[seqSubset[k]] > Cleft)
                        Cleft = firstAminoPos[seqSubset[k]];
                    if (lastAminoPos[seqSubset[k]] < Cright)
                        Cright = lastAminoPos[seqSubset[k]];
                }


                // Calculate the average number of different aminoacids

                const string residue_indices = "ARNDCQEGHILKMFPSTWYV";
                unsigned int count = 0;

                for (unsigned int p = Cleft; p <= Cright; p++)
                    for (unsigned int index = 0; index < 20; index++)
                        for (unsigned int k = 0; k < seqSubset.size(); k++)
                            if (seqSubset[k] == 0) {
                                if (ali.getTargetPos(p) == residue_indices[index]) {
                                    count++;
                                    break;
                                }
                            } else
                                if (ali.getTemplatePos(p, seqSubset[k] - 1) == residue_indices[index]) {
                                count++;
                                break;
                            }

                double F = (double) count / (double) (Cright - Cleft + 1);


                // Calculate the weight

                double neff = (1 / log(0.95)) * log(1 - (F / 20));
                unsigned int N = seqSubset.size();
                wMaster.push_back(neff / (double) N);
            }

        aliWeight.push_back(wMaster);


        // --------------------------------------------------
        // 4. Calculate weights for other sequences
        // --------------------------------------------------

        for (unsigned int j = 0; j < (numSeq - 1); j++) {
            vector<double> wOther;

            for (unsigned int i = 0; i < seqLen; i++)
                if (ali.getTemplatePos(i, j) == '-')
                    wOther.push_back(0.00);
                else {
                    // Calculate the subset of sequences

                    vector<unsigned int> seqSubset;

                    if (ali.getTargetPos(i) == ali.getTemplatePos(i, j))
                        seqSubset.push_back(0);
                    for (unsigned int k = 0; k < (numSeq - 1); k++)
                        if (ali.getTemplatePos(i, k) == ali.getTemplatePos(i, j))
                            seqSubset.push_back(k + 1);


                    // Calculate Cleft and Cright

                    unsigned int Cleft = firstAminoPos[seqSubset[0]];
                    unsigned int Cright = lastAminoPos[seqSubset[0]];

                    for (unsigned int k = 1; k < seqSubset.size(); k++) {
                        if (firstAminoPos[seqSubset[k]] > Cleft)
                            Cleft = firstAminoPos[seqSubset[k]];
                        if (lastAminoPos[seqSubset[k]] < Cright)
                            Cright = lastAminoPos[seqSubset[k]];
                    }


                    // Calculate the average number of different aminoacids

                    const string residue_indices = "ARNDCQEGHILKMFPSTWYV";
                    unsigned int count = 0;

                    for (unsigned int p = Cleft; p <= Cright; p++)
                        for (unsigned int index = 0; index < 20; index++)
                            for (unsigned int k = 0; k < seqSubset.size(); k++)
                                if (seqSubset[k] == 0) {
                                    if (ali.getTargetPos(p) == residue_indices[index]) {
                                        count++;
                                        break;
                                    }
                                } else
                                    if (ali.getTemplatePos(p, seqSubset[k] - 1) == residue_indices[index]) {
                                    count++;
                                    break;
                                }

                    double F = (double) count / (double) (Cright - Cleft + 1);


                    // Calculate the weight

                    double neff = (1 / log(0.95)) * log(1 - (F / 20));
                    unsigned int N = seqSubset.size();
                    wOther.push_back(neff / (double) N);
                }

            aliWeight.push_back(wOther);
        }




    }

    void
    PSICProfile::pCalculateRawFrequency(vector<double> &freq, double &freqGap,
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

    void
    PSICProfile::pConstructData(Alignment &ali) {
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

} // namespace
