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
// Description:     Calculate a frequency profile or PSSM using SeqDivergence
//                  weighting scheme. Some explanations can be found in:
//
//                  Guoli Wang, Roland L. Dunbrack jr.
//                  Scoring profile-to-profile sequence alignments.
//                  Institute for Cancer Research, Fox Chase Cancer Center,
//                  Philadelphia, Pennsylvania 19111, USA. March 16, 2004.
//
// -----------------x-----------------------------------------------------------

#include <SeqDivergenceProfile.h>

namespace Biopool {

    // CONSTRUCTORS:

    SeqDivergenceProfile::SeqDivergenceProfile() : Profile() {
    }

    SeqDivergenceProfile::SeqDivergenceProfile(const SeqDivergenceProfile &orig)
    : Profile(orig) {
        copy(orig);
    }

    SeqDivergenceProfile::~SeqDivergenceProfile() {
    }


    // OPERATORS:
    /**
     * @description
     * @param orig
     * @return 
     */
    SeqDivergenceProfile&
            SeqDivergenceProfile::operator =(const SeqDivergenceProfile &orig) {
        if (&orig != this)
            copy(orig);
        POSTCOND((orig == *this), exception);
        return *this;
    }


    // MODIFIERS:
    /**
     * @description
     * @param orig
     */
    void
    SeqDivergenceProfile::copy(const SeqDivergenceProfile &orig) {
        Profile::copy(orig);

        aliWeight.clear();
        for (unsigned int i = 0; i < orig.aliWeight.size(); i++)
            aliWeight.push_back(orig.aliWeight[i]);
    }
    /**
     * @description
     * @return 
     */
    SeqDivergenceProfile*
    SeqDivergenceProfile::newCopy() {
        SeqDivergenceProfile *tmp = new SeqDivergenceProfile(*this);
        return tmp;
    }


    // HELPERS:
    /**
     * @description
     * @param ali
     */
    void
    SeqDivergenceProfile::pCalculateWeight(Alignment &ali) {
        // --------------------------------------------------
        // 1. Load BLOSUM62 matrix
        // --------------------------------------------------

        string path = getenv("VICTOR_ROOT");
        if (path.length() < 3)
            cout << "Warning: environment variable VICTOR_ROOT is not set." << endl;
        string matrixFileName = path + "data/blosum62.dat";

        ifstream matrixFile(matrixFileName.c_str());
        if (!matrixFile)
            ERROR("Error opening substitution matrix file.", exception);

        SubMatrix sub(matrixFile);


        // --------------------------------------------------
        // 2. Calculate alignment scores for master sequence
        // --------------------------------------------------

        vector< vector<double> > A;
        vector<double> seqA;

        string seq1 = Alignment::getPureSequence(ali.getTarget());
        string seq2 = Alignment::getPureSequence(ali.getTarget());

        AGPFunction gf(12.00, 3.00);
        SequenceData ad(2, seq1, seq2);
        ScoringS2S ss(&sub, &ad, 0, 1.00);
        NWAlign nwAlign(&ad, &gf, &ss);

        seqA.push_back(nwAlign.getScore());

        for (unsigned int j = 1; j < numSeq; j++) {
            string seq2 = Alignment::getPureSequence(ali.getTemplate(j - 1));
            ad.setSequence(seq2, 2);

            ScoringS2S ss(&sub, &ad, 0, 1.00);
            NWAlign nwAlign(&ad, &gf, &ss);

            seqA.push_back(nwAlign.getScore());
        }

        A.push_back(seqA);


        // --------------------------------------------------
        // 3. Calculate alignment scores for other sequences
        // --------------------------------------------------

        for (unsigned int i = 1; i < numSeq; i++) {
            seqA.clear();

            string seq1 = Alignment::getPureSequence(ali.getTemplate(i - 1));
            ad.setSequence(seq1, 1);

            for (unsigned int j = i; j < numSeq; j++) {
                string seq2 = Alignment::getPureSequence(ali.getTemplate(j - 1));
                ad.setSequence(seq2, 2);

                ScoringS2S ss(&sub, &ad, 0, 1.00);
                NWAlign nwAlign(&ad, &gf, &ss);

                seqA.push_back(nwAlign.getScore());
            }

            A.push_back(seqA);
        }


        // --------------------------------------------------
        // 4. Calculate weights
        // --------------------------------------------------

        for (unsigned int i = 0; i < numSeq; i++) {
            double s = 0.00;

            for (unsigned int j = 0; j < numSeq; j++)
                if (i < j)
                    s += (max((A[i][j - i] / min(A[i][0], A[j][0])), 0.00) *
                        max((A[i][j - i] / min(A[i][0], A[j][0])), 0.00));
                else
                    s += (max((A[j][i - j] / min(A[i][0], A[j][0])), 0.00) *
                        max((A[j][i - j] / min(A[i][0], A[j][0])), 0.00));

            aliWeight.push_back(1 / (1 + s));
        }



    }
    /**
     * @description
     * @param freq
     * @param freqGap
     * @param ali
     * @param i
     */
    void
    SeqDivergenceProfile::pCalculateRawFrequency(vector<double> &freq, double &freqGap,
            Alignment &ali, unsigned int i) {
        if (aminoAcidOneLetterTranslator(ali.getTarget()[i]) != XXX)
            freq[aminoAcidOneLetterTranslator(ali.getTargetPos(i))] += aliWeight[0];
        else
            freqGap++;

        for (unsigned int j = 0; j < (numSeq - 1); j++)
            if (aminoAcidOneLetterTranslator(ali.getTemplatePos(i, j)) != XXX)
                freq[aminoAcidOneLetterTranslator(ali.getTemplatePos(i, j))] += aliWeight[j + 1];
            else
                freqGap++;
    }
    /**
     * @description
     * @param ali
     */
    void
    SeqDivergenceProfile::pConstructData(Alignment &ali) {
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
