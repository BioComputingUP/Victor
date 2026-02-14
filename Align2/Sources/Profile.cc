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
// Description:     Calculate a frequency profile or PSSM.
//
// -----------------x-----------------------------------------------------------

#include <Profile.h>
#include <ctime>
namespace Victor { namespace Align2{

    // CONSTRUCTORS:

    Profile::Profile() : profAliFrequency(), gapFreq(), seq(""), seqLen(0),
    numSeq(0), gap(false) {
    }

    Profile::Profile(const Profile &orig) {
        copy(orig);
    }

    Profile::~Profile() {
    }


    // OPERATORS:

    Profile&
            Profile::operator =(const Profile &orig) {
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
    Profile::copy(const Profile &orig) {
        profAliFrequency.clear();
        for (unsigned int i = 0; i < orig.profAliFrequency.size(); i++) {
            vector<double> tmp;
            for (unsigned int j = 0; j < orig.profAliFrequency[0].size(); j++)
                tmp.push_back(orig.profAliFrequency[i][j]);
            profAliFrequency.push_back(tmp);
        }

        gapFreq.clear();
        for (unsigned int i = 0; i < orig.gapFreq.size(); i++)
            gapFreq.push_back(orig.gapFreq[i]);

        seq = orig.seq;
        seqLen = orig.seqLen;
        numSeq = orig.numSeq;
        gap = orig.gap;
    }

    Profile*
    Profile::newCopy() {
        Profile *tmp = new Profile(*this);
        return tmp;
    }
    /**
     * 
     * @param ali
     */
    void
    Profile::setProfile(Alignment &ali) {
        struct tm* newtime;
        time_t t;
        pResetData();
        seqLen = ali.getTarget().size();
        numSeq = ali.size();
        time(&t);
        newtime = localtime(&t);
        cout << "ready for pConstructData " << newtime->tm_hour << "/" << newtime->tm_min << endl;
        pConstructData(ali);
    }
    /**
     * 
     * @param ali
     * @param is
     */
    void
    Profile::setProfile(Alignment &ali, istream &is) {
        pResetData();
        cout << "ali:  " << ali.getTarget();
        seqLen = ali.getTarget().size();
        numSeq = ali.size() - 2;

        if (!gap) {
            ali.purgeTargetInsertions();
            seqLen = ali.getTarget().size();
        }

        for (unsigned int i = 0; i < seqLen; i++) {
            double f;

            vector<double> tmp;
            for (unsigned int j = 0; j < 20; j++) {
                is >> f;
                tmp.push_back(f);
            }
            profAliFrequency.push_back(tmp);

            is >> f;
            gapFreq.push_back(f);
        }

        setSeq(ali.getTarget());
    }
    /**
     * 
     */
    void
    Profile::reverse() {
        vector< vector<double> > tmpPAF;
        for (unsigned int i = profAliFrequency.size(); i > 0; i--) {
            vector<double> tmp;
            for (unsigned int j = 0; j < profAliFrequency[i - 1].size(); j++)
                tmp.push_back(profAliFrequency[i - 1][j]);
            tmpPAF.push_back(tmp);
        }
        profAliFrequency.clear();
        profAliFrequency = tmpPAF;

        vector<double> tmpGF;
        for (unsigned int i = gapFreq.size(); i > 0; i--)
            tmpGF.push_back(gapFreq[i - 1]);
        gapFreq.clear();
        gapFreq = tmpGF;

        string tmpS = "";
        for (unsigned int i = seq.length(); i > 0; i--)
            tmpS.push_back(seq[i - 1]);
        seq = tmpS;
    }


    // HELPERS:
    /**
     * 
     * @param freq
     * @param freqGap
     * @param ali
     * @param i
     */
    void
    Profile::pCalculateRawFrequency(vector<double> &freq, double &freqGap,
            Alignment &ali, unsigned int i) {
        if (aminoAcidOneLetterTranslator(ali.getTarget()[i]) != XXX)
            freq[aminoAcidOneLetterTranslator(ali.getTarget()[i])]++;
        else
            freqGap++;

        for (unsigned int j = 0; j < (numSeq - 1); j++)
            if (aminoAcidOneLetterTranslator(ali.getTemplate(j)[i]) != XXX)
                freq[aminoAcidOneLetterTranslator(ali.getTemplate(j)[i])]++;
            else
                freqGap++;
    }
    /**
     * 
     * @param ali
     */
    void
    Profile::pConstructData(Alignment &ali) {
        if (!gap) {
            ali.purgeTargetInsertions();
            seqLen = ali.getTarget().size();
        }
        gapFreq.reserve(seqLen);
        for (unsigned int i = 0; i < seqLen; i++) {
            vector<double> tmp;
            for (AminoAcidCode i = ALA; i <= TYR; i++)
                tmp.push_back(0.00);
            profAliFrequency.push_back(tmp);
            gapFreq.push_back(0.00);
        }
        for (unsigned int i = 0; i < seqLen; i++) {
            pCalculateRawFrequency(profAliFrequency[i], gapFreq[i], ali, i);
            for (AminoAcidCode j = ALA; j <= TYR; j++)
                profAliFrequency[i][j] /= (numSeq - gapFreq[i]);
        }
        setSeq(ali.getTarget());
    }
    /**
     * 
     */
    void
    Profile::pResetData() {
        for (unsigned int i = 0; i < profAliFrequency.size(); i++)
            profAliFrequency[i].clear();
        profAliFrequency.clear();
        gapFreq.clear();
        seq = "";
        seqLen = 0;
        numSeq = 0;
    }

}} // namespace
