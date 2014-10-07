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
// Description:     Pairwise sequence and profile alignment originally based
//                  on the Java implementation from Peter Sestoft.
//                  http://www.dina.dk/~sestoft
//
// -----------------x-----------------------------------------------------------

#include <Align.h>

namespace Biopool {

    // CONSTRUCTORS:

    Align::Align(AlignmentData *ad, GapFunction *gf, ScoringScheme *ss) : ad(ad),
    gf(gf), ss(ss), F((ad->getSequence(1)).size() + 1),
    B((ad->getSequence(1)).size() + 1), n((ad->getSequence(1)).size()),
    m((ad->getSequence(2)).size()), res1Pos(), res2Pos() {
        //cout<<"building align objA\n";
        vector<double> frow(m + 1, 0);
        vector<Traceback> brow(m + 1);
        //cout<<"building align objB\n";
        for (unsigned int i = 0; i < F.size(); ++i) {
            F[i] = frow;
            B[i] = brow;
        }
        //cout<<"building align objC\n";
        setPenalties(0.98, 0.00);
        //cout<<"building align objD\n";
    }

    Align::Align(const Align &orig) {
        copy(orig);
    }

    Align::~Align() {
    }


    // OPERATORS:

    Align&
            Align::operator =(const Align &orig) {
        if (&orig != this)
            copy(orig);
        POSTCOND((orig == *this), exception);
        return *this;
    }


    // PREDICATES:

    vector<double>
    Align::getMultiMatchScore(unsigned int num) {
        vector<double> score;

        for (unsigned int i = 0; i < num; i++) {
            score.push_back(getScore());
            getMultiMatch();
        }

        return score;
    }

    vector<string>
    Align::getMatch() const {
        string res1, res2;
        res1Pos.clear();
        res2Pos.clear();
        Traceback tb = B0;
        int i = tb.i;
        int j = tb.j;

        tb = next(tb);
        while (((tb.i >= 0) || (tb.j >= 0)) && ((i != tb.i) || (j != tb.j))) {
            if (i == tb.i) {
                res1 = res1 + "-";
                res1Pos.push_back(INVALID_POS);
            } else {
                res1 = res1 + (ad->getSequence(1))[i - 1];
                res1Pos.push_back(i - 1);
            }

            if (j == tb.j) {
                res2 = res2 + "-";
                res2Pos.push_back(INVALID_POS);
            } else {
                res2 = res2 + (ad->getSequence(2))[j - 1];
                res2Pos.push_back(j - 1);
            }

            ad->calculateMatch(i, tb.i, j, tb.j);

            i = tb.i;
            j = tb.j;
            tb = next(tb);
        }

        reverse(res1.begin(), res1.end());
        reverse(res2.begin(), res2.end());
        reverse(res1Pos.begin(), res1Pos.end());
        reverse(res2Pos.begin(), res2Pos.end());
        vector<string> res(2);
        res[0] = res1;
        res[1] = res2;
        return res;
    }

    vector< vector<int> >
    Align::getMatchSubset() {
        getMatch(); // better: check if new computation is really necessary
        vector< vector<int> > resultVecs(2);
        resultVecs[0] = res1Pos;
        resultVecs[1] = res2Pos;
        return resultVecs;
    }

    vector<int>
    Align::shiftMatchSubset(vector<int> inputVector, int newStartPos) {
        unsigned int len = inputVector.size();
        vector<int> newInputVector;
        int gapCount = 0;

        for (unsigned int i = 0; i < len; i++)
            if (inputVector[i] == -1) {
                int gap = -1;
                gapCount++;
                newInputVector.push_back(gap);
            } else {
                int valueCurrent = inputVector[i];
                int newValue = valueCurrent + newStartPos;
                newInputVector.push_back(newValue);
            }

        return newInputVector;
    }

    void
    Align::outputMultiMatch(ostream &os, unsigned int num, bool fasta) {
        if (!fasta) {
            os << "Alignments:\n\n"
                    << "--xx--xx--xx--xx--xx--xx--xx--xx--xx--xx--xx--xx--xx--xx--xx--xx--xx--\n"
                    << endl;
        }

        for (unsigned int i = 0; i < num; i++) {
            if (!fasta)
                os << "Score = " << getScore() << "\t\t";
            getMultiMatch();
            ad->outputMatch(cout, fasta);
        }
    }

    vector<Alignment>
    Align::generateMultiMatch(unsigned int num) {
        vector<Alignment> va;

        for (unsigned int i = 0; i < num; i++) {
            double tmpScore = this->getScore();
            Alignment *tmp = new Alignment;
            getMultiMatch();
            *tmp = ad->generateMatch(tmpScore); // BEWARE!!!
            va.push_back(*tmp);
        }

        return va;
    }

    vector<double>
    Align::generateMultiMatchScore(unsigned int num) {
        vector<double> score;

        for (unsigned int i = 0; i < num; i++) {
            score.push_back(getScore());
            getMultiMatch();
        }

        return score;
    }

    void
    Align::doMatch(ostream &os) const {
        vector<string> match = getMatch();

        os << match[0] << "\n"
                << match[1] << "\n"
                << getScore() << "\n" << endl;
    }

    void
    Align::doMatchPlusHeader(ostream &os, string headerTarget,
            string headerTemplate) const {
        vector<string> match = getMatch();

        os << "> " << headerTarget << "\n"
                << match[0] << "\n\n"
                << "> " << headerTemplate << "\n"
                << match[1] << "\n\n"
                << "Score = " << getScore() << "\n" << endl;
    }


    // MODIFIERS:

    void
    Align::copy(const Align &orig) {
        ad = orig.ad->newCopy();
        gf = orig.gf->newCopy();
        ss = orig.ss->newCopy();

        F.clear();
        for (unsigned int i = 0; i < orig.F.size(); i++) {
            vector<double> *tmp = new vector<double>;
            for (unsigned int j = 0; j < orig.F[i].size(); j++)
                tmp->push_back(orig.F[i][j]);
            F.push_back(*tmp);
        }

        B.clear();
        for (unsigned int i = 0; i < orig.B.size(); i++) {
            vector<Traceback> *tmp = new vector<Traceback>;
            for (unsigned int j = 0; j < orig.B[i].size(); j++)
                tmp->push_back(orig.B[i][j]);
            B.push_back(*tmp);
        }

        B0 = orig.B0;
        n = orig.n;
        m = orig.m;

        res1Pos.clear();
        for (unsigned int i = 0; i < orig.res1Pos.size(); i++)
            res1Pos.push_back(orig.res1Pos[i]);

        res2Pos.clear();
        for (unsigned int i = 0; i < orig.res2Pos.size(); i++)
            res2Pos.push_back(orig.res2Pos[i]);

        penaltyMul = orig.penaltyMul;
        penaltyAdd = orig.penaltyAdd;
    }

    void
    Align::recalculateMatrix() {
        for (unsigned int i = 0; i < F.size(); i++)
            for (unsigned int j = 0; j < F[i].size(); j++)
                F[i][j] = -999;

        for (unsigned int i = 0; i < B.size(); i++)
            for (unsigned int j = 0; j < B[i].size(); j++) {
                B[i][j].i = -999;
                B[i][j].j = -999;
            }

        B0.i = -999;
        B0.j = -999;

        res1Pos.clear();
        res2Pos.clear();

        pCalculateMatrix(true);
    }

} // namespace
