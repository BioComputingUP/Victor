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
// Description:     Print alignment of two sequences considering also secondary
//                  structure.
//
// -----------------x-----------------------------------------------------------

#include <SecSequenceData.h>

namespace Victor { namespace Align2{

    // CONSTRUCTORS:
    /**
     * @description
     * @param n
     * @param seq1
     * @param seq2
     * @param sec1
     * @param sec2
     * @param _n1
     * @param _n2
     */
    SecSequenceData::SecSequenceData(int n, const string &seq1, const string &seq2,
            const string &sec1, const string &sec2, const string &_n1,
            const string &_n2) : AlignmentData(n, _n1, _n2), seq1(seq1), seq2(seq2),
    sec1(sec1), sec2(sec2) {
    }
    /**
     * @description
     * @param orig
     */
    SecSequenceData::SecSequenceData(const SecSequenceData &orig)
    : AlignmentData(orig) {
        copy(orig);
    }

    SecSequenceData::~SecSequenceData() {
    }


    // OPERATORS:
    /**
     * @description
     * @param orig
     * @return 
     */
    SecSequenceData&
            SecSequenceData::operator =(const SecSequenceData &orig) {
        if (&orig != this)
            copy(orig);
        POSTCOND((orig == *this), exception);
        return *this;
    }


    // PREDICATES:
    /**
     * @description
     * @param n
     * @return 
     */
    string
    SecSequenceData::getSequence(int n) {
        switch (n) {
            case 1:
                return seq1;
            case 2:
                return seq2;
            case 3:
                return sec1;
            default:
                return sec2;
        }
    }
    /**
     * @description
     * @param i
     * @param tbi
     * @param j
     * @param tbj
     */
    void
    SecSequenceData::calculateMatch(int i, int tbi, int j, int tbj) {
        string res1 = "", ret1 = "";
        string res2 = "", ret2 = "";

        if (i == tbi) {
            res1 = res1 + "-";
            ret1 = ret1 + "-";
        } else {
            res1 = res1 + seq1[i - 1];
            ret1 = ret1 + sec1[i - 1];
        }

        if (j == tbj) {
            res2 = res2 + "-";
            ret2 = ret2 + "-";
        } else {
            res2 = res2 + seq2[j - 1];
            ret2 = ret2 + sec2[j - 1];
        }

        add(res1, 0);
        add(ret1, 1);
        add(res2, 2);
        add(ret2, 3);
    }
    /**
     * @description
     */
    void
    SecSequenceData::getMatch() {
        for (int i = 0; i < n; i++) {
            string tmp = match[i];
            reverse(tmp.begin(), tmp.end());
            match[i] = tmp;
        }
    }
    /**
     * @description
     * @param os
     * @param fasta
     */
    void
    SecSequenceData::outputMatch(ostream &os, bool fasta) {
        string temp = "";
        unsigned int cons = 0;
        unsigned int sim = 0;

        for (unsigned int j = 0; j < match[0].length(); j++)
            if (match[0][j] == match[2][j]) {
                temp += match[0][j];
                cons++;
            } else
                if (similar(match[0][j], match[2][j])) {
                temp += '+';
                sim++;
            } else
                temp += " ";

        if (!fasta) {
            os << "Conservation = " << cons << " / " << cons + sim << " / "
                    << match[0].length() << "\t\t" << setw(5) << setprecision(3)
                    << (static_cast<double> (cons) / match[0].length() * 100) << "%  /  "
                    << setw(5) << setprecision(3)
                    << (static_cast<double> (cons + sim) / match[0].length() * 100) << "%\n\n"
                    << match[1] << "\n"
                    << match[0] << "\n"
                    << temp << "\n"
                    << match[2] << "\n"
                    << match[3] << "\n"
                    << "--xx--xx--xx--xx--xx--xx--xx--xx--xx--xx--xx--xx--xx--xx--xx--xx--xx--" << endl;
        } else {
            os << ">" << name1 << "\n"
                    << match[0] << "\n"
                    << ">" << name2 << "\n"
                    << match[2] << "\n"
                    << ">SecStr1\n"
                    << match[1] << "\n"
                    << ">SeqcStr2\n"
                    << match[3] << endl;
        }

        clear();
    }
    /**
     * @description
     * @param score
     * @return 
     */
    Alignment&
    SecSequenceData::generateMatch(double score) {
        Alignment *ali = new Alignment;
        ali->setTarget(match[0], name1);
        ali->setTemplate(match[2], name2, score);
        ali->setTemplate(match[1], "SecStr1");
        ali->setTemplate(match[3], "SecStr2");
        clear();
        return *ali;
    }


    // MODIFIERS:
    /**
     * @description
     * @param orig
     */
    void
    SecSequenceData::copy(const SecSequenceData &orig) {
        AlignmentData::copy(orig);
        seq1 = orig.seq1;
        seq2 = orig.seq2;
        sec1 = orig.sec1;
        sec2 = orig.sec2;
    }

    SecSequenceData*
    SecSequenceData::newCopy() {
        SecSequenceData *tmp = new SecSequenceData(*this);
        return tmp;
    }
    /**
     * @description
     * @param s
     * @param n
     */
    void
    SecSequenceData::setSequence(string s, int n) {
        switch (n) {
            case 1:
                seq1 = s;
                break;
            case 2:
                seq2 = s;
                break;
            case 3:
                sec1 = s;
                break;
            default:
                sec2 = s;
        }
    }

}} // namespace
