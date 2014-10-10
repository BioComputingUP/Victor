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


#include <ReverseScore.h>

namespace Victor { namespace Align2{

    // CONSTRUCTORS:
    /**
     * 
     * @param a
     */
    ReverseScore::ReverseScore(Align *a) {
        ali = a->newCopy();
        inv = a->newCopy();
        inv->getScoringScheme()->reverse();
    }

    ReverseScore::ReverseScore(const ReverseScore &orig) {
        copy(orig);
    }

    ReverseScore::~ReverseScore() {
    }


    // OPERATORS:
    /**
     * 
     * @param orig
     * @return 
     */
    ReverseScore&
            ReverseScore::operator =(const ReverseScore &orig) {
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
    ReverseScore::copy(const ReverseScore &orig) {
        ali = orig.ali;
        inv = orig.inv;
    }
    /**
     * 
     * @param forward
     * @param reverse
     * @param n
     * @return 
     */
    double
    ReverseScore::getZScore(double &forward, double &reverse, unsigned int n) {
        ali->recalculateMatrix();
        inv->recalculateMatrix();

        vector<double> score = inv->generateMultiMatchScore(n);

        //	cout << ">>>>>>> " << ali->getScore() << endl;
        //	for (unsigned int i = 0; i < score.size(); i++)
        //		cout << score[i] << endl;

        forward = ali->getScore();
        reverse = average(score);
        return ((forward - reverse) / (standardDeviation(score) != 0 ? standardDeviation(score) : 1));
    }

}} // namespace
