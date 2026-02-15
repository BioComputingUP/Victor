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
// Description:     Implement a standard substitution matrix.
//
// -----------------x-----------------------------------------------------------

#include <SubMatrix.h>

namespace Victor { namespace Align2{

    // CONSTRUCTORS:

    SubMatrix::SubMatrix() : Substitution() {
    }

    SubMatrix::SubMatrix(istream &is) : Substitution() {
        is >> *this;
        double res = residues.size()-2; //avoids '-'
        double avg1 = 0, avg2 = 0, avg3 = 0;
        for(int i=0;i<=res;i++){
        	for(int j=0;j<=i;j++){
        		double value = score[residues[i]][residues[j]];
        		avg1 += value;
        		if(residues[i] == residues[j]){
        			avg2 += value;
        		}
        		else{
        			avg3 += value;
        		}
    		}
        }
        minScore = score[0][0];
        for(int i=0;i<= res;i++){
			for(int j=0;j<=i;j++){
				double value = score[residues[i]][residues[j]];
				if(value < minScore){
					minScore = value;
				}
			}
        }
        averageScore = avg1 / ((res*res)/2);
        averageMatchScore = avg2 / res;
    	averageMismatchScore =  avg3 / (((float)(res*res-res))/2);
    }

    SubMatrix::SubMatrix(const SubMatrix &orig) {
        copy(orig);
    }

    SubMatrix::~SubMatrix() {
    }


    // OPERATORS:
    /**
     * 
     * @param orig
     * @return 
     */
    SubMatrix&
            SubMatrix::operator =(const SubMatrix &orig) {
        if (&orig != this)
            copy(orig);
        POSTCOND((orig == *this), exception);
        return *this;
    }
    /**
     * 
     * @param os
     * @param object
     * @return 
     */
    ostream&
    operator <<(ostream &os, const SubMatrix &object) {
        os << object.residues << endl;
        Substitution::pWriteDoubleVector(os, object.residuescores);
        return os;
    }
    /**
     * 
     * @param is
     * @param object
     * @return 
     */
    istream&
    operator >>(istream &is, SubMatrix &object) {
        is >> object.residues;
        Substitution::pReadDoubleVector(is, object.residuescores);
        object.buildscore(object.residues, object.residuescores);
        return is;
    }


    // MODIFIERS:
    /**
     * 
     * @param orig
     */
    void
    SubMatrix::copy(const SubMatrix &orig) {
        Substitution::copy(orig);

        residuescores.clear();
        for (unsigned int n = 0; n < orig.residuescores.size(); n++)
            residuescores.push_back(orig.residuescores[n]);

        residues = orig.residues;
        averageMatchScore = orig.averageMatchScore;
        averageMismatchScore = orig.averageMatchScore;
        averageScore = orig.averageScore;
        minScore = orig.minScore;
    }
    /**
     * 
     * @return 
     */
    SubMatrix*
    SubMatrix::newCopy() {
        SubMatrix *tmp = new SubMatrix(*this);
        return tmp;
    }

    double SubMatrix::getAverageMismatchScore(){
    	return -averageMismatchScore;
    }

    double SubMatrix::getAverageMatchScore(){
    	return averageMatchScore;
    }

    double SubMatrix::getAverageScore(){
    	return averageScore;
    }

    double SubMatrix::getMinScore(){
    	return minScore;
    }

}} // namespace
