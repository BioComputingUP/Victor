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


#ifndef __Substitution_H__
#define __Substitution_H__

#include <Debug.h>
#include <iostream>
#include <string>
#include <vector>

namespace Biopool {

    /** @brief    Base class for deriving substitution matrices.
     * 
     * @Description  
     * @This 
     **/
    class Substitution {
    public:

        // CONSTRUCTORS:

        /// Default constructor.
        Substitution();

        /// Copy constructor.
        Substitution(const Substitution &orig);

        /// Destructor.
        virtual ~Substitution();


        // OPERATORS:

        /// Assignment operator.
        Substitution& operator =(const Substitution &orig);

        /// Output operator.
        friend ostream& operator <<(ostream &os, const Substitution &object);

        /// Input operator.
        friend istream& operator >>(istream &is, Substitution &object);


        // PREDICATES:

        /// Dummy implementation.
        virtual string getResidues() const = 0;


        // MODIFIERS:

        /// Copy orig object to this object ("deep copy").
        virtual void copy(const Substitution &orig);

        /// Construct a new "deep copy" of this object.
        virtual Substitution* newCopy() = 0;

        /// Build scoring matrix from raw data.
        virtual void buildscore(const string &residues,
                const vector< vector<int> > &residuescores);


        // HELPERS:

        /// Helper function used to write a vector<vector> construct.

        /*template<class T> static void pWriteDoubleVector(ostream &os,
                vector<vector<T> > data);
         */
        static void pWriteDoubleVector(ostream &os, vector<vector<int> > data);
        /// Helper function used to read a vector<vector> construct.
        template<class T> static void pReadDoubleVector(istream &is,
                vector<vector<T> > &data);


        // ATTRIBUTES:

        vector< vector<int> > score; ///< Substitution score.


    protected:


    private:

    };

} // namespace

#endif
