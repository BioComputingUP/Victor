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


#ifndef __ThreadingInput_H__
#define __ThreadingInput_H__

#include <Substitution.h>
#include <iostream>
#include <string>

namespace Victor { namespace Align2{

    /** @brief    Implement I/O objects for handling threading files.
     * 
     * @Description  
     * @This 
     **/
    class ThreadingInput {
    public:

        // CONSTRUCTORS:

        /// Default constructor.
        ThreadingInput();

        /// istream constructor.
        ThreadingInput(istream &is);

        /// Copy constructor.
        ThreadingInput(const ThreadingInput &orig);

        /// Destructor.
        virtual ~ThreadingInput();


        // OPERATORS:

        /// Assignment operator.
        ThreadingInput& operator =(const ThreadingInput &orig);

        /// Output operator.
        friend ostream& operator <<(ostream &os, const ThreadingInput &object);

        /// Input operator.
        friend istream& operator >>(istream &is, ThreadingInput &object);


        // PREDICATES:

        /// Return the threading score.
        double score(int i, int j);

        /// Return the size of the object referred as the dimension of the matrix.
        virtual unsigned int size() const;


        // MODIFIERS:

        /// Copy orig object to this object ("deep copy").
        virtual void copy(const ThreadingInput &orig);

        /// Construct a new "deep copy" of this object.
        virtual ThreadingInput* newCopy();


        // HELPERS:

        /// Helper function used to write a vector<vector> construct.
        template<class T> static void pWriteDoubleVector(ostream &os,
                vector< vector<T> > data);

        /// Helper function used to read a vector<vector> construct.
        template<class T> static void pReadDoubleVector(istream &is,
                vector< vector<T> > &data);


    protected:


    private:

        // ATTRIBUTES:

        vector< vector<double> > residuescores; ///< Similarity scores.

    };

    // -----------------------------------------------------------------------------
    //                                ThreadingInput
    // -----------------------------------------------------------------------------

    // PREDICATES:

    inline double
    ThreadingInput::score(int i, int j) {
        return residuescores[j][i];
    }

    inline unsigned int
    ThreadingInput::size() const {
        return residuescores.size();
    }

}} // namespace

#endif
