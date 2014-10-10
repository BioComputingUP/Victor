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


#ifndef __PssmInput_H__
#define __PssmInput_H__

#include <Substitution.h>
#include <iostream>
#include <string>

namespace Victor { namespace Align2{

    /** @brief  Implement I/O objects for handling BLAST PSSM (Position
     *                  Specific Score Matrix).
     * 
     *   
     * @This 
     **/
    class PssmInput {
    public:

        // CONSTRUCTORS:

        /// Default constructor.
        PssmInput();

        /// istream constructor.
        PssmInput(istream &is);

        /// Copy constructor.
        PssmInput(const PssmInput &orig);

        /// Destructor.
        virtual ~PssmInput();


        // OPERATORS:

        /// Assignment operator.
        PssmInput& operator =(const PssmInput &orig);

        /// Output operator.
        friend ostream& operator <<(ostream &os, const PssmInput &object);

        /// Input operator.
        friend istream& operator >>(istream &is, PssmInput &object);


        // PREDICATES:

        /// Return the score of the aminoacid j in position i.
        double score(int i, int j);

        /// Return the size of the object referred as the dimension of the PSSM.
        virtual unsigned int size() const;


        // MODIFIERS:

        /// Copy orig object to this object ("deep copy").
        virtual void copy(const PssmInput &orig);

        /// Construct a new "deep copy" of this object.
        virtual PssmInput* newCopy();


        // HELPERS:

        /// Helper function used to write a vector<vector> construct.
        template<class T> static void pWriteDoubleVector(ostream &os,
                vector< vector<T> > data, vector<string> data1, vector<string> data2);

        /// Helper function used to read a vector<vector> construct.
        template<class T> static void pReadDoubleVector(istream &is,
                vector< vector<T> > &data, vector<string> &data1, vector<string> &data2);


    protected:


    private:

        // ATTRIBUTES:

        vector< vector<double> > residuescores; ///< PSSM scores.
        vector<string> allPosition; ///< Companion variable for I/O class.
        vector<string> allAa; ///< Companion variable for I/O class.

    };

    // -----------------------------------------------------------------------------
    //                                  PssmInput
    // -----------------------------------------------------------------------------

    // PREDICATES:

    inline double
    PssmInput::score(int i, int j) {
        return residuescores[i][j];
    }

    inline unsigned int
    PssmInput::size() const {
        return residuescores.size();
    }

}} // namespace

#endif
