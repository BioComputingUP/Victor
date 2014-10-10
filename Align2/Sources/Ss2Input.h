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


#ifndef __Ss2Input_H__
#define __Ss2Input_H__

#include <IoTools.h>
#include <Substitution.h>
#include <iostream>
#include <string>
#include <vector>

namespace Victor { namespace Align2{

    /** @brief   Implement I/O objects for handling PSI-PRED files.
     * 
     *   
     * @This 
     **/
    class Ss2Input {
    public:

        // CONSTRUCTORS:

        /// Default constructor.
        Ss2Input();

        /// istream constructor.
        Ss2Input(istream &is);

        /// Copy constructor.
        Ss2Input(const Ss2Input &orig);

        /// Destructor.
        virtual ~Ss2Input();


        // OPERATORS:

        /// Assignment operator.
        Ss2Input& operator =(const Ss2Input &orig);

        /// Output operator.
        friend ostream& operator <<(ostream &os, const Ss2Input &object);

        /// Input operator.
        friend istream& operator >>(istream &is, Ss2Input &object);


        // PREDICATES:

        /// Return the PSI-PRED score.
        double score(int i, int j);

        /// Return dimension of matrix.
        virtual unsigned int size() const;


        // MODIFIERS:

        /// Copy orig object to this object ("deep copy").
        virtual void copy(const Ss2Input& orig);

        /// Construct a new "deep copy" of this object.
        virtual Ss2Input* newCopy();


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
    //                                  Ss2Input
    // -----------------------------------------------------------------------------

    // PREDICATES:

    inline double
    Ss2Input::score(int i, int j) {
        return residuescores[i][j];
    }

    inline unsigned int
    Ss2Input::size() const {
        return residuescores.size();
    }

}} // namespace

#endif
