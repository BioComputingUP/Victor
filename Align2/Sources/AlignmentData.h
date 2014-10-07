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


#ifndef __AlignmentData_H__
#define __AlignmentData_H__

#include <Alignment.h>
#include <AlignmentBase.h>
#include <IoTools.h>
#include <algorithm>
#include <iostream>
#include <string>

namespace Biopool {

    /** @brief   Base class for printing alignments.
     * 
     * @Description  
     * @This 
     **/
    class AlignmentData {
    public:

        // CONSTRUCTORS:

        /// Default constructor.
        AlignmentData(int n, const string &name1 = "Seq1",
                const string &name2 = "Seq2");

        /// Copy constructor.
        AlignmentData(const AlignmentData &orig);

        /// Destructor.
        virtual ~AlignmentData();


        // OPERATORS:

        /// Assignment operator.
        AlignmentData& operator =(const AlignmentData &orig);


        // PREDICATES:

        /// Define if two residues are similar.
        virtual bool similar(char a, char b);

        /// Return the sequence at position n of the vector.
        virtual string getSequence(int n) = 0;

        /// Calculate single match positions.
        virtual void calculateMatch(int i, int tbi, int j, int tbj) = 0;

        /// Reverse the strings of the vector.
        virtual void getMatch() = 0;

        /// Control if the strings of the vector are similar and print them.
        virtual void outputMatch(ostream &os, bool fasta = false) = 0;

        /// Generate and return an ensemble of suboptimal alignments.
        virtual Alignment& generateMatch(double score = 0.00) = 0;


        // MODIFIERS:

        /// Copy orig object to this object ("deep copy").
        virtual void copy(const AlignmentData &orig);

        /// Construct a new "deep copy" of this object.
        virtual AlignmentData* newCopy() = 0;

        /// Insert a sequence at position n of the vector.
        virtual void add(string s, int n);

        /// Insert a void string in all positions of the vector.
        virtual void clear();


        // ATTRIBUTES:

        vector<string> match; ///< Matching alignment positions.
        int n; ///< Number of strings in the alignment.
        string name1; ///< Name of target sequence.
        string name2; ///< Name of template sequence.


    protected:


    private:

    };

    // -----------------------------------------------------------------------------
    //                                AlignmentData
    // -----------------------------------------------------------------------------

    // MODIFIERS:

    inline void
    AlignmentData::add(string s, int n) {
        match[n] += s;
    }

    inline void
    AlignmentData::clear() {
        for (int i = 0; i < n; i++)
            match[i] = "";
    }

} // namespace

#endif
