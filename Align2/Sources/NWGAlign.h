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

#ifndef __NWGAlign_H__
#define __NWGAlign_H__

#include <Align.h>

namespace Victor { namespace Align2{

    /** @brief  Implement Needleman-Wunsch-Gotoh global alignment.
     * 	@author Andrea Nasari (andrea.nasari\@studenti.unipd.it)
     *   

     **/
    class NWGAlign : public Align {
    public:

        // CONSTRUCTORS:

        /// Default constructor.
    	NWGAlign(AlignmentData *ad, GapFunction *gf, ScoringScheme *ss);

        /// Constructor with weighted alignment positions.
    	NWGAlign(AlignmentData *ad, GapFunction *gf, ScoringScheme *ss,
                const vector<unsigned int> &v1, const vector<unsigned int> &v2);

        /// Copy constructor.
    	NWGAlign(const NWGAlign &orig);

        /// Destructor.
        virtual ~NWGAlign();


        // OPERATORS:

        /// Assignment operator.
        NWGAlign& operator =(const NWGAlign &orig);


        // PREDICATES:

        /// Return two-element array containing an alignment with maximal score.
        virtual void getMultiMatch();

        ///Generate and return an ensemble of suboptimal alignments.
        virtual vector<Alignment> generateMultiMatch(unsigned int num = 1);

        // MODIFIERS:

        /// Copy orig object to this object ("deep copy").
        virtual void copy(const NWGAlign &orig);

        /// Construct a new "deep copy" of this object.
        virtual NWGAlign* newCopy();

        /// Return two-element array containing an alignment with maximal score.
        virtual vector<string> getMatch() const;

        /// Return alignment score.
         virtual double getScore() const;

        // HELPERS:

        /// Update/create matrix values.
        virtual void pCalculateMatrix(bool update = true);

        /// Update/create weighted matrix values.
        virtual void pCalculateMatrix(const vector<unsigned int> &v1,
                const vector<unsigned int> &v2, bool update = true);

    protected:

    private:
        //F from Align class as Gotoh diag score matrix.
        vector< vector<double> > F1; ///< Gotoh left score matrix.
        vector< vector<double> > F2; ///< Gotoh up score matrix.
        vector< vector<int> > D;  ///< Gotoh diag traceback matrix.
        vector< vector<int> > D1; ///< Gotoh left traceback matrix.
        vector< vector<int> > D2; ///< Gotoh up traceback matrix.
    };

}} // namespace

#endif
