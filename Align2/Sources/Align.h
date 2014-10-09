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


#ifndef __Align_H__
#define __Align_H__

#include <Alignment.h>
#include <AlignmentData.h>
#include <GapFunction.h>
#include <IoTools.h>
#include <ScoringScheme.h>
#include <Traceback.h>
#include <algorithm>
#include <iostream>
#include <limits>
#include <math.h>
#include <string>
#include <vector>

/// Numeric constant used for equality.
#define EQ_EPS   1E-8

/// Numeric approximation of equality for doubles.
#define EQUALS(a, b)   (fabs(a - b) < EQ_EPS ? true : false)

namespace Victor { namespace Align2{

    /** @brief  Pairwise sequence and profile alignment.
     * 
     * @Description   originally based
     *                  on the Java implementation from Peter Sestoft.
     *                  http://www.dina.dk/~sestoft
     * @This 
     **/
    class Align {
    public:

        /// Numeric value used to identify invalid alignment positions.

        enum {
            INVALID_POS = -1
        };


        // CONSTRUCTORS:

        /// Default constructor.
        Align(AlignmentData *ad, GapFunction *gf, ScoringScheme *ss);

        /// Copy constructor.
        Align(const Align &orig);

        /// Destructor.
        virtual ~Align();


        // OPERATORS:

        /// Assignment operator.
        Align& operator =(const Align &orig);


        // PREDICATES:

        /// Return AlignmentData pointer.
        AlignmentData* getAlignmentData();

        /// Return GapFunction pointer.
        GapFunction* getGapFunction();

        /// Return ScoringScheme pointer.
        ScoringScheme* getScoringScheme();

        /// Return next Traceback element.
        virtual Traceback next(const Traceback &tb) const;

        /// Return alignment score.
        virtual double getScore() const;

        /// Return alignment scores of an ensemble of suboptimal alignments.
        virtual vector<double> getMultiMatchScore(unsigned int num = 10);

        /// Return two-element array containing an alignment with maximal score.
        virtual vector<string> getMatch() const;

        /// Return two-element array containing an alignment with maximal score.
        virtual void getMultiMatch() = 0;

        /// Return subset corresponding to match.
        virtual vector< vector<int> > getMatchSubset();

        /// Return vector with positions shifted depending on new position.
        virtual vector<int> shiftMatchSubset(vector<int> inputVector,
                int newStartPos);

        
        virtual void outputMultiMatch(ostream &os, unsigned int num = 10,
                bool fasta = false);

        
        virtual vector<Alignment> generateMultiMatch(unsigned int num = 1);

        
        virtual vector<double> generateMultiMatchScore(unsigned int num = 10);

        /// Output of alignment result.
        virtual void doMatch(ostream &os) const;

        /// Output of alignment result (including headers of the two sequences).
        virtual void doMatchPlusHeader(ostream &os, string headerTarget,
                string headerTemplate) const;


        // MODIFIERS:

        /// 
        virtual void copy(const Align &orig);

        /// Construct a new "deep copy" of this object.
        virtual Align* newCopy() = 0;

        /// 
        void setPenalties(double mul, double add);

        ///
        void pModifyMatrix(int i, int j);

        /// 
        virtual void recalculateMatrix();


        // HELPERS:

        /// Update/create matrix values.
        virtual void pCalculateMatrix(bool update = true) = 0;

        /// Update/create weighted matrix values.
        virtual void pCalculateMatrix(const vector<unsigned int> &v1,
                const vector<unsigned int> &v2, bool update = true) = 0;


        // ATTRIBUTES:

        AlignmentData *ad; ///< Pointer to AlignmentData.
        GapFunction *gf; ///< Pointer to GapFunction.
        ScoringScheme *ss; ///< Pointer to ScoringScheme.
        vector< vector<double> > F; ///< Score matrix.
        vector< vector<Traceback> > B; ///< Traceback matrix.
        Traceback B0; ///< Starting point of the traceback.
        unsigned int n; ///< Length of target sequence.
        unsigned int m; ///< Length of template sequence.
        mutable vector<int> res1Pos; ///< Aligned positions for target sequence.
        mutable vector<int> res2Pos; ///< Aligned positions for template sequence.
        double penaltyMul; ///< Multiplicative penalty for suboptimal alignment.
        double penaltyAdd; ///< Additive penalty for suboptimal alignment.


    protected:


    private:

    };

    // -----------------------------------------------------------------------------
    //                                    Align
    // -----------------------------------------------------------------------------

    // PREDICATES:

    inline AlignmentData*
    Align::getAlignmentData() {
        return ad;
    }

    inline GapFunction*
    Align::getGapFunction() {
        return gf;
    }

    inline ScoringScheme*
    Align::getScoringScheme() {
        return ss;
    }

    inline Traceback
    Align::next(const Traceback& tb) const {
        if ((B.size() > 0) && (tb.i >= 0) && (tb.j >= 0) &&
                (tb.i < static_cast<int> (B.size())) &&
                (tb.j < static_cast<int> (B[tb.i].size())))
            return B[tb.i][tb.j];
        return Traceback::getInvalidTraceback();
    }

    inline double
    Align::getScore() const {
        return F[B0.i][B0.j];
    }


    // MODIFIERS:
    /**
     * @Description Set penalties for suboptimal alignments.
     * @param mul
     * @param add
     */
    
      inline void
    Align::setPenalties(double mul, double add) {
        penaltyMul = mul;
        penaltyAdd = add;
    }
    /**
     * @Description  Modify matrix during suboptimal alignment generation.
     * @param i
     * @param j
     */
    inline void
    Align::pModifyMatrix(int i, int j) {
        F[i][j] = penaltyMul * F[i][j] - penaltyAdd;
    }

}} // namespace

#endif
