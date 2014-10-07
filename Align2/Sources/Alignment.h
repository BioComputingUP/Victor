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

#ifndef __Alignment_H__
#define __Alignment_H__

#include <AlignmentBase.h>
#include <Debug.h>
#include <string>
#include <vector>

namespace Biopool {

    /** @brief  Implement a simple alignment type.
     * 
     * @Description  
     * @This 
     **/
    class Alignment : public AlignmentBase {
    public:

        // CONSTRUCTORS:

        /// Default constructor.
        Alignment();

        /// Copy constructor.
        Alignment(const Alignment &orig);

        /// Destructor.
        virtual ~Alignment();


        // OPERATORS:

        /// Assignment operator.
        Alignment& operator =(const Alignment &orig);


        // PREDICATES:

        /// Return score of template index.
        double getScore(unsigned int index = 0) const;

        /// Return E-value of template index.
        double getEvalue(unsigned int index = 0) const;

        /// Save as MSAF (Mannheim Sequence Alignment Format) like output.
        void saveMSAF(ostream &output) const;

        /// Save as FASTA like output.
        virtual void saveFasta(ostream &output) const;


        // MODIFIERS:

        /// Clear alignment data.
        virtual void clearAlignment();

        /// Swap templates index1 and index2.
        virtual void swapTemplate(unsigned int index1, unsigned int index2);

        /// Set template to t.
        virtual void setTemplate(string t, string tName = "template",
                double tScore = 0.0, double tEvalue = -1.0);

        /// Set score of template index.
        void setScore(double val, unsigned int index = 0);

        /// Set E-value of template index.
        void setEvalue(long double val, unsigned int index = 0);

        /// Clear template.
        virtual void clearTemplate();

        /// Remove all templates below index.
        virtual void cutTemplate(unsigned int index);

        /// Output of alignment including headers of the 2 sequences.
        virtual void doMatchPlusHeader(ostream &os, string headerTarget,
                int startTarget, int endTarget, string headerTemplate, int startTemplate,
                int endTemplate, string alignType, string gapOpen, string gapExtension,
                double pWeigth, double sWeigth, double tWeigth, string seqTarget,
                string seqTemplat) const;

        /// Read FASTA format input.
        void loadFasta(istream &input);

        /// Read output of CE program.
        void loadCEBody(istream &input);

        /// Read output of CE program.
        void loadCE(istream &input);

        /// Read map file and build a matrix 2*N .
        vector< vector<int> > loadMap(istream &is);

        /// Return the first aminoacids on the sequence that match the first on CE.
        vector<int> mapStructureCE2Sequence(vector< vector<int> > mapIn, istream &pdbstream,
                string startingChain, vector<int> shiftedSeqresEntries);

        /// Companion function to the previous: check if SEQRES and ATOM in pdb entry
        /// are the same otherwise it returns the correct enumeration based on ATOM.
        int getTrueMapFromPdb(istream &is, string startingChain, int structEntryMapCe);

        /// Read BLAST output produced with blast option -m 6 this stands for a
        /// convenient form of multiple sequence alignment.
        void loadBlastMode6(istream &input);
        void loadPsiBlastMode4(istream &input);
        /// Read BLAST output produced with blast option -m 6 this stands for a
        /// convenient form of multiple sequence alignment (safe version).
        void loadBlastMode6FullSeq(istream &input, string masterTarget);

        /// Read MSAF.
        void loadMSAF(istream &input);

        /// Copy orig object to this object ("deep copy").
        virtual void copy(const Alignment &orig);

        /// Combine two multiple sequence alignments of same target.
        virtual void addAlignment(const Alignment &other);

        /// For removing templates with identity lower than x (between 0 and 1);
        /// this compare only the first with the others.
        virtual void RemoveLowerSimple(double ID);

        /// For removing templates with identity lower than x (between 0 and 1);
        /// this compare all with all.
        virtual void RemoveLowerAll(double ID);

        /// For removing templates with identity upper than x (between 0 and 1);
        /// this compare only the first with the others.
        virtual void RemoveUpperSimple(double ID);

        /// For removing templates with identity upper than x (between 0 and 1);
        /// this compare all with all.
        virtual void RemoveUpperAll(double ID);


    protected:


    private:

        /// Read output of CE program.
        void loadCEHeader(istream &input);


        // ATTRIBUTES:

        vector<double> score; ///< Score, eg. (in bits) from BLAST.
        vector<long double> evalue; ///< Expectation value, eg. from BLAST.

    };

    // -----------------------------------------------------------------------------
    //                                  Alignment
    // -----------------------------------------------------------------------------

    // PREDICATES:

    inline double
    Alignment::getScore(unsigned int index) const {
        if (index >= score.size())
            ERROR("Invalid template requested.", exception);
        return score[index];
    }

    inline double
    Alignment::getEvalue(unsigned int index) const {
        if (index >= evalue.size())
            ERROR("Invalid template requested.", exception);
        return evalue[index];
    }


    // MODIFIERS:

    inline void
    Alignment::clearAlignment() {
        AlignmentBase::clearAlignment();
        score.clear();
        evalue.clear();
    }

    inline void
    Alignment::clearTemplate() {
        AlignmentBase::clearTemplate();
        score.clear();
        evalue.clear();
    }

} // namespace

#endif
