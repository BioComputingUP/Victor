
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

#ifndef __AlignmentBase_H__
#define __AlignmentBase_H__

#include <Debug.h>
#include <string>
#include <vector>

namespace Biopool {

    /** @brief      Abstract base class for all sorts of alignments.
     * 
     * @Description  
     * @This 
     **/
    class AlignmentBase {
    public:

        // CONSTRUCTORS:

        /// Default constructor.
        AlignmentBase();

        /// Copy constructor.
        AlignmentBase(const AlignmentBase &orig);

        /// Destructor.
        virtual ~AlignmentBase();


        // OPERATORS:

        /// Assignment operator.
        AlignmentBase& operator =(const AlignmentBase &orig);


        // PREDICATES:

        /// Return size of alignment.
        unsigned int size() const;

        /// Return length of alignment.
        unsigned int getLength();

        /// Return length of seq.
        unsigned int getSequenceLength(const string &seq);

        /// Return target name.
        virtual string getTargetName() const;

        /// Return target sequence.
        virtual string getTarget() const;

        /// Return target position p.
        virtual char getTargetPos(unsigned int p) const;

        /// Return target aa offset (only needed for alignment, default = 0).
        virtual int getTargetAminoAcidOffset() const;

        /// Return template index name.
        virtual string getTemplateName(unsigned int index = 0) const;

        /// Return template index sequence.
        virtual string getTemplate(unsigned int index = 0) const;

        /// Return template index position p.
        virtual char getTemplatePos(unsigned int p, unsigned int index = 0) const;

        /// Return template index aa offset (counting from zero).
        virtual int getTemplateAminoAcidOffset(unsigned int index = 0) const;

        /// Calculate pairwise identity between seq1 and seq2.
        virtual double calculatePairwiseIdentity(const string &seq1,
                const string &seq2);

        /// Calculate overall identity.
        virtual double calculateIdentity();

        /// Check for conservation. If index is 9999 check on @b all templates.
        virtual bool isConserved(unsigned int p, unsigned int index = 9999) const;

        /// Check for insertion at position p.
        virtual bool isInsertion(unsigned int p, unsigned int index = 0) const;

        /// Check for deletion at position p.
        virtual bool isDeletion(unsigned int p) const;

        /// Check for gap at position p.
        virtual bool isGap(unsigned int p, unsigned int index = 0) const;

        /// Return a vector of int containing the aligned aa in target and template.
        /// Each position of the vector contains the aa position in the target and
        /// in the template (-1 means gap).
        virtual vector< vector<int> > getMatchSubset();

        /// Return a new vector with positions shifted, depending on new position.
        virtual vector<int> shiftMatchSubset(vector<int> inputVector,
                int newStartPos);

        /// Companion class to the previous.
        virtual double matchPositionVector(vector<int> CeTarget,
                vector<int> CeTemplate, vector<int> seqTarget, vector<int> seqTemplate);

        /// Save single sequence in FASTA format.
        static void saveFasta(string t, string tName, ostream &output);

        /// Save as FASTA like output.
        virtual void saveFasta(ostream &output) const;

        /// Save single line in CLUSTAL format.
        static void saveClustal(string t, string tName, ostream &output,
                unsigned int from);

        /// Save as CLUSTAL like output.
        virtual void saveClustal(ostream &output) const;


        // MODIFIERS:

        /// Copy orig object to this object ("deep copy").
        virtual void copy(const AlignmentBase &orig);

        /// Construct a new "deep copy" of this object.
        virtual AlignmentBase* newCopy();

        /// Set target to t.
        virtual void setTarget(string t, string tName = "target");

        /// Set target residue to res.
        virtual void setTargetPos(unsigned int p, char res);

        /// Set target aa offset (only needed for alignment, default = 0).
        virtual void setTargetAminoAcidOffset(int orig);

        /// Set template to t.
        virtual void setTemplate(string t, string tName = "template");

        /// Set template index residue to res.
        virtual void setTemplatePos(unsigned int p, char res,
                unsigned int index = 0);

        /// Swap templates index1 and index2.
        virtual void swapTemplate(unsigned int index1, unsigned int index2);

        /// Set template index aa offset (counting from zero).
        virtual void setTemplateAminoAcidOffset(unsigned int index, int val);

        /// Insert character c in target and all templates at position p.
        void insertCharacter(unsigned int p, char c);

        /// Insert character '-' in target and all templates at position p.
        void insertDash(unsigned int p);

        /// Delete character (or '-') from target and all templates at position p.
        void deletePos(unsigned int p);

        /// Delete all gaps from target and all templates.
        void purgeTargetInsertions();

        /// Remove all templates below index.
        virtual void cutTemplate(unsigned int index);

        /// Clear template data.
        virtual void clearTemplate();

        /// Combine two multiple sequence alignments of same target.
        virtual void addAlignment(const AlignmentBase &other);

        /// Clear alignment data.
        virtual void clearAlignment();


        // HELPERS:

        /// Return sequence without '-' characters.
        static string getPureSequence(const string &s);

        /// Return position of index if '-' would not be there (counting from zero).
        static unsigned int getOrigPos(const string &s, unsigned int p);

        /// Return position of original index if '-' are now present.
        static unsigned int getNewPos(const string &s, unsigned int p);

        /// Return vector of words of a line of text.
        vector<string> getTokens(const string &text);

        /// Delete n-th character.
        string deleteChar(const string &s, unsigned int n);


    protected:

        // ATTRIBUTES:

        string targetName; ///< Target name.
        vector<string> seqTemplateName; ///< Template names.
        string target; ///< Target sequence.
        vector<string> seqTemplate; ///< Template sequences.
        int startAaTarget; ///< Start target aa offset.
        vector<int> startAaTemplates; ///< Start templates aa offsets.


    private:

    };

    // -----------------------------------------------------------------------------
    //                                AlignmentBase
    // -----------------------------------------------------------------------------

    // PREDICATES:

    inline unsigned int
    AlignmentBase::size() const {
        return seqTemplate.size();
    }

    inline unsigned int
    AlignmentBase::getLength() {
        return target.length();
    }

    inline string
    AlignmentBase::getTargetName() const {
        return targetName;
    }

    inline string
    AlignmentBase::getTarget() const {
        return target;
    }

    inline char
    AlignmentBase::getTargetPos(unsigned int p) const {
        if (p >= target.length())
            ERROR("AlignmentBase::getTargetPos() Invalid position requested.", exception);
        return target[p];
    }

    inline int
    AlignmentBase::getTargetAminoAcidOffset() const {
        return startAaTarget;
    }

    inline string
    AlignmentBase::getTemplateName(unsigned int index) const {
        if (index >= seqTemplateName.size())
            ERROR("AlignmentBase::getTemplateName() Invalid template requested.", exception);
        return seqTemplateName[index];
    }

    inline string
    AlignmentBase::getTemplate(unsigned int index) const {
        if (index >= seqTemplate.size())
            ERROR("AlignmentBase::getTemplate() Invalid template requested.", exception);
        return seqTemplate[index];
    }

    inline char
    AlignmentBase::getTemplatePos(unsigned int p, unsigned int index) const {
        if (index >= seqTemplate.size())
            ERROR("AlignmentBase::getTemplatePos() Invalid template requested.", exception);
        if (p >= seqTemplate[index].length())
            ERROR("AlignmentBase::getTemplatePos() Invalid position requested.", exception);
        return seqTemplate[index][p];
    }

    inline int
    AlignmentBase::getTemplateAminoAcidOffset(unsigned int index) const {
        return startAaTemplates[index];
    }

    inline void
    AlignmentBase::saveFasta(string t, string tName, ostream &output) {
        output << ">" << tName << "\n";

        for (unsigned int i = 0; i < t.length(); i++) {
            if ((i != 0) && ((i % 60) == 0))
                output << "\n";
            output << t[i];
        }

        output << "\n";
    }

    inline void
    AlignmentBase::saveClustal(string t, string tName, ostream &output,
            unsigned int from) {
        output << tName;

        for (int i = 0; i < (17 - static_cast<int> (tName.length())); ++i)
            output << " ";

        unsigned int max = ((from + 60) < t.length()) ? from + 60 : t.length();
        for (unsigned int i = from; i < max; i++)
            output << t[i];

        output << "\n";
    }


    // MODIFIERS:

    inline void
    AlignmentBase::setTarget(string t, string tName) {
        if (seqTemplate.size() > 0)
            if (t.length() != seqTemplate[0].length())
                ERROR("AlignmentBase::setTarget() Target length does not match template.", exception);
        targetName = tName;
        target = t;
    }

    inline void
    AlignmentBase::setTargetPos(unsigned int p, char res) {
        if (p >= target.length())
            ERROR("AlignmentBase::getTargetPos() Invalid position requested.", exception);
        target[p] = res;
    }

    inline void
    AlignmentBase::setTargetAminoAcidOffset(int orig) {
        startAaTarget = orig;
    }

    inline void
    AlignmentBase::setTemplatePos(unsigned int p, char res, unsigned int index) {
        if (index >= seqTemplate.size())
            ERROR("AlignmentBase::getTemplatePos() Invalid template requested.", expection);
        if (p >= seqTemplate[index].length())
            ERROR("AlignmentBase::getTemplatePos() Invalid position requested.", exception);
        seqTemplate[index][p] = res;
    }

    inline void
    AlignmentBase::setTemplateAminoAcidOffset(unsigned int index, int val) {
        PRECOND(index < startAaTemplates.size(), exception);
        startAaTemplates[index] = val;
    }

    inline void
    AlignmentBase::clearTemplate() {
        seqTemplate.clear();
        seqTemplateName.clear();
        startAaTemplates.clear();
    }

    inline void
    AlignmentBase::clearAlignment() {
        targetName = "";
        target = "";
        startAaTarget = 0;
        clearTemplate();
    }

} // namespace

#endif
