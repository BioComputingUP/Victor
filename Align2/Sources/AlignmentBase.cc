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
// -*- C++ -*-------x-----------------------------------------------------------
//
// Description:     Abstract base class for all sorts of alignments.
//
// -----------------x-----------------------------------------------------------

#include <AlignmentBase.h>
#include <IoTools.h>

namespace Biopool {

    // CONSTRUCTORS:

    AlignmentBase::AlignmentBase() : targetName(), target(), seqTemplateName(),
    seqTemplate(), startAaTarget(0) {
    }

    AlignmentBase::AlignmentBase(const AlignmentBase &orig) {
        copy(orig);
    }

    AlignmentBase::~AlignmentBase() {
    }


    // OPERATORS:
    /**
     * @description 
     * @param orig
     * @return 
     */
    AlignmentBase&
            AlignmentBase::operator =(const AlignmentBase &orig) {
        if (&orig != this)
            copy(orig);
        POSTCOND((orig == *this), exception);
        return *this;
    }


    // PREDICATES:
    /**
     * @description 
     * @param seq
     * @return 
     */
    unsigned int
    AlignmentBase::getSequenceLength(const string &seq) {
        return (getPureSequence(seq)).length();
    }
    /**
     * @description
     * @param seq1
     * @param seq2
     * @return 
     */
    double
    AlignmentBase::calculatePairwiseIdentity(const string &seq1, const string &seq2) {
        unsigned int maxL = (seq1.length() <= seq2.length()) ? seq1.length() : seq2.length();

        if (seq1.length() != seq2.length())
            cout << "Warning: sequence lengths do not match:\n"
                << "seq1 = " << seq1.length() << "\n"
            << "seq2 = " << seq2.length() << "\n";

        double tmp = 0;
        int countN = 0;

        for (unsigned int i = 0; i < maxL; i++)
            if ((seq1[i] != '-') && (seq1[i] == seq2[i]))
                tmp++;

        for (unsigned int i = 0; i < maxL; i++)
            if (seq2[i] == '-')
                countN++;


        return tmp / (seq1.length() - countN);
    }
    /**
     * @description
     * @return 
     */
    double
    AlignmentBase::calculateIdentity() {
        double tmp = 0;
        int countN = 0;

        for (unsigned int i = 0; i < target.length(); i++) {
            char tmpC = target[i];
            if (tmpC == '-')
                continue;

            bool add = true;
            for (unsigned int j = 0; j < seqTemplate.size(); j++)
                if (seqTemplate[j][i] == '-') {
                    countN++;
                    add = false;
                    break;
                } else
                    if (seqTemplate[j][i] != tmpC) {
                    add = false;
                    break;
                }

            if (add)
                tmp++;
        }

        return tmp / (getSequenceLength(target) - countN);
    }
    /**
     * @description
     * @param p
     * @param index
     * @return 
     */
    bool
    AlignmentBase::isConserved(unsigned int p, unsigned int index) const {
        if ((p > target.length()) || (index > this->size()))
            ERROR("Argument out of scope.", exception);

        if (index == 9999) {
            for (unsigned int i = 0; i < seqTemplate.size(); i++)
                if (target[p] != seqTemplate[i][p])
                    return false;
        } else
            if (target[p] != seqTemplate[index][p])
            return false;

        return true;
    }
    /**
     * @description
     * @param p
     * @param index
     * @return 
     */
    bool
    AlignmentBase::isInsertion(unsigned int p, unsigned int index) const {
        if ((p > target.length()) || (index >= seqTemplate.size()))
            ERROR("Argument out of scope.", exception);

        if (seqTemplate[index][p] != '-')
            return false;

        return true;
    }
    /**
     * @description
     * @param p
     * @return 
     */
    bool
    AlignmentBase::isDeletion(unsigned int p) const {
        if (p > target.length())
            ERROR("Argument out of scope.", exception);

        if (target[p] != '-')
            return false;

        return true;
    }
    /**
     * @description
     * @param p
     * @param index
     * @return 
     */
    bool
    AlignmentBase::isGap(unsigned int p, unsigned int index) const {
        if ((p > target.length()) || (index >= seqTemplate.size()))
            ERROR("Argument out of scope.", exception);

        if ((target[p] == '-') && (seqTemplate[index][p] == '-'))
            return true;

        return false;
    }
    /**
     * @description
     * @return 
     */
    vector< vector<int> >
    AlignmentBase::getMatchSubset() {
        string targetSequence = target;
        string templateSequence = seqTemplate[0];
        string stringGap = "A-"; // gap corresponds to position 1
        int countMatchTarget = -1; // I start from minus one because the first time it enters on the loop will give position zero meaning the start of the match
        int countMatchTemplate = -1;
        vector< vector<int> > resultVecs(2);

        for (unsigned int i = 0; i < templateSequence.size(); i++)
            if (targetSequence[i] == stringGap[1]) {
                int gap = -1;
                resultVecs[0].push_back(gap);
            } else {
                countMatchTarget++;
                resultVecs[0].push_back(countMatchTarget);
            }

        for (unsigned int i = 0; i < seqTemplate[0].size(); i++)
            if (templateSequence[i] == stringGap[1]) {
                int gap = -1;
                resultVecs[1].push_back(gap);
            } else {
                countMatchTemplate++;
                resultVecs[1].push_back(countMatchTemplate);
            }

        return resultVecs;
    }
    /**
     * @description
     * @param inputVector
     * @param newStartPos
     * @return 
     */
    vector<int>
    AlignmentBase::shiftMatchSubset(vector<int> inputVector, int newStartPos) {
        unsigned int len = inputVector.size();
        vector<int> newInputVector;
        int gapCount = 0;

        for (unsigned int i = 0; i < len; i++) // browse the whole vector
            if (inputVector[i] == -1) {
                int gap = -1;
                gapCount++;
                newInputVector.push_back(gap);
            } else {
                int valueCurrent = inputVector[i];
                int newValue = valueCurrent + newStartPos;
                newInputVector.push_back(newValue);
            }

        return newInputVector;
    }
    /**
     * @description
     * @param CeTarget
     * @param CeTemplate
     * @param seqTarget
     * @param seqTemplate
     * @return 
     */
    double
    AlignmentBase::matchPositionVector(vector<int> CeTarget, vector<int> CeTemplate,
            vector<int> seqTarget, vector<int> seqTemplate) {
        double overlap = 0.00;
        int counting = 0;
        int matchCount = 0;



        for (unsigned int i = 0; i < seqTarget.size(); i++) {
            int tmp;
            tmp = seqTarget[i];

            if ((tmp == -1) || (seqTemplate[i] == -1))
                continue;

            matchCount++;
            for (unsigned int j = 0; j < CeTarget.size(); j++)
                if (CeTarget[j] == tmp) {
                    if (CeTemplate[j] == seqTemplate[i])
                        counting++;
                    break;
                }
        }

        overlap = static_cast<double> (counting) / matchCount;
        return overlap;
    }
    /**
     * @description
     * @param output
     */
    void
    AlignmentBase::saveFasta(ostream &output) const {
        saveFasta(target, targetName, output);
        for (unsigned int j = 0; j < seqTemplate.size(); j++)
            saveFasta(seqTemplate[j], seqTemplateName[j], output);
    }
    /**
     * @description
     * @param output
     */
    void
    AlignmentBase::saveClustal(ostream &output) const {
        output << "CLUSTAL\n\n";

        for (unsigned int from = 0; from < target.length(); from += 60) {
            saveClustal(target, this->targetName, output, from);
            for (unsigned int j = 0; j < seqTemplate.size(); j++)
                saveClustal(seqTemplate[j], this->seqTemplateName[j], output, from);
            output << "\n";
        }
    }


    // MODIFIERS:
    /**
     * @description
     * @param orig
     */
    void
    AlignmentBase::copy(const AlignmentBase &orig) {
        targetName = orig.targetName;
        target = orig.target;
        startAaTarget = orig.startAaTarget;

        seqTemplateName.clear();
        seqTemplate.clear();
        startAaTemplates.clear();
        for (unsigned int i = 0; i < orig.seqTemplate.size(); i++) {
            seqTemplateName.push_back(orig.seqTemplateName[i]);
            seqTemplate.push_back(orig.seqTemplate[i]);
            startAaTemplates.push_back(orig.startAaTemplates[i]);
        }
    }
    /**
     * @description
     * @return 
     */
    AlignmentBase*
    AlignmentBase::newCopy() {
        AlignmentBase *tmp = new AlignmentBase(*this);
        return tmp;
    }
    /**
     * @description
     * @param t
     * @param tName
     */
    void    
    AlignmentBase::setTemplate(string t, string tName) {
        if (t.length() != target.length())
            ERROR("AlignmentBase::setTemplate() Template length does not match target.", exception);
        seqTemplateName.push_back(tName);
        seqTemplate.push_back(t);
        startAaTemplates.push_back(0);
    }
    /**
     * @description
     * @param index1
     * @param index2
     */
    void
    AlignmentBase::swapTemplate(unsigned int index1, unsigned int index2) {
        if ((index1 >= this->size()) || (index2 >= this->size()))
            ERROR("AlignmentBase::swapTarget() Index out of range.", exception);
        swap(seqTemplateName[index1], seqTemplateName[index2]);
        swap(seqTemplate[index1], seqTemplate[index2]);
        swap(startAaTemplates[index1], startAaTemplates[index2]);
    }
    /**
     * @description
     * @param p
     * @param c
     */
    void
    AlignmentBase::insertCharacter(unsigned int p, char c) {
        if (p <= target.size())
            target.insert(p, 1, c);

        for (unsigned int i = 0; i < seqTemplate.size(); ++i)
            if (p <= seqTemplate[i].size())
                seqTemplate[i].insert(p, 1, c);
    }
    /**
     * @description
     * @param p
     */
    void
    AlignmentBase::insertDash(unsigned int p) {
        insertCharacter(p, '-');
    }
    /**
     * @description
     * @param p
     */
    void
    AlignmentBase::deletePos(unsigned int p) {
        PRECOND((p < target.size()), exception);
        target = deleteChar(target, p);

        for (unsigned int i = 0; i < seqTemplate.size(); ++i) {
            ASSERT((p < seqTemplate[i].size()), exception);
            seqTemplate[i] = deleteChar(seqTemplate[i], p);
        }
    }
    /**
     * @description
     */
    void
    AlignmentBase::purgeTargetInsertions() {
        unsigned int i = 0;

        while (i < target.size())
            if ((target[i] == '-') || (target[i] == 'X')) {
                deletePos(i);
                i = 0;
            } else
                ++i;
    }
    /**
     * @description
     * @param index
     */
    void
    AlignmentBase::cutTemplate(unsigned int index) {
        // Removes all templates below index.
        if (index > seqTemplate.size())
            ERROR("Index out of scope.", exception);

        // Resize templates.
        seqTemplateName.resize(index);
        seqTemplate.resize(index);
        startAaTemplates.resize(index);

        // Cut empty positions.
        unsigned int i = 0;

        while (i < target.size())
            if (target[i] == '-') {
                bool gap = true;
                for (unsigned int j = 0; j < seqTemplate.size(); j++)
                    gap = gap && (seqTemplate[j][i] == '-');
                if (gap)
                    deletePos(i);
                else
                    i++;
            } else
                i++;
    }


    /// The algorithm is so long, because 6 cases are distinguished:
    /// - leading insertion for target A or target B
    /// - normal insertion for target A or target B
    /// - trailing insertion for target A or target B
    /// If ignoreInsertion is set, insertions of target A are ignored.
    /**
     * @description
     * @param orig
     */
    void
    AlignmentBase::addAlignment(const AlignmentBase &orig) {
        unsigned int i = 0;
        unsigned int origPos = 0;
        unsigned int newTemplatePos = 0;
        unsigned int leftGapCount1 = 0;
        unsigned int leftGapCount2 = 0;
        unsigned int rightGapCount1 = 0;
        unsigned int rightGapCount2 = 0;
        unsigned int leftGapCount2Orig = leftGapCount2;
        unsigned int rightGapCount2Orig = rightGapCount2;
        int diff = 0;
        AlignmentBase other = orig; // make safety copy
        string pure1 = getPureSequence(target); // sequence without '-'
        string pure2 = getPureSequence(other.target); // sequence without '-'

        if (target == other.target)
            goto COMBINE_END; // nothing to do if sequences have same insertions

        // Prepare sequence for adding (provisorically add 'Z' for undefined aa).
        if (orig.startAaTarget > startAaTarget)
            for (int i = 0; i < (orig.startAaTarget - startAaTarget); ++i)
                other.insertCharacter(0, '~');
        else
            if (orig.startAaTarget < startAaTarget) {
            for (int i = 0; i < (orig.startAaTarget - startAaTarget); ++i)
                insertCharacter(0, '~');
            startAaTarget = orig.startAaTarget;
        }

        // Count still missing letters and add at end (getPureSequence also counts '-').
        diff = static_cast<int> (getPureSequence(target).size()) -
                static_cast<int> (getPureSequence(other.target).size());
        if (diff < 0) // still not same size: add "X" at end
            for (int i = 0; i < abs(diff); ++i)
                insertCharacter(target.size(), '~');
        else
            if (diff > 0) // still not same size: add "X" at end
            for (int i = 0; i < abs(diff); ++i)
                other.insertCharacter(other.target.size(), '~');

#ifdef DEBUG_VERBOSE
        cout << "After preprocessing:\n"
                << target << "\n"
                << other.target << endl;
#endif


        // Count leading gaps.
        while ((leftGapCount1 < target.size()) && (target[leftGapCount1] == '-'))
            ++leftGapCount1;
        while ((leftGapCount2 < other.target.size()) && (other.target[leftGapCount2] == '-'))
            ++leftGapCount2;
        leftGapCount2Orig = leftGapCount2;

#ifdef DEBUG_VERBOSE
        cout << "Initial gap 1: " << leftGapCount1 << "\n"
                << "Initial gap 2: " << leftGapCount2 << endl;
#endif

        // Count trailing gaps.
        while ((rightGapCount1 < target.size()) && (target[target.size() - rightGapCount1 - 1] == '-'))
            ++rightGapCount1;
        while ((rightGapCount2 < other.target.size()) && (other.target[other.target.size() - rightGapCount2 - 1] == '-'))
            ++rightGapCount2;
        rightGapCount2Orig = rightGapCount2;

#ifdef DEBUG_VERBOSE
        cout << "Trailing gap 1: " << rightGapCount1 << "\n"
                << "Trailing gap 2: " << rightGapCount2 << endl;
#endif

        if (leftGapCount1 > leftGapCount2) {
            for (unsigned int k = 0; k < (leftGapCount1 - leftGapCount2); ++k) {
#ifdef DEBUG_VERBOSE
                cout << "Inserting gap into pos 0 of target 2: " << endl;
#endif

                other.insertDash(0);
            }
            leftGapCount2 = leftGapCount1;
        } else
            if (leftGapCount1 < leftGapCount2) {
            for (unsigned int k = 0; k < (leftGapCount2 - leftGapCount1); ++k) {
#ifdef DEBUG_VERBOSE
                cout << "Inserting gap into pos 0 of target 1: " << endl;
#endif

                insertDash(0);
            }
            leftGapCount1 = leftGapCount2;
        }

#ifdef DEBUG_VERBOSE
        cout << target << "\n"
                << other.target << endl;
#endif

        if (rightGapCount1 > rightGapCount2) {
            for (unsigned int k = 0; k < (rightGapCount1 - rightGapCount2); ++k) {
#ifdef DEBUG_VERBOSE
                cout << "Inserting gap into pos " << target.size() << " of target 2: " << endl;
#endif

                other.insertDash(other.target.size());
            }
            rightGapCount2 = rightGapCount1;
        } else
            if (rightGapCount1 < rightGapCount2) {
            for (unsigned int k = 0; k < (rightGapCount2 - rightGapCount1); ++k) {
#ifdef DEBUG_VERBOSE
                cout << "Inserting gap into pos " << target.size() << " of target 1: " << endl;
#endif

                insertDash(target.size());
            }
            rightGapCount1 = rightGapCount2;
        }

#ifdef DEBUG_VERBOSE
        cout << target << "\n" << other.target << endl;
#endif

        i = leftGapCount1; // now skip leading gaps
        while (i < (target.size() - rightGapCount1)) {
            if (target[i] == '-') {
                origPos = getOrigPos(target, i); // position of first non-dash char
                newTemplatePos = getNewPos(other.target, origPos); // compute position on other sequence

                unsigned int l = i + 1; // count length of insertion
                for (l = i + 1; (l < target.size()) && (target[l] == '-'); ++l) {
                }

                unsigned int l2 = l - i; // length of insertion
                unsigned int l3 = newTemplatePos;

#ifdef DEBUG_VERBOSE
                cout << "Insertion of length " << l2 << " at position " << i << " found in target 1 " << "\n"
                        << "Corresponding position on target 2: " << l3 << endl;
#endif

                for (l3 = newTemplatePos; (l3 < other.target.size()) && (other.target[l3] == '-'); ++l) {
                }

                unsigned int l4 = l3 - newTemplatePos; // length of insertion of template

#ifdef DEBUG_VERBOSE
                cout << "Found insertion of length " << l4 << " at position " << newTemplatePos << " of target 2" << endl;
#endif

                if (l4 < l2) // insert this many dashes
                {
                    unsigned int dl = l2 - l4;

                    for (unsigned int j = 0; j < dl; ++j) // insertion of other alignment is dl residues too short
                    {
#ifdef DEBUG_VERBOSE
                        cout << " inserting dash after position " << (newTemplatePos + 1) << " on target 2 " << endl;
#endif

                        other.insertDash(newTemplatePos + 1);
                    }

                    i += dl;
                } else
                    ++i;
            } else
                ++i;
        }

        i = leftGapCount2Orig; // start again from zero

#ifdef DEBUG_VERBOSE
        cout << "Other direction:" << endl;
#endif

        while (i < (orig.target.size() - rightGapCount2Orig)) // go through all '-' regions in target 2
        {
            if (orig.target[i] == '-') {
#ifdef DEBUG_VERBOSE
                cout << "character " << i << " of target 2\n"
                        << " dash found." << endl;
#endif

                origPos = getOrigPos(orig.target, i);
                newTemplatePos = getNewPos(target, origPos); // compute position on other sequence

                unsigned int l = i + 1; // count length of insertion
                for (l = i + 1; (l < orig.target.size()) && (orig.target[l] == '-'); ++l) {
                }

                unsigned int l2 = l - i; // length of insertion
                unsigned int l3 = newTemplatePos; // count length of insertion at other alignment

                for (l3 = newTemplatePos; (l3 < target.size())&& (target[l3] == '-'); ++l) {
                }

                unsigned int l4 = l3 - newTemplatePos; // length of insertion of template

#ifdef DEBUG_VERBOSE
                cout << "Found insertion of length " << l2 << " at position " << i << " of target 2\n"
                        << "Found corresponding insertion of length " << l4 << " at position " << newTemplatePos << " of target 1" << endl;
#endif

                if (l4 < l2) // insert this many dashes
                {
                    unsigned int dl = l2 - l4;

                    for (unsigned int j = 0; j < dl; ++j)// insertion of other alignement is dl residues too short
                    {
#ifdef DEBUG_VERBOSE
                        cout << "Inserting dash into template 1 at position " << (newTemplatePos + 1) << endl;
#endif

                        insertDash(newTemplatePos + 1);
                    }

                    i += dl;
                } else
                    ++i;
            } else
                ++i;
        }

COMBINE_END:

#ifdef DEBUG_VERBOSE
        cout << "Result:\n"
                << target << "\n"
                << other.target << endl;
#endif

        ASSERT((target.size() == other.target.size()), exception);

        // Append other templates.
        for (unsigned int i = 0; i < other.seqTemplate.size(); ++i) {
            seqTemplateName.push_back(other.seqTemplateName[i]);
            seqTemplate.push_back(other.seqTemplate[i]);
        }

        // Do clean up (replace '~' with '-').
        for (unsigned int i = 0; i < target.size(); ++i)
            if (target[i] == '~')
                target[i] = '-';

        for (unsigned int j = 0; j < seqTemplate.size(); ++j)
            for (unsigned int i = 0; i < seqTemplate[j].size(); ++i)
                if (seqTemplate[j][i] == '~')
                    seqTemplate[j][i] = '-';
    }


    // HELPERS:
    /**
     * @description
     * @param s
     * @return 
     */
    string
    AlignmentBase::getPureSequence(const string &s) {
        string result = "";

        for (unsigned int i = 0; i < s.size(); ++i)
            if (s[i] != '-')
                result.append(1, s[i]);

        return result;
    }


    /// If original index points to a dash, it returns the position of the first
    /// non-dash left character, again taking '-' not into account.
    /**
     * @description
     * @param s
     * @param p
     * @return 
     */
    unsigned int
    AlignmentBase::getOrigPos(const string &s, unsigned int p) {
        PRECOND((p < s.size()), exception);

        int i;

        for (i = static_cast<int> (p); i >= 0; --i) // move pointer to first non-dash character
            if (s[i] != '-')
                break;

        if (i < 0) {
            //		ERROR("No defined original position.", exception);
            cout << "getOrigPos.i = " << i << endl;
            return 0; // left bound dashes
        }

        unsigned int counter = 0;

        for (int j = 0; j < i; ++j) // count the dashes before this position
            if (s[j] == '-')
                ++counter;

        return i - counter; // position minus number of dashes counted
    }
    /**
     * @description
     * @param s
     * @param p
     * @return 
     */
    unsigned int
    AlignmentBase::getNewPos(const string &s, unsigned int p) {
        PRECOND((p < s.size()), exception);

        int count = -1; // count non-dash characters
        int pi = static_cast<int> (p);

        for (unsigned int i = 0; i < s.size(); ++i)
            if (s[i] != '-') {
                ++count;
                if (count == pi)
                    return i;
            }

        ERROR("Could not assign new index for sequence with insertions.", exception);

        return s.size(); // dummy
    }
    /**
     * @description
     * @param text
     * @return 
     */
    vector<string>
    AlignmentBase::getTokens(const string &text) {
        istringstream ist(text.c_str());
        char *charLine = new char[text.size() + 1]; // size of string
        vector<string> v;
        string s;

        while (!ist.eof()) {
            ist >> charLine;
            s = charLine; // assignment of c-strings to string
            //		DUMP(s);
            if (s != "")
                v.push_back(s);
        }

        delete[] charLine;

        return v;
    }
    /**
     * @description
     * @param s
     * @param n
     * @return 
     */
    string
    AlignmentBase::deleteChar(const string &s, unsigned int n) {
        string result = s.substr(0, n);

        if (n < (s.size() - 1))
            result = result + s.substr(n + 1, s.size());

        return result;
    }

} // namespace
