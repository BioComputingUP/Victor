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

//
// Description:     Implement a simple alignment type.
//
// -----------------x-----------------------------------------------------------

#include <Alignment.h>
#include <AlignmentBase.h>
#include <IoTools.h>
#include <stringtools.h>
#include <String2Number.h>



namespace Victor { namespace Align2{

    // CONSTRUCTORS:

    Alignment::Alignment() : AlignmentBase(), score(), evalue() {
    }

    Alignment::Alignment(const Alignment &orig) {
        this->copy(orig);
    }

    Alignment::~Alignment() {
    }


    // OPERATORS:

    Alignment&
            Alignment::operator =(const Alignment &orig) {
        if (&orig != this)
            copy(orig);
        POSTCOND((orig == *this), exception);
        return *this;
    }


    // PREDICATES:

    void
    sSaveMSAF(string t, string tName, int offset, ostream &output) {
        string tNameLocal;
        if (tName != "")
            tNameLocal = tName;
        else
            tName = "?"; // give default dummy name

        output << tName << " " << (offset + 1);
        for (int i = 0; i < 14 - static_cast<int> (tName.length()); ++i)
            output << " ";
        output << t << endl;
    }
/**
 * 
 * @param output
 */
    void
    Alignment::saveMSAF(ostream &output) const {
        sSaveMSAF(target, targetName, startAaTarget, output);
        for (unsigned int j = 0; j < seqTemplate.size(); j++)
            sSaveMSAF(seqTemplate[j], seqTemplateName[j], startAaTemplates[j],
                output);
    }
/**
 * 
 * @param output
 */
    void
    Alignment::saveFasta(ostream &output) const {
        AlignmentBase::saveFasta(target, targetName, output);

        for (unsigned int j = 0; j < seqTemplate.size(); j++) {
            string tmp = seqTemplateName[j] + " // " + dtosDEF(score[j]);
            AlignmentBase::saveFasta(seqTemplate[j], tmp, output);
        }
    }


    // MODIFIERS:
/**
 * 
 * @param index1
 * @param index2
 */
    void
    Alignment::swapTemplate(unsigned int index1, unsigned int index2) {
        AlignmentBase::swapTemplate(index1, index2);
        swap(score[index1], score[index2]);
        swap(evalue[index1], evalue[index2]);
    }

    void
    Alignment::setTemplate(string t, string tName, double tScore, double tEvalue) {
        AlignmentBase::setTemplate(t, tName);
        score.push_back(tScore);
        evalue.push_back(tEvalue);
    }

    void
    Alignment::setScore(double val, unsigned int index) {
        if (index >= score.size())
            ERROR("Index out of scope.", exception);
        score[index] = val;
    }

    void
    Alignment::setEvalue(long double val, unsigned int index) {
        if (index >= evalue.size())
            ERROR("Index out of scope.", exception);
        evalue[index] = val;
    }
/**
 * 
 * @param index
 */
    void
    Alignment::cutTemplate(unsigned int index) {
        AlignmentBase::cutTemplate(index);
        score.resize(index);
        evalue.resize(index);
    }
/**
 * Output of alignment including headers of the 2 sequences.
 * @param os
 * @param headerTarget
 * @param startTarget
 * @param endTarget
 * @param headerTemplate
 * @param startTemplate
 * @param endTemplate
 * @param alignType
 * @param gapOpen
 * @param gapExtension
 * @param pWeigth
 * @param sWeigth
 * @param tWeigth
 * @param seqTarget
 * @param seqTemplate
 */
    void
    Alignment::doMatchPlusHeader(ostream &os, string headerTarget, int startTarget,
            int endTarget, string headerTemplate, int startTemplate, int endTemplate,
            string alignType, string gapOpen, string gapExtension, double pWeigth,
            double sWeigth, double tWeigth, string seqTarget,
            string seqTemplate) const {
        os << "# " << getScore() << "\n" // score
                << "# " << alignType << "\n" // alignType
                << "# " << gapOpen << "\n" // gapOpen
                << "# " << gapExtension << "\n" // gExtension
                << "# " << pWeigth << "\n" // profile
                << "# " << sWeigth << "\n" // secondary structure
                << "# " << tWeigth << endl; // threading

        os << "> " << headerTarget << " " << startTarget << " " << endTarget << "\n"
                << seqTarget << "\n"
                << "> " << headerTemplate << " " << startTemplate << " " << endTemplate << "\n"
                << seqTemplate << endl;
    }


    /// Removes illegal characters from sequence.

    static
    string sRemoveIllegalChars(string tmp) {
        string t2 = "";

        for (unsigned int i = 0; i < tmp.length(); i++)
            if ((!isspace(tmp[i])) && (tmp[i] != '*'))
                t2 += tmp[i];

        return t2;
    }
/**
 * Read FASTA format input.
 * @param input
 */
    void
    Alignment::loadFasta(istream &input) {
        string tmp;
        tmp = readLine(input);
        if (tmp[0] != '>')
            ERROR("Unrecognized start character while loading FASTA alignment.", exception);

        targetName = tmp.substr(1, tmp.length());
        target = "";

        while (input) {
            tmp = readLine(input); 

            if (tmp[0] != '>')
                target += sRemoveIllegalChars(tmp);
            else
                break;
        }

        unsigned int count = 0;

        if (!input)
            ERROR("Abnormal input file end while loading FASTA alignment.", exception);

        while (input) {
            // NB: this version may not work properly if 'title' of 2nd, 3rd, etc.
            // sequence is in the format "> blah" rather than ">blah".
            seqTemplateName.push_back(tmp.substr(1, tmp.length()));
            seqTemplate.push_back("");
            score.push_back(0);
            evalue.push_back(-1);
            startAaTemplates.push_back(0);

            while (input) {
                tmp = readLine(input);

                // Warning: this 'if' serves to remove a bug which causes the last
                // line in a file to be read twice.
                if (!input)
                    break;

                if (tmp[0] != '>')
                    seqTemplate[count] += sRemoveIllegalChars(tmp);
                else
                    break;
            }

            count++;
        }

        // Trim length to match all entries
        for (unsigned int i = target.length(); i < seqTemplate[0].length(); i++)
            target += '-';
        for (unsigned int j = 0; j < seqTemplate.size(); j++)
            for (unsigned int i = seqTemplate[j].length(); i < target.length(); i++)
                seqTemplate[j] += '-';
    }

    void
    Alignment::loadCEHeader(istream &input) {
        PRINT_NAME;

        string tmp;
        vector<string> words;
        tmp = readLine(input);

        targetName = tmp;
        target = "";

        do // read header
        {
            tmp = readLine(input);
            words = getTokens(tmp);

            if ((words.size() > 0) && (words[0] == "Alignment"))
                break; // end of "header info"

            if ((words.size() > 2) && (words[0] == "Chain")) {
                if (words[1] == "1:")
                    targetName = words[2];
                else
                    if (words[1] == "2:")
                    seqTemplateName.push_back(words[2]); // name of template
            }
        } while (input);

        if (!input)
            ERROR("Abnormal input file end.", exception);
    }
/**
 * @param input
 */
    void
    Alignment::loadCEBody(istream &input) {
        PRECOND(input, exception);
        PRINT_NAME;

        string tmp;
        vector<string> words;

        if (targetName == "")
            targetName = "?";
        target = "";

        seqTemplate.push_back(""); // start with empty template sequence
        startAaTemplates.push_back(0); // yet undefined
        score.push_back(0);
        evalue.push_back(-1);
        if (seqTemplateName.size() < seqTemplate.size())
            seqTemplateName.push_back("?");

        bool startAaTargetRead = false;
        bool startAaTemplateRead = false;
        int counter = 1;

        do // read body
        {
            ++counter;
            tmp = readLine(input);
            words = getTokens(tmp);

            if (words.size() > 0) {
                if (((words[0] == "Starting") || (words[0] == "Alignment")) && (target.size() > 0))
                    break; // only quit loop if some chain has been read
            }

            if ((words.size() > 3) && (words[0] == "Chain")) {
                if (words[1] == "1:") {
                    if (!startAaTargetRead) {
                        startAaTarget = stoiDEF(words[2]) - 1; // read offset
                        startAaTargetRead = true;
                    }
                    target = target + words[3]; // add sequence
                } else
                    if (words[1] == "2:") {
                    if (!startAaTemplateRead) {
                        startAaTemplates[0] = stoiDEF(words[2]) - 1; // read offset
                        startAaTemplateRead = true;
                    }
                    seqTemplate[0] = seqTemplate[0] + words[3]; // add sequence
                }
            }

            if ((words.size() > 0) && (words[0] == "X2"))
                break; // end of body. start of transformation data
        } while (input);

        /* here should be the code for parsing the transformation data... :-) */

    }

    void
    Alignment::loadCE(istream &input) {
        PRINT_NAME;

        loadCEHeader(input);
        loadCEBody(input);
    }

    vector< vector<int> >
    Alignment::loadMap(istream &is) {
        vector< vector<int> > mapFile;
        //	mapFile.reserve(2);
        vector<int> str_pos;
        vector<int> internal_pos;

        string strmapPos;
        is >> strmapPos;
        string seqmapPos;
        is >> seqmapPos;

        int cou = 0;

        while (is) {
            cou++;
            //		cout << "We are at round " << cou << endl;

            int tmp = -1000;
            is >> tmp;
            if (tmp == -1000)
                break;
            //		cout << tmp << " " << endl;
            str_pos.push_back(tmp);

            int tmp2 = -1000;
            is >> tmp2;
            //		cout << tmp2 << " "<< endl;
            internal_pos.push_back(tmp2);
            //		cout << endl;
        }

        //	cout << "end of the cycle" << endl;
        mapFile.push_back(str_pos);
        mapFile.push_back(internal_pos);

        return mapFile;
    }

    vector<int>
    Alignment::mapStructureCE2Sequence(vector< vector<int> > mapIn,
            istream &pdbstream, string startingChain, vector<int> shiftedSeqresEntries) {
        vector<int> cleanPos;
        //	string strmapPos;
        //	is >> strmapPos;
        //	string seqmapPos;
        //	is >> seqmapPos;
        //	int startingSeqaa = -1;
        int truePosition = 0;
        for (unsigned int i = 0; i < shiftedSeqresEntries.size(); i++) {
            int goodPosition = -2;

            if (shiftedSeqresEntries[i] == -1) {
                goodPosition = -1;
                cleanPos.push_back(goodPosition);
                continue;
            }

            // Translate SEQRES position to ATOM position using Dunbrack tool. This
            // is done because input fasta files are based on ATOM entry while CE
            // alignments show starting position based on SEQRES entry.
            truePosition = getTrueMapFromPdb(pdbstream, startingChain,
                    shiftedSeqresEntries[i]);

            //		while (is) // very dummy implementation! It reads the file each time i got a guess
            for (unsigned int j = 0; j < mapIn[0].size(); j++) {
                int tmp;
                //			is >> tmp;
                tmp = mapIn[0][j];
                if (tmp == truePosition) // if first position match, return the second
                {
                    goodPosition = mapIn[1][j];
                    break; // it breaks the while loop
                }
                //			int tmp2;
                //			is >> tmp2;
            }

            // I push back the position: either the good one or -2 if truePosition
            // is not present on the map file.
            cleanPos.push_back(goodPosition);
        }

        return cleanPos;
    }

    int
    Alignment::getTrueMapFromPdb(istream &is, string startingChain,
            int structEntryMapCe) {
        string chain_ID;
        chain_ID = startingChain;
        if (chain_ID == "-")
            chain_ID = "0";
        string line;

        int newPos = structEntryMapCe;

        while (getline(is, line)) {
            string::iterator str_iter = line.begin();
            if (string(str_iter, str_iter + 6) == "SEQCRD" and string(str_iter + 7, str_iter + 8) == chain_ID and strip_and_return_int(string(str_iter + 20, str_iter + 24)) == structEntryMapCe) {
                newPos = strip_and_return_int(string(str_iter + 26, str_iter + 31));
                break;
            }
        }

        return newPos;
    }

    string verifySeq(string sequence) {
        unsigned int i = 0;
        while (i < sequence.length()) {
            if (sequence[i] == ' ')
                sequence[i] = '-';
            i++;
        }
        return sequence;

    }

    void
    Alignment::loadPsiBlastMode4(istream& input) {
        // clear existing info
        targetName = "";
        unsigned int StateResult = 0;
        target = "";
        unsigned int templateCount = 0;
        seqTemplate.clear();
        seqTemplateName.clear();
        int iniseq = 13;
        unsigned int pos;
        int state = 1; // state of parser
        while (input) {
            string line = readLine(input);

            vector<string> words = getTokens(line);
            if ((words.size() == 0) || (words[0] == "Gapped") || (words[0] == "Effective")) {
                continue; // skip empty lines
            }
            if ((words[0] == "Sequences")) {
                state = 0;
                continue;
            }
            if (words[0] == "Lambda") {
                line = readLine(input);
                continue;
            }
            if ((words[0].find("Results") == 0)&&(StateResult == 0)) {
                StateResult = 1;
            }
            switch (state) {
                case 0: // read scores & E-values
                    if (words[0].find("|") <= 10) {
                        if ((words[words.size() - 2] == "0.0") && (stodDEF(words[words.size() - 1]) == 0.0)) {
                            score.push_back(stodDEF(words[words.size() - 3]));

                            if (words[words.size() - 2][0] == 'e')
                                evalue.push_back(stodDEF('1' + words[words.size() - 2]));
                            else {
                                evalue.push_back(stodDEF(words[words.size() - 2]));
                            }

                        } else {
                            score.push_back(stodDEF(words[words.size() - 2]));

                            if (words[words.size() - 1][0] == 'e')
                                evalue.push_back(stodDEF('1' + words[words.size() - 1]));
                            else {
                                evalue.push_back(stodDEF(words[words.size() - 1]));
                            }
                        }
                    }

                case 1: // looking for start of alignment
                    if (words[0] == "Query=") { // name of query found
                        if (words.size() > 1) {
                            targetName = words[1]; // set proper name
                        } // do not change state of parser
                    } else if ((words[0] == "QUERY") || (words[0] == "Query_1")) {

                        iniseq = line.find(words[1], words[0].length()) + 5;
                        state = 2;
                        templateCount = 0;
                        if (targetName == "") {
                            targetName = "QUERY";
                        }
                        if (words.size() == 4) {
                            //target = words[2];
                            startAaTarget = stoiDEF(words[1]) - 1; // read offset
                            pos = line.find(words[3], iniseq - 1);
                            target += line.substr(iniseq, pos - iniseq - 2);
                        } else if (words.size() == 2) {
                            target = words[1]; // line consists only of "-" characters
                            startAaTarget = 0; // no offset
                        } else {
                            ERROR("Strange line case 1.", exception);
                        }
                    }
                    break;
                case 2: // reading body of alignment
                    if ((words[0] == "QUERY") || (words[0] == "Query_1")) {
                        state = 2;
                        iniseq = line.find(words[1], words[0].length()) + 5; // 5 because the quantity of space to use in the offset
                        templateCount = 0;
                        if (words.size() == 4) {
                            //target += words[2];
                            pos = line.find(words[3], iniseq - 1);
                            target += line.substr(iniseq, pos - iniseq - 2);
                        } else if (words.size() == 2) {
                            target += words[1]; // line consists only of "-" characters
                        } else {
                            ERROR("Strange line. case ", exception);
                        }
                    } else if ((words[0].find("Searching..") == 0) || ((words[0].find("Results") == 0) && (StateResult > 0))) {
                        // must be psi-blast! This was not the last psi-blast round. 
                        // discard the stuff read so far and continue searching
                        state = 1;
                        templateCount = 0;
                        target = "";
                        seqTemplate.clear(); // clear existing info
                        seqTemplateName.clear();
                        seqTemplate.resize(0);
                        score.clear();
                        evalue.clear();
                    } else if (words[0] == "Database:") {
                        return; // done!
                    }// reading next template:
                    else if ((words.size() <= 4) && (words.size() >= 2)) {
                        string seq;
                        // size == 3 is a strange case: 
                        // happens (I think) only, if template sequence is like :   
                        // identifier 0 ----------------------------


                        if ((words.size() == 4) || (words.size() == 3)) {
                            //seq = words[2];


                            if (words.size() == 4) {

                                pos = line.find(words[3], iniseq - 1);
                            } else {
                                pos = 60 + iniseq;
                            }

                            seq = line.substr(iniseq, pos - iniseq - 2);
                            seq = verifySeq(seq);
                            // cout<<"seq:"<<line.substr(19,pos-19-2)<<"---Pos: "<<pos<< "\n";

                        } else {

                            seq = words[1];
                        }

                        if (templateCount < seqTemplate.size()) {
                            seqTemplate[templateCount] += seq;
                        } else {

                            seqTemplate.push_back(seq);
                            seqTemplateName.push_back(words[0]);
                            score.push_back(0);
                            evalue.push_back(-1);
                            if ((words.size() == 4) || (words.size() == 3)) {
                                startAaTemplates.push_back(stoiDEF(words[1]) - 1);
                            } else { // no offset
                                startAaTemplates.push_back(0);
                            }
                        }

                        ++templateCount;
                    } else {
                        ERROR("Strange line.", exception);
                    }
                    break;
                default: ERROR("Internal error line 442.", exception);
            }
        }
    }

    void
    Alignment::loadBlastMode6(istream &input) {
        // Clear existing info
        targetName = "";
        target = "";
        unsigned int templateCount = 0;
        seqTemplate.clear();
        seqTemplateName.clear();

        int state = 1; // state of parser

        while (input) {
            string line = readLine(input);
            vector<string> words = getTokens(line);

            if (words.size() == 0)
                continue; // skip empty lines

            if (words[0] == "Sequences") {
                state = 0;
                continue;
            }

            switch (state) {
                case 0: // read scores & E-values
                    if (words[0].find("|") <= 10) {
                        score.push_back(stodDEF(words[words.size() - 2]));
                        if (words[words.size() - 1][0] == 'e')
                            evalue.push_back(stodDEF('1' + words[words.size() - 1]));
                        else
                            evalue.push_back(stodDEF(words[words.size() - 1]));
                    }

                case 1: // looking for start of alignment
                    if (words[0] == "Query=") // name of query found
                    {
                        if (words.size() > 1)
                            targetName = words[1]; // set proper name
                        // do not change state of parser
                    } else
                        if (words[0] == "QUERY") {
                        state = 2;
                        templateCount = 0;

                        if (targetName == "")
                            targetName = "QUERY";

                        if (words.size() == 4) {
                            target = words[2];
                            startAaTarget = stoiDEF(words[1]) - 1; // read offset
                        } else
                            if (words.size() == 2) {
                            target = words[1]; // line consists only of "-" characters
                            startAaTarget = 0; // no offset
                        } else
                            ERROR("Strange line.", exception);
                    }
                    break;

                case 2: // reading body of alignment
                    if (words[0] == "QUERY") {
                        state = 2;
                        templateCount = 0;
                        if (words.size() == 4)
                            target += words[2];
                        else
                            if (words.size() == 2)
                            target += words[1]; // line consists only of "-" characters
                        else
                            ERROR("Strange line.", exception);
                    } else
                        if (words[0].find("Searching..") == 0) {
                        // Must be PSI-BLAST! This was not the last psi-blast round.
                        // Discard the stuff read so far and continue searching.
                        state = 1;
                        templateCount = 0;
                        target = "";
                        seqTemplate.clear(); // clear existing info
                        seqTemplateName.clear();
                        score.clear();
                        evalue.clear();
                    } else
                        if (words[0] == "Database:")
                        return; // done!
                    else
                        if ((words.size() <= 4) && (words.size() >= 2)) {
                        // Reading next template

                        string seq;
                        // size == 3 is a strange case:
                        // happens (I think) only, if template sequence is like:
                        // identifier 0 ----------------------------

                        if ((words.size() == 4) || (words.size() == 3))
                            seq = words[2];
                        else
                            seq = words[1];

                        if (templateCount < seqTemplate.size())
                            seqTemplate[templateCount] += seq;
                        else {
                            seqTemplate.push_back(seq);
                            seqTemplateName.push_back(words[0]);
                            score.push_back(0);
                            evalue.push_back(-1);

                            if ((words.size() == 4) || (words.size() == 3))
                                startAaTemplates.push_back(stoiDEF(words[1]) - 1);
                            else // no offset
                                startAaTemplates.push_back(0);
                        }
                        ++templateCount;
                    } else
                        ERROR("Strange line.", exception);
                    break;

                default:
                    ERROR("Internal error line 442.", exception);
            }
        }
    }

    void
    Alignment::loadBlastMode6FullSeq(istream &input, string masterTarget) {
        // clear existing info
        targetName = "";
        target = "";
        unsigned int templateCount = 0;
        seqTemplate.clear();
        seqTemplateName.clear();
        int state = 1; // state of parser
        string addNterm = ""; // oscar
        string addNtermTemp = ""; // oscar
        bool controlNterm = false; // oscar
        int firstTime = 0;

        while (input) {
            string line = readLine(input);
            vector<string> words = getTokens(line);
            if (words.size() == 0)
                continue; // skip empty lines

            if (words[0] == "Sequences") {
                state = 0;
                continue;
            }

            switch (state) {
                case 0: // read scores & E-values
                    if (words[0].find("|") <= 10) {
                        score.push_back(stodDEF(words[words.size() - 2]));
                        if (words[words.size() - 1][0] == 'e')
                            evalue.push_back(stodDEF('1' + words[words.size() - 1]));
                        else
                            evalue.push_back(stodDEF(words[words.size() - 1]));
                    }

                case 1: // looking for start of alignment
                    if (words[0] == "Query=") { // name of query found
                        if (words.size() > 1)
                            targetName = words[1]; // set proper name
                    } else // do not change state of parser
                        if (words[0] == "QUERY") {
                        state = 2;
                        templateCount = 0;

                        if (targetName == "")
                            targetName = "QUERY";

                        if (words.size() == 4) {
                            target = words[2];
                            startAaTarget = stoiDEF(words[1]) - 1; // read offset
                        } else
                            if (words.size() == 2) {
                            target = words[1]; // line consists only of "-" characters
                            startAaTarget = 0; // no offset
                        } else
                            ERROR("Strange line.", exception);

                        firstTime = 0; // case 1 means we are at N terminus so put firstTime = 0 (default condition)
                    }

                    if (startAaTarget > 0) // truncated N terminus : add missing information to all sequences // oscar
                    {
                        unsigned int countNterminus = startAaTarget; // oscar

                        for (unsigned int i = 0; i < (countNterminus); i++) // oscar
                        {
                            addNterm += masterTarget[i]; // oscar
                            target = addNterm + target; // oscar
                            controlNterm = true; // unlock this also for templates
                            for (unsigned int i = 0; i < addNterm.size(); i++)
                                addNtermTemp += "-";
                        }
                    }
                    break;

                case 2: // reading body of alignment
                    if (words[0] == "QUERY") {
                        state = 2;
                        templateCount = 0;
                        firstTime++;
                        controlNterm = false; // reset to default
                        if (words.size() == 4)
                            target += words[2];
                        else
                            if (words.size() == 2)
                            target += words[1]; // line consists only of "-" characters
                        else
                            ERROR("Strange line.", exception);
                    } else
                        if (words[0].find("Searching..") == 0) {
                        // Must be PSI-BLAST! This was not the last psi-blast round.
                        // Discard the stuff read so far and continue searching.
                        state = 1;
                        templateCount = 0;
                        target = "";
                        seqTemplate.clear(); // clear existing info
                        seqTemplateName.clear();
                        score.clear();
                        evalue.clear();
                    } else
                        if (words[0] == "Database:") {
                        string pure1 = getPureSequence(target);

                        if (pure1.size() != masterTarget.size()) // since N terminus is fixed any difference should occur at C terminus
                        {
                            string addCterm = "";
                            int countCterm = (masterTarget.size() - pure1.size());
                            int cTermLeftPos = (masterTarget.size() - countCterm); // the adjustment is made in order to fit in the "for" loop
                            for (unsigned int i = cTermLeftPos; i < masterTarget.size(); i++)
                                addCterm += masterTarget[i];

                            target = target + addCterm; // does it really update the target??
                            string addCtermTemp = "";
                            for (unsigned int i = 0; i < addCterm.size(); i++)
                                addCtermTemp += "-";

                            for (unsigned int i = 0; i < seqTemplate.size(); i++)
                                seqTemplate[i] = seqTemplate[i] + addCtermTemp;
                        }

                        return; // done! previously i put all the stuff to handle C terminus
                    }// reading next template:
                    else
                        if ((words.size() <= 4) && (words.size() >= 2)) {
                        string seq;
                        // size == 3 is a strange case:
                        // happens (I think) only, if template sequence is like:
                        // identifier 0 ----------------------------

                        if ((words.size() == 4) || (words.size() == 3))
                            seq = words[2];
                        else
                            seq = words[1];

                        if (templateCount < seqTemplate.size())
                            seqTemplate[templateCount] += seq;
                        else {
                            if ((controlNterm) and (firstTime == 0)) { // oscar: check if we have N terminus to add
                                seq = addNtermTemp + seq;
                                //										seqTemplate[templateCount] += addNterm;  // do the job!!
                            }

                            seqTemplate.push_back(seq);
                            seqTemplateName.push_back(words[0]);
                            score.push_back(0);
                            evalue.push_back(-1);

                            if ((words.size() == 4) || (words.size() == 3))
                                startAaTemplates.push_back(stoiDEF(words[1]) - 1);
                            else // no offset
                                startAaTemplates.push_back(0);
                        }

                        ++templateCount;
                    } else
                        ERROR("Strange line.", exception);

                    break;

                default:
                    ERROR("Internal error line 442.", exception);
            }
        }
    }

    void
    Alignment::loadMSAF(istream &is) {
        seqTemplate.clear();
        seqTemplateName.clear();
        startAaTemplates.clear();
        int startTemp = 0;
        vector<string> words = getTokens(readLine(is));

        if (words.size() != 3)
            ERROR("Strange first line in loadMSAF.", exception);

        targetName = words[0];
        startAaTarget = stoiDEF(words[1]) - 1;
        target = words[2];

        while (is) {
            words = getTokens(readLine(is));
            if (words.size() != 3)
                return;
            string tmpName, tmpSeq;
            is >> tmpName >> startTemp >> tmpSeq;
            seqTemplateName.push_back(tmpName);
            seqTemplate.push_back(tmpSeq);
            score.push_back(0);
            evalue.push_back(-1);
            startAaTemplates.push_back(startTemp);
        }
    }

    void
    Alignment::copy(const Alignment &orig) {
        PRINT_NAME;

        AlignmentBase::copy(orig);

        score.clear();
        evalue.clear();

        for (unsigned int i = 0; i < orig.seqTemplate.size(); i++) {
            score.push_back(orig.score[i]);
            evalue.push_back(orig.evalue[i]);
        }
    }

    void
    Alignment::addAlignment(const Alignment &orig) {
        AlignmentBase::addAlignment(orig);

        for (unsigned int i = 0; i < orig.score.size(); i++) {
            score.push_back(orig.score[i]);
            evalue.push_back(orig.evalue[i]);
        }
    }

    void
    Alignment::RemoveLowerSimple(double ID) {
        unsigned int i = 0;

        while (i < seqTemplate.size()) {
            double x = calculatePairwiseIdentity(target, seqTemplate[i]);
            if (x < ID)
                seqTemplate.erase(seqTemplate.begin() + i);
            if (x >= ID)
                i++;
        }
    }

    void
    Alignment::RemoveLowerAll(double ID) {
        /* confront between the target and all templates */
        unsigned int i = 0;

        while (i < seqTemplate.size()) {
            double x = calculatePairwiseIdentity(target, seqTemplate[i]);
            if (x < ID)
                seqTemplate.erase(seqTemplate.begin() + i);
            if (x >= ID)
                i++;
        }

        /* confront between all templates */
        unsigned int k = 0;
        unsigned int j = k + 1;

        while (k < seqTemplate.size()) {
            while (j < seqTemplate.size()) {
                double x = calculatePairwiseIdentity(seqTemplate[k], seqTemplate[j]);
                if (x < ID)
                    seqTemplate.erase(seqTemplate.begin() + j);
                if (x >= ID)
                    j++;
            }

            k++;
            j = k + 1;
        }
    }

    void
    Alignment::RemoveUpperSimple(double ID) {
        unsigned int i = 0;

        while (i < seqTemplate.size()) {
            double x = calculatePairwiseIdentity(target, seqTemplate[i]);
            if (x > ID)
                seqTemplate.erase(seqTemplate.begin() + i);
            if (x <= ID)
                i++;
        }
    }

    void
    Alignment::RemoveUpperAll(double ID) {
        /* confront between the target and all templates */
        unsigned int i = 0;

        while (i < seqTemplate.size()) {
            double x = calculatePairwiseIdentity(target, seqTemplate[i]);
            if (x > ID)
                seqTemplate.erase(seqTemplate.begin() + i);
            if (x <= ID)
                i++;
        }

        /* confront between all templates */
        unsigned int k = 0;
        unsigned int j = k + 1;

        while (k < seqTemplate.size()) {
            while (j < seqTemplate.size()) {
                double x = calculatePairwiseIdentity(seqTemplate[k], seqTemplate[j]);
                if (x > ID)
                    seqTemplate.erase(seqTemplate.begin() + j);
                if (x <= ID)
                    j++;
            }

            k++;
            j = k + 1;
        }
    }

}} // namespace
