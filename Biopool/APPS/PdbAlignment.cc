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

/* 
 * Author: Stefano Cecconello
 *
 */

#include <cstdlib>
#include <SuperImpositor.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <iostream>
#include <GetArg.h>

using namespace std;
using namespace Victor;
using namespace Victor::Biopool;



//TODO:
//-SCRIVERE TEST
//-SCRIVERE la documentazione
//-SCRIVERE WIKI
//

void sShowHelp() {
    cout << "\n"
            << "   Options: \n"
            << "\t-i <filename>      \t\t Input file for PDB structure \n"
            << "\n"

            << "\t-r rmsd            \t\t run the rmsd algorithm and display the value \n"
            << "\t-m maxsub          \t\t run the maxsub algorithm and display the value \n"
            << "\t-g gdt             \t\t run the gdt algorithm and display the value \n"
            << "\t-t tmscore         \t\t run the tmscore algorithm and display the value \n"
            << "\n"

            << "\t-R rmsd            \t\t get in the output file the rmsd rototraslated proteins \n"
            << "\t-M maxsub          \t\t get in the output file the maxsub rototraslated proteins \n"
            << "\t-G gdt             \t\t get in the output file the gdt rototraslated proteins \n"
            << "\t-T tmscore         \t\t get in the output file the tmscore rototraslated proteins \n"
            << "\n"

            << "\t-x rotation method \t\t select the rotation method from: \n"
            << "\t                   \t\t -kabsch method \n"

            << "\n";
}

void savePdbOutput(vector <Spacer> spacers, string name);
string fromAlignmentToString(std::vector<std::pair<int, int> > range);
void saveAlignmentOutput(vector < std::vector<std::pair<int, int> > > align, string name);

int main(int nArgs, char* argv[]) {

    if (getArg("h", nArgs, argv)) {
        sShowHelp();
        return 1;
    };

    vector<string> inputFile;
    getArg("i", inputFile, nArgs, argv, "!");
    if (inputFile.size() != 2) {
        cout << "Error you need to get in input exactly 2 input file. Aborting. (-h for help)" << endl;
        return -1;
    }

    bool rmsd, maxsub, gdt, tmscore;
    bool rmsdOutput, maxsubOutput, gdtOutput, tmscoreOutput;

    rmsd = getArg("r", nArgs, argv);
    maxsub = getArg("m", nArgs, argv);
    gdt = getArg("g", nArgs, argv);
    tmscore = getArg("t", nArgs, argv);

    rmsdOutput = getArg("R", nArgs, argv);
    maxsubOutput = getArg("M", nArgs, argv);
    gdtOutput = getArg("G", nArgs, argv);
    tmscoreOutput = getArg("T", nArgs, argv);

    if (rmsd || maxsub || gdt || tmscore || rmsdOutput || maxsubOutput || gdtOutput || tmscoreOutput) {
        string rotationMethod;
        getArg("x", rotationMethod, nArgs, argv, "");

        ifstream inFile1(inputFile[0].c_str());
        ifstream inFile2(inputFile[1].c_str());

        // creates the PdbLoader objects
        PdbLoader pl1(inFile1, false, false, false, false, true);
        PdbLoader pl2(inFile2, false, false, false, false, true);

        Protein* prot1 = new Protein();
        Protein* prot2 = new Protein();

        pl1.setNoVerbose();
        pl2.setNoVerbose();
        pl1.setNoHAtoms();
        pl2.setNoHAtoms();

        prot1->load(pl1);
        prot2->load(pl2);


        SuperImpositor* superImpositor = new SuperImpositor(prot1, prot2, rotationMethod);

        if (rmsd) {
            superImpositor->calculateRMSD();
            double rmsd = superImpositor->getRmsdValue();
            if (!(rmsdOutput || maxsubOutput || gdtOutput || tmscoreOutput)) {
                vector <Spacer> spacers;
                Spacer newSet1 = superImpositor->getRMSDset1();
                Spacer newSet2 = superImpositor->getRMSDset2();

                spacers.push_back(newSet1);
                spacers.push_back(newSet2);
                savePdbOutput(spacers, "rmsd");
            }
            cout << "The rmsd value is': " << rmsd << "\n";
        }

        if (rmsdOutput) {
            superImpositor->calculateRMSD();
            vector <Spacer> spacers;
            Spacer newSet1 = superImpositor->getRMSDset1();
            Spacer newSet2 = superImpositor->getRMSDset2();

            spacers.push_back(newSet1);
            spacers.push_back(newSet2);
            savePdbOutput(spacers, "rmsd");
        }

        if (maxsub) {
            std::vector<std::pair<int, int> > vectorSet;

            for (unsigned int i = 0; i < superImpositor->getSet1()->size(); i++) {
                vectorSet.push_back(std::make_pair(i, i));
            }
            superImpositor->calculateMaxSub(3.5, vectorSet, 'y');
            double maxSub = superImpositor->getMaxsubValue();


            if (!(rmsdOutput || maxsubOutput || gdtOutput || tmscoreOutput)) {
                std::vector<std::pair<int, int> > range;
                range = superImpositor->getMaxsubAlignment();
                vector < std::vector<std::pair<int, int> > > align;
                align.push_back(range);
                saveAlignmentOutput(align, "maxsub");
                vector <Spacer> spacers;
                Spacer newSet1 = superImpositor->getMaxSubset1();
                Spacer newSet2 = superImpositor->getMaxSubset2();
                spacers.push_back(newSet1);
                spacers.push_back(newSet2);
                savePdbOutput(spacers, "maxsub");
            }
            cout << "The maxsub value is: " << maxSub << "\n";
        }

        if (maxsubOutput) {
            std::vector<std::pair<int, int> > vectorSet;

            for (unsigned int i = 0; i < superImpositor->getSet1()->size(); i++) {
                vectorSet.push_back(std::make_pair(i, i));
            }
            superImpositor->calculateMaxSub(3.5, vectorSet, 'y');

            std::vector<std::pair<int, int> > range;
            range = superImpositor->getMaxsubAlignment();
            vector < std::vector<std::pair<int, int> > > align;
            align.push_back(range);
            saveAlignmentOutput(align, "maxsub");


            vector <Spacer> spacers;
            Spacer newSet1 = superImpositor->getMaxSubset1();
            Spacer newSet2 = superImpositor->getMaxSubset2();
            spacers.push_back(newSet1);
            spacers.push_back(newSet2);
            savePdbOutput(spacers, "maxsub");
        }

        if (gdt) {
            std::vector<std::pair<int, int> > vectorSet;

            for (unsigned int i = 0; i < superImpositor->getSet1()->size(); i++) {
                vectorSet.push_back(std::make_pair(i, i));
            }
            superImpositor->calculateGdt(vectorSet);
            double gdt = superImpositor->getGdtValue();

            std::vector<std::pair<int, int> > range1;
            std::vector<std::pair<int, int> > range2;
            std::vector<std::pair<int, int> > range3;
            std::vector<std::pair<int, int> > range4;
            range1 = superImpositor->getGdtAlignment1();
            range2 = superImpositor->getGdtAlignment2();
            range3 = superImpositor->getGdtAlignment3();
            range4 = superImpositor->getGdtAlignment4();


            if (!(rmsdOutput || maxsubOutput || gdtOutput || tmscoreOutput)) {
                vector < std::vector<std::pair<int, int> > > align;
                align.push_back(range1);
                align.push_back(range2);
                align.push_back(range3);
                align.push_back(range4);
                saveAlignmentOutput(align, "gdt");
                vector <Spacer> spacers;
                Spacer newSet1 = superImpositor->getGdtset1_1();
                Spacer newSet2 = superImpositor->getGdtset1_2();

                Spacer newSet3 = superImpositor->getGdtset2_1();
                Spacer newSet4 = superImpositor->getGdtset2_2();

                Spacer newSet5 = superImpositor->getGdtset3_1();
                Spacer newSet6 = superImpositor->getGdtset3_2();

                Spacer newSet7 = superImpositor->getGdtset4_1();
                Spacer newSet8 = superImpositor->getGdtset4_2();

                spacers.push_back(newSet1);
                spacers.push_back(newSet2);
                spacers.push_back(newSet3);
                spacers.push_back(newSet4);
                spacers.push_back(newSet5);
                spacers.push_back(newSet6);
                spacers.push_back(newSet7);
                spacers.push_back(newSet8);
                savePdbOutput(spacers, "gdt");
            }
            cout << "The gdt value is: " << gdt << "\n";
        }

        if (gdtOutput) {
            std::vector<std::pair<int, int> > vectorSet;

            for (unsigned int i = 0; i < superImpositor->getSet1()->size(); i++) {
                vectorSet.push_back(std::make_pair(i, i));
            }
            superImpositor->calculateGdt(vectorSet);


            std::vector<std::pair<int, int> > range1;
            std::vector<std::pair<int, int> > range2;
            std::vector<std::pair<int, int> > range3;
            std::vector<std::pair<int, int> > range4;
            range1 = superImpositor->getGdtAlignment1();
            range2 = superImpositor->getGdtAlignment2();
            range3 = superImpositor->getGdtAlignment3();
            range4 = superImpositor->getGdtAlignment4();

            vector < std::vector<std::pair<int, int> > > align;
            align.push_back(range1);
            align.push_back(range2);
            align.push_back(range3);
            align.push_back(range4);
            saveAlignmentOutput(align, "gdt");


            vector <Spacer> spacers;
            Spacer newSet1 = superImpositor->getGdtset1_1();
            Spacer newSet2 = superImpositor->getGdtset1_2();

            Spacer newSet3 = superImpositor->getGdtset2_1();
            Spacer newSet4 = superImpositor->getGdtset2_2();

            Spacer newSet5 = superImpositor->getGdtset3_1();
            Spacer newSet6 = superImpositor->getGdtset3_2();

            Spacer newSet7 = superImpositor->getGdtset4_1();
            Spacer newSet8 = superImpositor->getGdtset4_2();

            spacers.push_back(newSet1);
            spacers.push_back(newSet2);
            spacers.push_back(newSet3);
            spacers.push_back(newSet4);
            spacers.push_back(newSet5);
            spacers.push_back(newSet6);
            spacers.push_back(newSet7);
            spacers.push_back(newSet8);
            savePdbOutput(spacers, "gdt");
        }

        if (tmscore) {
            std::vector<std::pair<int, int> > vectorSet;

            for (unsigned int i = 0; i < superImpositor->getSet1()->size(); i++) {
                vectorSet.push_back(std::make_pair(i, i));
            }
            superImpositor->calculateTMScore(vectorSet);
            double TMScore = superImpositor->getTMScoreValue();


            if (!(rmsdOutput || maxsubOutput || gdtOutput || tmscoreOutput)) {
                std::vector<std::pair<int, int> > range;
                range = superImpositor->getTMScoreAlignment();
                vector < std::vector<std::pair<int, int> > > align;
                align.push_back(range);
                saveAlignmentOutput(align, "TMScore");

                vector <Spacer> spacers;
                Spacer newSet1 = superImpositor->getTMScoreset1();
                Spacer newSet2 = superImpositor->getTMScoreset2();

                spacers.push_back(newSet1);
                spacers.push_back(newSet2);
                savePdbOutput(spacers, "TMScore");
            }
            cout << "The TMScore value is: " << TMScore << "\n";
        }

        if (tmscoreOutput) {
            std::vector<std::pair<int, int> > vectorSet;

            for (unsigned int i = 0; i < superImpositor->getSet1()->size(); i++) {
                vectorSet.push_back(std::make_pair(i, i));
            }
            superImpositor->calculateTMScore(vectorSet);

            std::vector<std::pair<int, int> > range;
            range = superImpositor->getTMScoreAlignment();
            vector < std::vector<std::pair<int, int> > > align;
            align.push_back(range);
            saveAlignmentOutput(align, "TMScore");


            vector <Spacer> spacers;
            Spacer newSet1 = superImpositor->getTMScoreset1();
            Spacer newSet2 = superImpositor->getTMScoreset2();

            spacers.push_back(newSet1);
            spacers.push_back(newSet2);
            savePdbOutput(spacers, "TMScore");
        }


    } else {
        cout << "No method select (-h for help)\n";
    }

    return 0;
}

void savePdbOutput(vector <Spacer> spacers, string name) {
    string proteineOUTPUT = name + "output.pdb";
    ofstream outFile(proteineOUTPUT.c_str());

    PdbSaver saveSet(outFile);
    for (unsigned int i = 0; i < spacers.size(); i++) {
        saveSet.saveSpacer(spacers[i]);
    }
}

void saveAlignmentOutput(vector < std::vector<std::pair<int, int> > > align, string name) {
    string proteineOUTPUT = name + "AlignmentOutput";
    ofstream outFile(proteineOUTPUT.c_str());

    for (unsigned int i = 0; i < align.size(); i++) {
        outFile << fromAlignmentToString(align[i]);
    }
}

string fromAlignmentToString(std::vector<std::pair<int, int> > range) {
    string String;

    string str1;
    string str2;
    for (unsigned int i = 0; i < range.size(); i++) {
        stringstream ss1;
        stringstream ss2;
        ss1 << range[i].first;
        ss2 << range[i].second;
        str1 = ss1.str();
        str2 = ss2.str();
        String = String + "(" + str1 + "," + str2 + ")\n";
    }
    String = String + "\n";
    return String;
}