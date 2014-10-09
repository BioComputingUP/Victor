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
/**
 * @Description LOop Building & Optimization This program generates an automatic model for a single indel
 */
#include <LoboTools.h>
#include <string>
#include <GetArg.h>
#include <LoopModel.h>
#include <LoopTableEntry.h>
#include <Spacer.h>
#include <PdbLoader.h>
#include <limits>
#include <RankAnalyzer.h>
#include <iostream>
#include <IoTools.h>
#include <PdbSaver.h>
#include <IntCoordConverter.h>
#include <VectorTransformation.h>



const unsigned int MAX_INDEL = 15;


using namespace Victor;

int main(int nArgs, char* argv[]) {
    if (getArg("h", nArgs, argv)) {
        showLoboAutoHelp();
        return 1;
    };
    if (getArg("-help", nArgs, argv)) {
        showMoreLoboAutoHelp();
        return 1;
    };

    // -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-
    // treat command-line arguments first:
    // -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-

    string inputFile, outputFile, sequenceFile, scwrlFile, scatterFile,
            tableFile, chainID;
    unsigned int pdbIndex1, pdbIndex2, pdbLength, num, num2, maxWrite;

    getArg("i", inputFile, nArgs, argv, "!");
    getArg("c", chainID, nArgs, argv, "!");
    getArg("o", outputFile, nArgs, argv, "test.pdb");
    getArg("-seq", sequenceFile, nArgs, argv, "!");
    getArg("-scwrl", scwrlFile, nArgs, argv, "!");
    getArg("s", pdbIndex1, nArgs, argv, 0);
    getArg("l", pdbLength, nArgs, argv, 0);
    getArg("-sol1", num, nArgs, argv, 60);
    getArg("-sol2", num2, nArgs, argv, 1);
    getArg("-maxWrite", maxWrite, nArgs, argv, 9999);
    getArg("-scatter", scatterFile, nArgs, argv, "!");
    getArg("-table", tableFile, nArgs, argv, "!");

    pdbIndex2 = pdbIndex1 + pdbLength;

    bool verbose = getArg("-verbose", nArgs, argv);
    bool test = getArg("-test", nArgs, argv);

    // variable defined in LoboTools.h:
    silent = getArg("-silent", nArgs, argv);

    LoopModel lm;

    if (scatterFile != "!") {
        ofstream scatFile(scatterFile.c_str());
        if (!scatFile)
            ERROR("Could not open file for writing. " + scatterFile, exception);
        lm.setScatterPlot(&scatFile);
    }

    if (tableFile != "!")
        lm.setTableFileName(tableFile);

    lm.setVerbose(2);

    double tmpEW = lm.getENDRMS_WEIGHT();
    getOption(tmpEW, "-endrmsWeigth", nArgs, argv);
    lm.setENDRMS_WEIGHT(tmpEW);
    treatSpecialOptions(nArgs, argv);

    getOption(LoopModel::OPT_MAX1, "-opt1", nArgs, argv);
    getOption(LoopModel::OPT_MAX2, "-opt2", nArgs, argv);
    getOption(LoopModel::OPT_NUM, "-optnum", nArgs, argv);
    vector<char> allCh;
    bool isDeletion = getOption("-del", nArgs, argv);
    bool isInsertion = getOption("-ins", nArgs, argv);
    if ((isInsertion && isDeletion) || (!isInsertion && !isDeletion))
        ERROR("Can only select either deletion or insertion not both.", exception);

    bool cluster = getOption("-cluster", nArgs, argv);
    bool norankRms = getOption("-norankRms", nArgs, argv);
    bool noFullModel = getOption("-noFullModel", nArgs, argv);
    bool noWrite = getOption("-noWrite", nArgs, argv);
    bool refine = getOption("-refine", nArgs, argv);
    bool optimize = getOption("-optimize", nArgs, argv);
    bool optall = getOption("-optall", nArgs, argv);
    bool withOxygen = getOption("-withOxygen", nArgs, argv);
    bool inter = getOption("-interpol", nArgs, argv);
    if (inter)
        lm.setInterpolated();

    if ((inputFile == "!") || (pdbIndex1 == 0) || (pdbIndex2 == 0)
            || ((isInsertion && (sequenceFile == "null")))) {
        cout << "Missing file specification. Aborting. (-h for help, --help "
                << "for long help)" << endl;
        return -1;
    }

    if (verbose)
        cout << "PDB-Index1 = " << pdbIndex1 << "\tPDB-Index2 = " << pdbIndex2
            << endl;

    // -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-
    // >>>>> END >>>>> treat command-line arguments
    // -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-


    if (!silent)
        cout << "Start, input: " << inputFile << " \n";
    Spacer *sp;

    ifstream inFile(inputFile.c_str());
    if (!inFile)
        ERROR("File not found.", exception);
    PdbLoader pl(inFile);
    pl.setNoHAtoms();
    allCh = pl.getAllChains();
    for (unsigned int i = 0; i < allCh.size(); i++)
        cout << "\t," << allCh[i] << ",";
    cout << "\n";

    /*check on validity of chain: 
    if user select a chain then check validity
     else select first valid one by default*/
    if (chainID != " ") {
        bool validChain = false;
        for (unsigned int i = 0; i < allCh.size(); i++) {
            if (allCh[i] == chainID[0]) {
                pl.setChain(chainID[0]);
                cout << "Loading chain " << chainID << "\n";
                validChain = true;
                break;
            }
        }
        if (!validChain) {
            cout << "Chain " << chainID << " is not available\n";
            return -1;
        }

    }


    pl.setPermissive();
    Protein prot;
    prot.load(pl);
    sp = prot.getSpacer(chainID[0]);

    cout << "Loaded protein." << endl;

    //find valid anchor residues:

    if (sp->isGap(pdbIndex1)) {

        unsigned int tmp = sp->getIndexFromPdbNumber(pdbIndex1) - 1;
        for (unsigned int i = 0; i < tmp; i++)
            if (!sp->isGap(pdbIndex1 - i)) {
                pdbIndex1 -= i;
                break;
            }

    }
    if (sp->isGap(pdbIndex2)) {
        unsigned int tmp = 0;
        if (sp->getIndexFromPdbNumber(pdbIndex2) < sp->sizeAmino())
            tmp = sp->sizeAmino() - sp->getIndexFromPdbNumber(pdbIndex2) - 1;

        for (unsigned int i = 0; i < tmp; i++) {
            if (!sp->isGap(pdbIndex2 + i)) {

                pdbIndex2 += i;

                break;
            }
        }
    }

    // if sequence file was passed, check if any AAs in the loop are missing:
    if (sequenceFile != "null") {
        // load seq file
        string seq = sGetSequence(sequenceFile, silent);

        // add missing AAs
        for (unsigned int i = (unsigned int) (pdbIndex1 + 1); i < (unsigned int) pdbIndex2; i++) {

            if (sp->isGap(i)) {
                string tmpType = oneLetter2ThreeLetter(seq[
                        i - sp->getStartOffset()]);

                if (sp->isGap(i - 1)) {
                    ERROR("Internal error adding missing amino acids.", exception);
                } else {
                    sp->insertAminoAfter(tmpType, sp->getIndexFromPdbNumber(i - 1));
                    sp->removeGap(i);
                }
            }
        }
    } else {// check that all loop AAs are valid: 

        for (unsigned int i = pdbIndex1 + 1; i < pdbIndex2; i++)
            if (sp->isGap(i))
                ERROR("Invalid loop residue specified without sequence template.",
                    exception);
    }

    // adapt indexes to internal AA count:
    unsigned int index1 = sp->getIndexFromPdbNumber(pdbIndex1);
    unsigned int index2 = sp->getIndexFromPdbNumber(pdbIndex2);
    unsigned int oldIndex1 = index1;

    if (verbose)
        cout << "Index1 = " << index1 << " Index2 = " << index2 << endl;

    if ((index1 >= sp->sizeAmino()) || (index2 >= sp->sizeAmino() - 1)
            || (index1 < 1) || (index2 < 1))
        ERROR("Invalid index encountered.", exception);

    if (sp->getAmino(index1)[C].distance(sp->getAmino(index1 + 1)[N]) == 0) {
        cout << "--> NB: Patching N-terminal anchor region.\n";
        index1--;
        if (index1 < 1)
            ERROR("Index is out of scope.", exception);
    };

    if (sp->getAmino(index2)[C].distance(sp->getAmino(oldIndex1)[C]) == 0) {
        cout << "--> NB: Patching C-terminal anchor region.\n";
        index2++;
        if (index2 >= sp->sizeAmino() - 1)
            ERROR("Index is out of scope.", exception);
    };

    // find ''optimal'' loop anchors:
    lm.defineLoopAnchors(*sp, index1, index2, isDeletion);

    treatProlineBug(index1, index2, *sp);

    if (index2 - index1 > MAX_INDEL) {
        // treat exceedingly large loops
        cout << "\nWARNING: Indel is too large (length = " << index2 - index1
                << ") - skipping...\n" << endl;
        exit(0);
    }

    IntCoordConverter icc;

    printSeqAndBfactors(index1, index2, *sp);
    fillLine();

    sp->setStateFromTorsionAngles();

    if (verbose)
        printSecStructureAndLoop(index1, index2, *sp);

    if (!test) {
        vector<string> typeVec;
        for (unsigned int i = index1 + 1; i < index2 + 1; i++)
            typeVec.push_back(sp->getAmino(i).getType());

        vector<Spacer> vsp;
        vsp = lm.createLoopModel(sp->getAmino(index1),
                sp->getAmino(index1 + 1)[N].getCoords(), sp->getAmino(index2),
                sp->getAmino(index2 + 1)[N].getCoords(), index1, index2, num,
                num2, typeVec);

        cout << "ranking...\n";

        if (refine) // refine all solutions
            lm.refineModel(*sp, index1, index2, vsp);

        if (optall) { // attempt local minimization of all solutions
            LoopModel::OPT_NUM = 1000;
            lm.optimizeModel(*sp, index1, index2, vsp);
        }

        if (!norankRms) {
            lm.rankRawScore(*sp, index1, index2, vsp, maxWrite);
            lm.doScatterPlot(*sp, index1, index2, vsp);
        }

        cout << "-----------------------------------------\n";
        cout << " Orig Loop energy = " << lm.calcOrigEnergy(*sp, index1, index2)
                << "  ( ?! ) \n";
        printResults(index1, index2, lm, *sp, vsp, withOxygen, maxWrite);

        if (cluster) {
            lm.clusterLoops(vsp);
            printResults(index1, index2, lm, *sp, vsp, withOxygen, maxWrite);
        }

        if (optimize) // attempt local minimization of top solution
            lm.optimizeModel(*sp, index1, index2, vsp);

        if (!noWrite)
            saveLoopEnsemble(index1, index2, lm, *sp, vsp, outputFile, maxWrite,
                noFullModel);

    } else {
        PdbSaver ps(cout);
        if (verbose)
            sp->save(ps);
    }

    // if SCWRL output file was requested, write out:
    if (scwrlFile != "null") {
        ofstream scwrlOut(scwrlFile.c_str());
        string tmp = lm.getSCWRLConservedSequence(*sp, index1, index2);
        scwrlOut << tmp << "\n";
    }

    fillLine();
    cout << endl;

    return 0;
}
