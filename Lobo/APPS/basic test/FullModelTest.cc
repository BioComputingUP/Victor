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
 * @Description This program generates all models for a single protein over a sliding window.
 */
#include <string>
#include <GetArg.h>
#include <LoopModel.h>
#include <LoopTableEntry.h>
#include <Spacer.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <IntCoordConverter.h>
#include <VectorTransformation.h>
#include <StatTools.h>
#include <XyzSaver.h>

#include <limits.h>    
using namespace Victor;

void sShowHelp() {
    LoopModel lm;
    cout << "Full Model Test\n"
            << "This program generates all models for a single protein "
            << " over a sliding window.\n"
            << " Options: \n"
            << "\t-i <filename> \t\t Input file\n"
            << "\t[-o <filename>] \t\t Output file basename\n"
            << "\t-s <offset> \t\t Size of *real* sliding window (default= 5)\n"
            << "\t[--sol1 <n-solutions>] \t Number of primary solutions"
            << " (default= " << LoopModel::MAX_ITER_SOL << ")\n"
            << "\t[--sol2 <n-solutions>] \t Number of secondary solutions"
            << " (default= 1)\n"
            << "\t[--maxWrite <n-solutions>] \t Number of solutions to write out"
            << " (default = sol1)\n"
            << "\t[--norank] \t\t Do not rank or filter results.\n"
            << "\t[--norankRms] \t\t Do not rank or filter results according "
            << "to end AA RMS.\n"
            << "\t[--fullRank] \t\t Write full ranking.\n"
            << "\t[--scwrl] \t\t Sequence output file for SCWRL (-s option)\n"
            << "\t[--endrmsWeigth] \t Choose weigthing of endrms"
            << " (default= " << lm.getENDRMS_WEIGHT() << ")\n"
            << "\t[--energyWeigth] \t Choose weigthing of energy"
            << " (default= " << LoopModel::ENERGY_WEIGTH << ")\n"
            << "\t[--secprefWeigth] \t Choose weigthing of secondary preference"
            << " (default= " << LoopModel::SECPREF_WEIGTH << ")\n"
            << "\t[--secprefTol] \t\t Choose tolerance for secondary preference"
            << " (default= " << LoopModel::SECPREF_TOL << ")\n"
            << "\t[--packingWeigth] \t Choose weigthing of packing density"
            << " (default= " << LoopModel::PACKING_WEIGTH << ")\n"
            << "\t[--weigthEP] \t\t Choose weigthing of EP compared to ED & EN"
            << " (default= " << LoopTableEntry::LAMBDA_EP << ")\n"
            << "\t[--maxSearch] \t\t Choose max fraction of tables to search "
            << "for best result (default= " << LoopTable::MAX_FACTOR << ")\n"
            << "\t[--vdwLimit] \t\t Choose threshold for VDW filter"
            << " (default= " << LoopModel::VDW_LIMIT << ")\n"
            << "\t[--energyLimit] \t Choose threshold for energy filter"
            << " (default= " << LoopModel::ENERGY_LIMIT << ")\n"
            << "\t[--simLimit] \t\t Choose threshold for similarity filter"
            << " (default= " << LoopModel::SIM_LIMIT << ")\n"
            << "\t[--useZScore]\t\t Use Z score ranking scheme instead of "
            << "raw scores.\n"
            << "\t[--withOxygen] \t\t Include Oxygen atoms in RMSD calculation\n"
            << "\t[--interpol] \t\t Use interpolated RAPDF energy\n"
            << "\t[--scatter <filename>] \t Write a scatter plot file.\n";
}

void sGetOption(double& param, char* optName, int nArgs, char* argv[], bool verbose = true) {
    double tmp = -100000.0;
    getArg(optName, tmp, nArgs, argv, -100000.0);
    if (tmp > -100000.0) {
        param = tmp;
        if (verbose)
            cout << optName << " = " << tmp << "\n";
    }
}

void sGetOption(unsigned int& param, char* optName, int nArgs, char* argv[], bool verbose = true) {
    unsigned int tmp = INT_MAX;
    getArg(optName, tmp, nArgs, argv, INT_MAX);
    if (tmp != INT_MAX) {
        param = tmp;
        if (verbose)
            cout << optName << " = " << tmp << "\n";
    }
}

void sGetOption(float& param, char* optName, int nArgs, char* argv[], bool verbose = true) {
    float tmp = -100000.0;
    getArg(optName, tmp, nArgs, argv, -100000.0);
    if (tmp > -100000.0) {
        param = tmp;
        if (verbose)
            cout << optName << " = " << tmp << "\n";
    }
}

int main(int nArgs, char* argv[]) {
    if (getArg("h", nArgs, argv)) {
        sShowHelp();
        return 1;
    };

    string inputFile, outputFile, scwrlFile, scatterFile, chainID;
    vector<char> allCh;
    unsigned int windowSize, num, num2, maxWrite;
    getArg("i", inputFile, nArgs, argv, "!");
    getArg("c", chainID, nArgs, argv, " ");
    getArg("o", outputFile, nArgs, argv, "null");
    getArg("s", windowSize, nArgs, argv, 5);
    windowSize++; // compensate counting scheme of LoopModel
    getArg("-sol1", num, nArgs, argv, LoopModel::MAX_ITER_SOL);
    getArg("-sol2", num2, nArgs, argv, 1);

    getArg("-scwrl", scwrlFile, nArgs, argv, "null");
    getArg("-scatter", scatterFile, nArgs, argv, "!");
    getArg("-maxWrite", maxWrite, nArgs, argv, 9999);

    ofstream scatFile(scatterFile.c_str());
    if (!scatFile)
        ERROR("File not found.", exception);

    LoopModel lm;

    lm.setVerbose();

    double tmpEW = lm.getENDRMS_WEIGHT();
    sGetOption(tmpEW, "-endrmsWeigth", nArgs, argv);
    lm.setENDRMS_WEIGHT(tmpEW);
    sGetOption(LoopModel::ENERGY_WEIGTH, "-energyWeigth", nArgs, argv);
    sGetOption(LoopModel::SECPREF_WEIGTH, "-secprefWeigth", nArgs, argv);
    sGetOption(LoopModel::SECPREF_TOL, "-secprefTol", nArgs, argv);
    sGetOption(LoopModel::PACKING_WEIGTH, "-packingWeigth", nArgs, argv);
    sGetOption(LoopTableEntry::LAMBDA_EP, "-weigthEP", nArgs, argv);
    sGetOption(LoopTableEntry::LAMBDA_ED, "-weigthED", nArgs, argv);
    sGetOption(LoopTableEntry::LAMBDA_EN, "-weigthEN", nArgs, argv);
    sGetOption(LoopTable::MAX_FACTOR, "-maxSearch", nArgs, argv);
    sGetOption(LoopModel::VDW_LIMIT, "-vdwLimit", nArgs, argv);
    sGetOption(LoopModel::ENERGY_LIMIT, "-energyLimit", nArgs, argv);
    sGetOption(LoopModel::SIM_LIMIT, "-simLimit", nArgs, argv);

    bool norank = getArg("-norank", nArgs, argv);
    if (norank)
        cout << "norank\n";

    bool norankRms = getArg("-norankRms", nArgs, argv);
    if (norankRms)
        cout << "norankRms\n";

    bool fullRank = getArg("-fullRank", nArgs, argv);
    if (fullRank)
        cout << "fullRank\n";

    bool useZScore = getArg("-useZScore", nArgs, argv);
    if (useZScore)
        cout << "useZScore\n";

    bool withOxygen = getArg("-withOxygen", nArgs, argv);
    if (withOxygen)
        cout << "withOxygen\n";

    bool inter = getArg("-interpol", nArgs, argv);
    if (inter) {
        cout << "interpol\n";
        lm.setInterpolated();
    }

    if (inputFile == "!") {
        cout << "Missing file specification. Aborting. (-h for help)" << endl;
        return -1;
    }

    cout << "Start\n";
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
    cout << "loaded.\n";

    IntCoordConverter icc;

    vector<double> counterMIX;
    vector<double> counterHELIX; // store single RMSDs
    vector<double> counterSTRAND;

    vector<double> min_counterMIX;
    vector<double> min_counterHELIX; // store single minimum RMSDs
    vector<double> min_counterSTRAND;

    for (unsigned int index1 = 2; index1 < sp->sizeAmino()-(windowSize + 2);
            index1++) {
        // Main iteration loop:
        unsigned int index2 = index1 + windowSize;

        cout << "Running from " << setw(3) << index1 + 1 << " to "
                << setw(3) << index2 + 1 << "\n";
        cout << "-----------------------------------------\n";

        vector<string> typeVec;

        for (unsigned int i = index1 + 1; i < index2 + 1; i++)
            typeVec.push_back(sp->getAmino(i).getType());

        vector<Spacer > vsp;
        vsp = lm.createLoopModel(sp->getAmino(index1),
                sp->getAmino(index1 + 1)[N].getCoords(), sp->getAmino(index2),
                sp->getAmino(index2 + 1)[N].getCoords(), index1, index2, num,
                num2, typeVec);




        cout << "-----------------------------------------\n";

        double tmp = 0.0;
        double tmpMin = 0.0;

        if (fullRank)
            for (unsigned int i = 0; i < vsp.size(); i++) {
                cout << setw(3) << i << "   ";
                lm.calculateRms(*sp, index1, index2, vsp[i], true,
                        withOxygen);
            } else
            if (vsp.size() > 0) {
            lm.calculateRms(*sp, index1, index2, vsp[0], true, withOxygen);
        }

        if (vsp.size() > 0) {
            tmp = lm.calculateRms2(*sp, index1, index2, vsp[0], false,
                    withOxygen);

            tmpMin = tmp;

            for (unsigned int i = 1; i < vsp.size(); i++) {
                double tmp2 = lm.calculateRms2(*sp, index1, index2, vsp[i],
                        false, withOxygen);
                if (tmp2 < tmpMin)
                    tmpMin = tmp2;
            }
        }

        // choose class:

        bool hel = true;
        for (unsigned int i = index1; i < index2; i++)
            if (sp->getAmino(i).getState() != HELIX) {
                hel = false;
                break;
            }

        bool str = false;
        if (!hel) {
            str = true;
            for (unsigned int i = index1; i < index2; i++)
                if (sp->getAmino(i).getState() != STRAND) {
                    str = false;
                    break;
                }
        }

        string segName = "."; // name of current offset block, 
        // used for output filenanme extension

        if (hel == true) {
            segName = segName + "HEL";
            cout << ">>> HELIX segment.\n";
            counterHELIX.push_back(tmp);
            min_counterHELIX.push_back(tmpMin);
        } else if (str == true) {
            segName = segName + "STR";
            cout << ">>> STRAND segment.\n";
            counterSTRAND.push_back(tmp);
            min_counterSTRAND.push_back(tmpMin);
        } else {
            segName = segName + "MIX";
            cout << ">>> MIXED segment.\n";
            counterMIX.push_back(tmp);
            min_counterMIX.push_back(tmpMin);
        }

        segName = segName + "_";

        if (index1 + 1 < 10)
            segName += "00";
        else if (index1 + 1 < 100)
            segName += "0";
        segName += itos(index1 + 1) + "_";

        if (index2 + 1 < 10)
            segName += "00";
        else if (index2 + 1 < 100)
            segName += "0";
        segName += itos(index2 + 1);

        // write output PDB files

        if (outputFile != "null") {
            unsigned int maxW = vsp.size() < maxWrite ? vsp.size() : maxWrite;
            for (unsigned int i = 0; i < maxW; i++) {
                string tmpStr = outputFile + segName + ".pdb.";

                Spacer *sp2; // set spacer for output
                sp2 = sp;

                if (i < 10)
                    tmpStr += "00";
                else if (i < 100)
                    tmpStr += "0";

                tmpStr += itos(i);

                ofstream outFile(tmpStr.c_str());
                if (!outFile)
                    ERROR("Could not create file.", exception);
                PdbSaver ps(outFile);

                lm.setStructure(*sp2, vsp[i], index1, index2);
                sp2->save(ps);
            };
        }

        // if SCWRL output file was requested, write out:
        if (scwrlFile != "null") {
            string tmp = scwrlFile + segName + ".seq"; // construct correct name
            ofstream scwrlOut(tmp.c_str());
            tmp = lm.getSCWRLConservedSequence(*sp, index1, index2);
            scwrlOut << tmp << "\n";
        }

        cout << "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n";

    }; // for

    // calculate statistics (average + standard deviation) of solutions:

    double avgHELIX = average(counterHELIX);
    double sdHELIX = standardDeviation(counterHELIX, avgHELIX);

    double min_avgHELIX = average(min_counterHELIX);
    double min_sdHELIX = standardDeviation(min_counterHELIX, min_avgHELIX);

    double avgSTRAND = average(counterSTRAND);
    double sdSTRAND = standardDeviation(counterSTRAND, avgSTRAND);

    double min_avgSTRAND = average(min_counterSTRAND);
    double min_sdSTRAND = standardDeviation(min_counterSTRAND, min_avgSTRAND);

    double avgMIX = average(counterMIX);
    double sdMIX = standardDeviation(counterMIX, avgMIX);

    double min_avgMIX = average(min_counterMIX);
    double min_sdMIX = standardDeviation(min_counterMIX, min_avgMIX);

    vector<double> counterALL;
    for (unsigned int i = 0; i < counterMIX.size(); i++)
        counterALL.push_back(counterMIX[i]);
    for (unsigned int i = 0; i < counterHELIX.size(); i++)
        counterALL.push_back(counterHELIX[i]);
    for (unsigned int i = 0; i < counterSTRAND.size(); i++)
        counterALL.push_back(counterSTRAND[i]);

    double avgALL = average(counterALL);
    double sdALL = standardDeviation(counterALL, avgALL);

    vector<double> min_counterALL;
    for (unsigned int i = 0; i < min_counterMIX.size(); i++)
        min_counterALL.push_back(min_counterMIX[i]);
    for (unsigned int i = 0; i < min_counterHELIX.size(); i++)
        min_counterALL.push_back(min_counterHELIX[i]);
    for (unsigned int i = 0; i < min_counterSTRAND.size(); i++)
        min_counterALL.push_back(min_counterSTRAND[i]);

    double min_avgALL = average(min_counterALL);
    double min_sdALL = standardDeviation(min_counterALL, min_avgALL);

    for (unsigned int i = 0; i < sp->sizeAmino(); i++) {
        if (sp->getAmino(i).getState() == HELIX)
            cout << "H";
        else if (sp->getAmino(i).getState() == STRAND)
            cout << "E";
        else
            cout << "-";
        if ((i + 1) % 60 == 0)
            cout << "\n";
    }
    cout << "\n";
    cout << "xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx\n";
    cout << "TOP 1 Performance\n";
    cout << "--x--x--x--x--x--x--x--x--x--\n";
    cout << "Class = HELIX:\n";
    cout << "Average RMSD = " << setw(5) << avgHELIX << "\t (sd = " << setw(5)
            << sdHELIX << " )\t # counts = " << setw(3) << counterHELIX.size()
            << "\n";
    cout << "Class = STRAND:\n";
    cout << "Average RMSD = " << setw(5) << avgSTRAND << "\t (sd = " << setw(5)
            << sdSTRAND << " )\t # counts = " << setw(3) << counterSTRAND.size()
            << "\n";
    cout << "Class = MIXTURE:\n";
    cout << "Average RMSD = " << setw(5) << avgMIX << "\t (sd = " << setw(5)
            << sdMIX << " )\t # counts = " << setw(3) << counterMIX.size() << "\n";
    cout << "Class = ALL:\n";
    cout << "Average RMSD = " << setw(5) << avgALL << "\t (sd = " << setw(5)
            << sdALL << " )\t # counts = " << setw(3)
            << counterALL.size()
            << "\n";
    cout << "xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx\n";
    cout << "Theoretical Best Performance\n";
    cout << "--x--x--x--x--x--x--x--x--x--\n";
    cout << "Class = HELIX:\n";
    cout << "Average RMSD = " << setw(5) << min_avgHELIX << "\t (sd = "
            << setw(5) << min_sdHELIX << " )\t # counts = " << setw(3)
            << min_counterHELIX.size() << "\n";
    cout << "Class = STRAND:\n";
    cout << "Average RMSD = " << setw(5) << min_avgSTRAND << "\t (sd = "
            << setw(5) << min_sdSTRAND << " )\t # counts = " << setw(3)
            << min_counterSTRAND.size() << "\n";
    cout << "Class = MIXTURE:\n";
    cout << "Average RMSD = " << setw(5) << min_avgMIX << "\t (sd = " << setw(5)
            << min_sdMIX << " )\t # counts = " << setw(3) << min_counterMIX.size()
            << "\n";
    cout << "Class = ALL:\n";
    cout << "Average RMSD = " << setw(5) << min_avgALL << "\t (sd = " << setw(5)
            << min_sdALL << " )\t # counts = " << setw(3)
            << min_counterALL.size()
            << "\n";
    cout << "xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx-xx\n";
    cout << endl;

    return 0;
}

