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
 * @Description This program allows to edit scatter plots from loop modeling
 */
#include <RankAnalyzer.h>
#include <GetArg.h>
#include <iostream>
#include <IoTools.h>

using namespace Biopool;

void sShowHelp() {
    cout << "Scat Edit\n"
            << "This program allows to edit scatter plots from loop modeling.\n"
            << " Options: \n"
            << "\t-i <filename> \t\t Input file\n"
            << "\t[-s <length>] \t\t Select only loops of this length (default "
            << "= all lengths)\n"
            << "\t[-c <value>] \t\t Maximum combined score cutoff value for"
            << " ''predictable'' loops, i.e.''good quality'' (default = 10000)\n"
            << "\t[-o <filename>] \t Output file (default = scat.plot)\n"
            << "\t[-d <weight>] \t\t Default weight (instead of 125)\n"
            << "\t-x <column> \t\t Select column for x axis in output\n"
            << "\t-y <column> \t\t Select column for y axis in output\n"
            << "\t\t NB: use code 99 to select optimized score I in -x & -y \n"
            << "\t\t     use code 999 to select optimized score II in -x & -y \n"
            << "\t\t     use code 88 to select optimized score III in -x & -y \n"
            << "\t [--noOpt1] \t do not optimize with scheme I\n"
            << "\t [--noOpt2] \t do not optimize with scheme II\n"
            << "\t [--noOpt3] \t do not optimize with scheme III\n"
            << "\t [--topfile <filename>] \t write top X stats to filename\n";
}

void sLine() {
    cout << "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n";
}

double sCorrelation(vector<double>& zs1, vector<double>& zs2) {
    double tmp = 0.0;
    for (unsigned int j = 0; j < zs1.size(); j++)
        tmp += zs1[j] * zs2[j];

    tmp /= zs1.size();
    return tmp;
}

int main(int nArgs, char* argv[]) {
    cout.setf(ios::fixed, ios::floatfield);

    if (getArg("h", nArgs, argv)) {
        sShowHelp();
        return 1;
    };

    string inputFile, outputFile, topFile;

    getArg("i", inputFile, nArgs, argv, "!");
    getArg("o", outputFile, nArgs, argv, "scat.plot");
    getArg("-topfile", topFile, nArgs, argv, "!");

    unsigned int select, colX, colY;
    getArg("s", select, nArgs, argv, 0);
    getArg("x", colX, nArgs, argv, 0);
    getArg("y", colY, nArgs, argv, 0);

    double maxScore;
    getArg("c", maxScore, nArgs, argv, 10000.0);

    double def;
    getArg("d", def, nArgs, argv, 125);

    if (inputFile == "!") {
        cout << "Missing file specification. Aborting. (-h for help)" << endl;
        return -1;
    }

    bool noOpt1 = getArg("-noOpt1", nArgs, argv);
    bool noOpt2 = getArg("-noOpt2", nArgs, argv);
    bool noOpt3 = getArg("-noOpt3", nArgs, argv);

    ifstream inFile(inputFile.c_str());
    if (!inFile)
        ERROR("File not found.", exception);

    RankAnalyzer ra;

    ra.load(inFile, select);
    ra.calcStatistics();

    // correlations:
    cout << "Correlations:\n";

    ra.printHeader();
    ra.printCorrelation("\nRMSD  ", 0);
    ra.printCorrelation("\nScore ", 1);
    ra.printCorrelation("\nendRms", 2);
    ra.printCorrelation("\nendRms", 3);
    ra.printCorrelation("\nsol_en", 4);
    cout << "\n";
    sLine();

    if (def > 0) {
        vector<double> data;

        for (unsigned int j = 0; j < ra.col[0].size(); j++)
            data.push_back(def * ra.col[2][j] + ra.col[3][j]);

        // average:
        double avgD = 0.0;
        for (unsigned int j = 0; j < data.size(); j++)
            avgD += data[j];
        avgD = avgD / data.size();

        //standard deviation:
        double sdD = 0.0;
        for (unsigned int j = 0; j < data.size(); j++)
            sdD += (data[j] - avgD) * (data[j] - avgD);
        sdD = sqrt(sdD / data.size());

        // do Z-score transform:
        vector<double> zSc;
        for (unsigned int j = 0; j < data.size(); j++)
            zSc.push_back((data[j] - avgD) / sdD);

        double cor = sCorrelation(ra.zScore[0], zSc);

        cout << "Correlation for default = " << def << "\t cor = " << cor
                << "\n";

        sLine();
    }

    cout << "Confidence threshold for Top X is " << maxScore << "\n";
    ra.printTopXResults(1, "!", maxScore);
    ra.printTopXResults(3, "!", maxScore);
    ra.printTopXResults(5, "!", maxScore);
    ra.printTopXResults(10, "!", maxScore);
    ra.printTopXResults(20, "!", maxScore);
    ra.printTopXResults(1000, "!", maxScore);
    sLine();

    if (topFile != "!")
        for (unsigned int i = 0; i <= 1000; i += 25) {
            ra.printTopXResults(1, topFile, i);
            ra.printTopXResults(3, topFile, i);
            ra.printTopXResults(5, topFile, i);
            ra.printTopXResults(10, topFile, i);
            ra.printTopXResults(20, topFile, i);
            ra.printTopXResults(1000, topFile, i);
        }

    // optimization:

    double best = sCorrelation(ra.zScore[0], ra.zScore[1]);
    double bestA = 125.0;

    vector<double> bestData;
    vector<double> bestRaw;

    if (!noOpt1) {

        cout << "Optimizing the scoring function:\n";
        cout << "Score = a * endRms + energy\n";

        for (double a = 1.0; a < 1500; a += 0.25) { // loop over parameter:

            vector<double> data;

            for (unsigned int j = 0; j < ra.col[0].size(); j++)
                data.push_back(a * ra.col[2][j] + ra.col[3][j]);

            // average:
            double avgD = 0.0;
            for (unsigned int j = 0; j < data.size(); j++)
                avgD += data[j];
            avgD = avgD / data.size();

            //standard deviation:
            double sdD = 0.0;
            for (unsigned int j = 0; j < data.size(); j++)
                sdD += (data[j] - avgD) * (data[j] - avgD);
            sdD = sqrt(sdD / data.size());

            // do Z-score transform:
            vector<double> zSc;
            for (unsigned int j = 0; j < data.size(); j++)
                zSc.push_back((data[j] - avgD) / sdD);

            double cor = sCorrelation(ra.zScore[0], zSc);

            // evaluate quality:
            if (cor > best) {
                best = cor;
                bestA = a;

                bestData.clear();
                for (unsigned int i = 0; i < zSc.size(); i++)
                    bestData.push_back(zSc[i]);

                bestRaw.clear();
                for (unsigned int i = 0; i < data.size(); i++)
                    bestRaw.push_back(data[i]);
            }

        }
        cout << "\n";
        sLine();
        cout << "---> Best a parameter is " << bestA << "\n";

        ra.printHeader();
        ra.printCorrelation("Best ", bestData);
        cout << "\n";
        sLine();
    }

    double best2 = sCorrelation(ra.zScore[0], ra.zScore[1]);
    double bestA2 = -1.0;

    vector<double> bestData2;
    vector<double> bestRaw2;

    if (!noOpt2) {
        cout << "Optimizing the scoring function II:\n";
        cout << "Score = a * (endRms+1)^2 + energy\n";

        for (double a = 0.125; a < 500; a += 0.125) {// loop over parameter:

            vector<double> data;

            for (unsigned int j = 0; j < ra.col[0].size(); j++)
                data.push_back(a * (1 + ra.col[2][j])*(1 + ra.col[2][j])
                    + ra.col[3][j]);

            // average:
            double avgD = 0.0;
            for (unsigned int j = 0; j < data.size(); j++)
                avgD += data[j];
            avgD = avgD / data.size();

            //standard deviation:
            double sdD = 0.0;
            for (unsigned int j = 0; j < data.size(); j++)
                sdD += (data[j] - avgD) * (data[j] - avgD);
            sdD = sqrt(sdD / data.size());

            // do Z-score transform:
            vector<double> zSc;
            for (unsigned int j = 0; j < data.size(); j++)
                zSc.push_back((data[j] - avgD) / sdD);

            double cor = sCorrelation(ra.zScore[0], zSc);

            // evaluate quality:
            if (cor > best2) {
                //  	      cout << "a= " << a << " corr= " << cor << "\t\t";
                best2 = cor;
                bestA2 = a;

                bestData2.clear();
                for (unsigned int i = 0; i < zSc.size(); i++)
                    bestData2.push_back(zSc[i]);

                bestRaw2.clear();
                for (unsigned int i = 0; i < data.size(); i++)
                    bestRaw2.push_back(data[i]);
            }
        }

        cout << "\n";
        sLine();

        cout << "---> Best a parameter II is " << bestA2 << "\n";

        ra.printHeader();
        ra.printCorrelation("Best ", bestData2);
        cout << "\n";
        sLine();
    }

    double best3 = sCorrelation(ra.zScore[0], ra.zScore[1]);
    double bestA3 = -1.0;
    double bestB3 = -1.0;

    vector<double> bestData3;
    vector<double> bestRaw3;

    if (!noOpt3) {
        cout << "Optimizing the scoring function III:\n";
        cout << "Score = a * endRms + energy + b * sol_en\n";

        for (double a = 1.0; a < bestA + 200; a += 1.0) // loop over parameter:
            for (double b = 0.0; b < 25; b += 0.25) {
                vector<double> data;

                for (unsigned int j = 0; j < ra.col[0].size(); j++)
                    data.push_back(a * ra.col[2][j] + ra.col[3][j]
                        + b * ra.col[4][j]);

                // average:
                double avgD = 0.0;
                for (unsigned int j = 0; j < data.size(); j++)
                    avgD += data[j];
                avgD = avgD / data.size();

                //standard deviation:
                double sdD = 0.0;
                for (unsigned int j = 0; j < data.size(); j++)
                    sdD += (data[j] - avgD) * (data[j] - avgD);
                sdD = sqrt(sdD / data.size());

                // do Z-score transform:
                vector<double> zSc;
                for (unsigned int j = 0; j < data.size(); j++)
                    zSc.push_back((data[j] - avgD) / sdD);

                double cor = sCorrelation(ra.zScore[0], zSc);

                // evaluate quality:
                if (cor > best3) {
                    best3 = cor;
                    bestA3 = a;
                    bestB3 = b;

                    bestData3.clear();
                    for (unsigned int i = 0; i < zSc.size(); i++)
                        bestData3.push_back(zSc[i]);

                    bestRaw3.clear();
                    for (unsigned int i = 0; i < data.size(); i++)
                        bestRaw3.push_back(data[i]);
                }
            }

        cout << "\n";
        sLine();
        cout << "---> Best a parameter III is " << bestA3 << " best b is "
                << bestB3 << "\n";

        ra.printHeader();
        ra.printCorrelation("Best ", bestData3);
        cout << "\n";
    }

    double best4 = sCorrelation(ra.zScore[0], ra.zScore[1]);
    double bestA4 = -1.0;
    double bestB4 = -1.0;

    vector<double> bestData4;
    vector<double> bestRaw4;

    if (!noOpt3) {
        cout << "Optimizing the scoring function IV:\n";
        cout << "Score = a * endRms + energy + b * prop\n";

        for (double a = 1.0; a < bestA + 200; a += 1.0) // loop over parameter:
            for (double b = -250.0; b < 250; b += 2.5) {
                vector<double> data;

                for (unsigned int j = 0; j < ra.col[0].size(); j++)
                    data.push_back(a * ra.col[2][j] + ra.col[3][j]
                        + b * ra.col[10][j]);

                // average:
                double avgD = 0.0;
                for (unsigned int j = 0; j < data.size(); j++)
                    avgD += data[j];
                avgD = avgD / data.size();

                //standard deviation:
                double sdD = 0.0;
                for (unsigned int j = 0; j < data.size(); j++)
                    sdD += (data[j] - avgD) * (data[j] - avgD);
                sdD = sqrt(sdD / data.size());

                // do Z-score transform:
                vector<double> zSc;
                for (unsigned int j = 0; j < data.size(); j++)
                    zSc.push_back((data[j] - avgD) / sdD);

                double cor = sCorrelation(ra.zScore[0], zSc);

                // evaluate quality:
                if (cor > best4) {
                    best4 = cor;
                    bestA4 = a;
                    bestB4 = b;

                    bestData4.clear();
                    for (unsigned int i = 0; i < zSc.size(); i++)
                        bestData4.push_back(zSc[i]);

                    bestRaw4.clear();
                    for (unsigned int i = 0; i < data.size(); i++)
                        bestRaw4.push_back(data[i]);
                }
            }

        cout << "\n";
        sLine();
        cout << "---> Best a parameter IV is " << bestA4 << " best b is "
                << bestB4 << "\n";

        ra.printHeader();
        ra.printCorrelation("Best ", bestData4);
        cout << "\n";
    }

    // write plot file (optional):

    if (colX * colY != 0) {
        ofstream out(outputFile.c_str());
        if (!out)
            ERROR("File not found.", exception);

        for (unsigned int i = 0; i < ra.col[0].size(); i++) {
            if (colX == 99)
                out << setw(5) << bestRaw[i] << "\t";
            else if (colX == 999)
                out << setw(5) << bestRaw2[i] << "\t";
            else if (colX == 88)
                out << setw(5) << bestRaw3[i] << "\t";
            else
                out << setw(5) << ra.col[colX - 1][i] << "\t";

            if (colY == 99)
                out << setw(5) << bestRaw[i] << "\n";
            else if (colY == 999)
                out << setw(5) << bestRaw2[i] << "\n";
            else if (colY == 88)
                out << setw(5) << bestRaw3[i] << "\n";
            else
                out << setw(5) << ra.col[colY - 1][i] << "\n";
        }
    }

    return 0;
}
