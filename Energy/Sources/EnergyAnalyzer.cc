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
 *@Class               EnergyAnalyzer
 *@Project        Victor
 *@Description   This class implements functions to analyze the ranking output of Energy.
 */
// Includes:
#include <EnergyAnalyzer.h>
#include <IoTools.h>
#include <StatTools.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;

static void sLine() {
    cout << "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n";
}

/**
 *@Description Calculates Correlation between two vectors
 *@param reference for the vectors of double(vector<double>&,vector<double>&)
 *@return corresponding value(static double)
 */
static double sCorrelation(vector<double>& zs1, vector<double>& zs2) {
    double tmp = 0.0;
    for (unsigned int j = 0; j < zs1.size(); j++)
        tmp += zs1[j] * zs2[j];

    tmp /= zs1.size();
    return tmp;
}

// CONSTRUCTORS/DESTRUCTOR:

/**
 *@Description  basic constructor
 */
EnergyAnalyzer::EnergyAnalyzer() {
    PRINT_NAME;
    for (unsigned int i = 0; i < MAX_COL; i++) {
        vector<double> tmp;
        col.push_back(tmp);
        vector<double> tmp2;
        zScore.push_back(tmp2);
        avg[i] = 0.0;
        sd[i] = 0.0;
    }

    for (unsigned int i = 0; i < MAX_DATA; i++) {
        vector<double> tmp;
        result.push_back(tmp);
        vector<double> tmp2;
        resScore.push_back(tmp);
    }

}

/**
 *@Description  constructor that uses another object copy
 *@param  reference to the original object to copy EnergyAnalyzer&
 */
EnergyAnalyzer::EnergyAnalyzer(const EnergyAnalyzer& orig) {
    PRINT_NAME;
    this->copy(orig);
}

/**
 *@Description  basic destructor 
 */
EnergyAnalyzer::~EnergyAnalyzer() {
    PRINT_NAME;
    for (unsigned int i = 0; i < MAX_COL; i++) {
        col[i].clear();
        zScore[i].clear();
    }
    for (unsigned int i = 0; i < MAX_DATA; i++) {
        result[i].clear();
        resScore[i].clear();
    }

}


// PREDICATES:

/**
 *@Description Prints the headers  
 *@param  none
 */
void EnergyAnalyzer::printHeader() {
    cout << "\t 1     \t 2     \t  3    \t  4    \t  5    \t 6    \n";
    cout << "\tGDT_TS \tScore  \tRapdf  \t  Solv \t HydB  \t Tors \n\n";
}

/**
 *@Description  Prints all the correlations
 *@param  none
 */
void EnergyAnalyzer::printAllCorrelations() {
    cout << "Correlations:\n";
    printHeader();
    printCorrelation((char*) "\nGDT_TS", 0);
    printCorrelation((char*) "\nScore ", 1);
    printCorrelation((char*) "\nRapdf ", 2);
    printCorrelation((char*) "\nSolv  ", 3);
    printCorrelation((char*) "\nHydB  ", 4);
    printCorrelation((char*) "\nTors  ", 5);
    cout << "\n";
    sLine();
}

/**
 *@Description  Prints a message and a specific correlation for a specific index of the zScore
 *@param  message to print (char * ), index of the Zscore (unsigned int) 
 */
void
EnergyAnalyzer::printCorrelation(char* intro, unsigned int index) {
    cout << intro;
    for (unsigned int i = 0; i < MAX_COL; i++)
        cout << " \t" << setw(5) << setprecision(3)
        << sCorrelation(zScore[index], zScore[i]);
    cout << "\n";
}

/**
 *@Description  Prints a message and the correlations between the values in the given vector and the Zscore indexes
 *@param  message to print (char * ), values for calculate the correlation(vector of doubles)  
 */
void EnergyAnalyzer::printCorrelation(char* intro, vector<double> data) {
    cout << intro;
    for (unsigned int i = 0; i < MAX_COL; i++)
        cout << " \t" << setw(5) << setprecision(3)
        << sCorrelation(data, zScore[i]);
    cout << "\n";
}

/**
 *@Description  Prints the results
 */
void EnergyAnalyzer::printTopXResults(unsigned int top, string topFile, double maxScore) {
    unsigned int maxIndex = 0;
    for (unsigned int i = 0; i < MAX_DATA; i++)
        if (result[i].size() == 0)
            break;
        else
            maxIndex++;

    vector<double> bestX;
    vector<double> bestXSc;
    unsigned int notCovered = 0;

    // find Top X answer for every single data entry:
    for (unsigned int i = 0; i < maxIndex; i++) {
        double tmp = result[i][0];
        double tmpSc = resScore[i][0];

        if (tmpSc > maxScore) {
            notCovered++;
            continue;
        }

        unsigned int maxTmp = (top < result[i].size() ? top
                : result[i].size());

        for (unsigned int j = 1; j < maxTmp; j++) {
            if (result[i][j] > tmp)
                tmp = result[i][j];

            if (j >= resScore[i].size())
                ERROR("Index out of scope.", exception);
            if (resScore[i][j] < tmpSc)
                tmpSc = resScore[i][j];
        }
        bestX.push_back(tmp);
        bestXSc.push_back(tmpSc);
    }

    // now calculate the statistics:
    double avgT = average(bestX);
    double avgSc = average(bestXSc);
    double sdT = standardDeviation(bestX, avgT);
    double sdSc = standardDeviation(bestXSc, avgSc);

    // if topFile is ''valid'', do file output, otherwise screen output:
    if (topFile != "!") {
        // do file I/O:
        ofstream outFile(topFile.c_str(), ios::app);
        if (!outFile)
            ERROR("Could not write file.", exception);

        outFile << setw(3) << maxScore << "\t" << setw(4)
                << top << "\t" << setw(5) << setprecision(3)
                << avgT << "\t" << setw(5) << setprecision(3)
                << sdT << "\t" << setw(5) << setprecision(3)
                << avgSc << "\t" << setw(5) << setprecision(3)
                << sdSc << "\t" << setw(5) << setprecision(3)
                << 100 * static_cast<double> (maxIndex - notCovered) / maxIndex
                << "\n";
    } else
        cout << "Top " << setw(4) << top << "\t results:   avg = "
            << setw(5) << setprecision(3)
        << avgT << "    sd = " << setw(5) << setprecision(2)
        << sdT << "\t score:   avg = "
            << setw(5) << setprecision(3)
        << avgSc << "   sd = " << setw(5) << setprecision(2)
        << sdSc << "\t  cov% = " << setw(5) << setprecision(3)
        << 100 * static_cast<double> (maxIndex - notCovered) / maxIndex
            << "\n";
}


// MODIFIERS:

/**
 *@Description  loads the information from the file
 *@param reference of the file(istream&)  ,not used variable(unsigned) , multiplies the value per 100(bool)
 *  *@return  changes are made internally(void)
 */
void EnergyAnalyzer::load(istream& inFile, unsigned int select, bool mult) {
    // result has to be cleared here

    result.clear();
    vector<double> tmpVec;
    result.push_back(tmpVec);
    unsigned int line = 0;
    int sum = 0;

    // read data from file:
    while (inFile) {
        eatComment(inFile);
        if (checkForKeyword(inFile, "DECOY")) {
            skipToNewLine(inFile);
            if (sum > 0)
                line++;
            continue;
        }

        sum++;

        if (!inFile)
            break;

        for (unsigned int i = 0; i < MAX_COL; i++) {
            double tmp = 0.0;
            inFile >> tmp;


            if (!inFile)
                break;

            if (i == 0) {
                if (mult)
                    tmp *= 100.0;
                result[line].push_back(tmp);
            } else if (i == 1)
                resScore[line].push_back(tmp);

            col[i].push_back(tmp);
        }
        skipToNewLine(inFile);

    }
}

/**
 *@Description calculates the statistics 
 *@param  none
 *@return   prints all the results (void)
 */
void EnergyAnalyzer::calcStatistics() {
    sLine();
    cout << "\t GDT_TS \t combScore \t Rapdf \t\t Solv  \t\t HydB  \t\t Tors \t"
            << " # entries = " << col[0].size() << "\n";

    // average:
    for (unsigned int i = 0; i < MAX_COL; i++)
        avg[i] = average(col[i]);

    cout << "Avg: ";
    for (unsigned int i = 0; i < MAX_COL; i++)
        cout << " \t " << setw(6) << avg[i];
    cout << "\n";

    //standard deviation:
    for (unsigned int i = 0; i < MAX_COL; i++)
        sd[i] = standardDeviation(col[i], avg[i]);

    cout << "SD: ";
    for (unsigned int i = 0; i < MAX_COL; i++)
        cout << " \t " << setw(6) << sd[i];
    cout << "\n";
    sLine();

    //Z-scores:
    for (unsigned int i = 0; i < MAX_COL; i++)
        for (unsigned int j = 0; j < col[i].size(); j++)
            zScore[i].push_back((col[i][j] - avg[i]) / sd[i]);
}

/**
 *@Description  Made a Linear optimization of the data, considering four sets of values
 *@param   four sets of values(double,double,double,double,double,double,double,double,double,double,double,double)
 *@return  vector containing the optimized data(vector <double>)
 */

vector<double> EnergyAnalyzer::linearOptimization(double a1, double a2, double a3, double b1,
        double b2, double b3, double c1, double c2, double c3, double d1, double d2,
        double d3) {
    double best = sCorrelation(zScore[0], zScore[1]);
    double bestA = 1.0;
    double bestB = 0.0;
    double bestC = 0.0;
    double bestD = 0.0;

    vector<double> bestData;
    vector<double> bestRaw;

    cout << "Optimizing the scoring function:\n\n";
    cout << "Score = a * rapdf + b * solv + c * hydB + d * tors\n";

    for (double d = d1; d <= d2; d += d3) // loop over parameter:
        for (double c = c1; c <= c2; c += c3) // loop over parameter:
            for (double b = b1; b <= b2; b += b3)
                for (double a = a1; a <= a2; a += a3) {
                    vector<double> data;

                    for (unsigned int j = 0; j < col[0].size(); j++)
                        data.push_back(a * col[2][j] + b * col[3][j]
                            + c * col[4][j] + d * col[5][j]);

                    vector<double> zSc = Zscore(data);

                    double cor = sCorrelation(zScore[0], zSc);

                    // evaluate quality:
                    if (cor < best) {
                        //  cout << "a= " << a << "   b= " << b << " c= " 
                        //       << c << " " << "  corr= " << cor << "\t\t";
                        best = cor;
                        bestA = a;
                        bestB = b;
                        bestC = c;
                        bestD = d;

                        bestData.clear();
                        for (unsigned int i = 0; i < zSc.size(); i++)
                            bestData.push_back(zSc[i]);

                        bestRaw.clear();
                        for (unsigned int i = 0; i < data.size(); i++)
                            bestRaw.push_back(data[i]);
                    }

                }

    cout << "\n--->\ta = " << bestA
            << "\t b = " << bestB << "\t c = " << bestC
            << "\t d = " << bestD << "\n\n";

    printHeader();
    printCorrelation((char*) "Best ", bestData);
    cout << "\n";
    sLine();

    return bestRaw;
}

/**
 *@Description  returns a random value
 *@param  none
 *@return  the corresponding value( double)
 */

inline double sGetRand() {
    double tmp = 0.0;
    for (unsigned int i = 0; i < 12; i++)
        tmp += static_cast<double> (rand()) / RAND_MAX;

    return tmp - 6;
}

/**
 *@Description  Returns three random values
 *@param  references for each random value(double&,double&,double&)
 *@return  changes are made internally (void)
 */
void sInit(double& a, double& b, double& c) {
    a = sGetRand() * 100 + 1000;
    b = sGetRand() * 10000;
    c = sGetRand() * -0.25 - 0.125;
}

/**
 *@Description  Returns three random values changed randomly
 *@param  references for each random value(double&,double&,double&)
 *@return  changes are made internally (void)
 */
void sMutate(double& a, double& b, double& c) {
    a *= 1 + (sGetRand() * 1);
    b *= 1 + (sGetRand() * 0.05);
    c *= 1 + (sGetRand() * 0.1);
}

double sCalcValue(double base, double a, double b, double c) {
    return a * (tanh(c * (base - b)) + 1) / 2;
}

/**
 *@Description  Made a non Linear optimization of the data, considering two max values
 *@param   two maximum values(unsigned int,unsigned int)
 *@return  vector containing the optimized data(vector <double>)
 */

vector<double>EnergyAnalyzer::nonLinearOptimization(unsigned int max1, unsigned int max2) {
    double best = 0.0;

    vector<double> bestA, bestB, bestC, a, b, c;
    for (unsigned int i = 0; i < 4; i++) {
        bestA.push_back(0.0);
        bestB.push_back(0.0);
        bestC.push_back(0.0);
        a.push_back(0.0);
        b.push_back(0.0);
        c.push_back(0.0);
    }

    vector<double> bestData;
    vector<double> bestRaw;

    cout << "Optimizing the non-linear scoring function:\n\n";
    cout << "Score = tanh( rapdf ) + tanh( solv ) + tanh( tors ) "
            << "+ tanh( hydB )\n";

    for (unsigned int num1 = 0; num1 <= max1; num1++) {// loop over parameter:
        for (unsigned int i = 0; i < 4; i++)
            sInit(a[i], b[i], c[i]);
        b[0] = -3000 + sGetRand() * 2700;
        b[1] = -1.5 + sGetRand() * 8.7;
        b[2] = -150 + sGetRand() * 85;
        b[3] = 20 + sGetRand() * 30;

        for (unsigned int num2 = 0; num2 <= max2; num2++) {
            vector<double> data;

            for (unsigned int i = 0; i < 4; i++)
                sMutate(a[i], b[i], c[i]);


            for (unsigned int j = 0; j < col[0].size(); j++)
                data.push_back(sCalcValue(-col[2][j], a[0], b[0], c[0])
                    + sCalcValue(-col[3][j], a[1], b[1], c[1])
                    + sCalcValue(-col[4][j], a[2], b[2], c[2])
                    + sCalcValue(-col[5][j], a[3], b[3], c[3]));

            vector<double> zSc = Zscore(data);

            double cor = sCorrelation(zScore[0], zSc);

            // evaluate quality:
            if (cor < best) {
                cout << "  corr= " << cor << "\t\t";
                best = cor;

                bestData.clear();
                for (unsigned int i = 0; i < zSc.size(); i++)
                    bestData.push_back(zSc[i]);

                bestRaw.clear();
                for (unsigned int i = 0; i < data.size(); i++)
                    bestRaw.push_back(data[i]);
            }
        }
    }


    cout << "\n";
    printHeader();
    printCorrelation((char*) "Best ", bestData);
    cout << "\n";
    sLine();

    return bestRaw;
}

/**
 *@Description  Copies the original Energy analyzer into another
 *@param reference to the original energy analyzer (EnergyAnalyzer&)
 *@return  the changes are made internally(void)
 */
void EnergyAnalyzer::copy(const EnergyAnalyzer& orig) {
    for (unsigned int i = 0; i < MAX_COL; i++) {
        col[i].clear();
        for (unsigned int j = 0; j < orig.col[i].size(); j++)
            col[i].push_back(orig.col[i][j]);

        zScore[i].clear();
        for (unsigned int j = 0; j < orig.col[i].size(); j++)
            zScore[i].push_back(orig.zScore[i][j]);

        avg[i] = orig.avg[i];
        sd[i] = orig.sd[i];
    }

    for (unsigned int i = 0; i < MAX_DATA; i++) {
        result[i].clear();
        resScore[i].clear();
        for (unsigned int j = 0; j < orig.result[i].size(); j++) {
            result[i].push_back(orig.result[i][j]);
            resScore[i].push_back(orig.resScore[i][j]);
        }
    }
}


// HELPERS:
