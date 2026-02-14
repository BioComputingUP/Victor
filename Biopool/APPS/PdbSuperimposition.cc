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
 * Author: Simone Piano
 *
 */

#include <GetArg.h>
#include <PdbLoader.h>
#include <Eigen/Core>
#include <StructurePdbSuperimposition.h>

using namespace Victor;using namespace Victor::Biopool;
using namespace Eigen;

void sShowHelp() {
    cout << "Pdb Superimposition -- calculate the superimposition scores (RMSD/MaxSub/GDT_TS/TM) of two pdb files\n"
        << "Options: \n"
        << "\t-i <filename> <filename>, \t Input PDB filename\n"
        << "\t[-r], \t Calculate only RMSD score\n"
        << "\t[-m], \t Calculate only MaxSub score\n"
        << "\t[-g], \t Calculate only GDT_TS score\n"
        << "\t[-t], \t Calculate only TM score\n"
        << "\t[-L <seed>], is a integer number for the initial seeds lenght of MaxSub\n"
	<< "\t\t used by MaxSub, GDT_TS, TM-Score (default is 4)\n"
	<< "\t[-d <threshold>], is a double number for the threshold of MaxSub (default is 3.5A)\n"
	<< "Usage: \n"
	<< "\tPdbSuperimposition -i <file_name> <file_name> [-<char_score>] [-L <integer_seed>] [-d <double_threshold>]\n";
}

int main(int nArgs, char* argv[]) {
    if (getArg("h", nArgs, argv)) {
        sShowHelp();
        return 1;
    }
    
    std::vector<string> inputFile;
    bool rmsd, all, maxsub, gdtts, tmscore;
    string seed, threshold;
    long int L; //seed
    double d; //threshold for MaxSub

    getArg("i", inputFile, nArgs, argv, "!");
    getArg("L", seed, nArgs, argv, "!");
    getArg("d", threshold, nArgs, argv, "!");
    rmsd = getArg("r", nArgs, argv);
    maxsub = getArg("m", nArgs, argv);
    gdtts = getArg("g", nArgs, argv);
    tmscore = getArg("t", nArgs, argv);

    // Set correct seed and threshold
    (seed == "!") ? L=4 : L = strtol(seed.c_str(), NULL, 10);
    (threshold == "!") ? d=3.5 : d=strtod(threshold.c_str(), NULL);

    if (nArgs == 1) {
        sShowHelp();
        return 1;
    }
    if (inputFile.size() != 2) {
        cout << "Missing file specification. (-h for help)" << endl;
        return 1;
    }
    if (L < 4) {
	cout << "L must be a valid integer number >= 4. (-h for help)" << endl;
        return 1;
    }
    if (d <= 0.0) {
        cout << "d must be a valid double number > 0.0. (-h for help)" << endl;
        return 1;
    }
    
    (!rmsd && !maxsub && !gdtts && !tmscore) ? all = true : all = false;

    ifstream inFile1(inputFile[0].c_str());
    ifstream inFile2(inputFile[1].c_str());
    
    PdbLoader pl1(inFile1);
    PdbLoader pl2(inFile2);
    
    // Set PdbLoader variables for first pdb file
    pl1.setNoHAtoms();
    
    // Set PdbLoader variables for second pdb file
    pl2.setNoHAtoms();
    	
    // Load the protein object
    Protein prot1, prot2;
    prot1.load(pl1);
    prot2.load(pl2);
    
    // chainsIdProt1 and chainsIdProt2 contains all IDs chain available
    vector<char> chainsIdProt1 = prot1.getAllChains();
    vector<char> chainsIdProt2 = prot2.getAllChains();
    
    // Get Spacer
    Spacer *sp1, *sp2;
    sp1 = prot1.getSpacer(chainsIdProt1[0]);
    sp2 = prot2.getSpacer(chainsIdProt2[0]);
    
    // Get Spacer coordinates in Eigen matrix format
    Eigen::MatrixX3d matrix1, matrix2;
    matrix1 = StructurePdbSuperimposition::fromSpacerToEigenMatrix(sp1);
    matrix2 = StructurePdbSuperimposition::fromSpacerToEigenMatrix(sp2);
    
    if(L > sp1->sizeAmino()) {
        cout << "L is larger than the number of points"  << endl;
        return -1;
    }
    
    // The number of AA of the two Spacer must be equal
    if (matrix1.rows() == matrix2.rows() && matrix1.cols() == matrix2.cols() && matrix1.rows() > 0) {        
        if (all || rmsd) {
            double rmsdScore = StructurePdbSuperimposition::RMSD(matrix1, matrix2, *sp1, *sp2);
            cout << "\nRMSD: " << rmsdScore << endl;
        }

	if (all || maxsub) {
	    double maxSubScore = StructurePdbSuperimposition::MaxSubScore(matrix1, matrix2, *sp1, *sp2, L, d, true);
	    cout << "MaxSub: " << maxSubScore << endl;
        }

	if (all || gdtts) { 
            double result[4];
            double gdtTsScore = StructurePdbSuperimposition::GDTTS(matrix1, matrix2, *sp1, *sp2, L, result);
            cout << "GDT_TS: " << gdtTsScore << endl;
	    cout << "\tGDT_TS (d<1): " << result[0] << endl;
	    cout << "\tGDT_TS (d<2): " << result[1] << endl;
	    cout << "\tGDT_TS (d<4): " << result[2] << endl;
	    cout << "\tGDT_TS (d<8): " << result[3] << endl;
	}

	if (all || tmscore) {
	    double tmScore = StructurePdbSuperimposition::TMScore(matrix1, matrix2, *sp1, *sp2, L);
	    cout << "TM-Score: " << tmScore << endl;
	}
    }
    else {
        cout << "The matrix has different size or null dimension" << endl;
	return -1;
    }
}
