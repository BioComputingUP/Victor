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

#include <iostream>
#include <iomanip>
#include <string>
#include <Protein.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <CifLoader.h>
#include <CifSaver.h>
#include <GetArg.h>

using namespace Victor;
using namespace Victor::Biopool;
using namespace std;

void sShowHelp() {
    cout << "CIF Reader / Writer \n"
	<< "Allows to read and write a CIF or a PDB file without modifying them.\n"
	<< " Options: \n"
	<< "\t-i <filename> \t\t Input CIF/PDB file\n"
	<< "\t-o <filename> \t\t Output CIF/PDB file\n"
	<< "\n";
}

int main(int argc, char** argv) {

    // --------------------------------------------------
    // 0. treat options
    // --------------------------------------------------

    if (getArg("h", argc, argv)) {
	sShowHelp();
	return 1;
    };

    string inputFile, outputFile;
    Protein prot;
    
    getArg("i", inputFile, argc, argv, "!");
    getArg("o", outputFile, argc, argv, "!");

    if ((inputFile == "!") || (outputFile == "!")) {
	cout << "Missing file specification. Aborting. (-h for help)" << endl;
	return -1;
    }

    // --------------------------------------------------
    // 1. read structure
    // --------------------------------------------------

    ifstream inFile(inputFile.c_str());
    if (!inFile) {
	ERROR("File does not exist.\n", exception);
    }
    
    if (inputFile.find("pdb") != string::npos) {
	cout << "Loading PDB file..." << endl;

	PdbLoader pl(inFile);
	prot.load(pl);

    } else if (inputFile.find("cif") != string::npos) {
	cout << "Loading CIF file..." << endl;

	CifLoader cl(inFile);
	prot.load(cl);

    } else {
	cout << "Uknown input file format. Aborting. (-h for help)" << endl;
	return -2;
    }
    
    inFile.close();

    // --------------------------------------------------
    // 2. write structure
    // --------------------------------------------------

    ofstream outFile(outputFile.c_str());
    if (!outFile) {
	ERROR("File not found.", exception);
    }
    
    if (outputFile.find("pdb") != string::npos) {
	cout << "Saving PDB file..." << endl;
	
	PdbSaver ps(outFile);
	prot.save(ps);
	
    } else if (outputFile.find("cif") != string::npos) {
	cout << "Saving CIF file..." << endl;
	
	CifSaver cs(outFile);
	prot.save(cs);
	
    } else {
	cout << "Uknown output file format. Aborting. (-h for help)" << endl;
	return -3;
    }

    outFile.close();
    
    return 0;
}

