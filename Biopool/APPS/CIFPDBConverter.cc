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
 */
#include <Protein.h>
#include <Loader.h>
#include <Saver.h>
#include <CIFLoader.h>
#include <CIFSaver.h>
#include <PDBLoader.h>
#include <PDBSaver.h>

#include <IoTools.h>
#include <GetArg.h>
#include <StringUtils.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <stdlib.h> 

using namespace std;
using namespace Victor;
using namespace Victor::Biopool;

void sShowHelp() {
    cout << "CIF PDB Converter $Revision: 1.0 $ -- converts a CIF file into a PDB file and vice versa." << endl
            << "Options:" << endl
            << "\t-i <filename> \t Input file" << endl
            << "\t-o <filename> \t Output file";
}

int main(int argc, char* argv[]) {

    if (getArg("h", argc, argv)) {
        sShowHelp();
        return 1;
    }

    string inputFile, outputFile;

    getArg("i", inputFile, argc, argv, "!");
    getArg("o", outputFile, argc, argv, "!");

    // Check input file
    if (inputFile == "!") {
        cout << "Missing input file specification. Aborting. (-h for help)" << endl;
        return -1;
    }
    ifstream inFile(inputFile.c_str());
    if (!inFile)
        ERROR("Input file not found.", exception);

    // Check output file
    if (outputFile == "!") {
        cout << "Missing output file specification. Aborting. (-h for help)" << endl;
        return -1;
    }

    inputFile = Trim(inputFile);
    outputFile = Trim(outputFile);
    
    string inExt = inputFile.substr((inputFile.size()-3),3);
    string outExt = outputFile.substr((outputFile.size()-3),3);
    
    transform(inExt.begin(), inExt.end(), inExt.begin(), ::toupper);
    transform(outExt.begin(), outExt.end(), outExt.begin(), ::toupper);
    
    if (inExt != "CIF" && inExt != "PDB")
    {
        cout << "Invalid input file format." << endl;
        return -1;
    }
    
    if (outExt != "CIF" && outExt != "PDB")
    {
        cout << "Invalid output file format." << endl;
        return -1;
    }
    
    if (inExt == outExt)
    {
        cout << "No conversion needed." << endl;
        return 0;
    }
    
    Loader* l;
    Saver* s;
    ofstream out(outputFile.c_str());
    if (!out)
    {
        ERROR("Output file not created.", exception);
    }
        
    
    if (inExt == "CIF")
    {
        l = new CIFLoader(inFile);
    }
    else
    {
        l = new PdbLoader(inFile);
    }

   // Load the protein object
    Protein prot;
    prot.load(*l);
    
    if (outExt == "CIF")
    {
        s = new CIFSaver(out);
    }
    else
    {
        s = new PdbSaver(out);
    }

    prot.save(*s);
    
    return 0;
}
