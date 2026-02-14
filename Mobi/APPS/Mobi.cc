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

    This file is part of Mobi extension created by
 	Luca Demo  -  University of Padova, department of Science
 */
/**
 */

#include <Protein.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <IoTools.h>
#include <GetArg.h>
#include <MobiProtein.h>
#include <TMScoreBin.h>
#include <MobiMethods.h>
#include <MobiUtils.h>

#include <String2Number.h>
#include <unistd.h>

using namespace Victor;
using namespace Victor::Biopool;
using namespace Victor::Mobi;

void sShowHelp() {
	cout << "Victor - Mobi NMR Mobility: 1.0\n"
		<< "Mobility defined accordingly with http://protein.bio.unipd.it/mobi/"
		<< "Options: \n"
		<< "\t-i <filename> \t Input PDB file (mandatory)\n"
		<< "\t-o <filename> \t Output mobility tracks to FASTA file (default stdout)\n"
		<< "\t-p <filename> \t Output ordered models in PDB file (default is no output)\n"
		<< "\t-s <filename> \t Output residue mobility \"score\" (default is no output)\n"
		<< "\t-c <id>       \t Chain identifier to read (default is first chain)\n"
		<< "\t-m <number>   \t NMR models numbers to read, comma separated list (default is all models)\n"
		<< "\t-t	\t TMscore binary path (default is ./TMscore)\n"
		<< "\t--d0 <number> \t d0 scaling factor for scaled distance. Default = 4\n"
		<< "\t--sd <number>  \t threshold for Avg Scaled Distance track. Default = 0.85\n"
		<< "\t--sdd <number> \t threshold for Scaled Distance Deviation track. Default = 0.09\n"
		<< "\t--phi <number> \t threshold for phi angles track. Default = 20\n"
		<< "\t--psi <number> \t threshold for psi angles track. Default = 20\n"
		<< "\t-v	\t verbose output\n"
		<< "\t-d	\t even more verbose output (lot of details)\n"
		<< "\t-h	\t shows this message\n";
}

int main(int argc, char* argv[]) {
	cout << "Victor Mobi Mobility Calculator" << endl;
	cout << "v. 0.1" << endl << endl;

	//If -h specified, show help message
    if (getArg("h", argc, argv)) {
        sShowHelp();
        return 1;
    }

    string inputFile, outputFile, outputPdb, outputScore, tmbin, chainID, modelList;
    string d0_s,sd_s,sdd_s,phi_s,psi_s;
    bool verbose, debug;
    vector<unsigned int> models;


    getArg("i", inputFile, argc, argv, "!");
    getArg("o", outputFile, argc, argv, "!");
    getArg("p", outputPdb, argc, argv, "!");
    getArg("s", outputScore, argc, argv, "!");
    getArg("t", tmbin, argc, argv, "!");
    getArg("c", chainID, argc, argv, "!");
    getArg("m", modelList, argc, argv, "!");
    getArg("-d0", d0_s, argc, argv, "!");
    getArg("-sd", sd_s, argc, argv, "!");
    getArg("-sdd", sdd_s, argc, argv, "!");
    getArg("-psi", psi_s, argc, argv, "!");
    getArg("-phi", phi_s, argc, argv, "!");
    verbose = getArg("v", argc, argv);
    debug = getArg("d", argc, argv);

    if (debug)
    	verbose = true;

    // Check input file
    if (inputFile == "!") {
        cout << "Missing input file specification. Aborting. (-h for help)" << endl;
        return -1;
    }
    ifstream inFile(inputFile.c_str());
    if (!inFile){
        ERROR("Input file not found.", exception);
    }else{
		if (verbose)
			cout << "Using PDB input file: " << inputFile << endl;
    }
    // check TM
    if (tmbin == "!") {
    	tmbin = "./TMscore";
    }
    if(access(tmbin.c_str(), X_OK) != 0){
    	ERROR("TM binary (" + tmbin + ") does not exist or not executable.", exception);
    }else{
    	if (verbose)
    		cout << "Using TMscore binary: " << tmbin << endl;
    }

    //Check mobi parameters
    double d0 = DEF_D0;
    if (d0_s != "!")
    	d0 = stodDEF(d0_s);
    double sd_th = DEF_SD_TH;
    if (sd_s != "!")
    	sd_th = stodDEF(sd_s);
    double sdsd_th = DEF_SDSD_TH;
    if (sdd_s != "!")
    	sdsd_th = stodDEF(sdd_s);
    double phi_th = DEF_PHI_TH;
    if (phi_s != "!")
    	phi_th = stodDEF(phi_s);
    double psi_th = DEF_PSI_TH;
	if (psi_s != "!")
		psi_th = stodDEF(psi_s);

    //MOBI
    PdbLoader pl(inFile);
    // Set PdbLoader variables
    if (!debug) {
        pl.setNoVerbose();
    }else{
    	pl.setVerbose();
    }

    // User selected chain
    if (chainID != "!") {
        if (chainID.size() > 1)
            ERROR("You can choose only 1 chain", error);
        pl.setChain(chainID[0]);
        pl.checkAndSetChain();
        if (verbose)
        	cout << "Selected chain: " << chainID[0] << endl;
    }// First chain
    else {
        if (pl.getAllChains().size() > 0){
			pl.setChain(pl.getAllChains()[0]);
			if (verbose)
				cout << "Selected chain: " << pl.getAllChains()[0] << endl;
			chainID = pl.getAllChains()[0];
        }else
        	if (verbose)
        		ERROR("No chains found. Quit...", exception);

    }

    // User selected models
    if (modelList != "!"){
    	models = sToVectorOfUIntDEF(translate(modelList,',',' '));
    		for (unsigned int i = 0; i < models.size(); i++)
    			if (models[i] > pl.getMaxModels() || models[i] < 1)
    				ERROR("You specified out of bound model number(s). Max for this input pdb is " + pl.getMaxModels(), exception);
    }else{
    	for (unsigned int i = 1; i <= pl.getMaxModels(); i++)
    		models.push_back(i);
    }
    if (verbose){
    	cout << "Selected models: ";
    	for (unsigned int i = 0; i < models.size(); i++)
    		cout << models[i] << " ";
    	cout << endl;
    }

    // Load the protein object
    MobiProtein* prot = new MobiProtein();
    if (debug)
    	prot->setVerbose(true);
    else
    	prot->setVerbose(false);
    if (verbose)
    	cout << "Loading models..." << endl;
    prot->load(pl,models);
    cout << endl;

    // TMscore: here we use the external binary
    TMScoreBin* tm = new TMScoreBin(tmbin,".");
    // Set MobiMethods Thresholds
    MobiMethods mm(d0, CA, sd_th, sdsd_th, psi_th, phi_th);
    if (verbose)
    	mm.verbosity(1);	//Outputs tracks on cout
    if (debug)
    	mm.verbosity(2);	//Output tracks and other calculation details

    //Mobility
    cout << "Mobility calculation..." << endl;
	vector<int> mobility = mm.mobiMobility(prot,tm);

	if (verbose){	//cout output
		MobiUtils::printFastaFile(mm);
		cout << endl;
	}
	if (outputFile != "!"){	//save fasta
		ofstream fastaOut(outputFile.c_str());
		MobiUtils::printFastaFile(mm,fastaOut);
		fastaOut.close();
		if (verbose)
			cout << "Tracks saved in FASTA file: " << outputFile << endl;
	}

	if (outputPdb != "!"){	//save Pdb
		ofstream pdbOut(outputPdb.c_str());
		MobiUtils::saveMobiPdb(prot,MobiUtils::sortModels(mm.getDistances()),mm.getSDMeans(), mm.getSDDevs(),pdbOut);
		pdbOut.close();
		if (verbose)
			cout << "Models saved in Pdb file: " << outputPdb << endl;
	}

	if (outputScore != "!"){	//Save Score
		ofstream scoreOut(outputScore.c_str());
		MobiUtils::saveScore(mm,chainID[0],scoreOut);
		scoreOut.close();
		if (verbose)
			cout << "Score saved in file: " << outputScore << endl;
	}



	//Close
	tm->cleanup();
	delete prot;
	delete tm;

	//Goodbye
	cout << endl << "Done." << endl;
    return 0;
}
