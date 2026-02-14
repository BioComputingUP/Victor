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
 * Author: Alberto Tilocca
 *
 */
#include <Spacer.h>
#include <Group.h>
#include <SideChain.h>
#include <vector3.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <XyzSaver.h>
#include <SeqSaver.h>
#include <IoTools.h>
#include <IntCoordConverter.h>
#include <IntSaver.h>
#include <GetArg.h>
#include "XyzLoader.h"
#include <cmath>
#include <Eigen/Geometry>

#include "CEFunctions.h"

using namespace Victor;using namespace Victor::Biopool;


void sShowHelp() {
    cout << "\n\nCE - Combinatiorial Extension -- calculate the protein structure alignment by incremental combinatorial extension of the optimal path\n"
        << "Options: \n"
        << "\t -i first pdb inputFile , \t \t(mandatory)\n"
        << "\t -j second pdb inputFile, \t \t(mandatory)\n"
        << "\t -m AFP size (defaul 8)\n"
        << "\t -c select the chain of the first pdb file   (default: first chain)\n"
        << "\t -d select the chain  of the second pdb file (default: first chain)\n"
        << "\t -l select the model of the first pdb file\n"
        << "\t -n select the model of the second pdb file\n"
        << "\t -o outputFile (deafult: chainID1_chainID2.pdb)\n"
        << "\t -v set verbose\n"
        << "\t -h this message\n"
	<< "Usage exemple: \n"
	<< "\t CEmain -i  input_1cpc.pdb -j  input_1col.pdb -c L -d A \n";
}

int main(int argc, char** argv) {

     if (getArg("h", argc, argv)) {
        sShowHelp();
        return 1;
    }

    string inputFile1; 
    string inputFile2; 
    string chainID1; 
    string chainID2; 
    int m;
    unsigned int modelNum1; 
    unsigned int modelNum2; 

    // Output
    string outputFile;
    string outputFile2="onlyModifiedChain.pdb";
   
     bool verbose; 

    // Input
    getArg("i", inputFile1, argc, argv, "!");
    getArg("j", inputFile2, argc, argv, "!");
    getArg("c", chainID1, argc, argv, "!");
    getArg("d", chainID2, argc, argv, "!");
    getArg("m", m, argc, argv, 8);
    getArg("l", modelNum1, argc, argv, 999);
    getArg("n", modelNum2, argc, argv, 999);
    verbose = getArg("v", argc, argv); 
    // Output
    string outputFileDefault="chains"+chainID1+"_"+chainID2+".pdb";
    getArg("o", outputFile, argc, argv, outputFileDefault);    

    // ---------------------------------------------------------------------- //


    if (inputFile1 == "!" || inputFile2 == "!" )  {
        cout << "Missing input files specification. Aborting. " << endl;
        return -1;
    }

    ifstream inFile1(inputFile1.c_str());
    ifstream inFile2(inputFile2.c_str());
    // Initialize PdbLoaders
    PdbLoader pl1(inFile1);
    PdbLoader pl2(inFile2);

    pl1.setModel(modelNum1);
    pl2.setModel(modelNum2);

    // Set chain(s) 1

    if (chainID1 != "!")
        pl1.setChain(chainID1[0]);
    else
        pl1.setChain(pl1.getAllChains()[0]);

    // Set chain(s) 2 
    if (chainID2 != "!")
        pl2.setChain(chainID2[0]);
    else
        pl2.setChain(pl2.getAllChains()[0]);

    if (!verbose){
        pl1.setNoVerbose();  // Set PdbLoader verbosity
        pl2.setNoVerbose();
    }

    if (verbose)
        cout << "PdbLoader set" << endl;

    // Create and load the protein object
    Protein* prot1 = new Protein ;
    Protein* prot2 = new Protein ;
    
    prot1->load(pl1);
    prot2->load(pl2);
    
    
    inFile1.close();
    inFile2.close();
    //Prendo gli spacer delle due catene

    Spacer* sp1 = prot1->getSpacer(chainID1[0]); // Get spacer
    Spacer* sp2 = prot2->getSpacer(chainID2[0]); // Get spacer
   
    int dimSpacer1=(int)sp1->sizeAmino();    
    int dimSpacer2=(int)sp2->sizeAmino();
    int combination=(m-1)*(m-2)/2; //to exclude the neighbors
    int smaller = ( dimSpacer1 < dimSpacer2 ) ? dimSpacer1 : dimSpacer2;
  
   //calculate the inner distances for each spacer.
    vector <vector<double> > distA(dimSpacer1);//inner distance matrix chain1
    vector <vector<double> > distB(dimSpacer2);//inner distance matrix chain1
 
     for (int i = 0; i < dimSpacer1; i++)
        distA[i].resize(dimSpacer1);
  
     for (int i = 0; i< dimSpacer2; i++)
        distB[i].resize(dimSpacer2);
    
   CEFunctions::innerDistance(sp1,distA);
   CEFunctions::innerDistance(sp2,distB);

   //Calculate the score Matrix
   vector<vector<double> > score(dimSpacer1);
   for (int i = 0; i< dimSpacer1; i++)
     score[i].resize(dimSpacer2);
  
   CEFunctions::createSMatrix(sp1,sp2, combination, m,score,distA,distB);
    
    
  //search of the path
    int MAXPATH=30;
    vector<double> scoreBuffer(MAXPATH);  
    vector<vector<afp> > tempPathBuffer(smaller);
    for (int i = 0; i< smaller; i++)
      tempPathBuffer[i].resize(smaller);
    vector<int> lenBuffer(MAXPATH);
    
    CEFunctions::CESearch(sp1,sp2,tempPathBuffer, scoreBuffer,score,smaller,combination,m,distA, distB,lenBuffer);
  
    //here I have the path
       
   // superimposition with kabash algorithm and take the path who have the lower RMSD
    
    int indexFinalPath = CEFunctions::finalSteps(sp1,sp2,tempPathBuffer,smaller,lenBuffer ,m,outputFile,outputFile2);
    
    cout << "scoreCE: "<<scoreBuffer[indexFinalPath]<<endl;

}
