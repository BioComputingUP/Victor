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

// --*- C++ -*------x-----------------------------------------------------------
//
//
// Description:     ClustalW allignemt.
//                  http://www.ebi.ac.uk/Tools/msa/clustalw2/
//
// -----------------x-----------------------------------------------------------


// Includes:
#include <SubMatrix.h>
#include <Alignment.h>
#include <AlignmentData.h>
#include <SequenceData.h>
#include <NWAlign.h>
#include <AGPFunction.h>
#include <ScoringS2S.h>

#include <Profile.h>
#include <HenikoffProfile.h>
#include <PSICProfile.h>
#include <SeqDivergenceProfile.h>

#include <NewickTree.h>
#include <PhyloSupport.h>
#include <ClustalW.h>

#include <sstream>
#include <string>
#include <vector>

// Global constants, typedefs, etc. (to avoid):
using namespace Victor::Align2;
using namespace Victor;
using namespace std;

namespace Victor { namespace Phylo{

	// CONSTRUCTORS:

	ClustalW::ClustalW(){
		score=0;
		tokenSize=0;
	}

	ClustalW::ClustalW(NewickTree gT) {

		cout<<"constructor"<<endl;
		guideTree=gT;
        progressiveAlign();
	}

	/**
	 *@Description Basic destructor
	 */
	ClustalW::~ClustalW() {
	}

    // OPERATORS:
	/**
	 *@Description Create progressive Align strating from a guide tree.
	 *
	 *@return
	 */
	void ClustalW::progressiveAlign(){

		cout<<" Progressive alignament Start "<<endl;
		vector<string> tmpV(2);
		int tokenSize=PhyloSupport::tokenSize;


		vector <iNode*> nodeTree(guideTree.getNumberOfLeaf());
		for(unsigned int i=0;i<guideTree.getNumberOfLeaf();i++){
			nodeTree[i]=guideTree.getLeafInPosition(i).getRoot();
		}

		if(tokenSize<1)
		{
			tokenSize=nodeTree[0]->seq.size()/10;
			cout<<"Use Default token Size "<<tokenSize<<endl;
		}

		setTokenSize(tokenSize);

		vector <string> seqV(1);
		vector <double> weigthV(1);
		vector <string> nameV(1);
		vector <string> nameFinal(guideTree.getNumberOfLeaf());
		cout<<"percent to the end:"<<endl;
		unsigned int j=1;
		while(nodeTree.size()>1){
			cout<<""<<((j)*100/(guideTree.getNumberOfLeaf())*100)/100<<"% "<<endl;
			for(unsigned int i=0;i<nodeTree.size();i++){//for each node in nodeTree
				if(nodeTree.size()==2 && i==0) {
					if(nodeTree[0]->allignSeq[0].size()<nodeTree[1]->allignSeq[0].size()){
						tmpV=PhyloSupport::AlingMultiSvsMultiS(nodeTree[0]->allignSeq,nodeTree[1]->allignSeq,nodeTree[0]->weigthV,nodeTree[1]->weigthV,false,tokenSize);
						//insert name
						nodeTree[0]->nameV=PhyloSupport::mergeStringVector(nodeTree[0]->nameV,nodeTree[1]->nameV);
						//insert weigth
						nodeTree[0]->weigthV=PhyloSupport::mergeDoubleVector(nodeTree[0]->weigthV,nodeTree[1]->weigthV);
					}
					else {
						tmpV=PhyloSupport::AlingMultiSvsMultiS(nodeTree[1]->allignSeq,nodeTree[0]->allignSeq,nodeTree[1]->weigthV,nodeTree[0]->weigthV,false,tokenSize);
						//insert name
						nodeTree[0]->nameV=PhyloSupport::mergeStringVector(nodeTree[1]->nameV,nodeTree[0]->nameV);
						//insert weigth
						nodeTree[0]->weigthV=PhyloSupport::mergeDoubleVector(nodeTree[1]->weigthV,nodeTree[0]->weigthV);
					}
					nodeTree[1]->ClustalW=true;
					nodeTree[0]->allignSeq=tmpV;

				}
				else if(j==1){//populate nodeTree
					seqV[0]=nodeTree[i]->seq;
					weigthV[0]=nodeTree[i]->weigth;
					nameV[0]=nodeTree[i]->name;
					nodeTree[i]->nameV=nameV;
					nodeTree[i]->allignSeq=seqV;
					nodeTree[i]->weigthV=weigthV;
				}
				else if(j==2 && nodeTree[i]->parent->numberOfChildLeaf==2){//first lap of align only seq VS seq
					tmpV=PhyloSupport::AlingSvsS(nodeTree[i]->seq,nodeTree[i]->seq);
					weigthV[0]=nodeTree[i]->weigth;
					nameV[0]=nodeTree[i]->name;
					if(weigthV.size()==1)
						weigthV.push_back(nodeTree[i+1]->weigth);
					else
						weigthV[1]=nodeTree[i+1]->weigth;
					if(nameV.size()==1)
						nameV.push_back(nodeTree[i+1]->name);
					else
						nameV[1]=nodeTree[i+1]->name;

					nodeTree.erase(nodeTree.begin()+i+1);
					nodeTree[i]->parent->allignSeq=tmpV;
					nodeTree[i]->parent->weigthV=weigthV;
					nodeTree[i]->parent->nameV=nameV;
					nodeTree[i]=nodeTree[i]->parent;
				}
				else if(j>2 && nodeTree[i]->parent->numberOfChildLeaf==(int)j && !nodeTree[i]->ClustalW && i<nodeTree.size()){

					if(nodeTree[i]->isLeft){
						if(nodeTree[i]->allignSeq[0].size()<nodeTree[i]->parent->right->allignSeq[0].size()) {
							tmpV=PhyloSupport::AlingMultiSvsMultiS(nodeTree[i]->allignSeq,nodeTree[i]->parent->right->allignSeq,nodeTree[i]->weigthV,nodeTree[i]->parent->right->weigthV,false,tokenSize);
							//insert name
							nodeTree[i]->parent->nameV=PhyloSupport::mergeStringVector(nodeTree[i]->nameV,nodeTree[i]->parent->right->nameV);
							//insert weigth
							nodeTree[i]->parent->weigthV=PhyloSupport::mergeDoubleVector(nodeTree[i]->weigthV,nodeTree[i]->parent->right->weigthV);
						}
						else {
							tmpV=PhyloSupport::AlingMultiSvsMultiS(nodeTree[i]->parent->right->allignSeq,nodeTree[i]->allignSeq,nodeTree[i]->parent->right->weigthV,nodeTree[i]->weigthV,false,tokenSize);
							//insert name
							nodeTree[i]->parent->nameV=PhyloSupport::mergeStringVector(nodeTree[i]->parent->right->nameV,nodeTree[i]->nameV);
							//insert weigth
							nodeTree[i]->parent->weigthV=PhyloSupport::mergeDoubleVector(nodeTree[i]->parent->right->weigthV,nodeTree[i]->weigthV);
						}
						nodeTree[i]->parent->right->ClustalW=true;
					}
					else{//is right child
						if(nodeTree[i]->allignSeq[0].size()<nodeTree[i]->parent->left->allignSeq[0].size()) {
							tmpV=PhyloSupport::AlingMultiSvsMultiS(nodeTree[i]->allignSeq,nodeTree[i]->parent->left->allignSeq,nodeTree[i]->weigthV,nodeTree[i]->parent->right->weigthV,false,tokenSize);
							//insert name
							nodeTree[i]->parent->nameV=PhyloSupport::mergeStringVector(nodeTree[i]->nameV,nodeTree[i]->parent->left->nameV);
							//insert weigth
							nodeTree[i]->parent->weigthV=PhyloSupport::mergeDoubleVector(nodeTree[i]->weigthV,nodeTree[i]->parent->left->weigthV);
						}
						else {
							tmpV=PhyloSupport::AlingMultiSvsMultiS(nodeTree[i]->parent->left->allignSeq,nodeTree[i]->allignSeq,nodeTree[i]->parent->right->weigthV,nodeTree[i]->weigthV,false,tokenSize);
							//insert name
							nodeTree[i]->parent->nameV=PhyloSupport::mergeStringVector(nodeTree[i]->parent->left->nameV,nodeTree[i]->nameV);
							//insert weigth
							nodeTree[i]->parent->weigthV=PhyloSupport::mergeDoubleVector(nodeTree[i]->parent->left->weigthV,nodeTree[i]->weigthV);
						}
						nodeTree[i]->parent->left->ClustalW=true;
					}
					nodeTree[i]=nodeTree[i]->parent;
					nodeTree[i]->allignSeq=tmpV;
				}
				else if(nodeTree[i]->ClustalW && i<nodeTree.size()){
					nodeTree.erase(nodeTree.begin()+i);
				}

			}
			j++;
		}
		cout<<"...complete!"<<endl;

		string outString=ClustalW::printClustalWFromat(nodeTree[0]->allignSeq,nodeTree[0]->nameV);

		cout<<endl<<outString<<endl;

		//write clustalW
	    ofstream outFile("out.clustalw");
	    if(!outFile) {
	        cout<<"Error During creation File out.clustalw";
	    }
	   	outFile<<outString<<endl;

	    setScore(scoreClustalW(nodeTree[0]->allignSeq,nodeTree[0]->weigthV));
	    cout<<endl<<"score of ClustalW "<<getScore()<<" Token Size for Graph in MultiAling Change TokenSize for best align command --t. Actualy "<<tokenSize<<endl;

	    outFile.close();
        cout<<"Creation File out.clustalw Complete"<<endl;

		cout<<"endl clustaW"<<endl;

	}


    // HELPERS:
	/**
	 *@Description Returns String which represents the correct rappresentation in ClustalW format with no name of seq
	 *@param vector<string> seq Vector of sequence.
	 *@return string txt Correct rappresetation in CLustalW format(50 char by lane)
	 */
	string ClustalW::printClustalWFromat(vector<string> seq){
		vector <string> names(seq.size());
		for(unsigned int i=0; i<names.size();i++)
			names[i]="seq_"+PhyloSupport::intToString(i);
		return ClustalW::printClustalWFromat(seq,names);
	}

	/**
	 *@Description Returns String which represents the correct rappresentation in ClustalW format with correct name of seq
	 *@param vector<string> seq Vector of sequence.
	 *@return string txt Correct rappresetation in CLustalW format(50 char by lane)
	 */
	string ClustalW::printClustalWFromat(vector<string> seq, vector <string> names){
		string txt="CLUSTALW \n\n";
		string tmp="";
		unsigned int j=0;
		unsigned int index=0;
		unsigned int dim=50;
		vector <unsigned int> countSimb(seq.size(),0);
		bool onlyGap=true;
		while(j<seq[0].size()){
			for(unsigned int i=0;i<seq.size();i++){
				txt+=names[i]+" ";
				for(j=dim*index;j<=dim+dim*index && j<seq[i].size();j++){
					txt+=seq[i][j];
					if(seq[i][j]!='-') {
						countSimb[i]++;
						onlyGap=false;
					}
				}
				if(!onlyGap)
					txt+=" "+PhyloSupport::intToString(countSimb[i]);
				else
					txt+=" ";
				onlyGap=true;
				txt+="\n";
			}
			txt+="\n\n";
			index++;
		}
		return txt;
	}


	/**
	 *@Description Calculate A sort of simalarity score on ClustalW align consider how much each column of char are similar no weigth consider
	 *@param vector<string> seq Vector of sequence.
	 *@return double score the score of similarity.
	 */
	double ClustalW::scoreClustalW(vector<string> allignSeq){
		//matrix config
    	string path = getenv("VICTOR_ROOT");
		if (path.length() < 3)
			cout << "Warning: environment variable VICTOR_ROOT is not set." << endl;

		string dataPath = path + "data/";

		//Default matrix
		string matrixFileName=PhyloSupport::matrix;
		matrixFileName = dataPath + matrixFileName;
		ifstream matrixFile(matrixFileName.c_str());
		if (!matrixFile)
			ERROR("Error opening substitution matrix file.", exception);
		SubMatrix sub(matrixFile);
		//end matrix config
		double score=0;
		for(unsigned int x=0;x<allignSeq.size();x++){//for each vector do score
			for(unsigned int j=0;j<allignSeq[x].size();j++){//num of char in one string
				for(unsigned int i=0; i<allignSeq.size();i++){//num of string in the vector
				if(i!=x)
					score+=sub.score[allignSeq[x][j]][allignSeq[i][j]];
				}
			}
		}
		return score;
	}
	/**
	 *@Description Calculate A sort of simalarity score on ClustalW align consider how much each column of char are similar no weigth consider
	 *@param vector<string> seq Vector of sequence.
	 *@param vector<double> w vector of weigth from guide tree
	 *@return double score the score of similarity.
	 */
	double ClustalW::scoreClustalW(vector<string> allignSeq, vector<double> w){
		//matrix config
    	string path = getenv("VICTOR_ROOT");
		if (path.length() < 3)
			cout << "Warning: environment variable VICTOR_ROOT is not set." << endl;

		string dataPath = path + "data/";

		//Default matrix
		string matrixFileName=PhyloSupport::matrix;
		matrixFileName = dataPath + matrixFileName;
		ifstream matrixFile(matrixFileName.c_str());
		if (!matrixFile)
			ERROR("Error opening substitution matrix file.", exception);
		SubMatrix sub(matrixFile);
		//end matrix config
		double score=0;
		for(unsigned int x=0;x<allignSeq.size();x++){//for each vector do score
			for(unsigned int j=0;j<allignSeq[x].size();j++){//num of char in one string
				for(unsigned int i=0; i<allignSeq.size();i++){//num of string in the vector
				if(i!=x)
					score+=sub.score[allignSeq[x][j]][allignSeq[i][j]]*w[x]*w[i];
				}
			}
		}
		return score;
	}

    // PREDICATES:
    void ClustalW::setTokenSize(unsigned int t){
    	tokenSize=t;
    }
    void ClustalW::setScore(double s){
    	score=s;
    }

    // PREDICATES:
    unsigned int ClustalW::getTokenSize(){
    	return tokenSize;
    }
    double ClustalW::getScore(){
    	return score;
    }

}} // namespace
