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


// Includes:
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include "PhyloSupport.h"
#include <SubMatrix.h>
#include <Alignment.h>
#include <AlignmentData.h>
#include <SequenceData.h>
#include <NWAlign.h>
#include <AGPFunction.h>
#include <ScoringS2S.h>
#include <IoTools.h>
#include <math.h>
#include <SeqNodeGraph.h>

// Global constants, typedefs, etc. (to avoid):
using namespace Victor::Align2;
using namespace Victor;
using std::string;
using std::vector;
using std::ifstream;
using std::cout;
using std::endl;


namespace Victor { namespace Phylo{


    /** @brief
     *
     **/
	// CONSTRUCTORS/DESTRUCTOR:

	/**
	 *@Description Basic contructor
	 */
	PhyloSupport::PhyloSupport(){

	}


	/**
	 *@Description Basic destructor
	 */
	PhyloSupport::~PhyloSupport() {

	}

    // ATTRIBUTES:
	double PhyloSupport::openGapPenalty=12;
	double PhyloSupport::extensionGapPenalty=3;
	double PhyloSupport::cSeq=0.8;
	double PhyloSupport::downs=999.9;
	double PhyloSupport::downa=999.9;
	double PhyloSupport::ups=999.9;
	double PhyloSupport::upa=999.9;
	int PhyloSupport::tokenSize=-1;
	unsigned int PhyloSupport::weightingScheme=0;
	string PhyloSupport::matrix="blosum62.dat";
	bool PhyloSupport::verbose=false;


	vector<Alignment> PhyloSupport::calcAlignmentV(Alignment *aliSec, vector<vector<double> > &distance , bool ktuples,bool verbose){
		string seq1Name, seq2Name, seq1, seq2, sec1, sec2;

		Alignment newAli;
		string path = getenv("VICTOR_ROOT");
		if (path.length() < 3)
			cout << "Warning: environment variable VICTOR_ROOT is not set." << endl;

		string dataPath = path + "data/";

		//n=4? no idea. Default in APP/subali.cc
		int n=4;

		AlignmentData *ad;

		//Default matrix
		string matrixFileName=PhyloSupport::matrix;
		matrixFileName = dataPath + matrixFileName;
		ifstream matrixFile(matrixFileName.c_str());
		if (!matrixFile)
			ERROR("Error opening substitution matrix file.", exception);
		SubMatrix sub(matrixFile);

		//Default gap function

		GapFunction *gf = new AGPFunction(openGapPenalty, extensionGapPenalty);


		//Default Structure
		 Structure *str;
		 str = 0;

		//Default ScoringScheme
		ScoringScheme *ss;

		//Global Align
		Align *a;
		if(verbose)
			cout << "Start Suboptimal Needleman-Wunsch alignments:" << endl;

		double suboptPenaltyMul=1;
		double suboptPenaltyAdd=1;

		unsigned int suboptNum=1;

		vector<Alignment> a2;

		Alignment a3;

		vector <Alignment> alignV(aliSec->size()+1);//vector of align

		double score=0;
		//size not count target. so +1
		for(unsigned int index=0;index<aliSec->size()+1;index++){
			if(verbose)
				cout<<" index="<<index;
			if(index==aliSec->size()){
				seq1 = Alignment::getPureSequence(aliSec->getTarget());
				seq1Name = aliSec->getTargetName();
			}
			else{
				seq1 = Alignment::getPureSequence(aliSec->getTemplate(index));
				seq1Name = aliSec->getTemplateName(index);
			}
			for(unsigned int j=0;j<aliSec->size()+1;j++){
				if(verbose)
					cout<<" j="<<index;
				if(j==aliSec->size()){
					seq2 = Alignment::getPureSequence(aliSec->getTarget());
					seq2Name = aliSec->getTargetName();
				}else{
					seq2 = Alignment::getPureSequence(aliSec->getTemplate(j));
					seq2Name = aliSec->getTemplateName(j);
				}

				ad = new SequenceData(n, seq1, seq2, seq1Name, seq2Name);

				ss = new ScoringS2S(&sub, ad, str, cSeq);

				a = new NWAlign(ad, gf, ss);

				a->setPenalties(suboptPenaltyMul, suboptPenaltyAdd);
				a2 = a->generateMultiMatch(suboptNum);
				if (a2.size() == 0)
					ERROR("No output alignments generated.", exception);

				a2[0].cutTemplate(1);
				if(ktuples==true)
				{
					score=PhyloSupport::distanceCalcTwoSeqktuples(a2[0].getTarget(),a2[0].getTemplate(0));
				}
				else
					score=PhyloSupport::distanceCalcTwoSeq(a2[0].getTarget(),a2[0].getTemplate(0));
				distance[index][j]=score;
				if(verbose)
					cout<<"SCORE="<<distance[index][j]<<" ";
				if(j==0)
					newAli = a2[0];
				else
					newAli.addAlignment(a2[0]);
			}//end for j
			if(verbose){
				cout<<" seq1= "<<seq1Name<<endl;
				cout<<" seq2= "<<seq2Name<<endl;
			}
			alignV[index]=newAli;
		}//end for index

		if(verbose)
			cout<<"NWAlign ok"<<endl;
		return alignV;
	}


	vector<string> PhyloSupport::AlingSvsS(string seq1,string seq2,bool verbose){

		string path = getenv("VICTOR_ROOT");
		if (path.length() < 3)
			cout << "Warning: environment variable VICTOR_ROOT is not set." << endl;

		string dataPath = path + "data/";

		//n=4? no idea. Default in APP/subali.cc
		int n=4;

		AlignmentData *ad;

		//Default matrix
		string matrixFileName=PhyloSupport::matrix;
		matrixFileName = dataPath + matrixFileName;
		ifstream matrixFile(matrixFileName.c_str());
		if (!matrixFile)
			ERROR("Error opening substitution matrix file.", exception);
		SubMatrix sub(matrixFile);

		//Default gap function

		GapFunction *gf = new AGPFunction(openGapPenalty, extensionGapPenalty);


		//Default Structure
		 Structure *str;
		 str = 0;

		//Default ScoringScheme
		ScoringScheme *ss;

		//Global Align
		Align *a;
		if(verbose)
			cout << "Start Suboptimal Needleman-Wunsch alignments:" << endl;

		double suboptPenaltyMul=1;
		double suboptPenaltyAdd=1;

		unsigned int suboptNum=1;

		vector<Alignment> a2;


		ad = new SequenceData(n, seq1, seq2, "seq1", "seq2");

		ss = new ScoringS2S(&sub, ad, str, cSeq);

		a = new NWAlign(ad, gf, ss);

		a->setPenalties(suboptPenaltyMul, suboptPenaltyAdd);
		a2 = a->generateMultiMatch(suboptNum);
		if (a2.size() == 0)
			ERROR("No output alignments generated.", exception);
		vector<string> result(2);

		result[0]=a2[0].getTarget();
		result[1]=a2[0].getTemplate();
		if(verbose){
			cout<<"result "<<a2.size()<<" result[1] "<<result[1]<<endl<<" result[0] "<<result[0]<<endl;
			cout<<"NWAlign SvsS Oexit-------------"<<endl;
		}
		return result;
}


    vector<string> PhyloSupport::AlingMultiSvsMultiS(vector <string> seq1,vector <string> seq2,vector <double> vWeigth1,vector <double> vWeigth2,bool verbose,int tokenSize){

    	string gap="-";
    	string tokenSizeGap="";
    	string tokenSizeGap2="";
    	int tokenSize2=tokenSize;

    	//insert gap at end for string with same size//not need
    	while(seq1[0].size()>seq2[0].size()){
        	for(unsigned int i=0; i<seq2.size();i++){
        		seq2[i]+=gap;
        	}
    	}
    	while(seq1[0].size()<seq2[0].size()){
        	for(unsigned int i=0; i<seq1.size();i++){
        		seq1[i]+=gap;
        	}
    	}

    	for(int i=0;i<tokenSize;i++){
    		tokenSizeGap+=gap;
    	}
    	for(int i=0;i<tokenSize2;i++){
    		tokenSizeGap2+=gap;
    	}

    	if(verbose)
    		cout<<"AlingMultiSvsMultiS"<<endl;
 		while(seq1[0].size()%tokenSize!=0){
 			for(unsigned int i=0;i<seq1.size();i++){
 				seq1[i]+="-";
 			}
 		}

 		while(seq2[0].size()%tokenSize2!=0){
 			for(unsigned int i=0;i<seq2.size();i++){
 				seq2[i]+="-";
 			}
 		}

    	vector <SeqNodeGraph*> tokenS1(seq1[0].size()/tokenSize);
    	vector <SeqNodeGraph*> tokenS2(seq2[0].size()/tokenSize2);

    	for(unsigned int i=0; i<seq1[0].size()/tokenSize;i++)
    	{//for all char
    		vector<string> tempSV(seq1.size());
    		for(unsigned int j=0; j<seq1.size();j++)
    		{//for all seq[j]
    			tempSV[j]=seq1[j].substr(i*tokenSize,tokenSize);
    		}

    		tokenS1[i]= new SeqNodeGraph(i,i-1+tokenSize,seq1.size(),tempSV,seq2[0].size()/tokenSize);

    	}

    	for(unsigned int i=0; i<seq2[0].size()/tokenSize2;i++)
    	{
    		vector<string> tempSV(seq2.size());
    		for(unsigned int j=0; j<seq2.size();j++)
    		{
    			tempSV[j]=seq2[j].substr(i*tokenSize2,tokenSize2);
    		}

    		tokenS2[i]= new SeqNodeGraph(i,i-1+tokenSize2,seq2.size(),tempSV,seq1[0].size()/tokenSize2);
    	}



    	for(unsigned int i=0; i<tokenS1.size();i++){
    		SeqNodeGraph::setNodeWithWeigth(tokenS1[i],tokenS2,vWeigth1,vWeigth2);
		}



    	vector <int> edgeFor1(tokenS1.size()+tokenS2.size(),-1);
    	vector <string> finalS(seq1.size()+seq2.size(),"");
    	unsigned int count1=0;
    	unsigned int count2=0;
    	while(count1<tokenS1.size() && count2<tokenS2.size()){

    		if(count1==0) {
        		edgeFor1[count1]=tokenS1[count1]->returnBestEdge();
        		count2=edgeFor1[count1];
        		count1++;
    		}
    		else {
    			edgeFor1[count1]=tokenS1[count1]->returnBestEdgeAfterIndex(count2);
    			count2=edgeFor1[count1];
    			count1++;
    		}

    	}


    	//insert in final vector of string
    	for(unsigned int i=0; i<tokenS1.size();i++){
			if(edgeFor1[i]!=-1 ){
				for(unsigned x=0;x<seq1.size();x++){
					if(i==0){
						for(int j=0; j<edgeFor1[i];j++){
							finalS[x]+=tokenSizeGap2;
						}//end for j
						finalS[x]+=tokenS1[i]->getTokenSeq(x);
					}//if
					else{
						for(int j=0; j<edgeFor1[i]-edgeFor1[i-1]-1;j++){
							finalS[x]+=tokenSizeGap2;
						}
						finalS[x]+=tokenS1[i]->getTokenSeq(x);
					}
				}//end for x
			}//end if
			else{
				for(unsigned x=0;x<seq1.size();x++){
					finalS[x]+=tokenS1[i]->getTokenSeq(x);
				}
			}
    	}//ed for i

    	for(unsigned int i=0; i<seq2.size();i++){
    		finalS[i+seq1.size()]+=seq2[i];
    	}

    	//insert gap at end
    	while(finalS[0].size()>finalS[seq1.size()].size()){
        	for(unsigned int i=0; i<seq2.size();i++){
        		finalS[i+seq1.size()]+=gap;
        	}
    	}
    	while(finalS[0].size()<finalS[seq1.size()].size()){
        	for(unsigned int i=0; i<seq1.size();i++){
        		finalS[i]+=gap;
        	}
    	}

    	//delete useless gap
    	bool uselessGap=true;
    	while(uselessGap){
    		uselessGap=true;
    		for(unsigned int i=1; i<finalS.size();i++)
    		{
    			if(finalS[i][finalS[i].size()-1]!=finalS[i-1][finalS[i-1].size()-1] || finalS[i][finalS[i].size()-1]!='-'){
    				uselessGap=false;
    			}
    		}

    		if(uselessGap){
				for(unsigned int i=0; i<finalS.size();i++)
				{
					finalS[i]=finalS[i].substr(0,finalS[i].size()-1);
				}
    		}
    	}
    	return finalS;
    }


    /**
     *  Calculate the distance for two sequence. distance = 1-identity(PID)
     *  Percent sequence identity (PID)=Calculates the identity as the number of identical positions,
     *  it gives a positive value when the two sequences have a gap in the same position,
     *  divided by the length of the alignment.
     *
     */
   double PhyloSupport::distanceCalcTwoSeq(string seq1,string seq2){
        if (seq1.length() != seq2.length())
            cout << "Warning: sequence lengths do not match:\n"
                << "seq1 = " << seq1.length() << "\n"
            << "seq2 = " << seq2.length() << "\n";

        double tmp = 0;
        for (unsigned int i = 0; i < seq1.length(); i++){
            if ( seq1[i] == seq2[i])
                tmp++;
        }

        return  (1- (tmp / seq1.length()));
    }


   /**
    *  Calculate the distance for two sequence by ktuples formula.
    *  return number of k-tuples (number of residues identical) minus
	*	constant penality for each gap.
    *
    */
  double PhyloSupport::distanceCalcTwoSeqktuples(string seq1,string seq2){
       if (seq1.length() != seq2.length())
           cout << "Warning: sequence lengths do not match:\n"
               << "seq1 = " << seq1.length() << "\n"
           << "seq2 = " << seq2.length() << "\n";

       double penality = 0;
       double d=0;
       for (unsigned int i = 0; i < seq1.length(); i++){
    	   if(seq1[i]=='-' || seq2[i]=='-')
    		   penality+=1;
    	   else if( seq1[i] == seq2[i] )
               d++;
       }
       penality=(penality/seq1.length())*0.9*d;

       return 1-((d-penality)/seq1.length());
   }


  	  string PhyloSupport::insertGapPosition(string seq, int position){
  		  string part1="";
  		  string part2="";
  		  part1=seq.substr(0,position);
  		  part2=seq.substr(position);
  		  seq=part1+"-"+part2;
  		  return seq;
  	  }



	std::vector<std::string> &PhyloSupport::split(const std::string &s, char delim, std::vector<std::string> &elems) {
		    std::stringstream ss(s);
		    std::string item;
		    while (std::getline(ss, item, delim)) {
		        elems.push_back(item);
		    }
		    return elems;
		}

		std::vector<std::string> PhyloSupport::split(const std::string &s, char delim) {
		    std::vector<std::string> elems;
		    split(s, delim, elems);
		    return elems;
		}


		double PhyloSupport::calcDivR(vector<double> vDist){
			double tmpR=0;
			for(unsigned int i=0;i<vDist.size();i++){
				tmpR+=vDist[i];
			}
			return tmpR;
		}
	    /**
	     * Create a string of matrix value;
	     *
	     */
		void PhyloSupport::printMatrix(vector<vector<double> > &distance){
			string sMatrix="";
			//std::setprecision(3)//for minus precision
			for(unsigned int i=0; i<distance.size();i++){
				for(unsigned int j=0; j<distance.size();j++){
					cout<<distance[i][j]<<" ";
				}
				cout<<endl;
			}
		}


		string PhyloSupport::intToString ( int num )
		  {
		     ostringstream ss;
		     ss << num;
		     return ss.str();
		  }

		vector<string> PhyloSupport::mergeStringVector( vector<string> v1,vector<string> v2 ){
			vector<string> v(v1.size()+v2.size());
			for(unsigned int i=0; i<v1.size();i++){
				v[i]=v1[i];
			}
			for(unsigned int i=0; i<v2.size();i++){
				v[i+v1.size()]=v2[i];
			}
			return v;

		}
		vector<double> PhyloSupport::mergeDoubleVector( vector<double> v1,vector<double> v2 ){
			vector<double> v(v1.size()+v2.size());
			for(unsigned int i=0; i<v1.size();i++){
				v[i]=v1[i];
			}
			for(unsigned int i=0; i<v2.size();i++){
				v[i+v1.size()]=v2[i];
			}
			return v;

		}


}} // namespace
