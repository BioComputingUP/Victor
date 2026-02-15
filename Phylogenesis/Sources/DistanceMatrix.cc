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

#include <DistanceMatrix.h>
#include <AGPFunction.h>
#include <SequenceData.h>
#include <ScoringS2S.h>
#include <NWGAlign.h>
#include <SWAlign.h>
#include <FSAlign.h>
#include <IoTools.h>
#include <sstream>
#include <string>
#include <iomanip>

namespace Victor { namespace Phylogenesis{

	DistanceMatrix::~DistanceMatrix(){}

	DistanceMatrix::DistanceMatrix(std::string fileName){

		std::ifstream inputFile(fileName.c_str());
		if (!inputFile)
			ERROR("Error opening distance matrix file.", exception);

		distanceMatrix.push_back(vector<double>());
		identityMatrix.push_back(vector<double>());
		int index_i = 0;
		int index_j = 1;

		std::string line;
		while (std::getline(inputFile, line))
		{
			std::istringstream iss(line);
			int a, b;
			double c;
			if (!(iss >> a >> b >> c))
				ERROR("Error reading distance matrix file.", exception);

			if(a == index_i && b == index_j)
			{
				distanceMatrix[a].push_back(c);
				identityMatrix[a].push_back(1-c);
				index_j++;
				continue;
			}

			if(a == index_i+1 && b == index_i+2)
			{
				distanceMatrix.push_back(vector<double>());
				identityMatrix.push_back(vector<double>());
				distanceMatrix[a].push_back(c);
				identityMatrix[a].push_back(1-c);
				index_i = a;
				index_j = b + 1;
				continue;
			}

			ERROR("Error reading input distance matrix file, invalid index.", exception);
		}
	}

	DistanceMatrix::DistanceMatrix(vector<Sequence*> sequences, AlignAlgorithm alignAlgorithm, SubMatrix* sub, double gop, double gep) :
			alignAlgorithm(alignAlgorithm), sequences(sequences), sub(sub) {
		   gf = new AGPFunction(gop, gep);
	}

	void DistanceMatrix::init(){

		if(distanceMatrix.size() > 0)
			return;

		cout<<"Start of Pairwise alignments\nAligning...\n"<<endl;

		for(unsigned int i=0;i<sequences.size()-1;i++){
			for(unsigned int j=i+1;j<sequences.size();j++){
				vector<string> result = pairwiseSequenceAlign(sequences[i]->sequence, sequences[j]->sequence);
				double distanceScore = calculateDistance(result[0], result[1]);
				appendDistance(i, distanceScore);
				double identityScore = calculateIdentity(result[0], result[1]);
				appendIdentity(i, identityScore);
				//Output
				cout<<fixed;
				cout<<setprecision(2);
				cout<<"Sequences ("<<setw(2)<< sequences[i]->id+1<<":"<<setw(2)<<sequences[j]->id+1<<") Aligned. Identity Score: "<<setw(6)<<identityScore*100<<" Score: "<<(1.0 - distanceScore)*100<<endl;
			}
		}
	}

	double DistanceMatrix::calculateDistance(string alignedS1, string alignedS2){
		double matches = 0;

		for(unsigned int i=0;i<alignedS1.size();i++){
			if(alignedS1[i] == '-' || alignedS2[i] == '-'){
				continue;
			}
			if(alignedS1[i] == alignedS2[i]){
				matches++;
			}
		}

		double sim = matches / min(AlignmentBase::getPureSequence(alignedS1).size(),AlignmentBase::getPureSequence(alignedS2).size());
		return 1 - sim;
	}

	double DistanceMatrix::calculateIdentity(string alignedS1, string alignedS2){
		double ungappedPos = 0;
		double matches = 0;

		for(unsigned int i=0;i<alignedS1.size();i++){
			if(alignedS1[i] == '-' || alignedS2[i] == '-'){
				continue;
			}
			ungappedPos++;
			if(alignedS1[i] == alignedS2[i]){
				matches ++;
			}
		}

		double sim = matches / ungappedPos;
		return sim;
	}

	void DistanceMatrix::appendDistance(int i, double value){
		if(i >= (int)distanceMatrix.size())
			distanceMatrix.push_back(vector<double>());

		distanceMatrix[i].push_back(value);
	}

	void DistanceMatrix::appendIdentity(int i, double value){
		if(i >= (int)identityMatrix.size())
			identityMatrix.push_back(vector<double>());

		identityMatrix[i].push_back(value);
	}

	double DistanceMatrix::at(int i, int j){
		int ii = std::min(i,j);
		int jj = std::max(i,j);

		return distanceMatrix[ii][jj-ii-1];
	}

	double DistanceMatrix::identityAt(int i, int j){
		int ii = std::min(i,j);
		int jj = std::max(i,j);

		return identityMatrix[ii][jj-ii-1];
	}

	vector<string> DistanceMatrix::pairwiseSequenceAlign(string s1, string s2){
		SequenceData *ad = new SequenceData(2, s1, s2);
		ScoringS2S *ss = new ScoringS2S(sub, ad, 0, 1.00);
		Align* a;
		switch(alignAlgorithm) {
			case global:
				a =  new NWGAlign(ad, gf, ss);
				break;
			case local:
				a = new SWAlign(ad, gf, ss);
				break;
			default:
				a = new FSAlign(ad, gf, ss);
				break;
		}

		vector<string> res = a->getMatch();
		delete a;
		return res;
	}

}}
