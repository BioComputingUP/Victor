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

#include <ClustalW.h>
#include <SingleCluster.h>
#include <iostream>
#include <algorithm>
#include <map>
#include <ClusterUtils.h>
#include <CWScoringScheme.h>
#include <CWGapFunction.h>
#include <MsaNWG.h>
#include <limits>
#include <cmath>
#include <iomanip>

namespace Victor { namespace Phylogenesis{

ClustalW::ClustalW(DistanceMatrix* dm,  HierarchicalClusterBuilder* clusterBuilder, vector<Sequence*> sequences, double cutoff, double gop, double gep, bool pam, bool verbose)
 : dm(dm), clusterBuilder(clusterBuilder), sequences(sequences), cutoff(cutoff), gop(gop), gep(gep), pam(pam), verbose(verbose){
	subMatrixFactory = new SubMatrixFactory(true);
}

MultipleAlignment* ClustalW::performMultipleSequenceAlignment(){

	// Stage 1: all pairs of sequences are aligned separately in order
	// to calculate a distance matrix giving the divergence of each pair
	// of sequences;

	dm->init();
	bool divSeqs = selectMostDivergentSequences();

	// Stage 2: guide tree is calculated from the distance matrix;

	guideTree = clusterBuilder->performClustering();
	if(!guideTree->isRooted()){
		// Midpoint root for NJ clustering
		guideTree = Cluster::midpointRooting(guideTree, sequences.size());
	}

	calculatesWeights(guideTree, 0.0);
	double max = getMaxWeight(guideTree, 0.0);
	normalizesWeights(guideTree, max);
	//print(guideTree);

	// Stage 3: the sequences are progressively aligned according
	// to the branching order in the guide tree;
	// The basic procedure at this stage is to use a series of pairwise
	// alignments to align larger and larger groups of sequences,
	// following the branching order in the guide tree.

	cout<<"Start of Multiple Alignment\n\nAligning..."<<endl;
	vector<Sequence*> result = alignCluster(guideTree);

	if(divSeqs){
		cout<<"\nAligning most divergent sequences"<<endl;
		result = alignMostDivergentSequences(result);
	}

	return new MultipleAlignment(result);
}

bool ClustalW::selectMostDivergentSequences(){
	bool divSeqs = false;
	for(unsigned int i=0;i<sequences.size();i++){
		double avg = 0; //numeric_limits<double>::max();
		for(unsigned int j=0;j<sequences.size();j++){
			if(i != j)
				avg += dm->identityAt(sequences[i]->id,sequences[j]->id);
		}
		avg = avg/(sequences.size()-1);
		if(avg < cutoff){
			sequences[i]->divergent = true;
			cout<<"Warning seq "<<sequences[i]->id<<" is divergent (avg identity score "<<avg<<")"<<endl;
			divSeqs = true;
		}
	}
	return divSeqs;
}

vector<Sequence*> ClustalW::alignMostDivergentSequences(vector<Sequence*> msa){

	vector<Sequence*> div = vector<Sequence*>();
	for(unsigned int i=0;i<sequences.size();i++){
		if(sequences[i]->divergent){
			div.push_back(sequences[i]);
		}
	}

	while(!div.empty()){

		int it = 0;
		vector<Sequence*> lessDiv = vector<Sequence*>();
		lessDiv.push_back(div[0]);

		double pcid = getPercentIdentity(lessDiv, msa);
		for(unsigned int i=1;i<div.size();i++){
			vector<Sequence*> tmp = vector<Sequence*>();
			tmp.push_back(div[i]);
			double currentPcid = getPercentIdentity(tmp, msa);
			if(currentPcid > pcid){
				it = i;
				pcid = currentPcid;
				lessDiv = tmp;
			}
		}
		if(msa.size() > 0){
			msa = multipleSequenceAlignment(msa, lessDiv);
		}
		else{
			msa = lessDiv;
		}

		div.erase(div.begin()+it);
	}

	return msa;
}

vector<Sequence*> ClustalW::alignCluster(Cluster* c){
	if(c->isSingle()){
		SingleCluster* sc = static_cast<SingleCluster*>(c);
		vector<Sequence*> res = vector<Sequence*>();
		if(!sc->sequence->divergent){
			res.push_back(sc->sequence);
		}
		return  res;
	}

	vector<Sequence*> leftMsa = alignCluster(c->branches[0]);
	vector<Sequence*> rightMsa = alignCluster(c->branches[1]);
	return multipleSequenceAlignment(leftMsa, rightMsa);
}

vector<Sequence*> ClustalW::multipleSequenceAlignment(const vector<Sequence*>& left, const vector<Sequence*>& right){

	//For divergent seqs
	if(left.size() == 0 && right.size() == 0){
		cout<<"Delayed"<<endl;
		return vector<Sequence*>();
	}
	if(left.size() == 0){
		cout<<"Delayed"<<endl;
		return right;
	}
	if(right.size() == 0){
		cout<<"Delayed"<<endl;
		return left;
	}

	//Avg length
	double leftLength = Sequence::GetAvgLength(left);
	double rightLength = Sequence::GetAvgLength(right);

	//Pcid
	double pcid = getPercentIdentity(left, right);

	//Matrix and scale factor
	SubMatrix* sub = NULL;
	double scale;
	getSubMatrix(pcid, sub, scale);

	//Initial gop & gep
	double initialGop = (gop + log(min(leftLength,rightLength))) * sub->getAverageMismatchScore() * scale;
	double initialGep = gep * (1.0 + std::abs(log(leftLength/rightLength)));

	CWScoringScheme* ss = new CWScoringScheme(sub, left, right);
	CWGapFunction *gf = new CWGapFunction(initialGop, initialGep, left, right);
	MsaNWG* nwmsa = new MsaNWG(ss, gf);
	std::vector<Sequence*> result = nwmsa->getMatch(left, right);

	cout<<"Sequences: "<<setw(3)<<left.size() + right.size() <<"  Score: "<<nwmsa->getScore()<<endl;
	if(verbose){
		cout<<"- aligned seqs ["<<left[0]->id;
		for(unsigned int i=1;i<left.size();i++)
			cout<<","<<left[i]->id;
		cout<<"] with ["<<right[0]->id;
		for(unsigned int i=1;i<right.size();i++)
			cout<<","<<right[i]->id;
		cout<<"]\n  initial gop "<<initialGop<<"\n  initial gep "<<initialGep<<endl;
	}

	delete nwmsa;
	delete gf;
	delete ss;

	return result;
}

double ClustalW::getPercentIdentity(const vector<Sequence*>& left, const vector<Sequence*>& right){
	double totalDistance = 0;
	for(unsigned int i=0;i<left.size();i++){
		for(unsigned int j=0;j<right.size();j++){
			totalDistance += (dm->identityAt(left[i]->id, right[j]->id));
		}
	}
	double pcid = totalDistance/(left.size()*right.size());
	return pcid*100;
}

void ClustalW::getSubMatrix(double distance, SubMatrix*& sub, double& scale){
	if(pam){
		//Pam series
		scale=0.75;
		if(distance > 80.0){
			sub = subMatrixFactory->getSubMatrix(PAM20);
		}
		else if(distance > 60.0){
			sub = subMatrixFactory->getSubMatrix(PAM60);
		}
		else if(distance > 40.0){
			sub = subMatrixFactory->getSubMatrix(PAM120);
		}
		else{
			sub = subMatrixFactory->getSubMatrix(PAM350);
		}
	}
	else{
		//Blosum series
		scale = 0.75;
		if(distance > 80.0){
			sub = subMatrixFactory->getSubMatrix(BLOSUM80);
		}
		else if(distance > 60.0){
			sub = subMatrixFactory->getSubMatrix(BLOSUM62);
		}
		else if (distance > 40.0)
		{
			sub = subMatrixFactory->getSubMatrix(BLOSUM45);
		}
		else if(distance > 30.0){
			scale = 0.5;
			sub = subMatrixFactory->getSubMatrix(BLOSUM45);
		}
		else if (distance > 20.0)
		{
			scale=0.6;
			sub = subMatrixFactory->getSubMatrix(BLOSUM45);
		}
		else
		{
			scale=0.6;
			sub = subMatrixFactory->getSubMatrix(BLOSUM30);
		}
	}
	scale *= 0.5;
}

Cluster* ClustalW::getGuideTree(){
	return guideTree;
}

}}
