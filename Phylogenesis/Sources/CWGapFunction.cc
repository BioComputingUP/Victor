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

#include <CWGapFunction.h>
#include <iostream>
#include <iomanip>

namespace Victor { namespace Phylogenesis {

// Utility functions

std::map<char, double> init(){
	std::map<char, double> m;
    m['A'] = 1.13;
    m['C'] = 1.13;
    m['D'] = 0.96;
    m['E'] = 1.31;
    m['F'] = 1.20;
    m['G'] = 0.61;
    m['H'] = 1.00;
    m['I'] = 1.32;
    m['K'] = 0.96;
    m['L'] = 1.21;
    m['M'] = 1.29;
    m['N'] = 0.63;
    m['P'] = 0.74;
    m['Q'] = 1.07;
    m['R'] = 0.72;
    m['S'] = 0.76;
    m['T'] = 0.89;
    m['V'] = 1.25;
    m['Y'] = 1.00;
    m['W'] = 1.23;
    return m;
}

static std::map<char, double> PascarellaAndArgosSpecificGapModificationFactors = init();

int countGaps(const std::vector<Sequence*>& seqs, int pos){
	int tot = 0;
	for(unsigned int i = 0; i<seqs.size();i++){
		if(seqs[i]->isGap(pos)){
			tot++;
		}
	}
	return tot;
}

CWGapFunction::CWGapFunction(double initialGop, double initiapGep, const std::vector<Sequence*>& vertical, const std::vector<Sequence*>& horizontal) :
		initialGop(initialGop), initialGep(initiapGep) {
	verticalGops.reserve(vertical.size());
	verticalGeps.reserve(vertical.size());
	horizontalGops.reserve(horizontal.size());
	horizontalGeps.reserve(horizontal.size());
	initGapValues(vertical, verticalGops, verticalGeps);
	initGapValues(horizontal, horizontalGops, horizontalGeps);
}

CWGapFunction::~CWGapFunction() {}

void CWGapFunction::initGapValues(const std::vector<Sequence*>& seqs, std::vector<double>& gops, std::vector<double>& geps){

	//Init position-specific gap penalties
	unsigned int totalPos = seqs[0]->sequence.size();

	for(unsigned int POS=0; POS<totalPos; POS++){

		int withGap = 0;
		int withoutGap = 0;
		for(unsigned int i=0;i<seqs.size();i++){
			if(seqs[i]->isGap(POS)){
				withGap++;
			}
			else{
				withoutGap++;
			}
		}

		//1) Lowered gap penalties at existing gaps

		geps.push_back(withGap > 0 ? (initialGep / 2.0) : initialGep);
		if(withGap>0){
			gops.push_back(initialGop*0.3*((double)withoutGap/seqs.size()));
			continue;
		}

		//2) If a position does not have any gaps but is within 8 residues of an existing gap, the GOP is increased
		// 	 by GOP -> GOP*(2+((8-distance from gap)*2)/8)

		bool closeGap = false;
		int distance = 0;
		for(unsigned int i=1;i<=8 && !closeGap ;i++){
			distance++;
			if(countGaps(seqs,POS+i) > 0 || countGaps(seqs,POS-i) > 0){
				closeGap = true;
			}
		}
		if(closeGap){
			gops.push_back(initialGop*(2.0+((8.0-(double)distance)*2.0)/8.0));
			continue;
		}

		//3) Reduced gap penalties in hydrophilic stretches.
		//   Any run of 5 hydrophilic residues is considered to be a hydrophilic stretch.

		bool foundHydrophilicStretch = false;
		for(unsigned int i =0; i< seqs.size() && !foundHydrophilicStretch ;i++){
			for(int c = (int)POS-4; c <= (int)POS && !foundHydrophilicStretch; c++){
				foundHydrophilicStretch = seqs[i]->isHydrophilicStretch(c);
			}
		}
		if(foundHydrophilicStretch){
			gops.push_back(initialGop* 2.0/3.0);
			continue;
		}

		//4) Residue specific penalties
		//   If there is no hydrophilic stretch and the position does not contain any gaps, then the GOP is multiplied by
		//	 Pascarella and Argos modification factors.  If there is a mixture of residues at a position, the multiplication
		//   factor is the average of all the contributions from each sequence.

		double avg = 0;
		for(unsigned int i =0; i< seqs.size() ;i++){
			avg += PascarellaAndArgosSpecificGapModificationFactors[seqs[i]->sequence[POS]];
		}
		avg = avg / seqs.size();
		gops.push_back(initialGop * avg);
	}

	#ifdef DEBUG
	for(unsigned int POS=0; POS<totalPos; POS++){
		std::cout<<POS<<"\t";
	}
	std::cout<<std::endl;
	for(unsigned int i =0; i< seqs.size(); i++){
		std::string seq = seqs[i]->sequence;
		for(unsigned int j=0;j<seq.size();j++){
			std::cout<<seq[j]<<"\t";
		}
		std::cout<<std::endl;
	}
	for(unsigned int POS=0; POS<totalPos; POS++){
		std::cout<<std::setprecision(2)<<std::fixed <<gops[POS]<<"\t";
	}
	std::cout<<std::endl;
	for(unsigned int POS=0; POS<totalPos; POS++){
		std::cout<<std::setprecision(2)<<std::fixed <<geps[POS]<<"\t";
	}
	std::cout<<std::endl<<std::endl;
	#endif
}

double CWGapFunction::getHorizontalOpenPenalty(int pos){
	return horizontalGops[pos-1];
}

double CWGapFunction::getHorizontalExtensionPenalty(int pos){
	return horizontalGeps[pos-1];
}

double CWGapFunction::getVerticalOpenPenalty(int pos){
	return verticalGops[pos-1];
}

double CWGapFunction::getVerticalExtensionPenalty(int pos){
	return verticalGeps[pos-1];
}

}}
