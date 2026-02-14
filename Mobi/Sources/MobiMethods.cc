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
/*!
 *  \author    Luca Demo
 *  \date      2015
 */
#include <Atom.h>
#include <Spacer.h>
#include <PdbLoader.h>
#include <MobiMethods.h>
#include <MobiUtils.h>
#include <TMScoreBin.h>

using namespace Victor::Biopool;
using namespace Victor::Mobi;

double const Victor::Mobi::DEF_D0 = 4;
AtomCode const Victor::Mobi::DEF_ATOM = CA;
double const Victor::Mobi::DEF_PHI_TH = 20;
double const Victor::Mobi::DEF_PSI_TH = 20;
double const Victor::Mobi::DEF_SD_TH = .85;
double const Victor::Mobi::DEF_SDSD_TH = .09;
int const Victor::Mobi::DEF_VERBOSE = 0;

template<typename T>
/**
 * Helper method for verbose output
 */
void printVector(vector<T> const &v, string pre = ""){
	if (pre != "")
		cout << pre << endl;
	for (unsigned int i = 0; i < v.size(); i++)
		cout << v[i];
	cout << endl;
}

void MobiMethods::distances(MobiProtein* protein, TMScoreBinder* tm, VectorCollection<double>& scaledDist, VectorCollection<double>& dist, MobiMethods &mm){
	scaledDist.clear();
	dist.clear();
	MobiProtein* imposed;
	vector<double> tmpDist;
	//Perform superimposition for every model pair
	for (unsigned int i = 0; i < protein->size(); i++)
		for (unsigned int j = i+1; j < protein->size(); j++){
			imposed = tm->TMScore(*protein,i,j);
			tmpDist = scaledDistance(imposed->getModel(0),protein->getModel(j), mm.getAtom(), mm.getD0());
			scaledDist.addValue(i,j, tmpDist);
			tmpDist = distance(imposed->getModel(0),protein->getModel(j), mm.getAtom());
			dist.addValue(i,j, tmpDist);
			delete imposed;
		}
}

vector<double> MobiMethods::distance(Spacer* mod1, Spacer* mod2, AtomCode atom){
	if (mod2->size() != mod1->size())
		ERROR("Reference protein spacer has a different number of Aminos",exception);
	vector<double> dist(mod1->size());
	for (unsigned int i = 0; i < mod2->size(); i++){	//foreach Amino
		if (mod2->getAmino(i).getType() != mod1->getAmino(i).getType())
			ERROR("Aminos order in the two sequences is not the same",exception);
		Atom& mAtom = mod1->getAmino(i).getAtom(atom);
		Atom& refAtom = mod2->getAmino(i).getAtom(atom);

		dist[i] = sqrt((mAtom.getCoords() - refAtom.getCoords()).square());
	}
	return dist;
}


vector<double> MobiMethods::scaledDistance(Spacer* mod1, Spacer* mod2, AtomCode atom, double d0){
	if (mod2->size() != mod1->size())
		ERROR("Reference protein spacer has a different number of Aminos",exception);
	vector<double> sd(mod1->size());
	for (unsigned int i = 0; i < mod2->size(); i++){	//foreach Amino
		if (mod2->getAmino(i).getType() != mod1->getAmino(i).getType())
			ERROR("Aminos order in the two sequences is not the same",exception);
		Atom& mAtom = mod1->getAmino(i).getAtom(atom);
		Atom& refAtom = mod2->getAmino(i).getAtom(atom);

		double dist = sqrt((mAtom.getCoords() - refAtom.getCoords()).square());
		dist = mAtom.distance(refAtom);
		sd[i] = 1/(1 + pow(dist / d0, 2.0));
	}
	return sd;
}


void MobiMethods::phis(MobiProtein* protein, VectorCollection<double>& phis){
	phis.clear();
	vector<double> modelPhis(protein->getModel(0)->sizeAmino());
	for (unsigned int i = 0; i < protein->size(); i++){
		for (unsigned int j = 0; j < protein->getModel(i)->sizeAmino(); j++){
			AminoAcid& aa = protein->getModel(i)->getAmino(j);
			modelPhis[j] = aa.getPhi();
		}
		phis.addValue(i,modelPhis);
	}
}

void MobiMethods::psis(MobiProtein* protein, VectorCollection<double>& psis){
	psis.clear();
	vector<double> modelPsis(protein->getModel(0)->sizeAmino());
	for (unsigned int i = 0; i < protein->size(); i++){
		for (unsigned int j = 0; j < protein->getModel(i)->sizeAmino(); j++){
			AminoAcid& aa = protein->getModel(i)->getAmino(j);
			modelPsis[j] = aa.getPsi();
		}
		psis.addValue(i,modelPsis);
	}
}

vector<int> MobiMethods::secondaryMobi(MobiProtein* protein){
	//Check for errors in protein
	unsigned int const protLen = protein->size();
	if (protLen < 1)
		ERROR("Empty protein object, no models",exception);
	unsigned int const modLen = protein->getModel(0)->sizeAmino();
	if (modLen < 1)
		ERROR("Empty spacer or incorrect size",exception);
	vector<int> dssp(modLen,0);

	//Build DSSP matrix for this protein
	//1st index = Rows = models; 2nd index = Columns = Residues
	vector<vector<char> > protDSSP(protLen);
	//char protDSSP[protLen][modLen];		//not ISO C++
	for (unsigned int i = 0; i < protLen; i++){	//foreach i model
		//get DSSP for i-model
		protDSSP[i] = vector<char>(modLen,0);
		vector<set<char> > modelDSSP = protein->getModel(i)->getDSSP();
		if (modelDSSP.size() != modLen)
			ERROR("DSSP length differs from expected",exception);
		for (unsigned int j = 0; j < modLen; j++)	//foreach j amino in i-model
			if (modelDSSP[j].size() > 0)
				protDSSP[i][j] = (*(modelDSSP[j].begin()));
			else
				protDSSP[i][j] = 'C';
	}
	//Build DSSP mobility
	for (unsigned int j = 0; j < modLen; j++){	//foreach residue
		//assign Coil (2) value
		for (unsigned int i = 0; i < protLen; i++) //foreach model
			if (protDSSP[i][j] == 'C' || protDSSP[i][j] == 'S')
				dssp[j] = 2;
			else{
				dssp[j] = 0;
				break;
			}
		//assign Mobile (1) value
		if (dssp[j] != 2)	//not coil
			for (unsigned int i = 0; i < protLen; i++) //foreach model
				if (protDSSP[i][j] != protDSSP[0][j]){
					dssp[j] = 1;
					break;
				}
	}
	return dssp;
}


vector<double> MobiMethods::mobilityIndex(VectorCollection<double>& distances){
	vector<double> mi(distances.vectorsSize(),0);
	for (typename std::map<int,std::vector<double> >::const_iterator it = distances.begin(); it != distances.end(); ++it)	//foreach distance record
		for (unsigned int i = 0; i < distances.vectorsSize(); i++)	//foreach residue
			mi[i] += (pow(it->second[i],2) / distances.size());	//cumulate the 2pow of distance
	for (unsigned int i = 0; i < mi.size(); i++)
		mi[i] = sqrt(mi[i]);
	return mi;
}



vector<int> MobiMethods::mobiMobility(MobiProtein* protein, TMScoreBinder* tm){
	done = false;
	//Scaled Distances
    distances(protein,tm,scaledDist,dist,*this);

    //Set sequence
    for (unsigned int i = 0; i < protein->getModel(0)->sizeAmino(); i++)
    	sequence.push_back(protein->getModel(0)->getAmino(i).getType1L());

    //Scaled Distances Mean
    SDMeans = scaledDist.mean();
    SDMeanMobility= vector<int>(SDMeans.size());
    for (unsigned int i = 0; i < SDMeans.size(); i++)
    	SDMeanMobility[i] = SDMeans[i] < sd_th ? 1 : 0;
    //Scaled Distances deviations
    SDDevs = scaledDist.stdDev();
    SDDevsMobility = vector<int>(SDDevs.size());
	for (unsigned int i = 0; i < SDDevs.size(); i++)
		SDDevsMobility[i] = SDDevs[i] > sdsd_th ? 1 : 0;
	//Psi angles
	psis(protein,psiAngles);
	psiDev = psiAngles.stdDev();
	PsiMobility = vector<int>(psiDev.size());
	for (unsigned int i = 0; i < psiDev.size(); i++)
		PsiMobility[i] = psiDev[i] > psi_th ? 1 : 0;
	//Phi angles
	phis(protein,phiAngles);
	phiDev = phiAngles.stdDev();
	PhiMobility = vector<int>(phiDev.size());
	for (unsigned int i = 0; i < phiDev.size(); i++)
		PhiMobility[i] = phiDev[i] > phi_th ? 1 : 0;
	//Secondary mobility
	secMobility = secondaryMobi(protein);

	//Calculate Mobi Mobility
	mobiMob = SDFilters();


	done = true;
	return mobiMob;
}


vector<int> MobiMethods::SDFilters(vector<int> const &sd, vector<int> const &sdsd,
		vector<int> const &dssp, vector<int> const &phis, vector<int> const &psis, MobiMethods &mm){
	vector<int> sdm = sd;
	if (mm.verbosity() > 1)
		printVector(sdm,"Initial SD");
	//DSSP non-mobile => 1->0 in SD
	for (unsigned int i = 0; i < sdm.size(); i++)
		if (dssp[i] == 0)
			sdm[i] = 0;
	if (mm.verbosity() > 1)
		printVector(sdm,"After DSSP");

	//Patterns in SD
	//1011 -> 1111
	for (unsigned int i = 0; i < sdm.size() - 3; i++)
		if (sdm[i] == 1 && sdm[i+1] == 0 && sdm[i+2] == 1 && sdm[i+3] == 1){
			sdm[i+1] = 1;
			i+=2; //jump to next possible pattern (+3 considering ++ of for cycle)
			if (mm.verbosity() > 1)
				printVector(sdm, "1011->1111");
		}
	//1101 -> 1111
	for (unsigned int i = 0; i < sdm.size() - 3; i++)
		if (sdm[i] == 1 && sdm[i+1] == 1 && sdm[i+2] == 0 && sdm[i+3] == 1){
			sdm[i+2] = 1;
			i+=2; //jump to next possible pattern (+3 considering ++ of for cycle)
			if (mm.verbosity() > 1)
				printVector(sdm, "1101->1111");
		}
	//10011 -> 11111
	for (unsigned int i = 1; i < sdm.size() - 4; i++)
		if (sdm[i] == 1 && sdm[i+1] == 0 && sdm[i+2] == 0 && sdm[i+3] == 1 && sdm[i+4] == 1){
			sdm[i+1] = 1;
			sdm[i+2] = 1;
			i+=3; //jump to next possible pattern (+4 considering ++ of for cycle)
			if (mm.verbosity() > 1)
				printVector(sdm,"10011->11111");
		}
	//11001 -> 11111
	for (unsigned int i = 1; i < sdm.size() - 4; i++)
		if (sdm[i] == 1 && sdm[i+1] == 1 && sdm[i+2] == 0 && sdm[i+3] == 0 && sdm[i+4] == 1){
			sdm[i+2] = 1;
			sdm[i+3] = 1;
			i+=3; //jump to next possible pattern (+4 considering ++ of for cycle)
			if (mm.verbosity() > 1)
				printVector(sdm,"10011->11111");
		}
	//01010 -> 00000
	for (unsigned int i = 1; i < sdm.size() - 4; i++)
		if (sdm[i] == 0 && sdm[i+1] == 1 && sdm[i+2] == 0 && sdm[i+3] == 1 && sdm[i+4] == 0){
			sdm[i+1] = 0;
			sdm[i+3] = 0;
			i+=3; //jump to next possible pattern (+4 considering ++ of for cycle)
			if (mm.verbosity() > 1)
				printVector(sdm,"10011->11111");
		}
	//00100 -> 00000
	for (unsigned int i = 1; i < sdm.size() - 4; i++)
		if (sdm[i] == 0 && sdm[i+1] == 0 && sdm[i+2] == 1 && sdm[i+3] == 0 && sdm[i+4] == 0){
			sdm[i+2] = 0;
			i+=2; //jump to next possible pattern (+3 considering ++ of for cycle)
			if (mm.verbosity() > 1)
				printVector(sdm,"10011->11111");
		}
	//001100 -> 000000
	for (unsigned int i = 1; i < sdm.size() - 5; i++)
		if (sdm[i] == 0 && sdm[i+1] == 0 && sdm[i+2] == 1 && sdm[i+3] == 1 && sdm[i+4] == 0 && sdm[i+5] == 0){
			sdm[i+2] = 0;
			sdm[i+3] = 0;
			i+=3; //jump to next possible pattern (+4 considering ++ of for cycle)
			if (mm.verbosity() > 1)
				printVector(sdm,"10011->11111");
		}

	//Final patterns
	//110->111
	for (unsigned int i = 0; i < sdm.size() - 2; i++)
		if (sdm[i] == 1 && sdm[i+1] == 1 && sdm[i+2] == 0)
			if (phis[i+2] == 1 && psis[i+2] == 1 && sdsd[i+2] == 1 && psis[i+1] == 1){
				sdm[i+2] = 1;
				if (mm.verbosity() > 1)
					printVector(sdm, "110->111");
				i += 2; //jump to next possible pattern (+3 considering ++ of for cycle)
			}
	//011->111
	for (unsigned int i = 0; i < sdm.size() - 2; i++)
		if (sdm[i] == 0 && sdm[i+1] == 1 && sdm[i+2] == 1)
			if (phis[i] == 1 && psis[i] == 1 && sdsd[i] == 1 && phis[i+1] == 1){
				sdm[i] = 1;
				if (mm.verbosity() > 1)
					printVector(sdm, "011->111");
				i += 2; //jump to next possible pattern (+3 considering ++ of for cycle)
			}

	return sdm;
}
