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
#ifndef MOBI_SOURCES_MOBIMETHODS_H_
#define MOBI_SOURCES_MOBIMETHODS_H_

#include <Atom.h>
#include <AtomCode.h>
#include <Spacer.h>
#include <VectorCollection.h>
#include <TMScoreBin.h>

using namespace Victor::Biopool;
//using namespace Victor::Mobi;

namespace Victor{ namespace Mobi{

extern const double DEF_D0;
extern const AtomCode DEF_ATOM;
extern const double DEF_PHI_TH;
extern const double DEF_PSI_TH;
extern const double DEF_SD_TH;
extern const double DEF_SDSD_TH;

extern const int DEF_VERBOSE;

/**
 * This class provides methods to calculate mobility.\n
 * Can be instantiated to provide settings to the static methods and execute Mobi calculation.\n
 * Default constructor provides default mobi parameters.
 */
class MobiMethods{
public:

	/**
	 * Constructor with parameters.
	 * @param _d0  d0 normalization value for scaled distance
	 * @param _atom  atom to use for distance calculations
	 * @param _sd_th  Scaled distance mean threshold
	 * @param _sdsd_th  Scaled distance deviation threshold
	 * @param _psi_th  Psi angles deviation threshold
	 * @param _phi_th  Phi angles deviation threshold
	 */
	MobiMethods(double _d0 = DEF_D0, AtomCode _atom = DEF_ATOM, double _sd_th = DEF_SD_TH,
			double _sdsd_th = DEF_SDSD_TH, double _psi_th = DEF_PSI_TH, double _phi_th = DEF_PHI_TH) :
			 d0(_d0), atom(_atom), sd_th(_sd_th), sdsd_th(_sdsd_th), psi_th(_psi_th), phi_th(_phi_th), verbose(DEF_VERBOSE), done(false){}

	/**
	 * Set verbosity level
	 * @param v  verbosity to set
	 */
	void verbosity(int v){
		this->verbose = v;
	}

	/**
	 * Get verbosity level
	 * @return verbosity
	 */
	int verbosity(){
		return this->verbose;
	}

	/**
	 * Calculates the Scaled and "simple" distances given the ProteinModel object containing all the models.
	 * For all possible pair of models, Superimposition is performed and distances are calculated.
	 * @param protein  pointer to protein object containing all the models
	 * @param tm pointer to TMScore binder
	 * @param scaledDist destination of scaled distances
	 * @param dist destination of distances
	 * @param mm reference to MobiMethod object for mobi parameters
	 */
	static void distances(MobiProtein* protein, TMScoreBinder* tm, VectorCollection<double>& scaledDist, VectorCollection<double>& dist, MobiMethods &mm);


	/**
	 * Given a protein containing models, populate a VectorCollection of phi angles
	 * @param protein pointer to the protein to elaborate
	 * @param phis output VectorCollection
	 */
	static void phis(MobiProtein* protein, VectorCollection<double>& phis);

	/**
	 * Given a protein containing models, populate a VectorCollection of psi angles
	 * @param protein pointer to the protein to elaborate
	 * @param psis output VectorCollection
	 */
	static void psis(MobiProtein* protein, VectorCollection<double>& psis);

	/**
	 * Given a protein containing models, build DSSP mobility estimation.\n The returning
	 * vector contains values according to Mobi Methods guidelines:\n
	 * 0 = non mobile\n
	 * 1 = mobile\n
	 * 2 = coil
	 * @param protein reference to the protein
	 * @return a vector<int> containing the secondary mobility estimation (0,1,2 values)
	 */
	static vector<int> secondaryMobi(MobiProtein* protein);


	/**
	 * Returns mobility indexes per residue. Static method
	 * @param distances vectorCollection ("simple" distances)
	 * @return vector<double> with mobility indexes
	 */
	static vector<double> mobilityIndex(VectorCollection<double>& distances);

	/**
	 * Return mobility indexes per residue
	 * @return vector<double> with mobility indexes
	 */
	vector<double> mobilityIndex(){
		return MobiMethods::mobilityIndex(dist);
	}

	/**
	 * Performs mobility calculation
	 * @param protein pointer to the protein to process
	 * @param tm pointer to TM imposer
	 * @return a vector<int> of mobility values (0,1)
	 */
	vector<int> mobiMobility(MobiProtein* protein, TMScoreBinder* tm);

	/**
	 * Applies Mobi filters to Scaled Distance mobility track to obtain final mobility values
	 * @param sdm vector<int> of scaled distance mobility track
	 * @param sdsd vector<int> scaled distance std dev mobility track
	 * @param sec vector<int> secondary structure mobility track
	 * @param phis vector<int> phi angles mobility track
	 * @param psis vector<int> psi angles mobility track
	 * @param settings reference to MobiMethods object for mobi parameters
	 * @return vector<int> of filtered mobility ("all" mobility)
	 */
	static vector<int> SDFilters(vector<int> const &sdm, vector<int> const &sdsd,
			vector<int> const &sec, vector<int> const &phis, vector<int> const &psis, MobiMethods &settings);

	/**
	 * Get D0 value
	 */
	double getD0(){	return this->d0;}

	/**
	 * Get Atom used for distance calculation
	 */
	AtomCode getAtom(){return this->atom;}

	/**
	 * Get Average scale distance threshold
	 */
	double getSDTh(){return this->sd_th;}

	/**
	 * get Scaled distance deviance threshold
	 */
	double getSDSDTh(){return this->sdsd_th;}

	/**
	 * Get Phi angles deviation threshold
	 */
	double getPhiTh(){return this->phi_th;}

	/**
	 * Get Psi angles deviation threshold
	 */
	double getPsiTh(){return this->psi_th;}

	/**
	 * get average scaled distance
	 */
	const vector<double> getSDMeans(){return SDMeans;}

	/**
	 * get average scaled distance mobility track
	 */
	const vector<int> getSDMeanMobility(){return SDMeanMobility;}

	/**
	 * get scaled distance deviation vector
	 */
	const vector<double> getSDDevs(){return SDDevs;}

	/**
	 * get scaled distance deviation mobility track
	 */
	const vector<int> getSDDevsMobility(){return SDDevsMobility;}

	/**
	 * get psi angles deviations
	 */
	const vector<double> getPsiDevs(){return psiDev;}

	/**
	 * get psi angles mobility track
	 */
	const vector<int> getPsiMobility(){return PsiMobility;}

	/**
	 * get phi angles deviations
	 */
	const vector<double> getPhiDevs(){return phiDev;}

	/**
	 * get phi mobility track
	 */
	const vector<int> getPhiMobility(){return PhiMobility;}

	/**
	 * get secondary structure mobility
	 */
	const vector<int> getSecMobility(){return secMobility;}

	/**
	 * get mobi mobility "all" track
	 */
	const vector<int> getMobiMobility(){return mobiMob;}

	/**
	 * get sequence
	 */
	const vector<char> getSequence(){return sequence;}

	/**
	 * get distance collection
	 */
	const VectorCollection<double>& getDistances(){return dist;}

	/**
	 * get scaled distance collection
	 */
	const VectorCollection<double>& getScaledDistances(){return scaledDist;}


	/**
	 * get phi angles
	 */
	const VectorCollection<double>& getPhiAngles(){return phiAngles;}

	/**
	 * get psi angles
	 */
	const VectorCollection<double>& getPsiAngles(){return psiAngles;}

	/**
	 * Return true if mobi calculation has been performed
	 * @return true if mobi mobility is calculated, false otherwise
	 */
	bool isDone(){
		return done;
	}

protected:


	/**
	 * Calculate the "simple" distance given two models.
	 * @param mod1 the first model (Spacer pointer)
	 * @param mod2 the second model (Spacer pointer)
	 * @param atom (AtomCode) atom to base distance calcutation on, default atom if not specified (CA)
	 */
	static vector<double> distance(Spacer* mod1, Spacer* mod2, AtomCode atom = DEF_ATOM);


	/**
	 * Calculate the scaled distance given two models
	 * @param mod1 the first model (Spacer pointer)
	 * @param mod2 the second model (Spacer pointer)
	 * @param atom atom to base distance calculation on, default atom if not specified (CA)
	 * @param d0  d0 normalization value, default value if not specified (4)
	 */
	static vector<double> scaledDistance(Spacer* mod1, Spacer* mod2, AtomCode atom = DEF_ATOM, double d0 = DEF_D0);

private:

	/**
	 * Applies Mobi filters to Scaled Distance mobility track to obtain final mobility values.
	 * Private method as shortcut in mobiMobility method.
	 * @return vector<int> of filtered mobility ("all" mobility)
	 */
	vector<int> SDFilters(){
		return SDFilters(getSDMeanMobility(), getSDDevsMobility(), getSecMobility(), getPhiMobility(),
				getPsiMobility(), *this);
	}




	double d0;
	double sdsd_th;
	double sd_th;
	double psi_th;
	double phi_th;
	AtomCode atom;
	int verbose;

	vector<char> sequence;
	vector<double> SDMeans;
	vector<int> SDMeanMobility;
	vector<double> SDDevs;
	vector<int> SDDevsMobility;
	vector<double> psiDev;
	vector<int> PsiMobility;
	vector<double> phiDev;
	vector<int> PhiMobility;
	vector<int> secMobility;
	vector<int> mobiMob;

	VectorCollection<double> scaledDist;	//Scaled distances
	VectorCollection<double> dist;			//"Simple" distances
	VectorCollection<double> phiAngles;
	VectorCollection<double> psiAngles;

	bool done;
};

}}	//Namespaces


#endif /* MOBI_SOURCES_MOBIMETHODS_H_ */
