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

#include <vector>
#include <Debug.h>
#include <Spacer.h>
#include <MobiMethods.h>
#include <MobiProtein.h>
#include <PdbSaver.h>
#include <vector3.h>


#ifndef MOBI_SOURCES_UTILS_H_
#define MOBI_SOURCES_UTILS_H_

namespace Victor{ namespace Mobi{

/**
 *	@brief Provides utilities for Mobi operations
 */
class MobiUtils{


public:

	/**
	 * Given a vector of numeric values, the function order it. It also switches values in another vector according to
	 * order operations on the first vector.
	 * @param vals reference to the first vector, to be ordered
	 * @param ord reference to the second vector
	 * @param left left bound to order
	 * @param right right bound to order
	 */
	static void sortVector(vector<double> &vals, vector<int> &ord, int left, int right) {
		int i = left, j = right;
		double tmpVal;
		int tmpOrd;
		double pivot = vals[(left + right) / 2];
		/* partition */
		while (i <= j) {
			while (vals[i] < pivot)
				i++;
			while (vals[j] > pivot)
				j--;
			if (i <= j) {
				tmpVal = vals[i];
				tmpOrd = ord[i];
				vals[i] = vals[j];
				ord[i] = ord[j];
				vals[j] = tmpVal;
				ord[j] = tmpOrd;
				i++;
				j--;
			}
		};
		/* recursion */
		if (left < j)
			sortVector(vals, ord, left, j);
		if (i < right)
			sortVector(vals, ord, i, right);
	}


	/**
	 * Sort models from the most representative to the least one.\n
	 * Models are searched using in the protein provided by applying superimposition and calculatin all the RMSD between model pairs.
	 * The most representative model has the lowest average RMSD.
	 * @param protein pointer to protein
	 * @param tm pointer to the TM imposer object
	 * @return int vector containing the ordered model IDs
	 */
	static vector<int> sortModels(MobiProtein* protein, TMScoreBinder* tm){
		vector<MobiProtein*> models(protein->size());
		vector<double> distance;

		if (protein->size() < 2)
			ERROR("Cannot calculate average model of a protein with less of 2 models!",exception);

		VectorCollection<double> distances;
		VectorCollection<double> scDist;
		MobiMethods mm;
		MobiMethods::distances(protein, tm, scDist, distances, mm);

		//Calculate RMSDs
		vector<int> modNames(models.size());
		vector<double> modRmsd(models.size(),0.0);


		for (unsigned int i=0; i < models.size(); i++){
			VectorCollection<double> d2 = distances.getValuesByModel(i);
			modNames[i]=i;
			modRmsd[i] = d2.RMSD();
		}

		//Order the models
		sortVector(modRmsd,modNames,0,modNames.size()-1);


//		for (unsigned int i = 0; i < modNames.size(); i++)
//			cout << modNames[i]+1 << " => " << modRmsd[i] << endl;
		return modNames;
	}

	/**
	 * Sort models from the most representative to the least one.\n
	 * Models are searched using the provided VectorCollection of distances between the models.
	 * RMSD is calculated and the most representative model has the lowest average RMSD.
	 * @param distances distances
	 * @return vector containing the ordered model IDs
	 */
	static vector<int> sortModels(const VectorCollection<double>& distances){
		vector<int> models = distances.getModels();
		//Calculate RMSDs
		vector<int> modNames;
		vector<double> modRmsd;

		for (unsigned int i=0; i < models.size(); i++){
			VectorCollection<double> d2 = distances.getValuesByModel(models[i]);
			modNames.push_back(models[i]);
			modRmsd.push_back(d2.RMSD());
		}

		//Order the models
		sortVector(modRmsd,modNames,0,modNames.size()-1);
		return modNames;
	}

	/**
	 * Prints the fasta output on the specified stream. If output stream is cout, then print additional informations.
	 * @param mobi reference to MobiMethods with mobility tracks
	 * @param output reference to output stream, default is cout
	 */
	static void printFastaFile(MobiMethods& mobi, ostream& output = cout){
		if (!mobi.isDone())
			ERROR("Cannot save score of a non-done MobiMethod object",exception);
		if (output == cout){	//if output is cout, add more informations
			output << "Victor Mobi Mobility results." << endl;
			output << "d0 scaling value \t\td0 = " << mobi.getD0() << endl;
			output << "Avg Scaled Distance \t\tTh = " << mobi.getSDTh() << endl;
			output << "Standard deviation distance \tTh = " << mobi.getSDSDTh() << endl;
			output << "Phi Angle \t\t\tTh = " << mobi.getPhiTh() << endl;
			output << "Psi Angle \t\t\tTh = " << mobi.getPsiTh() << endl;
		}

		output << "> Sequence" << endl;
		for (unsigned int i = 0; i < mobi.getSequence().size(); i++)
			output << mobi.getSequence()[i];
		output << endl << "> Average scale distance 0" << endl;
		for (unsigned int i = 0; i < mobi.getSDMeanMobility().size(); i++)
			output << (mobi.getSDMeanMobility()[i] == 1 ? "M" : ".");
		output << endl << "> Standard deviation distance 1" << endl;
		for (unsigned int i = 0; i < mobi.getSDDevsMobility().size(); i++)
			output << (mobi.getSDDevsMobility()[i] == 1 ? "M" : ".");
		output << endl << "> Phi angle 2" << endl;
		for (unsigned int i = 0; i < mobi.getPhiMobility().size(); i++)
			output << (mobi.getPhiMobility()[i] == 1 ? "M" : ".");
		output << endl << "> Psi angle 3" << endl;
		for (unsigned int i = 0; i < mobi.getPsiMobility().size(); i++)
			output << (mobi.getPsiMobility()[i] == 1 ? "M" : ".");
		output << endl << "> Secondary structure 4" << endl;
		for (unsigned int i = 0; i < mobi.getSecMobility().size(); i++)
			output << (mobi.getSecMobility()[i] == 1 ? "M" : (mobi.getSecMobility()[i] == 2 ? "c" : "."));
		output << endl << "> All 5"<< endl;
		for (unsigned int i = 0; i < mobi.getMobiMobility().size(); i++)
			output << (mobi.getMobiMobility()[i] == 1 ? "M" : ".");
	}

	/**
	 * Saves a spacer in PDB format. Also writes MODEL and ENDMDL.
	 * Victor PdbSaver does not offer enough functionalities to do this, so this method exists...
	 * "Write residues" and "write secondary" options are ignored.
	 *@param prot pointer to protein containing models
	 *@param order vector<int> ordered models IDs
	 *@param avgDist average scaled distance vector for this protein's models
	 *@param distDev scaled distance standard deviation vector for this protein's models
	 *@param output output stream
	 *@return void
	 */
	static void saveMobiPdb(MobiProtein* prot, vector<int> order, const vector<double>& avgDist, const vector<double>& distDev, ostream& output = cout){
		PdbSaver pdbSaver(output);
		pdbSaver.setWriteAtomOnly();

		output << "REMARK   1 Victor Mobi Mobility PDB output\n";
		output << "REMARK   1 Occupancy field has been replaced by scaled distance deviations\n";
		output << "REMARK   1 T-Factor field has been replaced by average scaled distance\n";
		for (unsigned int m = 0; m < order.size(); m++){	//Foreach model
			output << "MODEL     " << prot->getModelsID()[order[m]] << "\n";
			int aminoOffset = 0; //prot->getModel(m)->getStartOffset();
			for (unsigned int r = 0; r < prot->getModel(order[m])->sizeAmino(); r++){	//Foreach residue
				//Amino offset
				aminoOffset++;
				while ((prot->getModel(order[m])->isGap(aminoOffset)) && (aminoOffset < prot->getModel(order[m])->maxPdbNumber()))
					aminoOffset++;
				AminoAcid& aa = prot->getModel(order[m])->getAmino(r);

				for (unsigned int a = 0; a < aa.size(); a++){	//Foreach Atom in AminoAcid
					Atom& atom = aa[a];
					string aName = aa[a].getType();
					if (!isdigit(aName[0]) && (aName.size() < 4))	//Fix atom name string length to 4 chars
						aName = ' ' + aName;
					while (aName.size() < 4)
						aName += ' ';
					//write
					output << "ATOM" << setw(7) << atom.getNumber() << " " << aName	<< " "
						<< aa.getType() << " " << prot->getChainLetter(0) << setw(4) << aminoOffset << "    "
						<< setw(8) << fixed << setprecision(3) << atom.getCoords().x
						<< setw(8) << fixed << setprecision(3) << atom.getCoords().y
						<< setw(8) << fixed << setprecision(3) << atom.getCoords().z
						<< setw(6) << fixed << setprecision(2) << avgDist[a]
						<< setw(6) << fixed << setprecision(2) << distDev[a]*100
						<< "           " <<  (isdigit(aName[1]) ? aName[2] : aName[1]) << "\n";
				}

			}
			output << "ENDMDL" << "\n";
		}
	}


	/**
	 * Saves mobility scores.
	 * @param score MobiMethod object with calculated mobility
	 * @param chain char value for chain name
	 * @param output output stream, default is cout
	 */
	static void saveScore(MobiMethods& score, char chain, ostream& output = cout){
		if (!score.isDone())
			ERROR("Cannot save score of a non-done MobiMethod object",exception);
		vector<char> seq = score.getSequence();
		vector<double> avg = score.getSDMeans();
		vector<double> dev = score.getSDDevs();
		vector<double> phi = score.getPhiAngles().mean();
		vector<double> psi = score.getPsiAngles().mean();
		vector<double> mi = score.mobilityIndex();
		vector<int> mobi = score.getMobiMobility();
		if (output == cout)
			output << "Mobi Mobility Score" << endl;
		for (unsigned int i = 0; i < seq.size(); i++){
			output << i+1 << "\t" << chain << "\t" << seq[i] << "\t" << (mobi[i] == 1 ? "M" : "-") << "\t"
			<< setw(6) << fixed << setprecision(2) << avg[i] << "\t"
			<< setw(6) << fixed << setprecision(2) << dev[i] << "\t"
			<< setw(6) << fixed << setprecision(2) << phi[i] << "\t"
			<< setw(6) << fixed << setprecision(2) << psi[i] << "\t"
			<< setw(6) << fixed << setprecision(2) << mi[i] << endl;
		}
	}
};
}}
#endif
