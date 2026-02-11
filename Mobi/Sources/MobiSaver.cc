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

/**
 * @file MobiSaver.h
 * @author Riccardo Zanella
 * @date 21 lug 2015
 * @version 0.1
 */

//Includes:
#include "MobiSaver.h"

using namespace Victor;
using namespace Victor::Biopool;
using namespace Victor::Mobi;

const string OUT = "ResultMobi.fasta";

/**
 *   Default Constructor
 * 	 @param  ProteinModels& , object PreteinModels
 * 	 @param  string output, path file output
 * 	 @param  bool verbose , verbose
 * 	 @param double boundSD, default = 0.85
 * 	 @param double boundStandD, default = 0.09
 * 	 @param double anglePHI, default = 20
 * 	 @param double anglePSI, default = 20
 */
MobiSaver::MobiSaver(ProteinModels prot, string output, ofstream& stream,
		bool verbose, double boundSD, double boundStandD, double anglePHI,
		double anglePSI) :
		PdbSaver(stream), verbose(verbose), ScalD(boundSD), StandD(boundStandD), angPHI(
				anglePHI), angPSI(anglePSI) {
	// TODO Auto-generated constructor stub

	//if file does not exist then create and insert aminoacid sequence in head otherwise don't create file
	out = output + OUT;

	if (access(out.c_str(), F_OK) != 0) {

		if (verbose)
			cout << "\n ### Create fasta file in output ###" << endl;

		ofstream fout;			//stream in output
		vector<char>::iterator walk = prot.sequence.begin();

		fout.open(out.c_str());
		if (!fout) {
			ERROR("Error to create file!", error);
		} else {
			fout << "> sequence" << endl;
			while (walk != prot.sequence.end()) {
				fout << *(walk);
				walk++;
			}

			fout.close();

			if (verbose)
				cout << "file success created!" << endl;
		}

	}

}

/**
 *   Default Deconstructor
 *
 */
MobiSaver::~MobiSaver() {
	// TODO Auto-generated destructor stub
}

// PREDICATES:

/**
 * save mobility from average scale distance vector
 * @param vector <double>, vector of average scale distance
 * @return void
 */
void MobiSaver::mob_eveScalD(vector<double> everageDistance) {

	if (everageDistance.size() == 0)
		ERROR("Vector everage distance empty!", exception);

	ofstream fout(out.c_str(), ofstream::app);

	if (!fout)
		ERROR("Could not open file for writing.", error);

	if (verbose)
		cout << "\n ### Save mobility from average distance ###" << endl;

	fout << endl;
	fout << "> Average scale distance 0" << endl;

	for (vector<double>::iterator walk = everageDistance.begin();
			walk != everageDistance.end(); walk++) {
		if (*(walk) < ScalD)
			fout << "M";
		else
			fout << ".";
	}

	fout.close();

	if (verbose)
		cout << "Mobility from average scaled distance saved!" << endl;

}

/**
 * save mobility from standard deviation distance vector
 * @param vector <double>, vector of scale distance
 * @return void
 */
void MobiSaver::mob_stanD(vector<double> Scale_distance) {

	if (Scale_distance.size() == 0)
		ERROR("Vector standard deviation empty!", exception);

	ofstream fout(out.c_str(), ofstream::app);

	if (!fout)
		ERROR("Could not open file for writing.", error);

	if (verbose)
		cout << "\n ### Mobility from standard deviation saved on file ###"
				<< endl;

	fout << endl;
	fout << "> Standard deviation distance 1" << endl;

	for (vector<double>::iterator walk = Scale_distance.begin();
			walk != Scale_distance.end(); walk++) {
		if (*(walk) > StandD)
			fout << "M";
		else
			fout << ".";
	}

	fout.close();

	if (verbose)
		cout << "Mobility from standard deviation saved!" << endl;

}

/**
 * save mobility from vector of phi angle
 * @param vector <double>, vector of scale distance
 * @return void
 */
void MobiSaver::mob_aPHI(vector<double> angle_PHI) {

	if (angle_PHI.size() == 0)
		ERROR("Vector angle PHI empty!", exception);

	ofstream fout(out.c_str(), ofstream::app);

	if (!fout)
		ERROR("Could not open file for writing.", error);

	if (verbose)
		cout << "\n ### Mobility from angle phi saved on file ###" << endl;

	fout << endl;
	fout << "> Phi angle 2" << endl;

	for (vector<double>::iterator walk = angle_PHI.begin();
			walk != angle_PHI.end(); walk++) {
		if (*(walk) > angPHI)
			fout << "M";
		else
			fout << ".";
	}

	fout.close();

	if (verbose)
		cout << "Mobility from angle_PHI saved!" << endl;

}

/**
 * save mobility from vector of psi angle
 * @param vector <double>,  from vector of psi angle
 * @return void
 */
void MobiSaver::mob_aPSI(vector<double> angle_PSI) {

	if (angle_PSI.size() == 0)
		ERROR("Vector angle PSI empty!", exception);

	ofstream fout(out.c_str(), ofstream::app);

	if (!fout)
		ERROR("Could not open file for writing.", error);

	if (verbose)
		cout << "\n ### Mobility from angle phi saved on file ###" << endl;

	fout << endl;
	fout << "> Psi angle 3" << endl;

	for (vector<double>::iterator walk = angle_PSI.begin();
			walk != angle_PSI.end(); walk++) {
		if (*(walk) > angPSI)
			fout << "M";
		else
			fout << ".";
	}

	fout.close();

	if (verbose)
		cout << "Mobility from angle_PSI saved!" << endl;

}

/**
 * save mobility structure
 * @param vector<char>,  from vector of mobility structure
 * @return void
 */
void MobiSaver::mob_SecS(vector<char> Mob_SecStructure) {

	if (Mob_SecStructure.size() == 0)
		ERROR("Vector Mob_SecStructure empty!", exception);

	ofstream fout(out.c_str(), ofstream::app);

	if (!fout)
		ERROR("Could not open file for writing.", error);

	if (verbose)
		cout << "\n ### Mobility from secondary structures saved on file ###"
				<< endl;

	fout << endl;
	fout << "> Secondary structure 4" << endl;

	for (vector<char>::iterator walk = Mob_SecStructure.begin();
			walk != Mob_SecStructure.end(); walk++)
		fout << *(walk);

	fout.close();

	if (verbose)
		cout << "Mobility from secondary structure saved!" << endl;

}

/**
 * save mobility from average scale distance filtered by secondary structures
 * @param vector <double>,  vector average scaled distance
 * @param vector <char>,  vector secondary structures
 * @return void
 */
void MobiSaver::mob_eveScalD_filtSecS(vector<double> everageDistance,
		vector<char> mobSecS) {

	if (everageDistance.size() == 0 || mobSecS.size() == 0)
		ERROR("Exception: Vectors must be full.", exception);

	ofstream fout(out.c_str(), ofstream::app);

	if (!fout)
		ERROR("Could not open file for writing.", error);

	if (verbose)
		cout
				<< "\n ### Save mobility from average Scale distance filtered by Secondary Structure...  ###"
				<< endl;

	fout << endl;
	fout << "> Average scale distance filtered with SecStr 0" << endl;

	vector<double>::iterator everage = everageDistance.begin();
	vector<char>::iterator secStruct = mobSecS.begin();

	//If average is fix stay fix, if average is mobile and DSSP is mobile then change
	//to fix
	while (everage != everageDistance.end() || secStruct != mobSecS.end()) {

		if (*(everage) < ScalD && *(secStruct) != '.')
			fout << "M";
		else
			fout << ".";

		everage++;
		secStruct++;
	}

	fout.close();

	if (verbose)
		cout
				<< "Mobility from average scaled  distance filtered by Secondary Structure saved!"
				<< endl;

}

/**
 * save mobility from average scale distance filtered by angles and scaled distance
 * @param vector <double>,  vector average scaled distance
 * @param vector <double>, angle phi
 * @param vector <double>, angle psi
 * @param vector <double>, Scaled Distance
 * @return void
 */
void MobiSaver::mob_eveScalD_filteredByPHI_PSI_standD(
		vector<double> everageDistance, vector<double> angle_PHI,
		vector<double> angle_PSI, vector<double> Scale_distance) {

	vector<string> vettoreSupporto;

	if (everageDistance.size() == 0 || angle_PHI.size() == 0
			|| angle_PSI.size() == 0 || Scale_distance.size() == 0)
		ERROR("One vectors is empty!", exception);

	ofstream fout(out.c_str(), ofstream::app);

	if (!fout)
		ERROR("Could not open file for writing.", error);

	if (verbose)
		cout
				<< "\n ### Save mobility from average Scale distance filtered by PHI PSI Standard deviation...  ###"
				<< endl;

	fout << endl;
	fout << "> Average scale distance filtered with PHI PSI StdD 0" << endl;

	const string filter1 = "MM.";
	const string filter2 = ".MM";

	for (vector<double>::iterator walk = everageDistance.begin();
			walk != everageDistance.end(); walk++) {
		if (*(walk) < ScalD)
			vettoreSupporto.push_back("M");
		else
			vettoreSupporto.push_back(".");
	}

	string test;

	vector<string>::iterator mob = vettoreSupporto.begin();
	vector<double>::iterator phi = angle_PHI.begin();
	vector<double>::iterator psi = angle_PSI.begin();
	vector<double>::iterator stdD = Scale_distance.begin();

	int i;
	bool flag;
	vector<string>::iterator j;
	vector<double>::iterator h, s, d;

	while (mob != vettoreSupporto.end()) {

		i = 0;
		flag = true;
		test.clear();
		j = mob;
		h = phi;
		s = psi;
		d = stdD;

		while (j != vettoreSupporto.end() && i < 3 && flag) {

			test = test + (*(j));

			if (test.compare(filter1) == 0) {
				if (*(h) > angPHI && *(s) > angPSI && *(d) > StandD
						&& (*(s - 1) > angPSI))
					flag = false;

			} else if (test.compare(filter2) == 0) {

				if (*(h - 2) > angPHI && *(s - 2) > angPSI && *(d - 2) > StandD
						&& (*(h - 1) > angPHI))
					flag = false;

			}

			j++;
			h++;
			s++;
			d++;
			i++;
		}

		if (!flag) {

			fout << "MMM";
			mob = j;
			phi = h;
			psi = s;
			stdD = d;

		} else {
			fout << *(mob);
			mob++;
		}

	}

	fout.close();

	if (verbose)
		cout
				<< "Mobility from average Scale distance filtered by PHI PSI Standard deviation whith mask saved!"
				<< endl;

}

/**
 * save mobility from standard deviation with mask
 * @param vector <double>,  vector scaled distance
 * @return void
 */
void MobiSaver::mob_stanD_withMask(vector<double> Scale_distance) {

	vector<string> vettoreSupporto;
	map<string, string> map;

	if (Scale_distance.size() == 0)
		ERROR("Vector standard deviation empty!", exception);

	ofstream fout(out.c_str(), ofstream::app);

	if (!fout)
		ERROR("Could not open file for writing.", error);

	if (verbose)
		cout << "\n ### Save on file standard deviation with Mask ###" << endl;

	fout << endl;
	fout << "> Standard deviation distance with mask 1" << endl;

	const pair<string, string> product1("M.MM", "MMMM");
	const pair<string, string> product2("MM.M", "MMMM");
	const pair<string, string> product3("M..MM", "MMMMM");
	const pair<string, string> product4("MM..M", "MMMMM");
	const pair<string, string> product5(".M.M.", ".....");
	const pair<string, string> product6("..M..", ".....");
	const pair<string, string> product7("..MM..", "......");

	map.insert(product1);
	map.insert(product2);
	map.insert(product3);
	map.insert(product4);
	map.insert(product5);
	map.insert(product6);
	map.insert(product7);

	for (vector<double>::iterator walk = Scale_distance.begin();
			walk != Scale_distance.end(); walk++) {
		if (*(walk) > StandD)
			vettoreSupporto.push_back("M");
		else
			vettoreSupporto.push_back(".");
	}

	string test;

	vector<string>::iterator mob = vettoreSupporto.begin();
	std::map<string, string>::iterator iter = map.begin();

	int i;
	bool flag;
	vector<string>::iterator j;

	while (mob != vettoreSupporto.end()) {

		i = 0;
		flag = true;
		test.clear();
		j = mob;

		while (j != vettoreSupporto.end() && i < 6 && flag) {

			test = test + (*(j));
			iter = map.find(test);
			if (iter != map.end()) {
				test = iter->second;
				flag = false;
			}

			j++;
			i++;
		}

		if (!flag) {
			fout << test;
			mob = j;
		} else {
			fout << *(mob);
			mob++;
		}

	}

	fout.close();

	if (verbose)
		cout << "Mobility from standard deviation whith mask saved!" << endl;

}

/**
 * save all type of mobility
 * @param vector <double>,  vector average scaled distance
 * @param vector <double>,  vector Scaled Distance
 * @param vector <double>, angle phi
 * @param vector <double>, angle psi
 * @param vector <char> , secondary structures
 *
 * @return void
 */
void MobiSaver::allMobility(vector<double> everageDistance,
		vector<double> Scale_distance, vector<double> angle_PHI,
		vector<double> angle_PSI, vector<char> Mob_SecStructure) {

	mob_eveScalD(everageDistance);
	mob_eveScalD_filtSecS(everageDistance, Mob_SecStructure);
	mob_eveScalD_filteredByPHI_PSI_standD(everageDistance, angle_PHI, angle_PSI,
			Scale_distance);
	mob_stanD(Scale_distance);
	mob_stanD_withMask(Scale_distance);
	mob_aPHI(angle_PHI);
	mob_aPSI(angle_PSI);
	mob_SecS(Mob_SecStructure);

}

/**
 * Method that shadow save Protein of PdbSaver, this method save in .pdb format
 * all  original models in ProteinModels with column standard deviation modified
 * and column B factor modified with average scaled distance.
 * @param ProteinModels , object ProteinModels
 * @param vector <double>,  vector average scaled distance
 * @param vector <double>,  vector standard deviation
 *  @return void
 */

void MobiSaver::saveProtein(ProteinModels prot,
		vector<double> everageDistance, vector<double> Scale_distance) {

	if (everageDistance.size() == 0 || Scale_distance.size() == 0)
		ERROR("One vector is empty!", exception);

	vector<double>::iterator walk;
	vector<double>::iterator walk2;

	if (verbose)
		cout << "\n ### Save Pdb file with column b factor modify ###" << endl;

	output << "REMARK   1  template_A_NMRMOV.pdb" << "\n"
			<< "REMARK   2 B-factors changed to averaged Scaled Distance x 100 \n"
			<< "REMARK   3 after TMScore super imposition with d0 = 4.0; Occupancy changed to SD \n";

	int g;
	for (unsigned int i = 0; i < prot.original_models.size(); i++) {
		if (prot.original_models[i].size() > 0) {

			aminoOffset = 1;
			g = i + 1;
			output << "MODEL        " << g << endl;
			//method of class Component. It checks how deep is the spacer

			walk = everageDistance.begin();
			walk2 = Scale_distance.begin();
			//saving is one ammino at a time
			for (unsigned int j = 0; j < prot.original_models[i].sizeAmino();
					j++) {
				AminoAcid gr = prot.original_models[i].getAmino(j);
				gr.sync();

				for (unsigned int y = 0; y < gr.size(); y++) {
					string atName = gr[y].getType();

					if (atName == "OXT") // cosmetics: OXT has to be output after
						continue; // the sidechain and therefore goes in saveSpacer

					// Added variable for correcting atom type H (last column in PDBs)
					char atomOneLetter;
					if (!isdigit(atName[0])) {
						atomOneLetter = atName[0];
					} else {
						atomOneLetter = atName[1];
					}

					// Added control for size by Damiano Piovesan
					// example HG12
					if (!isdigit(atName[0]) && (atName.size() < 4))
						atName = ' ' + atName;
					while (atName.size() < 4)
						atName += ' ';

					output << "ATOM" << setw(7) << gr[y].getNumber() << " "
							<< atName << " " << gr.getType() << " " << chain
							<< setw(4) << aminoOffset << "    " << setw(8)
							<< setprecision(3) << gr[y].getCoords().x << setw(8)
							<< setprecision(3) << gr[y].getCoords().y << setw(8)
							<< setprecision(3) << gr[y].getCoords().z << "   "
							<< fixed << setprecision(2)<< *(walk2)
							<< setw(6) << fixed << setprecision(2)
							<< (*(walk) * 100) << "           " << atomOneLetter << "\n";

					atomOffset = gr[y].getNumber() + 1;

				}
				walk++;
				walk2++;
				aminoOffset++;
			}

			output << "ENDMODL\n";
		}

	}

	endFile();

	if (verbose)
		cout << "Pdb file saved!" << endl;

}
