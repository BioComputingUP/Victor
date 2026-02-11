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


#ifndef MOBI_SOURCES_MOBISAVER_H_
#define MOBI_SOURCES_MOBISAVER_H_

//Includes:
#include <PdbSaver.h>
#include <ProteinModels.h>
#include <vector>
#include <iterator>
#include <iostream>
#include <fstream>
#include <unistd.h>

using namespace Victor::Biopool;

namespace Victor {
namespace Mobi {

/**@brief Implements a MobiSaver.
 * This class inherit PdbSaver class. This class allow to Saver in output file all type of vector generated from
 * Standard Deviation and SecondaryStructure classes to have the mobility.
 * You can save Mobility also filtered by same type of mask.
 */

class MobiSaver : public PdbSaver{
public:

	// CONSTRUCTORS/DESTRUCTOR:
	MobiSaver(ProteinModels __protein,string __output, ofstream& __stream, bool __verbose = false, double __boundSD = 0.85, double __boundStandD = 0.09, double __anglePHI = 20, double __anglePSI = 20);
	virtual ~MobiSaver();

	// PREDICATES:

	void mob_eveScalD(vector <double> __everageDistance);
	void mob_eveScalD_filtSecS(vector <double> __everageDistance, vector<char> __mobSecS);
	void mob_eveScalD_filteredByPHI_PSI_standD(vector <double> __everageDistance, vector <double> __angle_PHI, vector <double> __angle_PSI, vector <double> __Scale_distance);
	void mob_stanD(vector <double> __Scale_distance);
	void mob_stanD_withMask(vector <double> __Scale_distance);
	void mob_aPHI(vector <double> __angle_PHI);
	void mob_aPSI(vector <double> __angle_PSI);
	void mob_SecS(vector <char> __Mob_SecStructure);
	void allMobility(vector<double> __everageDistance,
			vector<double> __Scale_distance, vector<double> __angle_PHI,
			vector<double> __angle_PSI, vector<char> __Mob_SecStructure);
	virtual void saveProtein(ProteinModels __protein, vector<double> __everageDistance, vector<double> __StandardDeviation);

private:
	bool verbose;
	string out;
	double ScalD, StandD, angPHI, angPSI;
};
}
}

#endif /* MOBI_SOURCES_MOBISAVER_H_ */
