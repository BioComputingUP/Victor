/* This file is part of Victor
 Victor is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 Victor is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.
 You should have received a copy of the GNU General Public License
 along with Victor. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file ScaleDistance.h
 * @author Riccardo Zanella
 * @date 11 lug 2015
 * @version 0.1
 */

#ifndef MOBI_SOURCES_STANDARDDEVIATION_H_
#define MOBI_SOURCES_STANDARDDEVIATION_H_

#include<ProteinModels.h>
#include<AtomCode.h>
#include <string.h>
#include <math.h>

using namespace Victor::Biopool;

namespace Victor {
namespace Mobi {

class StandardDeviation {

	/**
	 * @brief Extends Protein class with functionalities related to manipulation of NMR model and models comparations.
	 * Scale Distance metric is provided.
	 */

public:

	// CONSTRUCTORS/DESTRUCTOR:
	StandardDeviation(const ProteinModels& _orig, bool _verbose = false);

	virtual ~StandardDeviation();

	vector<double> get_everage_distance();
	vector<double> get_standard_deviation();
	vector<double> get_StandarDev_angle_PHI();
	vector<double> get_StandarDev_angle_PSI();

private:
	void getCaAtom(Spacer* s, bool flag);

private:
	bool verbose;
	vector<Spacer> original_models;
	vector<Spacer> models;
	vector<Atom> CaVector1, CaVector2;
	vector<vector<double> > dist_from_Ca_atoms;
	vector<double> dist_everage;
	vector<double> ScD;
	vector<double> angle_PHI;
	vector<double> angle_PSI;

};
}
}

#endif /* MOBI_SOURCES_STANDARDDEVIATION_H_ */
