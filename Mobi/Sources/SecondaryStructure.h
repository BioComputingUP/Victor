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
 * @file SecondaryStructure.h
 * @author Riccardo Zanella
 * @date 21 lug 2015
 * @version 0.1
 */

// Includes:
#ifndef MOBI_SOURCES_SECONDARYSTRUCTURE_H_
#define MOBI_SOURCES_SECONDARYSTRUCTURE_H_

#include<ProteinModels.h>

using namespace Victor::Biopool;

namespace Victor {
namespace Mobi {

/**@brief Implements a Secondary Structure.
 * Includes methods to obtain mobility from the analyze of Secondary Structure.
 */

class SecondaryStructure {
public:

	// CONSTRUCTORS/DESTRUCTOR:

	SecondaryStructure(const ProteinModels& _orig, bool _verbose = false);
	virtual ~SecondaryStructure();

	// PREDICATES:
	vector< vector<char> >getSecStructFromModels();
	vector<char> getMobilitySecondaryStructure();

private:
	bool verbose;
	vector<Spacer> models;
	vector<char> mobility;

};
}
}
#endif /* MOBI_SOURCES_SECONDARYSTRUCTURE_H_ */
