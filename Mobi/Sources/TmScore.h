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
 * @file TmScore.h
 * @author Riccardo Zanella
 * @date Lug 2015
 * @version 0.1
 */

#ifndef MOBI_SOURCES_TMSCOREBIN_H_
#define MOBI_SOURCES_TMSCOREBIN_H_

//Include:
#include <iostream>
#include <string>
#include <ProteinModels.h>

using namespace Victor::Biopool;
using namespace Victor::Mobi;
using namespace std;

extern const string TMTMP_OUT;

namespace Victor {

namespace Mobi {

/**
 * @brief TMScore functionalities through external binary.
 */

class TmScore {

	// CONSTRUCTORS/DESTRUCTOR:
public:

	TmScore(string _binary = "TMScore", string output = ".", bool _verbose =
			false);

	// PREDICATES:
	Protein* TmImpose(string modelFile, string nativeFile);

	// MODIFIERS:
	/**
		 You can chose if you want verbose output from this class
		 @param  none
		 @return void
		 */
	void setVerbose(bool v) {
		verbose = v;
	}
	/**
	 * Default deconstructor
	 */
	~TmScore();

private:
	string binary;
	const string outputTmScore;
	bool verbose;
};
}
}
#endif /* MOBI_SOURCES_TMSCOREBIN_H_ */
