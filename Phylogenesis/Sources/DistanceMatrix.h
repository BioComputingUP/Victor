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

#ifndef __DistanceMatrix_H__
#define __DistanceMatrix_H__

#include <Sequence.h>
#include <GapFunction.h>
#include <SubMatrix.h>
#include <string>
#include <vector>

using namespace Victor::Align2;

namespace Victor { namespace Phylogenesis{

	/// Align algorithm type
	enum AlignAlgorithm { global = 0, local = 1, freeshift = 2 };

	/** @brief    Implements a triangular matrix useful to store distances between sequences without wasting memory.
	 *			  Access and update operations are performed in constant time.
	 *	@author   Andrea Nasari (andrea.nasari\@studenti.unipd.it)
	 *
	 **/
	class DistanceMatrix{
	public:

		DistanceMatrix(){};

		/// Constructor that init the matrix from file. The file must be a sequence indexes of in the correct order, followed by the value (see samples/distance.dat).
		DistanceMatrix(std::string fileName);

		/// Constructor.
		DistanceMatrix(vector<Sequence*> sequences, AlignAlgorithm alignAlgorithm, SubMatrix* sub, double gop, double gep);

		/// Destructor.
		virtual ~DistanceMatrix();

		/// Init the matrix from sequences calculating the pairwise alignments.
		void init();

		/// Add a new value to the matrix, insertion order must be respected.
		void appendDistance(int i, double value);

		/// Add a new value to the matrix, insertion order must be respected.
		void appendIdentity(int i, double value);

		/// Distance value access element
		double at(int i, int j);

		/// Identity value access element
		double identityAt(int i, int j);

	protected:

		//// ClustalW distance score from two optimal aligned sequences
		virtual double calculateDistance(string alignedS1, string alignedS2);

		//// ClustalW identity score from two optimal aligned sequences
		virtual double calculateIdentity(string alignedS1, string alignedS2);

	private:

		// ATTRIBUTES:

		std::vector<std::vector<double> > distanceMatrix;
		std::vector<std::vector<double> > identityMatrix;

		GapFunction* gf;
		SubMatrix* sub;
		AlignAlgorithm alignAlgorithm;
		vector<Sequence*> sequences;

		/// Align 2 sequences
		vector<string> pairwiseSequenceAlign(string s1, string s2);
	};

}}

#endif
