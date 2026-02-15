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

#ifndef __Clustal_W__
#define __Clustal_W__

#include <Cluster.h>
#include <DistanceMatrix.h>
#include <list>
#include <map>
#include <HierarchicalClusterBuilder.h>
#include <SubMatrixFactory.h>
#include <MultipleAlignment.h>

namespace Victor { namespace Phylogenesis{

	/** @brief Multiple Sequence Alignment with ClustalW heuristic.
	* 		   (Thompson, Higgins & Gibson, 1994).
	*   @author Andrea Nasari (andrea.nasari\@studenti.unipd.it)
	*
	**/
	class ClustalW
	{
	public:
		/// Constructor.
		ClustalW(DistanceMatrix* dm, HierarchicalClusterBuilder*  clusterBuilder, vector<Sequence*> sequences,
				double cutoff, double gop, double gep, bool pam, bool verbose);

		/// Perform a Multiple Sequence Alignment.
		MultipleAlignment* performMultipleSequenceAlignment();

		/// Returns the guide tree
		Cluster* getGuideTree();

	private:

		// ATTRIBUTES:

		vector<Sequence*> sequences;
		double cutoff, gop, gep;
		bool pam, verbose;
		Cluster* guideTree;
		DistanceMatrix* dm;
		HierarchicalClusterBuilder*  clusterBuilder;
		SubMatrixFactory* subMatrixFactory;

		// FUNCTIONS:

		vector<Sequence*> multipleSequenceAlignment(const vector<Sequence*>& left, const vector<Sequence*>& right);
		vector<Sequence*> alignCluster(Cluster* c);
		vector<Sequence*> alignMostDivergentSequences(vector<Sequence*> msa);
		double getPercentIdentity(const vector<Sequence*>& left, const vector<Sequence*>& right);
		void getSubMatrix(double distance, SubMatrix*& sub, double& scale);
		bool selectMostDivergentSequences();
	};

}}

#endif
