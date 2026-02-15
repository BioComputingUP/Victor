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

#ifndef __NewickTreeBuilder_H__
#define __NewickTreeBuilder_H__

#include <Cluster.h>
#include <PhylogeneticTreeExporter.h>
#include <string>

namespace Victor { namespace Phylogenesis{

	/** @brief   Clustering tree exporter in the Newick tree format.
	 *  @author  Andrea Nasari (andrea.nasari\@studenti.unipd.it)
	 *
	 **/
	class NewickPhylogeneticTreeExporter : public PhylogeneticTreeExporter
	{
	public:
		/// Constructor.
		NewickPhylogeneticTreeExporter();

		/// Default constructor.
		NewickPhylogeneticTreeExporter(std::string fileName);

		/// Export the clustering tree in the Newick tree format
		void save(Cluster* tree);

		/// Export the clustering tree in the Newick tree format to output
		void save(Cluster* tree, ostream &output);

	private:

		/// Recursive call to build the clustering tree
		std::string buildRec(Cluster* cluster, double fatherDistance);

		/// Output filename
		std::string fileName;
	};

}}

#endif
