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

#include <IoTools.h>
#include <SingleCluster.h>
#include <NewickPhylogeneticTreeExporter.h>
#include <sstream>
#include <fstream>
#include <iostream>

namespace Victor { namespace Phylogenesis{

	//g++ bug? http://stackoverflow.com/questions/19122574/to-string-isnt-a-member-of-std
	template <typename T>
	string to_string(T val)
	{
		stringstream stream;
		stream << val;
		return stream.str();
	}

	NewickPhylogeneticTreeExporter::NewickPhylogeneticTreeExporter(){
	}

	NewickPhylogeneticTreeExporter::NewickPhylogeneticTreeExporter(std::string fileName) : fileName(fileName) {
	}

	std::string NewickPhylogeneticTreeExporter::buildRec(Cluster* cluster, double distanceFromFather) {

		std::stringstream result;

		// leaf node A:dist
		if(cluster->isSingle())
		{
			string name = static_cast<SingleCluster*>(cluster)->sequence->name;
			name = name.substr(0,name.find(' '));
			return name + ":" + to_string(distanceFromFather);
		}

		// internal node (..,..,..):distanceFromFather - current
		result << "(";
		for(unsigned int i = 0; i < cluster->branches.size(); i++)
		{
			if(i>0)
				result << ",";
			result << buildRec(cluster->branches[i], cluster->lengths[i]);
		}
		result << ")";

		// root node
		if(distanceFromFather == 0)
		{
			result<<";";
			return result.str();
		}

		result<<":"<< to_string(distanceFromFather);
		return result.str();
	}

	void NewickPhylogeneticTreeExporter::save(Cluster* tree){
		std::string result = buildRec(tree, 0.0);

		if(!fileName.empty())
		{
			std::cout<<"\nWriting newick file "<<fileName<<std::endl;
			std::ofstream  output;
			output.open(fileName.c_str());
			output<<result;
			output.close();
		}
		else{
			std::cout<<std::endl<<result<<std::endl;
		}
	}

	void NewickPhylogeneticTreeExporter::save(Cluster* tree, ostream &output){
		std::string result = buildRec(tree, 0.0);
		output << result;
	}

}}
