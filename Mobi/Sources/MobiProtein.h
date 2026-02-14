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
#ifndef MOBI_SOURCES_PROTEINMODEL_
#define MOBI_SOURCES_PROTEINMODEL_

#include<Protein.h>
#include<PdbLoader.h>
#include<AtomCode.h>

using namespace Victor::Biopool;

namespace Victor { namespace Mobi {
	/**
	 * @brief Extends Protein class with functionalities related to manipulation of NMR model and models comparations.
	 */
	class MobiProtein : public Protein{
	public:
		MobiProtein() : Protein(), verbose(false){modelsID = vector<unsigned int>(0);};

		/**
		 * Load Pdb data into Protein object. The specified models are loaded.
		 * There is no direct relation between models id in pdb file and model
		 * index inside protein object. Index inside protein depends on the order
		 * in which the models are loaded.
		 * @param pl reference to the PdbLoader to load from
		 * @param models models to load (accordind to pdb file model names), vector<unsigned int>
		 */
		vector<unsigned int> load(PdbLoader& pl, vector<unsigned int> models);

		/**
		 * Load Pdb data into Protein object. All models are loaded.
		 * There is no direct relation between models id in pdb file and model
		 * index inside protein object. Index inside protein depends on the order
		 * in which the models are loaded.
		 * @param pl (PdbLoader&) reference to the PdbLoader to load from
		 */
		vector<unsigned int> load(PdbLoader& pl);

		/**
		 * Get the specified model in this Protein.
		 * @param model (unsigned int) the model internal id (not original Pdb model name) to get
		 * @return the model
		 */
		Spacer* getModel(unsigned int model);

		/**
		 * Set verbosity to this class methods
		 * @param v verbosity true/false
		 */
		void setVerbose(bool v){
			verbose = v;
		}

		/**
		 * Returns all original (as PDB) loaded modell ids.
		 * @return vector containing original Pdb ids of loaded models. Vector index matches internal id.
		 */
		const vector<unsigned int> getModelsID(){ return modelsID;}

		/**
		 * Returns the the model by its original pdb ID
		 * @param model (unsigned int) the model pdb id
		 * @return the model
		 */
		Spacer* getModelByPdbID(unsigned int model);

	private:
		bool verbose;
		vector<unsigned int> modelsID;
	};
}}


#endif /* MOBI_SOURCES_PROTEINMODEL_ */
