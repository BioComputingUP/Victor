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

#ifndef MOBI_SOURCES_TMSCOREBINDER_H_
#define MOBI_SOURCES_TMSCOREBINDER_H_

#include<Spacer.h>
#include<MobiProtein.h>

using namespace Victor::Biopool;
using namespace Victor::Mobi;

namespace Victor { namespace Mobi {
	/**
	 * @brief base class that represent a binder to TMScore functionalities.\n
	 * Class that extends TMScoreBinder must provide reimplementation for virtual methods and can be used in Mobi as Interfaces to TMScore imposition
	 */
	class TMScoreBinder{
	public:

		/**
		 * @brief Given a MobiProtein superimpose two models contained in it.
		 * The superimposed (rotated/traslated) model is then loaded in the heap and pointer is returned.
		 * @param prot (MobiProtein&) reference to ProteinModel object
		 * @param model1 (unsigned int) model#1 name in ProteinModel object
		 * @param model2 (unsigned int) model#2 name in ProteinModel object
		 * @return MobiProtein reference containing the super imposed model
		 */
		virtual MobiProtein* TMScore(MobiProtein& prot, unsigned int model1, unsigned int model2){
			ERROR("Unimplemented Method TMImpose(ProteinModel&, uint, uint, ProteinModel**)!",exception);
			return 0;
		}

		/**
		 * @brief Superimposition of two models
		 * The superimposed (rotated/traslated) model is then loaded in the heap and pointer is returned.
		 * @param prot1 (MobiProtein&) reference to the fist ProteinModel
		 * @param model1 (unsigned int) model#1 name in the first ProteinModel object
		 * @param prot2 (MobiProtein&) reference to the second ProteinModel
		 * @param model2 (unsigned int) model#2 name in the second ProteinModel object
		 * @return MobiProtein reference containing the super imposed model
		 */
		virtual MobiProtein* TMScore(MobiProtein& prot1, unsigned int model1, MobiProtein& prot2, unsigned int model2){
			ERROR("Unimplemented Method TMImpose (ProteinModel&, uint, ProteinModel&, uint, ProteinModel**)!",exception);
			return 0;
		}

		/**
		 * Default destructor
		 */
		virtual ~TMScoreBinder(){}
	};

}}



#endif /* MOBI_SOURCES_TMSCOREBINDER_H_ */
