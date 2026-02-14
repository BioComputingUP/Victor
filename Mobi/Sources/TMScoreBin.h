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
#ifndef MOBI_SOURCES_TMSCOREBIN_H_
#define MOBI_SOURCES_TMSCOREBIN_H_


#include <iostream>
#include <string>
#include <MobiProtein.h>
#include <TMScoreBinder.h>

using namespace Victor::Biopool;
using namespace Victor::Mobi;
using namespace std;


extern const std::string TMTMP_IN1;
extern const std::string TMTMP_IN2;
extern const std::string TMTMP_OUT;



namespace Victor { namespace Mobi {
	/**
	 * @brief TMScore functionalities through external binary.
	 */
	class TMScoreBin: public TMScoreBinder{
	public:

		/**
		 * @brief Setup.
		 * @param _binary (string) full path to binary TMScore file, must have execution permission
		 * @param _tmp (string) full path to temp dir, must have write permission
		 * @param _verbose (bool) verbosity
		 */
		TMScoreBin(std::string _binary = "TMscore", std::string _tmp = ".", bool _verbose = false) :
			binary(_binary),
			tmp(_tmp.substr(_tmp.length()-1,1) == "/" ? _tmp : _tmp + "/"),
			verbose(_verbose)
		{};

		/**
		 * @brief Given two pdb files containing each a model of the same protein,
		 * call TMScore binary to superimpose the first over the second. The superimposed (rotated/traslated)
		 * model is then loaded in a MobiProtein and a pointer is returned.
		 * @param model1 full path to model#1 file
		 * @param model2 full path to model#2 file
		 * @return MobiProtein reference containing the super imposed model
		 */
		virtual MobiProtein* TMScore(string model1, string model2);

		/**
		 * @brief Given a ProteinModel call TMScore binary to superimpose two models contained in it.
		 * The superimposed (rotated/traslated)model is then loaded in a MobiProtein and a pointer is returned.
		 * @param prot reference to MobiProtein object
		 * @param model1 model#1 index in MobiProtein object
		 * @param model2 model#2 index in MobiProtein object
		 * @return MobiProtein reference containing the super imposed model
		 */
		virtual MobiProtein* TMScore(MobiProtein& prot, unsigned int model1, unsigned int model2);


		/**
		 * @brief Superimposition of two models through call to TMScore binary
		 * The superimposed (rotated/traslated)model is then loaded in a MobiProtein and a pointer is returned.
		 * @param prot1 reference to the fist ProteinModel
		 * @param model1 model#1 name in the first ProteinModel object
		 * @param prot2  reference to the second ProteinModel
		 * @param model2 model#2 name in the second ProteinModel object
		 * @return MobiProtein reference containing the super imposed model
		 */
		virtual MobiProtein* TMScore(MobiProtein& prot1, unsigned int model1, MobiProtein& prot2, unsigned int model2);


		/**
		 * Set verbosity.
		 * Verbosity is applied on PdbLoader instances.
		 * @param v new verbosity option
		 */
		void setVerbose(bool v){
			verbose = v;
		}

		/**
		 * Delete temporary files
		 */
		void cleanup();

		/**
		 * Default deconstructor
		 */
		virtual ~TMScoreBin(){};


	private:
		std::string binary;
		std::string tmp;
		bool verbose;
	};

}}



#endif /* MOBI_SOURCES_TMSCOREBIN_H_ */
