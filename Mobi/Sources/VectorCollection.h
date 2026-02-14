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

#ifndef MOBI_SOURCES_VECTORCOLLECTION_H_
#define MOBI_SOURCES_VECTORCOLLECTION_H_

#include <Spacer.h>
#include <AtomCode.h>
#include <Atom.h>
#include <AtomCode.h>
#include <utility>
#include <map>
#include <cmath>
using namespace Victor::Biopool;

//extern double psi;

namespace Victor { namespace Mobi {

/**
 * @brief Class used to manage collections of Mobility results.
 * We use this class to process vectors of values (usually doubles) resulting from mobility operations,
 * like scaled distance values or angles. This class provides also methods to calculate means and std deviations
 * for corresponding (same position) values in different vectors.\n
 * Internally, an id is assigned to each vector. You can specify your own id or use the provided id translator
 * so that, given two (ordered) id values (for example 2 models ids), a unique internal id is assigned to the id pair,
 * then you can:\n
 * - get a list of all the models who partecipate in this collection\n
 * - get vectors related to a model given the model id\n
 * - get the vector related to a model pair\n
 * \b Attention : all vectors must have the same dimension, the first value added to the collection
 * determines the length of all vectors. You can clear the collection to reset this value.
 */
template <class V>
class VectorCollection{

public:
	/**
	 * Default constructor
	 */
	VectorCollection(){
		this->results = std::map<int,std::vector<V> >();
		this->safePairID = true;
	}

	/**
	 * Clear the collection deleting all elements in it
	 */
	void clear(){
		results.clear();
		safePairID = true;
	}

	/**
	 * Add values (for example scaled distance vector) to this collection
	 * @param id (int) values id
	 * @param result (vector<double>&) vector of values
	 */
	void addValue(int id, const std::vector<V>& result){
		if (size() > 0)
			if (result.size() != vectorsSize())
				ERROR("Trying to add a result of non compatible size",exception);
		results.insert(std::make_pair(id,result));
		safePairID = false;
	}

	/**
	 * Add values (for example scaled distance vector) to this collection
	 * @param m1 (int) first model id
	 * @param m2 (int) second model id
	 * @param result (vector<double>&) vector of values
	 */
	void addValue(int m1, int m2, const std::vector<V>& result){
		if (size() > 0)
			if (result.size() != vectorsSize())
				ERROR("Trying to add a result of non compatible size",exception);
		results.insert(std::make_pair(VectorCollection::id(m1,m2),result));
	}

	/**
	 * Return values given the id
	 * @param id (int) id to search
	 * @return (vector<double>&) values, if found
	 */
	vector<V> getValue(int id){
		typename std::map<int,std::vector<V> >::const_iterator it = this->results.find(id);
		if (it != this->results.end())
			return it->second;
		else
			ERROR("Unable to find values with given id",exception);
			return vector<V>();	//never reached, just to hide warnings
	}

	/**
	 * Returns a list of all the model partecipating in this collection
	 * @return (vector<int>) models
	 */
	vector<int> getModels() const{
		if (!safePairID)
			cout << "[VectorCollection] getModels() Warning: vectors has been added with CUSTOM ids!";
		vector<int> models;
		for (typename std::map<int,std::vector<V> >::const_iterator it = this->results.begin(); it != this->results.end(); ++it){
			if (std::find(models.begin(), models.end(), VectorCollection::modelsFromID(it->first)[0]) == models.end())
				models.push_back(VectorCollection::modelsFromID(it->first)[0]);
			if (std::find(models.begin(), models.end(), VectorCollection::modelsFromID(it->first)[1]) == models.end())
				models.push_back(VectorCollection::modelsFromID(it->first)[1]);
		}
		return models;
	}

	/**
	 * Returns all records relative to the specified model
	 * @param model (int) the model to search
	 * @return (VectorCollection) subcollection relative to the model
	 */
	VectorCollection<V> getValuesByModel(int model) const{
		if (!safePairID)
			cout << "[VectorCollection] getValuesByModel() Warning: vectors has been added with CUSTOM ids!";
		VectorCollection outVC;
		typename std::map<int,std::vector<V> >::const_iterator it = this->results.begin();

		while (it != this->results.end()){
			if (it->first / delimiter == model)
				outVC.addValue(it->first, it->second);
			else
				if (it->first % delimiter == model)
					outVC.addValue(it->first, it->second);
			it++;
		}
		return outVC;
	}

	/**
	 * Returns the length of vectors in this collection
	 * @return (unsigned int) length
	 */
	unsigned int vectorsSize() const{
		if (size() < 1)
			ERROR("Cannot get model lenght, since there are no models in this SDResult object",exception);
		return this->results.begin()->second.size();
	}

	/**
	 * Returns the size of this collection
	 * @return (unsigned int) size
	 */
	unsigned int size() const{
		return this->results.size();
	}

	/**
	 * Read-Only Iterator to this collection, begin value
	 * @return (map::<int,std::vector<double> >) iterator
	 */
	typename std::map<int,std::vector<V> >::const_iterator begin() const{
		return this->results.begin();
	}

	/**
	 * Read-Only Iterator to this collection, end value
	 * @return (map::<int,std::vector<double> >) iterator
	 */
	typename std::map<int,std::vector<V> >::const_iterator end() const{
		return this->results.end();
	}

	/**
	 * Calculates the means of values in same position in the vectors
	 * @return (vector<double>) means
	 */
	vector<V> mean() const{
		vector<V> mean = vector<V>(this->vectorsSize(),0.0);
		typename std::map<int,std::vector<V> >::const_iterator it;
		for (it = this->results.begin(); it != this->results.end(); ++it)
			for (unsigned int a = 0; a < this->vectorsSize(); a++)
				mean[a] += it->second[a];
		for (unsigned int a = 0; a < this->vectorsSize(); a++)
			mean[a] /= this->size();
		return mean;
	}

	/**
	 * Calculates the standard deviation of values in same position in the vectors
	 * @return (vector<double>) standard deviations
	 */
	vector<V> stdDev() const{
		vector<V> mean = this->mean();
		vector<V> sd = vector<V>(this->vectorsSize(),0.0);
		typename std::map<int,std::vector<V> >::const_iterator it;
		for (it = this->results.begin(); it != this->results.end(); ++it)	//foreach distance record
			for (unsigned int a = 0; a < this->vectorsSize(); a++)	//foreach residue
				sd[a] += pow(it->second[a] - mean[a],2);	//cumulate the 2pow of distance minus mean
		for (unsigned int a = 0; a < this->vectorsSize(); a++)
			sd[a] = sqrt(sd[a] / this->size());
		return sd;
	}

	/**
	 * Calculate rmsd for all distance records inside this VectorCollection and return their average
	 * @return (double) mean rmsd
	 */
	V RMSD() const{
		V rmsd = 0;
		V singleRmsd;
		typename std::map<int,std::vector<V> >::const_iterator it;
		for (it = this->results.begin(); it != this->results.end(); ++it){	//foreach distance record
			singleRmsd = 0;
			for (unsigned int i = 0; i < this->vectorsSize(); i++)	//foreach residue
				singleRmsd += pow(it->second[i],2);	//cumulate the 2pow of distance
			singleRmsd = sqrt(singleRmsd / this->vectorsSize());
			rmsd += singleRmsd;
		}
		return (rmsd / results.size());
	}

	/**
	 * Short output... debug purpouses
	 */
	void print(){
		typename std::map<int,std::vector<V> >::iterator it;
		for (it = this->results.begin(); it != this->results.end(); ++it)
			cout << it->first << " => [" << it->second[0] << "," << it->second[1] << ", ...]" << endl;
	}

	/**
	 * multiplying value to separate two models, it also represent the max number of model acceptable
	 */
	static const int delimiter = 10000;

	/**
	 * build model pair id from models ids
	 * @param m1 first model id
	 * @param m2 second model id
	 * @return pair id
	 */
	static int id(int m1, int m2){
		return m1*delimiter + m2;
	}

	/**
	 * return two models id given the pair id
	 * @param id pair id
	 * @return 2 sized vector with the 2 model ids
	 */
	static vector<int> modelsFromID(int id){
		vector<int> models(2);
		models[0] = id / delimiter;
		models[1] = id % delimiter;
		return models;
	}


private:
	bool safePairID;

protected:
	/**
	 * Map container for the collection
	 */
	std::map<int,std::vector<V> > results;
};



}}	//namespaces

#endif /* MOBI_SOURCES_VECTORCOLLECTION_H_ */
