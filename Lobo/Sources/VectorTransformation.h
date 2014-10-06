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
/** 
 *@Class:              VectorTransformation 
*/
#ifndef _VECTORTRANSFORMATION_H_
#define _VECTORTRANSFORMATION_H_

// Includes:
#include <Debug.h>
#include <vector>
#include <vector3.h>
#include <matrix3.h>
#include <IntCoordConverter.h>

namespace Biopool {
  
// Global constants, typedefs, etc. (to avoid):

  /**@brief  This class allows to store transformation steps for transforming a 
*     vector v into v' according to the series of steps performed earlier.
 * 
*@Description  
 * */
class VectorTransformation{
public: 

// CONSTRUCTORS/DESTRUCTOR:
  VectorTransformation();
  VectorTransformation(const VectorTransformation& orig);
  virtual ~VectorTransformation();

// PREDICATES:
  vgVector3<float> transform(vgVector3<float> orig);

// MODIFIERS:
  void addAlignVectors(vgVector3<float> v1, vgVector3<float> v2);
  void addRot(vgMatrix3<float> rm);
  void addTrans(vgVector3<float> t);
  void clear();

  void copy(const VectorTransformation& orig);
  
// OPERATORS:
  VectorTransformation& operator=(const VectorTransformation& orig);

protected:

private:

// HELPERS: 
  void addNewElem();

// ATTRIBUTES:
  vector<vgMatrix3<float> > rot;
  vector<vgVector3<float> > trans;

};

// ---------------------------------------------------------------------------
//                                  LoopModel
// -----------------x-------------------x-------------------x-----------------

// MODIFIERS:
/**
 *@Description Adds the new aligned vector by rotation
*@param  the two vectors to align(vgVector3<float> ,vgVector3<float> )
 *@return  changes are made internally(void)
*/
inline void VectorTransformation::addAlignVectors(vgVector3<float> v1, 
vgVector3<float> v2){
  vgMatrix3<float> tmpRefMat(1);
  alignVectors(v1, v2, tmpRefMat);
  addRot(tmpRefMat);
}


// HELPERS: 
/**
 *@Description Adds new element in the rotation matrix and in the translation vector.
*@param  none
 *@return  changes are made internally(void)
*/
inline void VectorTransformation::addNewElem(){
  vgMatrix3<float> tmpM(1);
  rot.push_back(tmpM);
  vgVector3<float> tmpV(0,0,0);
  trans.push_back(tmpV);
}

/**
 *@Description Converts the vector of float into a vector of doubles 
*@param  the original vector of floats(vgVector3<float>)
 *@return   vector of double(vgVector3<double>)
*/
inline vgVector3<double> convert( vgVector3<float> v){
  vgVector3<double> _v;
  _v.x = v.x;
  _v.y = v.y;
  _v.z = v.z;
  return _v;
}
/**
 *@Description Converts the vector of doubles  into a vector of  float
*@param  the original vector of double(vgVector3<double>)
 *@return   vector of float(vgVector3<float>)
*/
inline vgVector3<float> convert( vgVector3<double> v){
  vgVector3<float> _v;
  _v.x = v.x;
  _v.y = v.y;
  _v.z = v.z;
  return _v;
}
/**
 *@Description Converts the matrix of doubles  into a matrix of  float
*@param  the original matrix of double(vgMatrix3<double>)
 *@return   matrix of float(vgMatrix3<float>)
*/
inline vgMatrix3<float> convert( vgMatrix3<double> v){
  vgMatrix3<float> _v;
  for (unsigned int i = 0; i < 9; i++)
    _v[i] = v[i];
  return _v;
}
/**
 *@Description Converts the matrix of float into a matrix of doubles 
*@param  the original matrix of floats(vgMatrix3<float>)
 *@return   matrix of double(vgMatrix3<double>)
*/
inline vgMatrix3<double> convert( vgMatrix3<float> v){
  vgMatrix3<double> _v;
  for (unsigned int i = 0; i < 9; i++)
    _v[i] = v[i];
  return _v;
}

} // namespace

#endif //_VECTORTRANSFORMATION_H_

