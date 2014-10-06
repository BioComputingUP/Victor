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
@Description  matrix with 3x3 entries

*/

#ifndef _TYPES_MATRIX3_H_
#define _TYPES_MATRIX3_H_

#include <vglStd.h>
#include <vector3.h>

template <class T>
/** @brief class implements matrix methods
* 
* @Description Includes methods that allow to get and set , calculate, etc.
 **/
class vgMatrix3
{
public:
  vgMatrix3() {}
  vgMatrix3(T d) { setDiagonal(d); }
  vgMatrix3(T d1, T d2, T d3) { setDiagonal(d1,d2,d3); }
  vgMatrix3(const vgMatrix3<T>& other);
  vgMatrix3<T>& operator = (const vgMatrix3<T>& other);

  void setDiagonal(T d1, T d2, T d3);
  void setDiagonal(T d) { setDiagonal(d,d,d); }

  T& operator [] (ptrdiff_t index) { return ((T*)this)[index]; } 
  const T& operator [] (ptrdiff_t index) const { return ((T*)this)[index]; } 

  vgMatrix3<T> operator * (T scale) const;
  vgVector3<T> operator * (const vgVector3<T>& v) const;
  vgMatrix3<T> operator * (const vgMatrix3<T>& other) const;
  vgMatrix3<T> operator + (const vgMatrix3<T>& other) const;
  vgMatrix3<T> operator - (const vgMatrix3<T>& other) const;
  vgMatrix3<T>& operator *= (T scale);
  vgMatrix3<T>& operator *= (const vgMatrix3<T>& other);
  vgMatrix3<T>& operator += (const vgMatrix3<T>& other);
  vgMatrix3<T>& operator -= (const vgMatrix3<T>& other);
  vgMatrix3<T> operator - () const;

  bool operator == (const vgMatrix3<T>& other) const;
  bool operator != (const vgMatrix3<T>& other) const { return !(*this == other); }

  vgMatrix3<T> inverse() const;
  vgMatrix3<T> transpose() const;

  vgVector3<T> x;		// 1st row vector
  vgVector3<T> y;		// 2nd row vector
  vgVector3<T> z;		// 3rd row vector

  static vgMatrix3<T> createRotationMatrix(vgVector3<T> axis, vg_ieee64 angle);
  static vgMatrix3<T> createScaleMatrix(vgVector3<T> s) { return vgMatrix3<T>(s.x, s.y, s.z); }
};

/* helper functions */

template <class T>
inline
vgMatrix3<T>
vgCreateScaleMatrix3(const vgVector3<T>& s)
{
  return vgMatrix3<T>::createScaleMatrix(s);
}

template <class T>
inline
vgMatrix3<T>
vgCreateRotationMatrix3(vgVector3<T> n, vg_ieee64 alpha)
{
  return vgMatrix3<T>::createRotationMatrix(n, alpha);
}

template <class T>
inline 
vgVector3<T>
operator * (const vgVector3<T>& v,const vgMatrix3<T>& m)
{
  return vgVector3<T>(v.x*m.x.x + v.y*m.y.x + v.z*m.z.x,
		      v.y*m.x.y + v.y*m.y.y + v.z*m.z.y,
		      v.x*m.x.z + v.y*m.y.z + v.z*m.z.z);
}

#ifdef __GNUC__			
typedef vgMatrix3<vg_ieee32> __vgMatrix3_vg_ieee32; // for inline optimization
typedef vgMatrix3<vg_ieee64> __vgMatrix3_vg_ieee64;
#endif

#endif /* _TYPES_MATRIX3_H_ */
