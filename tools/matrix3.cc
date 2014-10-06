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

static const char* __rcsid__ = "@(#) $Id: matrix3.cc,v 1.5 2008-05-09 13:58:50 biocomp Exp $"; static const char __use_rcsid__ = (&__use_rcsid__ - __rcsid__);


#include <math.h>
#include <matrix3.h>
#include <vglSwap.h>
#include <vglMath.h>
using namespace Biopool;


template <class T> vgMatrix3<T>& vgMatrix3<T>:: operator = (const vgMatrix3<T>& other)
{
  x = other.x;
  y = other.y;
  z = other.z;
  return *this;
}

template <class T> vgMatrix3<T>:: vgMatrix3(const vgMatrix3<T>& other) {
  x = other.x;
  y = other.y;
  z = other.z;
}

template <class T> void vgMatrix3<T>:: setDiagonal(T d1, T d2, T d3)
{
  x.x = d1; x.y = (T)0; x.z = (T)0; 
  y.x = (T)0; y.y = d2; y.z = (T)0; 
  z.x = (T)0; z.y = (T)0; z.z = d3; 
}

template <class T> vgMatrix3<T> vgMatrix3<T>::operator * (T scale) const
{
  vgMatrix3<T> result;
  result.x = x*scale;
  result.y = y*scale;
  result.z = z*scale;
  return result;
}

template <class T> vgVector3<T> vgMatrix3<T>::operator * (const vgVector3<T>& v) const
{
  vgVector3<T> result;
  result.x = x*v;
  result.y = y*v;
  result.z = z*v;
  return result;
}

template <class T> vgMatrix3<T> vgMatrix3<T>:: operator * (const vgMatrix3<T>& other) const
{
  vgMatrix3<T> result;

  result.x.x = x.x*other.x.x + x.y*other.y.x + x.z*other.z.x;
  result.x.y = x.x*other.x.y + x.y*other.y.y + x.z*other.z.y;
  result.x.z = x.x*other.x.z + x.y*other.y.z + x.z*other.z.z;

  result.y.x = y.x*other.x.x + y.y*other.y.x + y.z*other.z.x;
  result.y.y = y.x*other.x.y + y.y*other.y.y + y.z*other.z.y;
  result.y.z = y.x*other.x.z + y.y*other.y.z + y.z*other.z.z;

  result.z.x = z.x*other.x.x + z.y*other.y.x + z.z*other.z.x;
  result.z.y = z.x*other.x.y + z.y*other.y.y + z.z*other.z.y;
  result.z.z = z.x*other.x.z + z.y*other.y.z + z.z*other.z.z;

  return result;
}

template <class T> vgMatrix3<T> vgMatrix3<T>:: operator + (const vgMatrix3<T>& other) const
{
  vgMatrix3<T> result;
  result.x = x+other.x;
  result.y = y+other.y;
  result.z = z+other.z;
  return result;
}

template <class T> vgMatrix3<T> vgMatrix3<T>::operator - (const vgMatrix3<T>& other) const
{
  vgMatrix3<T> result;
  result.x = x-other.x;
  result.y = y-other.y;
  result.z = z-other.z;
  return result;
}

template <class T> vgMatrix3<T>& vgMatrix3<T>:: operator *= (T scale)
{
  x *= scale;
  y *= scale;
  z *= scale;
  return *this;
}

template <class T> vgMatrix3<T>& vgMatrix3<T>:: operator *= (const vgMatrix3<T>& other)
{
  vgMatrix3<T> result;

  result.x.x = x.x*other.x.x + x.y*other.y.x + x.z*other.z.x;
  result.x.y = x.x*other.x.y + x.y*other.y.y + x.z*other.z.y;
  result.x.z = x.x*other.x.z + x.y*other.y.z + x.z*other.z.z;

  result.y.x = y.x*other.x.x + y.y*other.y.x + y.z*other.z.x;
  result.y.y = y.x*other.x.y + y.y*other.y.y + y.z*other.z.y;
  result.y.z = y.x*other.x.z + y.y*other.y.z + y.z*other.z.z;

  result.z.x = z.x*other.x.x + z.y*other.y.x + z.z*other.z.x;
  result.z.y = z.x*other.x.y + z.y*other.y.y + z.z*other.z.y;
  result.z.z = z.x*other.x.z + z.y*other.y.z + z.z*other.z.z;

  *this = result;
  return *this;
}

template <class T> vgMatrix3<T>& vgMatrix3<T>:: operator += (const vgMatrix3<T>& other)
{
  x += other.x;
  y += other.y;
  z += other.z;
  return *this;
}

template <class T> vgMatrix3<T>& vgMatrix3<T>:: operator -= (const vgMatrix3<T>& other)
{
  x -= other.x;
  y -= other.y;
  z -= other.z;
  return *this;
}

template <class T> vgMatrix3<T> vgMatrix3<T>:: operator - () const
{
  vgMatrix3<T> result;
  result.x = -x;
  result.y = -y;
  result.z = -z;
  return result;
}

template <class T> bool vgMatrix3<T>:: operator == (const vgMatrix3<T>& other) const
{
  T* values = (T*)this;
  T* otherValues = (T*)(&other);
  ptrdiff_t k;
  for (k=0; k<9; k++) {
    if (values[k] != otherValues[k]) return false;
  }
  return true;
}

template <class T> vgMatrix3<T> vgMatrix3<T>:: inverse() const
{
  vgMatrix3<T> result((T)1);
  vgMatrix3<T> src = *this;
  // inverse (gauss algorithm)
  T pivot, newPivot;
  pivot = vg_abs(src.x.x);
  if ((newPivot = vg_abs(src.y.x)) > pivot) { pivot = newPivot; swap(src.y,src.x); swap(result.y,result.x); }
  if ((newPivot = vg_abs(src.z.x)) > pivot) { pivot = newPivot; swap(src.z,src.x); swap(result.z,result.x); }
  // now: x.x is biggest element in first row
  if (pivot == 0) return result;
  pivot = 1/src.x.x;
  result.x *= pivot; src.x *= pivot;
  result.y -= result.x * src.y.x; src.y -= src.x * src.y.x;
  result.z -= result.x * src.z.x; src.z -= src.x * src.z.x;
  // now: first row is (1, 0, 0)
  pivot = vg_abs(src.y.y);
  if ((newPivot = vg_abs(src.z.y)) > pivot) { pivot = newPivot; swap(src.z,src.y); swap(result.z,result.y); }
  // now: y.y is biggest element in second row
  if (pivot == 0) return result;
  pivot = 1/src.y.y;
  result.y *= pivot; src.y *= pivot;
  result.z -= result.y * src.z.y; src.z -= src.y * src.z.y;
  // now: second row is (?, 1, 0)
  if (src.z.z == 0) return result;
  pivot = 1/src.z.z;
  result.z *= pivot; src.z *= pivot;
  // now: third row is (?, ?, 1)
  result.y -= result.z * src.y.z; src.y -= src.z * src.y.z; 
  result.x -= result.z * src.x.z; src.x -= src.z * src.x.z; 
  // now: third row is (0, 0, 1, 0)
  result.x -= result.y * src.x.y; 
  // now: second row is (0, 1, 0, 0)
  return result;
}

template <class T> vgMatrix3<T> vgMatrix3<T>:: transpose() const
{
  vgMatrix3<T> result;
  result.x.x = x.x;
  result.x.y = y.x;
  result.x.z = z.x;

  result.y.x = x.y;
  result.y.y = y.y;
  result.y.z = z.y;

  result.z.x = x.z;
  result.z.y = y.z;
  result.z.z = z.z;
  
  return result;
}


template <class T> vgMatrix3<T> vgMatrix3<T>:: createRotationMatrix(vgVector3<T> n, vg_ieee64 alpha)
{
  vgMatrix3<T> R((T)1);
  vgMatrix3<T> S;
  vgMatrix3<T> U;

  if ((n.x == 0) && (n.y == 0) && (n.z == 0)) 
    return R;

  n.normalize();

  vgVector3<T> nSin = n * (T)sin(alpha);

  S.x.x = (T)0;    S.x.y = -nSin.z; S.x.z = nSin.y; 
  S.y.x = nSin.z;  S.y.y = (T)0;    S.y.z = -nSin.x;
  S.z.x = -nSin.y; S.z.y = nSin.x;  S.z.z = (T)0;      

  U.x.x = n.x*n.x; U.x.y = n.x*n.y; U.x.z = n.x*n.z;
  U.y.x = n.y*n.x; U.y.y = n.y*n.y; U.y.z = n.y*n.z;
  U.z.x = n.z*n.x; U.z.y = n.z*n.y; U.z.z = n.z*n.z;

  R -= U;

  R *= (T)cos(alpha); 
  R += U;
  R += S;

  return R;
}

#if !defined __GNUC__		// gcc has problems with global initializers (!?)
const char* vgClassInfo<vgMatrix3<vg_ieee32> >::name = "vgMatrix3<vg_ieee32>";
const size_t vgClassInfo<vgMatrix3<vg_ieee32> >::size = sizeof(vgMatrix3<vg_ieee32>);
const bool vgClassInfo<vgMatrix3<vg_ieee32> >::numeric = false;
const bool vgClassInfo<vgMatrix3<vg_ieee32> >::ordered = false;
const vgMatrix3<vg_ieee32> vgClassInfo<vgMatrix3<vg_ieee32> >::min = (vg_ieee32)0;
const vgMatrix3<vg_ieee32> vgClassInfo<vgMatrix3<vg_ieee32> >::max = (vg_ieee32)0;
const char* vgClassInfo<vgMatrix3<vg_ieee64> >::name = "vgMatrix3<vg_ieee64>";
const size_t vgClassInfo<vgMatrix3<vg_ieee64> >::size = sizeof(vgMatrix3<vg_ieee64>);
const bool vgClassInfo<vgMatrix3<vg_ieee64> >::numeric = false;
const bool vgClassInfo<vgMatrix3<vg_ieee64> >::ordered = false;
const vgMatrix3<vg_ieee64> vgClassInfo<vgMatrix3<vg_ieee64> >::min = (vg_ieee64)0;
const vgMatrix3<vg_ieee64> vgClassInfo<vgMatrix3<vg_ieee64> >::max = (vg_ieee64)0;
#endif

#if defined __GNUC__		// for inlining of template class methods (!?)
typedef vgVector3<vg_ieee32> __vgVector3_vg_ieee32;
typedef vgVector3<vg_ieee64> __vgVector3_vg_ieee64;

void (*__swap_vgVector3_vg_ieee32_vgVector3_vg_ieee32)(vgVector3<vg_ieee32>&, vgVector3<vg_ieee32>&) = &swap;
void (*__swap_vgVector3_vg_ieee64_vgVector3_vg_ieee64)(vgVector3<vg_ieee64>&, vgVector3<vg_ieee64>&) = &swap;
#endif

#if defined __GNUC__
template class vgMatrix3<vg_ieee32>;
template class vgMatrix3<vg_ieee64>;
#endif

#if defined __MSCC__
template vgMatrix3<vg_ieee32>; 
template vgMatrix3<vg_ieee64>; 
#endif

#if defined __DECCXX 
#pragma define_template vgMatrix3<vg_ieee32>
#pragma define_template vgMatrix3<vg_ieee64>
#endif

#if defined __SGICC__
#pragma instantiate vgMatrix3<vg_ieee32>
#pragma instantiate vgMatrix3<vg_ieee64>
#endif


