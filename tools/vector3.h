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
@Description  vector with 3 components
 */

#ifndef _TYPES_VECTOR3_H_
#define _TYPES_VECTOR3_H_
#include <vglStd.h>
#include <vglMath.h>

template <class T>
/** @brief class implements a simple atom type.
 * 
 * @Description Includes methods that allow to get and set type, bind, unbind,coordinates , code,  etc.NB Angles are in degrees.
 **/
class vgVector3 {
public:

    vgVector3() {
    };

    vgVector3(T x, T y, T z = (T) 0) {
        this->x = x;
        this->y = y;
        this->z = z;
    };

    vgVector3(const vgVector3<T>& other) {
        x = other.x;
        y = other.y;
        z = other.z;
    };

    vgVector3<T>& operator=(const vgVector3<T>& other) {
        x = other.x;
        y = other.y;
        z = other.z;
        return *this;
    };

    T sum() const {
        return x + y + z;
    };

    T square() const {
        return x * x + y * y + z*z;
    };

    vg_ieee64 length() const {
        return sqrt((vg_ieee64) square());
    };
    inline vgVector3<T>& normalize();
    inline bool operator ==(const vgVector3<T>& other) const;
    inline bool operator !=(const vgVector3<T>& other) const;
    inline bool operator<(const vgVector3<T>& other) const;
    inline bool operator <=(const vgVector3<T>& other) const;
    inline bool operator>(const vgVector3<T>& other) const;
    inline bool operator >=(const vgVector3<T>& other) const;
    inline bool operator<(T value) const;
    inline bool operator <=(T value) const;
    inline bool operator>(T value) const;
    inline bool operator >=(T value) const;

    T& operator [] (ptrdiff_t index) {
        return ((T*)this)[index];
    };
    const T& operator [] (ptrdiff_t index) const {
        return ((T*)this)[index];
    };

    inline vgVector3<T> operator +(const vgVector3<T>& other) const;
    inline vgVector3<T> operator -(const vgVector3<T>& other) const;
    inline vgVector3<T> operator *(T faktor) const;
    inline T operator *(const vgVector3<T>& other) const;
    inline vgVector3<T> operator /(T divisor) const;
    inline vgVector3<T>& operator +=(const vgVector3<T>& other);
    inline vgVector3<T>& operator -=(const vgVector3<T>& other);
    inline vgVector3<T>& operator *=(T faktor);
    inline vgVector3<T>& operator /=(T divisor);
    inline vgVector3<T> operator -() const;

    inline vgVector3<T> cross(const vgVector3<T>& other) const; // cross product

    T x; // 1st component
    T y; // 2nd component
    T z; // 3rd component

};

/* inline functions */

template <class T> inline vgVector3<T>& vgVector3<T>::normalize() {
    vg_ieee64 l = length();
    if (l != 0.0) {
        l = 1.0 / l;
        x = (T) (x * l);
        y = (T) (y * l);
        z = (T) (z * l);
    }
    return *this;
}

template <class T> inline bool vgVector3<T>::operator ==(const vgVector3<T> & other) const {
    return ((x == other.x) && (y == other.y) && (z == other.z));
}

template <class T>inline bool vgVector3<T>::operator !=(const vgVector3<T>& other) const {
    return ((x != other.x) || (y != other.y) || (z != other.z));
}

template <class T> inline bool vgVector3<T>::operator<(const vgVector3<T>& other) const {
    return (this->square() < other.square());
}

template <class T> inline bool vgVector3<T>::operator <=(const vgVector3<T>& other) const {
    return (this->square() <= other.square());
}

template <class T> inline bool vgVector3<T>::operator>(const vgVector3<T>& other) const {
    return (this->square() > other.square());
}

template <class T> inline bool vgVector3<T>::operator >=(const vgVector3<T>& other) const {
    return (this->square() >= other.square());
}

template <class T> inline bool vgVector3<T>::operator<(T length) const {
    return (this->square() < length * length);
}

template <class T> inline bool vgVector3<T>::operator <=(T length) const {
    return (this->square() <= length * length);
}

template <class T> inline bool vgVector3<T>::operator>(T length) const {
    return (this->square() > length * length);
}

template <class T> inline bool vgVector3<T>::operator >=(T length) const {
    return (this->square() >= length * length);
}

template <class T> inline vgVector3<T> vgVector3<T>::operator +(const vgVector3<T>& other) const {
    return vgVector3<T>(x + other.x, y + other.y, z + other.z);
}

template <class T> inline vgVector3<T> vgVector3<T>::operator -(const vgVector3<T>& other) const {
    return vgVector3<T>(x - other.x, y - other.y, z - other.z);
}

template <class T> inline vgVector3<T> vgVector3<T>::operator *(T factor) const {
    return vgVector3<T>(factor*x, factor*y, factor * z);
}

template <class T> inline T vgVector3<T>::operator *(const vgVector3<T>& other) const {
    return x * other.x + y * other.y + z * other.z;
}

template <class T> inline vgVector3<T> vgVector3<T>::operator /(T divisor) const {
    return vgVector3<T>(x / divisor, y / divisor, z / divisor);
}

template <class T> inline vgVector3<T>& vgVector3<T>::operator +=(const vgVector3<T>& other) {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
}

template <class T> inline vgVector3<T>& vgVector3<T>::operator -=(const vgVector3<T>& other) {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
}

template <class T> inline vgVector3<T>& vgVector3<T>::operator *=(T factor) {
    x *= factor;
    y *= factor;
    z *= factor;
    return *this;
}

template <class T> inline vgVector3<T>& vgVector3<T>::operator /=(T divisor) {
    x /= divisor;
    y /= divisor;
    z /= divisor;
    return *this;
}

template <class T> inline vgVector3<T> vgVector3<T>::operator -() const {
    return vgVector3<T>(-x, -y, -z);
}

template <class T> inline vgVector3<T> operator *(T factor, const vgVector3<T>& other) {
    return vgVector3<T>(factor * other.x, factor * other.y, factor * other.z);
}

template <class T> inline vgVector3<T> vgVector3<T>::
cross(const vgVector3<T>& other) const {
    vgVector3<T> result;
    result.x = y * other.z - z * other.y;
    result.y = z * other.x - x * other.z;
    result.z = x * other.y - y * other.x;
    return result;
}

#ifdef __GNUC__   
typedef vgVector3<vg_ieee32> __vgVector3_vg_ieee32; // for inline optimization
typedef vgVector3<vg_ieee64> __vgVector3_vg_ieee64;
#endif

#endif /* _TYPES_VECTOR3_H_ */
