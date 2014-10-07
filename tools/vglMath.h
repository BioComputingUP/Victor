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
@Description  mathematical tool functions

 */

#ifndef _TOOLS_MATH_H_
#define _TOOLS_MATH_H_

#include <std.h>
#include <math.h>		// provides fabs(double) and sometimes abs(double) 

inline vg_int8 vg_abs(vg_int8 x) {
    return (x > 0) ? x : -x;
}

inline vg_uint8 vg_abs(vg_uint8 x) {
    return x;
}

inline vg_int16 vg_abs(vg_int16 x) {
    return (x > 0) ? x : -x;
}

inline vg_uint16 vg_abs(vg_uint16 x) {
    return x;
}

inline vg_int32 vg_abs(vg_int32 x) {
    return (x > 0) ? x : -x;
}

inline vg_uint32 vg_abs(vg_uint32 x) {
    return x;
}

inline vg_ieee32 vg_abs(vg_ieee32 x) {
    return (x > 0) ? x : -x;
}

inline vg_ieee64 vg_abs(vg_ieee64 x) {
    return (x > 0) ? x : -x;
}

#if !defined M_PI
#define M_PI                3.14159265358979323846
#endif
#if !defined M_PI_2
#define M_PI_2              1.57079632679489661923
#endif
#if !defined M_PI_4
#define M_PI_4              0.78539816339744830962
#endif

#endif /* _TOOLS_MATH_H_ */
