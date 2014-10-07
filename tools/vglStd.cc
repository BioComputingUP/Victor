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
@Description  system independent data types

 */

static const char* __rcsid__ = "@(#) $Id: vglStd.cc,v 1.5 2008-05-09 13:58:51 biocomp Exp $";
static const char __use_rcsid__ = (&__use_rcsid__ - __rcsid__);

#include <vglStd.h>

// I modified this include because for non WIN32 platforms,
// I also need "float.h" and "values.h" has to be replaced
// by "limits.h"
#include <limits.h>
#include <float.h>
#define MAXFLOAT FLT_MAX
#define MAXDOUBLE DBL_MAX

template<> const char* vgClassInfo<vg_int8>::name = "vg_int8";
template<> const size_t vgClassInfo<vg_int8>::size = sizeof (vg_int8);
template<> const bool vgClassInfo<vg_int8>::numeric = true;
template<> const bool vgClassInfo<vg_int8>::ordered = true;
template<> const vg_int8 vgClassInfo<vg_int8>::min = -0x80;
template<> const vg_int8 vgClassInfo<vg_int8>::max = 0x7f;

template<> const char* vgClassInfo<vg_uint8>::name = "vg_uint8";
template<> const size_t vgClassInfo<vg_uint8>::size = sizeof (vg_uint8);
template<> const bool vgClassInfo<vg_uint8>::numeric = true;
template<> const bool vgClassInfo<vg_uint8>::ordered = true;
template<> const vg_uint8 vgClassInfo<vg_uint8>::min = 0;
template<> const vg_uint8 vgClassInfo<vg_uint8>::max = 0xff;

template<> const char* vgClassInfo<vg_int16>::name = "vg_int16";
template<> const size_t vgClassInfo<vg_int16>::size = sizeof (vg_int16);
template<> const bool vgClassInfo<vg_int16>::numeric = true;
template<> const bool vgClassInfo<vg_int16>::ordered = true;
template<> const vg_int16 vgClassInfo<vg_int16>::min = -0x8000;
template<> const vg_int16 vgClassInfo<vg_int16>::max = 0x7fff;

template<> const char* vgClassInfo<vg_uint16>::name = "vg_uint16";
template<> const size_t vgClassInfo<vg_uint16>::size = sizeof (vg_uint16);
template<> const bool vgClassInfo<vg_uint16>::numeric = true;
template<> const bool vgClassInfo<vg_uint16>::ordered = true;
template<> const vg_uint16 vgClassInfo<vg_uint16>::min = 0;
template<> const vg_uint16 vgClassInfo<vg_uint16>::max = 0xffff;

template<> const char* vgClassInfo<vg_int32>::name = "vg_int32";
template<> const size_t vgClassInfo<vg_int32>::size = sizeof (vg_int32);
template<> const bool vgClassInfo<vg_int32>::numeric = true;
template<> const bool vgClassInfo<vg_int32>::ordered = true;
template<> const vg_int32 vgClassInfo<vg_int32>::min = (vg_int32) 0x80000000;
template<> const vg_int32 vgClassInfo<vg_int32>::max = 0x7fffffff;

template<> const char* vgClassInfo<vg_uint32>::name = "vg_uint32";
template<> const size_t vgClassInfo<vg_uint32>::size = sizeof (vg_uint32);
template<> const bool vgClassInfo<vg_uint32>::numeric = true;
template<> const bool vgClassInfo<vg_uint32>::ordered = true;
template<> const vg_uint32 vgClassInfo<vg_uint32>::min = 0;
template<> const vg_uint32 vgClassInfo<vg_uint32>::max = 0xffffffff;

template<> const char* vgClassInfo<vg_ieee32>::name = "vg_ieee32";
template<> const size_t vgClassInfo<vg_ieee32>::size = sizeof (vg_ieee32);
template<> const bool vgClassInfo<vg_ieee32>::numeric = true;
template<> const bool vgClassInfo<vg_ieee32>::ordered = true;
template<> const vg_ieee32 vgClassInfo<vg_ieee32>::min = -MAXFLOAT;
template<> const vg_ieee32 vgClassInfo<vg_ieee32>::max = MAXFLOAT;

template<> const char* vgClassInfo<vg_ieee64>::name = "vg_ieee64";
template<> const size_t vgClassInfo<vg_ieee64>::size = sizeof (vg_ieee64);
template<> const bool vgClassInfo<vg_ieee64>::numeric = true;
template<> const bool vgClassInfo<vg_ieee64>::ordered = true;
template<> const vg_ieee64 vgClassInfo<vg_ieee64>::min = -MAXDOUBLE;
template<> const vg_ieee64 vgClassInfo<vg_ieee64>::max = MAXDOUBLE;




