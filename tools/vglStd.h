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

#ifndef _TYPES_STD_H_
#define _TYPES_STD_H_

#define OLD_IOSTREAMS		// all systems seem to have "oldstyle" IO streams

#include <stddef.h>		// for "ptrdiff_t" and others
#include <vglEndian.h>	// for little/big endian defines
#include <vglCompilers.h>	// for compiler dependent settings

typedef signed char vg_int8; // signed integer with exactly 8 bits
typedef unsigned char vg_uint8; // unsigned integer with exactly 8 bits
typedef signed short vg_int16; // signed integer with exactly 16 bits
typedef unsigned short vg_uint16; // unsigned integer with exactly 16 bits
typedef signed int vg_int32; // signed integer with exactly 32 bits
typedef unsigned int vg_uint32; // unsigned integer with exactly 32 bits
typedef float vg_ieee32; // IEEE floating point with exactly 32 bits
typedef double vg_ieee64; // IEEE floating point with exactly 64 bits

extern const char* const vglVersion;

// simple structure to provide some information about classes.
// (for internal use only)

template <class T>
class vgClassInfo {
public:
    static const char* name; // name of the class
    static const size_t size; // size in bytes
    static const bool numeric; // true if +-*/ defined
    static const bool ordered; // true if comparison defined
    static const T min; // representable minimum
    static const T max; // representable maximum
};

#endif /* _TYPES_STD_H_ */
