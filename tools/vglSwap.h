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
@Description functions for value exchange and byte swapping

 */

#ifndef _TOOLS_SWAP_H_
#define _TOOLS_SWAP_H_

#include <std.h>

namespace Victor {

    /* exchange values */

    template <class T> inline void rotate(T& a, T& b, T&c) {
        T tmp = a;
        a = b;
        b = c;
        c = tmp;
    }

    template <class T> inline void swap(T& a, T& b) {
        T tmp = a;
        a = b;
        b = tmp;
    }

    /* little/big endian conversion */

    inline void swap(vg_uint16& a) {
        a = (a << 8) | (a >> 8);
    }

    inline void swap(vg_int16& a) {
        a = (vg_int16) ((vg_uint16) a << 8) | ((vg_uint16) a >> 8);
    }

    inline void swap(vg_uint32& a) {
        vg_uint16* p = (vg_uint16*) (&a);
        *p = (*p << 8) | (*p >> 8);
        p++;
        *p = (*p << 8) | (*p >> 8);
        a = (a >> 16) | (a << 16);
    }

    inline void swap(vg_int32& a) {
        vg_uint16* p = (vg_uint16*) (&a);
        *p = (*p << 8) | (*p >> 8);
        p++;
        *p = (*p << 8) | (*p >> 8);
        a = (vg_int32) (((vg_uint32) a >> 16) | ((vg_uint32) a << 16));
    }

    template <class T> inline void swap(T& a) {
        switch (sizeof (T)) {
            case 2:
                swap(*(vg_uint16*) & a);
                break;
            case 4:
                swap(*(vg_uint32*) & a);
                break;
            default:
                break;
        }
    }

} // end-namespace

#endif /* _TOOLS_SWAP_H_ */
