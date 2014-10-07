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
@Description  compiler dependent definitions

 */

#ifndef _TYPES_COMPILERS_H_
#define _TYPES_COMPILERS_H_

// SGI CC seems to have no own "marker"
#if (defined _OS_IRIX_ || defined _OS_IRIX64_) && !defined __GNUC__ && !defined __SGICC__
#define __SGICC__
#endif

// HP CC seems to have no own "marker"
#if defined _OS_HPUX_ && !defined __GNUC__ && !defined __HPCC__
#define __HPCC__
#endif

// Microsoft Visual C++ compiler
#if defined _MSC_VER && !defined __MSCC__
#define __MSCC__
#endif

// HP CC doesn't support "signed" and "volatile" declarations
#ifdef __HPCC__
#define signed
#define volatile
#endif

// only GCC and SGI CC (in a higher version!) seems to know "bool"
#if !defined __GNUC__ && !defined __SGICC__
// SGI CC in a newer version knows "_BOOL" as marker for bool type
#if !defined __SGICC__ || !defined _BOOL
#define HAS_NO_BOOL
#endif
#endif

#if defined HAS_NO_BOOL
#define bool vg_boolean_t
#define true vg_true_value
#define false vg_false_value

typedef int bool;
static const bool true = 1;
static const bool false = 0;

#undef HAS_NO_BOOL
#endif

#endif /* _TYPES_COMPILERS_H_ */
