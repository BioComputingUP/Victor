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
@Description  system independent byte order definitions

*/

#ifndef _TYPES_ENDIAN_H_
#define _TYPES_ENDIAN_H_

#if defined _OS_Linux_
#include <endian.h>
#endif
#if defined _OS_IRIX_ || defined _OS_IRIX64_
#include <sys/endian.h>
#endif
#if defined _OS_OSF1_
#include <machine/endian.h>
#endif
#if defined _OS_AIX_
#include <sys/machine.h>
#endif
#if defined __APPLE__
#include <sys/param.h>
//#else
//#include <i386/endian.h>
#endif

#if !defined BIG_ENDIAN
#define BIG_ENDIAN 4321
#endif

#if !defined LITTLE_ENDIAN
#define LITTLE_ENDIAN 1234
#endif

#if !defined BYTE_ORDER
#if defined __i386 || defined WIN32 || defined __x86_64
#define BYTE_ORDER LITTLE_ENDIAN
#endif
#if defined __sparc__ || defined __mips__
#define BYTE_ORDER BIG_ENDIAN
#endif
#endif

#if !defined BYTE_ORDER || !defined LITTLE_ENDIAN || !defined BIG_ENDIAN
#error FATAL compile error: byte order not defined
#endif

#endif /* _TYPES_ENDIAN_H_ */
