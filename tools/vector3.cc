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

static const char* __rcsid__ = "@(#) $Id: vector3.cc,v 1.4 2008-05-09 13:58:51 biocomp Exp $"; static const char __use_rcsid__ = (&__use_rcsid__ - __rcsid__);


#include <vector3.h>

#if !defined __GNUC__		// gcc has problems with global initializers (!?)
const char* vgClassInfo<vgVector3<vg_ieee32> >::name = "vgVector3<vg_ieee32>";
const size_t vgClassInfo<vgVector3<vg_ieee32> >::size = sizeof(vgVector3<vg_ieee32>);
const bool vgClassInfo<vgVector3<vg_ieee32> >::numeric = true;
const bool vgClassInfo<vgVector3<vg_ieee32> >::ordered = false;
const vgVector3<vg_ieee32> vgClassInfo<vgVector3<vg_ieee32> >::min = vgVector3<vg_ieee32>(0.0f,0.0f,0.0f);
const vgVector3<vg_ieee32> vgClassInfo<vgVector3<vg_ieee32> >::max = vgVector3<vg_ieee32>(0.0f,0.0f,0.0f);
const char* vgClassInfo<vgVector3<vg_ieee64> >::name = "vgVector3<vg_ieee64>";
const size_t vgClassInfo<vgVector3<vg_ieee64> >::size = sizeof(vgVector3<vg_ieee64>);
const bool vgClassInfo<vgVector3<vg_ieee64> >::numeric = true;
const bool vgClassInfo<vgVector3<vg_ieee64> >::ordered = false;
const vgVector3<vg_ieee64> vgClassInfo<vgVector3<vg_ieee64> >::min = vgVector3<vg_ieee64>(0.0,0.0,0.0);
const vgVector3<vg_ieee64> vgClassInfo<vgVector3<vg_ieee64> >::max = vgVector3<vg_ieee64>(0.0,0.0,0.0);
#endif

#if defined __GNUC__		// for inlining of template class methods (!?)
typedef vgVector3<vg_ieee32> __vgVector3_vg_ieee32;
typedef vgVector3<vg_ieee64> __vgVector3_vg_ieee64;
#endif

#if defined __GNUC__
template class vgVector3<vg_ieee32>;
template class vgVector3<vg_ieee64>;
#endif

#if defined __MSCC__
// template vgVector3<vg_ieee32>; <- not neccessary: classinfo above gives instantiation!
// template vgVector3<vg_ieee64>; <- not neccessary: classinfo above gives instantiation!
#endif

#if defined __DECCXX 
#pragma define_template vgVector3<vg_ieee32>
#pragma define_template vgVector3<vg_ieee64>
#endif

#if defined __SGICC__
#pragma instantiate vgVector3<vg_ieee32>
#pragma instantiate vgVector3<vg_ieee64>
#endif
