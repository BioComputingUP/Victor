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

#ifndef __MsaGapFunction_h___
#define __MsaGapFunction_h___

namespace Victor { namespace Phylogenesis {

	/** @brief Base class for Multiple Sequence Alignment (MSA) gap functions.
	 *  @author   Andrea Nasari (andrea.nasari\@studenti.unipd.it)
	 *
	 **/
	class MsaGapFunction{
	public:

		/// Default constructor.
		MsaGapFunction(){}

		/// Destructor.
		virtual ~MsaGapFunction(){}

		/// Return open gap penalty for horizontal seq at position p.
		virtual double getHorizontalOpenPenalty(int pos) = 0;

		/// Return extension gap penalty for  horizontal seq at position p.
		virtual double getHorizontalExtensionPenalty(int pos) = 0;

		/// Return open gap penalty for vertical seq at position p.
		virtual double getVerticalOpenPenalty(int pos) = 0;

		/// Return extension gap penalty for vertical seq at position p.
		virtual double getVerticalExtensionPenalty(int pos) = 0;
	};

}}

#endif
