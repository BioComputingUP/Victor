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

#ifndef __SingleCluster_H__
#define __SingleCluster_H__

#include <Cluster.h>
#include <string>
#include "Sequence.h"

namespace Victor { namespace Phylogenesis{

	/** @brief    It represents a leaf of the clustering tree containing a Sequence object.
	 *  @author   Andrea Nasari (andrea.nasari\@studenti.unipd.it)
	 *
	 **/
	class SingleCluster : public Cluster{
	public:
		/// Constructor.
		SingleCluster(Sequence* sequence);

		/// Returns if is a single cluster (true)
		virtual bool isSingle();

		// ATTRIBUTES:
		Sequence* sequence;
	};

	inline SingleCluster::SingleCluster(Sequence* sequence) : sequence(sequence)
	{
		distance = 0.0;
	}

	inline bool SingleCluster::isSingle()
	{
		return true;
	}

}}

#endif
