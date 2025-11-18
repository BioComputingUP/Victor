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

#ifndef _VICTOR_PHYLO_PHYLO_
#define _VICTOR_PHYLO_PHYLO_


namespace Victor {
/** @brief Molecular PHYLOgenenesys package.
 * 
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 */
namespace Phylo {}
}

#include <AminoAcid.h>
#include <Sequence.h>
#include <SubstitutionMatrix.h>

#include <DistanceMatrix.h>
#include <DistanceMatrixBuilder.h>
#include <IdentityPercentage.h>
#include <LevenshteinDistance.h>
#include <FengDoolittleDistance.h>

#include <PhylogeneticTree.h>
#include <RootedTree.h>
#include <UnrootedTree.h>
#include <NewickParser.h>

#include <PhylogeneticAlgorithm.h>
#include <ClusteringAlgorithm.h>
#include <UPGMA.h>
#include <FitchMargoliash.h>
#include <NJ.h>

#include <MultipleAlignment.h>

#include <MultipleAlignmentAlgorithm.h>
#include <FengDoolittle.h>
#include <ClustalW.h>

#include <Visitor.h>

#endif  // _VICTOR_PHYLO_PHYLO_