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
 *@Class:             LigandSet 
 *@Project Name:      Victor
 *@ATTENTION:         This class is *NOT* finished yet.
 *                     There are many methods which have only a head
 *                     declaration and no code inside.
 *                     2012-Francesco Lovo: insert, delete, getLigand, load and save
 *                     methods are completed
 */
//Includes:
#include <LigandSet.h>
#include <Debug.h>
#include <IntCoordConverter.h>
#include <limits.h>
#include <float.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;

// CONSTRUCTORS/DESTRUCTOR:

LigandSet::LigandSet() : Polymer(1, 1), startOffset(1), gaps()//, startAtomOffset(0)
{
    PRINT_NAME;
}

LigandSet::LigandSet(const LigandSet& orig) {
    PRINT_NAME;
    this->copy(orig);
}

LigandSet::~LigandSet() {
    PRINT_NAME;
}

void
LigandSet::copy(const LigandSet& orig) {
    PRINT_NAME;
    Polymer::copy(orig);
}

LigandSet*
LigandSet::clone() {
    LigandSet* tmp = new LigandSet;
    tmp->copy(*this);
    return tmp;
}

//PREDICATES

/**
 *@Description  Predicate used to determine if the PDB ligand number is valid (false)  or a gap (true).
 *@param  int
 *@return  bool
 */
bool LigandSet::isGap(int index) {

    if (index < startOffset)
        return true;

    for (unsigned int i = 0; i < gaps.size(); i++) {
        if (index == gaps[i])
            return true;
    }

    if (index >= static_cast<int> (sizeLigand() + startOffset + gaps.size()))
        return true;
    else
        return false;
}

//MODIFIERS:

/**
 *@Description  Modifier used to add a gap at the PDB ligand number.
 *@param  int
 *@return  void
 */
void
LigandSet::addGap(int index) {
    if (index < startOffset) {
        startOffset = index;
        return;
    }

    if ((gaps.size() == 0) || (index > gaps[gaps.size() - 1]))
        gaps.push_back(index);
    else
        for (unsigned int i = gaps.size() - 1; i >= 0; i--)
            if (index == gaps[i])
                break;
            else
                if (index < gaps[i]) {
                gaps.insert(gaps.begin() + i, 1, index);
                break;
            }
}

// OPERATORS:

LigandSet&
        LigandSet::operator=(const LigandSet& orig) {
    PRINT_NAME;
    if (&orig != this)
        copy(orig);
    return *this;
}

void
LigandSet::insertComponent(Component* g) {
    if (g->hasSuperior())
        ERROR("Ligand does have a superior", exception);
    if (g->getClassName() != "Ligand")
        ERROR("LigandSet can only insert a Ligand", exception);
    Polymer::insertComponent(g);
    setModified();
}


// HELPERS:

void
LigandSet::resetBoundaries() {
    vgVector3<double> tmpV(DBL_MAX - 1, DBL_MAX - 1, DBL_MAX - 1);
    lowerBound = tmpV;
    upperBound = -tmpV;
    unsigned int sizeL = sizeLigand();
    for (unsigned int i = 0; i < sizeL; i++)
        for (unsigned int j = 0; j < 3; j++) {
            if (getLigand(i).getLowerBound()[j] < lowerBound[j])
                lowerBound[j] = getLigand(i).getLowerBound()[j];
            if (getLigand(i).getUpperBound()[j] > upperBound[j])
                upperBound[j] = getLigand(i).getUpperBound()[j];
        }
}




