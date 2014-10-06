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
*@Class:             Protein
*@Base Class(es):    Polymer
*@Derived Class(es): -
*@Containing:        Spacer, Ligandset
*@Project Name:      Victor
*@Description:       This class is a container of Polymers, each one
*                     storing a Spacer and (eventually) a LigandSet for each
*                     valid chain in the PDB input file.
*  Warning:           Effective construction of Polymer with these characteristics
*                     is made by load method
*/

// Includes:
#include <Protein.h>
#include <iostream>
using namespace std;
using namespace Biopool;

// CONSTRUCTORS/DESTRUCTOR:
Protein::Protein() : Polymer(1,1) 
{ PRINT_NAME; }

Protein::Protein(const Protein& orig)           
{
  PRINT_NAME;
  this->copy(orig);  
}

Protein::~Protein()
{ PRINT_NAME;  } 

// PREDICATES:
Spacer*                     //return a pointer to the Spacer with the correct chainID
Protein::getSpacer(char c)
{
    unsigned int n = getChainNum(c);
    return getSpacer(n);
}   
/**
 *@Description 
 *@param 
 */
Polymer&
Protein::getPolymer(unsigned int n)
{
  if ( n > components.size() - 1 )
    ERROR("Index out of range",exception);
  return *(dynamic_cast<Polymer*> (components[n]));
}
/**
 *@Description 
 *@param 
 */
Spacer*
Protein::getSpacer(unsigned int n)
{
    if ( n > sizeProtein() - 1 )
      ERROR("Index out of range",exception);
    Polymer& p=getPolymer(n);
    return &(dynamic_cast<Spacer&>(p[0]));
}
/**
 *@Description 
 *@param 
 */
LigandSet*   //return a pointer to the ligandSet (if exists) with the correct chainID
Protein::getLigandSet(char c)
{
    unsigned int n = getChainNum(c);
    return getLigandSet(n);
}   

LigandSet*
Protein::getLigandSet(unsigned int n)
{
    if ( n > sizeProtein() - 1 )
      ERROR("Index out of range",exception);
    Polymer& p = getPolymer(n);
    if (p.size()<2)
        return NULL;     //No Ligands for this chain
    return &(dynamic_cast<LigandSet&>(p[1]));
}

unsigned int
Protein::getChainNum (char c)                                           
{
    for (unsigned int i=0;i<chains.size();i++)
        if (chains[i]==c)
            return i;
    ERROR("Chain not found",exception);
}

char
Protein::getChainLetter (unsigned int i)                                
{
    if (i>chains.size())
        ERROR("Index out of range",exception);
    return chains[i];
}

// MODIFIERS:
void
Protein::removeComponent(Component* c)
{
  Polymer::removeComponent(c);
  setModified();
}
void
Protein::deleteComponent(Component* c)
{
  Polymer::deleteComponent(c);
  setModified();
}
void                 //not finished because Component::copy is not finished too
Protein::copy(const Protein& orig)
{
  PRINT_NAME; 
  chains = orig.chains;
  Polymer::copy(orig);      
}

Protein*                                                          
Protein::clone()
{
    Protein* tmp = new Protein;
    tmp->copy(*this);
    return tmp;
}
/**
 *@Description Insert a Polymer in the Protein at the back side.
*               It is related to a single chain, and must
*               contain a Spacer and (eventually) a LigandSet
 *@param component reference
 *@return void 
 */ void 
Protein::insertComponent(Component* p) 
{  
  if ( p->hasSuperior() )
    ERROR("Component does have a superior",exception);
  if ( p->getClassName() != "Polymer")
    ERROR("Protein can insert directly only Polymer",exception);
   Polymer::insertComponent(p);
} 

// OPERATORS:
      /**
 *@Description 
 *@param 
 */
Protein& 
Protein::operator=(const Protein& orig)
{
  PRINT_NAME;
  if (&orig != this)
    copy(orig);
  return *this;
}
// HELPERS:
/**
 *@Description 
 *@param 
 */
void 
Protein::printComponents() 
{
    ERROR("Not finished yet",exception);
}
