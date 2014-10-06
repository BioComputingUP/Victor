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

// Includes:
#include <Monomer.h>
#include <Debug.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;

// CONSTRUCTORS/DESTRUCTOR:
Monomer::Monomer(unsigned int mI, unsigned int mO) : Component(mI, mO) 
{ PRINT_NAME; }

Monomer::Monomer(const Monomer& orig)
{
  PRINT_NAME;
  this->copy(orig);  
}

Monomer::~Monomer()
{ 
  PRINT_NAME;
  if ( hasSuperior() )
    {
      getSuperior().removeComponent(this);
      setSuperior(NULL);
    }
} 


// PREDICATES:

// MODIFIERS:

void 
Monomer::copy(const Monomer& orig)
{
  PRINT_NAME; 
  

}

Component* 
Monomer::clone()
{
  return new Monomer;
}

// OPERATORS:

// HELPERS:

 
