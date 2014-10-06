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


#ifndef _MONOMER_H_
#define _MONOMER_H_
 
// Includes:
#include <Component.h>
#include <Debug.h>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
/**@brief Base class for components without elements  
 * 
*@Description  
*@This 
 * */
class Monomer : public Component
{
public: 

// CONSTRUCTORS/DESTRUCTOR:
  Monomer(unsigned int mI = 1, unsigned int mO = 1); 
  Monomer(const Monomer& orig);
  ~Monomer();
  
// PREDICATES:
  virtual string getClassName() const
    { return "Monomer"; };
  virtual vgVector3<double> getTrans() const
    { ERROR(" Monomer::getTrans is not a viable method.", exception); };
  virtual vgMatrix3<double> getRot() const
    { ERROR("Monomer::getRot is not a viable method.", exception); };

  virtual void save(Saver& s)
    { ERROR("Monomer::save is not a viable method.", exception); };

// MODIFIERS:
  void insertComponent(Component* c)
    { ERROR("Monomer::insertComponent is not a viable method.", exception); };
 
  void removeComponent(Component* c)
    { ERROR("Monomer::removeComponent is not a viable method.", exception); }; 

  void deleteComponent(Component* c)
    { ERROR("Monomer::deleteComponent is not a viable method.", exception); }; 

  void copy(const Monomer& orig);
  virtual Component* clone();

  virtual void load(Loader& l)
    { ERROR("Monomer::load is not a viable method.", exception); };
  virtual void setTrans(vgVector3<double> t)
    { ERROR("Monomer::setTrans is not a viable method.", exception); };
  virtual void addTrans(vgVector3<double> t)
    { ERROR("Monomer::addTrans is not a viable method.", exception); };
  virtual void setRot(vgMatrix3<double> r)
    { ERROR("Monomer::setRot is not a viable method.", exception); };
  virtual void addRot(vgMatrix3<double> r)
    { ERROR("Monomer::addRot is not a viable method.", exception); };
  virtual void sync()   // synchronize coords with structure
    { ERROR("Monomer::sync is not a viable method.", exception); };

  virtual void acceptCalculator(EnergyVisitor* v)
    { ERROR("Monomer::acceptCalculator is not a viable method.", exception); };
  virtual void acceptOptimizer(OptimizationVisitor* v)
    { ERROR("Monomer::acceptOptimizer is not a viable method.", exception); };

// OPERATORS:

protected:
  
// HELPERS:
  virtual void resetBoundaries()
    { ERROR("Monomer::resetBoundaries is not a viable method.", exception); };

private:
  
};


// ---------------------------------------------------------------------------
//                                    Monomer
// -----------------x-------------------x-------------------x-----------------

// PREDICATES:

// MODIFIERS:

// OPERATORS:

} // namespace
#endif //_MONOMER_H_
