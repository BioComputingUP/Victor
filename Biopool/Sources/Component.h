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
 

#ifndef _COMPONENT_H_
#define _COMPONENT_H_
 
// Includes:
#include <vector>
#include <string>
#include <Bond.h>
#include <Loader.h>
#include <Saver.h>
#include <vector3.h>
#include <Visitor.h>

// Global constants, typedefs, etc. (to avoid):


namespace Biopool {
/**@brief Base class for composite structures. 
 * 
*@Description Implementing the Composite pattern.
*@This 
 **/
class Component : public Bond
{
public: 

// CONSTRUCTORS/DESTRUCTOR:
  Component(unsigned int mI = 2, unsigned int mO = 2);
  Component(const Component& orig);
  virtual ~Component();

// PREDICATES:
  unsigned int size() const;
  virtual string getClassName() const = 0;
  unsigned int getDepth() const;

  virtual vgVector3<double> getLowerBound(double dist = 0.0);
  virtual vgVector3<double> getUpperBound(double dist = 0.0);
  bool collides(Component& other, double dist = 0.0);  
  // is other within dist of this' bounding box?

  virtual vgVector3<double> getTrans() const = 0;
  virtual vgMatrix3<double> getRot() const = 0;

  Component& getSuperior();
  bool hasSuperior() const;
  bool inSync(); // coords in-sync with changes?

  virtual void save(Saver& s) = 0;  // data saver

// MODIFIERS:

  virtual void connectIn(Component* c, unsigned int offset = 1);
  virtual void connectOut(Component* c, unsigned int offset = 1);
  virtual Component* unconnectIn();
  virtual Component* unconnectOut();

  virtual void insertComponent(Component* c) = 0;
  virtual void removeComponent(Component* c) = 0;
  virtual void deleteComponent(Component* c) = 0;

  void copy(const Component& orig);
  virtual Component* clone() = 0;
  virtual void load(Loader& l) = 0;  // data loader

  virtual void setTrans(vgVector3<double> t) = 0;
  virtual void addTrans(vgVector3<double> t) = 0;
  virtual void setRot(vgMatrix3<double> r) = 0;
  virtual void addRot(vgMatrix3<double> r) = 0;
  virtual void sync() = 0; // synchronize coords with structure
  virtual void setModified();

  void setSuperior(Component* c);  
  // ought to be 'protected' but the compiler will not handle 
  // it correctly, if it is ?!

  virtual void acceptCalculator(EnergyVisitor* v) = 0;
  virtual void acceptOptimizer(OptimizationVisitor* v) = 0;

// OPERATORS:
  Component& operator=(const Component& orig);

protected:

// HELPERS:
  virtual void resetBoundaries() = 0;

  void setLowerBound(vgVector3<double> _lb);
  void setUpperBound(vgVector3<double> _ub);

// ATTRIBUTES:

  Component* superior;
  vector<Component*> components;
  vgVector3<double> lowerBound;  // lower bound for bounding box
  vgVector3<double> upperBound;  // upper bound for bounding box
  bool modified;                 // component modified?
  
private:
  
};


// ---------------------------------------------------------------------------
//                                    Component
// -----------------x-------------------x-------------------x-----------------

// PREDICATES:
inline unsigned int 
Component::size() const
{
  return components.size();
}

inline unsigned int 
Component::getDepth() const
{
  if (superior != NULL)
    return superior->getDepth() + 1;
  else
    return 0;
}

inline bool 
Component::collides(Component& other, double dist)
{
  return ( (lowerBound.x <= other.upperBound.x+dist) 
	   || (lowerBound.y <= other.upperBound.y+dist)
	   || (lowerBound.z <= other.upperBound.z+dist)
	   || (upperBound.x >= other.lowerBound.x-dist) 
	   || (upperBound.y >= other.lowerBound.y-dist) 
	   || (upperBound.z >= other.lowerBound.z-dist) );  
}  

inline Component& 
Component::getSuperior()
{
  if (superior != NULL)
    return *superior;
  else
    {
      DEBUG_MSG("Component::getSuperior() Warning: no superior found.");
      return *this;
    }
}

inline bool
Component::hasSuperior() const
{
  return (superior != NULL);
}

inline bool 
Component::inSync()
{
  return (!modified && ( !hasSuperior() || (getSuperior().inSync()))) ;
}

// MODIFIERS:

inline void Component::setModified()
{
  if (modified)
    return;
  modified = true;
  if (hasSuperior())
    getSuperior().setModified();
}

// OPERATORS:

// HELPERS:

inline void 
Component::setSuperior(Component* c)
{
  superior = c;
}

inline void 
Component::setLowerBound(vgVector3<double> _lb)
{
  lowerBound = _lb;
}

inline void 
Component::setUpperBound(vgVector3<double> _ub)
{
  upperBound = _ub;
}

} // namespace
#endif //_COMPONENT_H_

