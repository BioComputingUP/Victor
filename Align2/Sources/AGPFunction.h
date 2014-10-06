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
#ifndef __AGPFunction_H__
#define __AGPFunction_H__

#include <GapFunction.h>

namespace Biopool
{
/** @brief Implement AGP (Affine Gap Penalty) function.
 * 
* @Description  
* @This 
 **/
class AGPFunction : public GapFunction
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	AGPFunction() : o(12.00), e(3.00)
	{ }

	/// Constructor assigning o and e.
	AGPFunction(double o, double e) : o(o), e(e)
	{ }

	/// Copy constructor.
	AGPFunction(const AGPFunction &orig) : GapFunction(orig)
	{
		copy(orig);
	}

	/// Destructor.
	virtual ~AGPFunction()
	{ }


// OPERATORS:

	/// Assignment operator.
	AGPFunction& operator = (const AGPFunction &orig);


// PREDICATES:

	/// Return open gap penalty for template position p.
	virtual double getOpenPenalty(int p);

	/// Return extension gap penalty for template position p.
	virtual double getExtensionPenalty(int p);


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const AGPFunction &orig);

	/// Construct a new "deep copy" of this object.
	virtual AGPFunction* newCopy();

	/// Set open gap penalty.
	virtual void setOpenPenalty(double pen);

	/// Set extension gap penalty.
	virtual void setExtensionPenalty(double pen);


protected:


private:

// ATTRIBUTES:

	double o;    ///< Open gap penalty.
	double e;    ///< Extension gap penalty.

};

// -----------------------------------------------------------------------------
//                                 AGPFunction
// -----------------------------------------------------------------------------

// OPERATORS:

inline AGPFunction&
AGPFunction::operator = (const AGPFunction &orig)
{
	if (&orig != this)
		copy(orig);
	POSTCOND((orig == *this), exception);
	return *this;
}


// PREDICATES:

inline double
AGPFunction::getOpenPenalty(int p)
{
	return o;
}


inline double
AGPFunction::getExtensionPenalty(int p)
{
	return e;
}


// MODIFIERS:

inline void
AGPFunction::copy(const AGPFunction &orig)
{
	GapFunction::copy(orig);
	o = orig.o;
	e = orig.e;
}


inline AGPFunction*
AGPFunction::newCopy()
{
	AGPFunction *tmp = new AGPFunction(*this);
	return tmp;
}


inline void
AGPFunction::setOpenPenalty(double pen)
{
	o = pen;
}


inline void
AGPFunction::setExtensionPenalty(double pen)
{
	e = pen;
}

} // namespace

#endif
