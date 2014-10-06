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


#ifndef __Blosum_H__
#define __Blosum_H__

#include <iostream>
#include <string>
#include <Substitution.h>

/** @brief  This class implemenents a standard substitution matrix
 * 
* @Description  
* @This 
 **/
class Blosum : public Substitution 
{
public:
    
// CONSTRUCTORS:
  /// Default constructor.
  Blosum();
  /// Copy constructor.
  Blosum(const Blosum& orig);
  /// Constructor from istream.
  Blosum(istream& is);
  /// Destructor.
  virtual ~Blosum();
  
// OPERATORS:
  /// Assignment operator.
  Blosum& operator = (const Blosum& orig);
  /// Output operator.
  friend ostream& operator << (ostream& os, const Blosum& object);
  /// Input operator.
  friend istream& operator >> (istream& is, Blosum& object);

// PREDICATES:
  /// Implementation of abstract class method.
  virtual string getResidues() const { return residues; }
  /// Return dimension of matrix.
  virtual unsigned int size() const { return residuescores.size(); }

protected:
  
// MODIFIERS:
  /// Copies orig object to this object ("deep copy").
  virtual void copy(const Blosum& orig);

private:

// ATTRIBUTES:
  vector<vector<int> > residuescores; ///< Similarity scores.
  string residues;                    ///< Alphabet of allowed characters.

};

#endif /* __Blosum_H__ */

