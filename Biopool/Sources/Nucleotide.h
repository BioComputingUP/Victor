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


#ifndef _NUCLEOTIDE_H_
#define _NUCLEOTIDE_H_

//unuseful comment
// Includes:
#include <Group.h>
#include <IntCoordConverter.h>
#include <Visitor.h>
#include <NucleotideCode.h>


// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
/**@brief class implements a simple nucleotide.
 * 
*@Description Includes methods that allow to get and set angles, connections, bonds, side chain info, etc.
*@This NB Angles are in degrees.
 * */
class Nucleotide : public Group
{
public: 

  // CONSTRUCTORS/DESTRUCTOR:
  Nucleotide();
  Nucleotide(const Nucleotide& orig);
  virtual ~Nucleotide(); 
  
  // PREDICATES:
  virtual string getClassName() const
    { return "Nucleotide"; }
  virtual unsigned int getCode() const;
   unsigned int size() const;

   unsigned int sizeNucleotide() const;
 
  
  bool isMember(const AtomCode& ac) const;

  void save(Saver& s);   // data saver

  // MODIFIERS:
  

  void copy(const Nucleotide& orig);
  void setType(string _name);
  
  void setBonds(double NToCaLen, double CaToCLen, double CToOLen, double atCaToCAng, double CaToOAng); 
  bool setBondsFromPdbCode(bool connect, Nucleotide* prev, bool permissive);
  
  virtual void sync(); // synchronize coords with structure
  virtual void setModified();

  void load(Loader& l);  // data loader

  virtual Component* clone();
   
  // OPERATORS:
  bool operator==(const Nucleotide& other) const;
  bool operator!=(const Nucleotide& other) const;
  Nucleotide& operator=(const Nucleotide& orig);
  Atom& operator[](unsigned int n);
  const Atom& operator[](unsigned int n) const;
  Atom& operator[](const AtomCode& ac);
  const Atom& operator[](const AtomCode& ac) const;

protected:

private:

  // ATTRIBUTES:

  NucleotideCode type;	         // Nucleotide type
  // vector atoms (inherited from Group) contains the backbone
  IntCoordConverter icc; // should be static, but compiler won't accept it?
};

// ---------------------------------------------------------------------------
//                                    Nucleotide
// -----------------x-------------------x-------------------x-----------------

// PREDICATES:
inline unsigned int 
Nucleotide::getCode() const
{
  return nucleotideThreeLetterTranslator(id);
}


inline unsigned int 
Nucleotide::size() const 
{
  return Group::size();
}

inline unsigned int 
 Nucleotide::sizeNucleotide() const
 {
   return Group::size();
 }

inline bool 
Nucleotide::isMember(const AtomCode& ac) const
{
  return (pGetAtom(ac) != NULL);
}

inline void 
Nucleotide::save(Saver& s)
{
  s.saveNucleotide(*this);
}


// MODIFIERS:

inline void 
Nucleotide::setType(string _name)
{ 
  type = nucleotideThreeLetterTranslator(_name);
  id.setName(_name);
}


inline void 
Nucleotide::load(Loader& l)
{
  l.loadNucleotide(*this);
  resetBoundaries();
}

inline void 
Nucleotide::setModified()
{
  Group::setModified();
}




// OPERATORS:

inline bool 
Nucleotide::operator==(const Nucleotide& other) const
{
  return (dynamic_cast<const Identity*>(this)) == 
    dynamic_cast<const Identity*>(&other);
}
  

inline bool 
Nucleotide::operator!=(const Nucleotide& other) const
{
  return (dynamic_cast<const Identity*>(this)) != 
    dynamic_cast<const Identity*>(&other);
}


inline Atom& 
Nucleotide::operator[](unsigned int n)
{
  PRECOND(n < this->size(), exception);
  return Group::getAtom(n);
}

inline const Atom& 
Nucleotide::operator[](unsigned int n) const
{
  PRECOND(n < this->size(), exception);
  return Group::getAtom(n);
}

inline Atom& 
Nucleotide::operator[](const AtomCode& ac)
{
  Atom* a = pGetAtom(ac);
  INVARIANT(a != NULL, exception);
  return *a;
}

inline const Atom& 
Nucleotide::operator[](const AtomCode& ac) const
{
  Atom* a = pGetAtom(ac);
  INVARIANT(a != NULL, exception);
  return *a;
}
}
#endif //_NUCLEOTIDE_H_
