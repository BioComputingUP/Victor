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


#ifndef _ATOM_H_
#define _ATOM_H_

// Includes:
#include <AtomCode.h>
#include <vector3.h>
#include <matrix3.h>
#include <SimpleBond.h>

namespace Victor { namespace Biopool { 


    // Global constants, typedefs, etc. (to avoid):

    class Group;

    /**@brief Implements a simple atom type.
     * 
     *  Includes methods that allow you to get and set type, bind, unbind, coordinates , code,  etc.
     * NB: Angles are in degrees.
     **/
    class Atom : public SimpleBond {
    public:

        // CONSTRUCTORS/DESTRUCTOR:
        Atom(unsigned int mI = 1, unsigned int mO = 4);
        Atom(const Atom& orig);
        virtual ~Atom();

        // PREDICATES:
        AtomCode getCode() const;
        unsigned long getNumber() const;

        vgVector3<double> getCoords();

        double getBFac() {
            return Bfac;
        }
        
        double getOccupancy(){
            return Occupancy;
        }

        double distance(Atom& other);

        bool hasSuperior();
        Group& getSuperior();
        const Group& getSuperior() const;

        vgVector3<double> getTrans() const;
        vgMatrix3<double> getRot() const;
        virtual bool inSync(); // coords in-sync with changes?

        // MODIFIERS:
        void clear();
        void copy(const Atom& orig);
        void bindStructure(Atom& before, bool connect);
        // sets this relative to before if connect == true, otherwise like bindIn
        virtual void bindIn(SimpleBond& c);
        virtual void bindOut(SimpleBond& c);
        virtual void unbindIn(SimpleBond& c);
        virtual void unbindOut(SimpleBond& c);
        void setType(string _name);
        void setCode(AtomCode ac);
        void setNumber(unsigned long _number);

        void setCoords(double _x, double _y, double _z);
        void setCoords(vgVector3<double> c);

        void setBFac(double _b) {
            Bfac = _b;
        }
        
        void setOccupancy(double _o){
            Occupancy = _o;
        }

        void setTrans(vgVector3<double> t);
        void addTrans(vgVector3<double> t);
        void setRot(vgMatrix3<double> r);
        void addRot(vgMatrix3<double> r);
        virtual void sync(); // synchronize coords with structure

        virtual const Atom& getInBond(unsigned n) const;
        virtual Atom& getInBond(unsigned n);
        virtual const Atom& getOutBond(unsigned n) const;
        virtual Atom& getOutBond(unsigned n);

        void setSuperior(Group* gr);
        void setModified();
        void setUnModified(); //used in Qmean: this flag cause problem there.

        // OPERATORS:
        Atom& operator=(const Atom& orig);

    protected:

    private:

        // HELPERS: 
        bool isNotFirstAtomInStructure();
        void propagateRotation(); // pass on rotation matrix to following atoms


        // ATTRIBUTES:
        Group* superior; // structure to which atom belongs,
        AtomCode type; // Atom type
        vgVector3<double> coords; // xyz-Coords

        double Bfac; // B-factor
        double Occupancy; // Occupancy

        vgVector3<double> trans; // relative translation
        vgMatrix3<double> rot; // relative rotation
        bool modified; // --""--  modified?  
    };



    // ---------------------------------------------------------------------------
    //                                    Atom
    // -----------------x-------------------x-------------------x-----------------

    // PREDICATES:

    inline AtomCode
    Atom::getCode() const {
        return type;
    }

    inline unsigned long
    Atom::getNumber() const {
        return id;
    }

    inline vgVector3<double>
    Atom::getCoords() {
        if (!inSync())
            sync();
        return coords;
    }

    inline Group&
    Atom::getSuperior() {
        PRECOND(superior != NULL, exception);
        return *superior;
    }

    inline const Group&
    Atom::getSuperior() const {
        PRECOND(superior != NULL, exception);
        return *superior;
    }

    inline bool
    Atom::hasSuperior() {
        return (superior != NULL);
    }

    inline vgVector3<double>
    Atom::getTrans() const {
        return trans;
    }

    inline vgMatrix3<double>
    Atom::getRot() const {
        return rot;
    }

    inline bool
    Atom::inSync() {
        return (!modified);
    }

    // MODIFIERS:

    inline void
    Atom::bindStructure(Atom& before, bool connect) {
        if (connect)
            this->setTrans(this->getCoords() - before.getCoords());
        this->bindIn(before);
    }

    inline void
    Atom::bindIn(SimpleBond& c) {
        SimpleBond::bindIn(c);
        setModified();
    }

    inline void
    Atom::bindOut(SimpleBond& c) {
        SimpleBond::bindOut(c);
        setModified();
    }

    inline void
    Atom::unbindIn(SimpleBond& c) {
        SimpleBond::unbindIn(c);
        setModified();
    }

    inline void
    Atom::unbindOut(SimpleBond& c) {
        SimpleBond::unbindOut(c);
        setModified();
    }

    inline void
    Atom::setType(string _name) {
        id.setName(_name);
        type = AtomTranslator(_name);
    }

    inline void
    Atom::setCode(AtomCode ac) {
        type = ac;
        id.setName(AtomTranslator(ac));
    }

    inline void
    Atom::setNumber(unsigned long _number) {
        id.setNumber(_number);
    }

    inline void
    Atom::setTrans(vgVector3<double> t) {
        trans = t;
        setModified();
    }

    inline void
    Atom::addTrans(vgVector3<double> t) {
        trans += t;
        setModified();
    }

    inline void
    Atom::setRot(vgMatrix3<double> r) {
        rot = r;
        setModified();
    }

    inline void
    Atom::addRot(vgMatrix3<double> r) {
        rot = r * rot;
        setModified();
    }

    inline const Atom&
    Atom::getInBond(unsigned n) const {
        return dynamic_cast<const Atom&> (SimpleBond::getInBond(n));
    }

    inline Atom&
    Atom::getInBond(unsigned n) {
        return dynamic_cast<Atom&> (SimpleBond::getInBond(n));
    }

    inline const Atom&
    Atom::getOutBond(unsigned n) const {
        return dynamic_cast<const Atom&> (SimpleBond::getOutBond(n));
    }

    inline Atom&
    Atom::getOutBond(unsigned n) {
        return dynamic_cast<Atom&> (SimpleBond::getOutBond(n));
    }

    inline void
    Atom::setSuperior(Group* gr) {
        this->superior = gr;
    }

    // OPERATORS:


}} //namespace
#endif //_ATOM_H_
