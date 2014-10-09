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

#ifndef _GROUP_H_
#define _GROUP_H_

// Includes:
#include <vector>
#include <Atom.h>
#include <Monomer.h>
#include <matrix3.h>
#include <Visitor.h>

namespace Victor { namespace Biopool { 
    
    /**@brief This class implements a simple chemical group
     * 
     * */
    
    // Global constants, typedefs, etc. (to avoid):

    class Group : public Monomer {
    public:

        // CONSTRUCTORS/DESTRUCTOR:
        Group(unsigned int mI = 1, unsigned int mO = 4);
        Group(const Group& orig);
        virtual ~Group();

        virtual string getClassName() const {
            return "Group";
        }

        // PREDICATES:
        virtual unsigned int getCode() const;
        virtual unsigned int size() const;

        Atom& getAtom(unsigned int n);
        const Atom& getAtom(unsigned int n) const;

        virtual void save(Saver& s); // data saver
        vector<Atom>& giveAtoms();

        vgVector3<double> getTrans() const;
        vgMatrix3<double> getRot() const;

        virtual bool isMember(const AtomCode& ac) const;

        // MODIFIERS:
        void copy(const Group& orig);
        void setAtom(unsigned int n, Atom& at);
        void addAtom(Atom& a);
        void removeAtom(Atom& a);

        virtual void load(Loader& l); // data loader

        void setTrans(vgVector3<double> t);
        void addTrans(vgVector3<double> t);
        void setRot(vgMatrix3<double> r);
        void addRot(vgMatrix3<double> r);
        virtual void sync(); // synchronize coords with structure
        virtual void setModified();

        virtual const Group& getInBond(unsigned int n) const;
        virtual Group& getInBond(unsigned int n);
        virtual const Group& getOutBond(unsigned int n) const;
        virtual Group& getOutBond(unsigned int n);

        virtual void acceptCalculator(EnergyVisitor* v);
        virtual void acceptOptimizer(OptimizationVisitor* v);

        // OPERATORS:

        Group& operator=(const Group& orig);
        virtual Atom& operator[](unsigned int n);
        virtual const Atom& operator[](unsigned int n) const;
        virtual Atom& operator[](const AtomCode& ac);
        virtual const Atom& operator[](const AtomCode& ac) const;

    protected:

        // HELPERS:
        void resetBoundaries();
        Atom* pGetAtom(const AtomCode& ac) const;

    private:

        // HELPERS:
        unsigned int pGetIndex(const unsigned long cmp) const;

        // ATTRIBUTES:
        vector<Atom> atoms;
        vgVector3<double> trans; // group translation
        vgMatrix3<double> rot; // group rotation

    };

    // ---------------------------------------------------------------------------
    //                                    Group
    // -----------------x-------------------x-------------------x-----------------

    // PREDICATES:

    inline unsigned int
    Group::getCode() const {
        ERROR("Group::getCode() undefined.", exception);
        return 0;
    }

    inline unsigned int
    Group::size() const {
        return atoms.size();
    }

    inline void
    Group::save(Saver& s) {
        s.saveGroup(*this);
    }

    inline vector<Atom>&
    Group::giveAtoms() {
        return atoms;
    }

    inline Atom&
    Group::getAtom(unsigned int n) {
        PRECOND(n < atoms.size(), exception);
        return atoms[n];
    }

    inline const Atom&
    Group::getAtom(unsigned int n) const {
        PRECOND(n < atoms.size(), exception);
        return atoms[n];
    }

    inline vgVector3<double>
    Group::getTrans() const {
        return trans;
    }

    inline vgMatrix3<double>
    Group::getRot() const {
        return rot;
    }

    inline bool
    Group::isMember(const AtomCode& ac) const {
        return (pGetAtom(ac) != NULL);
    }

    // MODIFIERS:

    inline void
    Group::setAtom(unsigned int n, Atom& at) {
        PRECOND(n < atoms.size(), exception);
        atoms[n] = at;
        atoms[n].setSuperior(this);
        Group::setModified();
    }

    inline void
    Group::load(Loader& l) {
        l.loadGroup(*this);
        resetBoundaries();
    }

    inline void
    Group::setTrans(vgVector3<double> t) {
        trans = t;
        Group::setModified();
    }

    inline void
    Group::addTrans(vgVector3<double> t) {
        trans += t;
        Group::setModified();
    }

    inline void
    Group::setRot(vgMatrix3<double> r) {
#ifndef NDEBUG
#endif
        rot = r;
        Group::setModified();
    }

    inline void
    Group::addRot(vgMatrix3<double> r) {
#ifndef NDEBUG
#endif
        rot = r * rot;
        Group::setModified();
    }

    inline const Group&
    Group::getInBond(unsigned int n) const {
        return dynamic_cast<const Group&> (Bond::getInBond(n));
    }

    inline Group&
    Group::getInBond(unsigned int n) {
        return dynamic_cast<Group&> (Bond::getInBond(n));
    }

    inline const Group&
    Group::getOutBond(unsigned int n) const {
        return dynamic_cast<const Group&> (Bond::getOutBond(n));
    }

    inline Group&
    Group::getOutBond(unsigned int n) {
        return dynamic_cast<Group&> (Bond::getOutBond(n));
    }

    inline void
    Group::acceptCalculator(EnergyVisitor* v) {
        v->PrepareGroup(*this);
    }

    inline void
    Group::acceptOptimizer(OptimizationVisitor* v) {
        v->PrepareGroup(*this);
    }


    // OPERATORS:

    inline Atom&
    Group::operator[](unsigned int n) {
        PRECOND(n < atoms.size(), exception);
        return atoms[n];
    }

    inline const Atom&
    Group::operator[](unsigned int n) const {
        PRECOND(n < atoms.size(), exception);
        return atoms[n];
    }

    inline Atom&
    Group::operator[](const AtomCode& ac) {
        Atom* a = pGetAtom(ac);
        INVARIANT(a != NULL, exception);
        if (a == NULL)
            ERROR("Inexistent atom requested.", exception);
        return *a;
    }

    inline const Atom&
    Group::operator[](const AtomCode& ac) const {
        Atom* a = pGetAtom(ac);
        INVARIANT(a != NULL, exception);
        if (a == NULL)
            ERROR("Inexistent atom requested.", exception);
        return *a;
    }

    // HELPERS:

}} //namespace
#endif //_GROUP_H_
