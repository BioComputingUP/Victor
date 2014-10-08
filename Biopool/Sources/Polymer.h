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


#ifndef _POLYMER_H_
#define _POLYMER_H_

// Includes:
#include <Component.h>
#include <Debug.h>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {

    /**@brief Implements methods to manage a polymer created by components.
     * 
     *@Description Includes methods that allow to insert and remove components.
     */
    class Polymer : public Component {
    public:

        // CONSTRUCTORS/DESTRUCTOR:
        Polymer(unsigned int mI = 2, unsigned int mO = 2);
        Polymer(const Polymer& orig);
        virtual ~Polymer();

        // PREDICATES:
        //char getChainID() { return chainID; }

        virtual string getClassName() const {
            return "Polymer";
        }

        virtual vgVector3<double> getTrans() const {
            ERROR("Polymer::getTrans is not a viable method.", exception);
        }

        virtual vgMatrix3<double> getRot() const {
            ERROR("Polymer::getRot is not a viable method.", exception);
        }

        virtual void save(Saver& s) {
            ERROR("Polymer::save is not a viable method.", exception);
        }

        // MODIFIERS:
        //void setChainID(char _chainID){ chainID = _chainID; }
        virtual void insertComponent(Component* c);
        virtual void removeComponent(Component* c);
        //---------------------------------
        virtual void removeComponentFromIndex(unsigned int i);
        //---------------------------------
        virtual void deleteComponent(Component* c);

        void copy(const Polymer& orig);
        virtual Component* clone();

        virtual void load(Loader& l) {
            ERROR("Polymer::load is not a viable method.", exception);
        }

        virtual void setTrans(vgVector3<double> t) {
            ERROR("Polymer::setTrans is not a viable method.", exception);
        }

        virtual void addTrans(vgVector3<double> t) {
            ERROR("Polymer::addTrans is not a viable method.", exception);
        }

        virtual void setRot(vgMatrix3<double> r) {
            ERROR("Polymer::setRot is not a viable method.", exception);
        }

        virtual void addRot(vgMatrix3<double> r) {
            ERROR("Polymer::addRot is not a viable method.", exception);
        }

        virtual void sync() // synchronize coords with structure
        {
            ERROR("Polymer::sync is not a viable method.", exception);
        }

        virtual void acceptCalculator(EnergyVisitor* v) {
            ERROR("Polymer::acceptCalculator is not a viable method.", exception);
        }

        virtual void acceptOptimizer(OptimizationVisitor* v) {
            ERROR("Polymer::acceptOptimizer is not a viable method.", exception);
        }

        // OPERATORS:
        virtual Component& operator[](unsigned int n);
        virtual const Component& operator[](unsigned int n) const;

    protected:

        // HELPERS:

        virtual void resetBoundaries() {
            ERROR("Polymer::resetBoundaries is not a viable method.", exception);
        }
        // ATTRIBUTES


    private:

    };


    // ---------------------------------------------------------------------------
    //                                    Polymer
    // -----------------x-------------------x-------------------x-----------------

    // PREDICATES:

    // MODIFIERS:

    // OPERATORS:
    inline Component&
    Polymer::operator[](unsigned int n) {
        PRECOND(n < size(), exception);
        return *components[n];
    }

    inline const Component&
    Polymer::operator[](unsigned int n) const {
        PRECOND(n < size(), exception);
        return *components[n];
    }

} // namespace
#endif //_POLYMER_H_


