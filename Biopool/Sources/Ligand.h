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


#ifndef _Ligand_H_
#define _Ligand_H_


// Includes:
#include <Bond.h>
#include <Polymer.h>
#include <Group.h>
#include <Visitor.h>
#include <AminoAcidCode.h>
#include <NucleotideCode.h>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {

    /**@brief Implements methods to verify the ligand properties
     * 
     * */
 
    class Ligand : public Group {
    public:

        // CONSTRUCTORS/DESTRUCTOR:
        Ligand();
        Ligand(const Ligand& orig);
        virtual ~Ligand();

        virtual string getClassName() const {
            return "Ligand";
        }

        // PREDICATES:
        bool isMetalCompound(); 
        bool isCommonMetal(); 

        bool isSimpleMetalIon(); //so I introduced thoose to distinguish between 
        bool isWater(); //water, metal ions and other cofactors
        bool isCofactor();
        virtual void save(Saver& s); // data saver


        // MODIFIERS:
        void copy(const Ligand& orig);

        virtual void load(Loader& l); // data loader

        // OPERATORS:
        Ligand& operator=(const Ligand& orig);
        virtual Atom& operator[](unsigned int n);
        virtual const Atom& operator[](unsigned int n) const;


    protected:
        // HELPERS:

        // ATTRIBUTES:

    private:
    };

    // ---------------------------------------------------------------------------
    //                                    Ligand
    // -----------------x-------------------x-------------------x-----------------

    // PREDICATES:

    inline void
    Ligand::save(Saver& s) {
        s.saveLigand(*this);
    }

    // MODIFIERS:

    inline void
    Ligand::load(Loader& l) {
        l.loadLigand(*this);
        resetBoundaries();
    }


    // OPERATORS:

    inline Atom&
    Ligand::operator[](unsigned int n) {
        return Group::operator[](n);
    }

    inline const Atom&
    Ligand::operator[](unsigned int n) const {
        return Group::operator[](n);
    }


} // namespace
#endif //_Ligand_H_

