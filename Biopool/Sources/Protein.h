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


#ifndef PROTEIN_H
#define	PROTEIN_H

// Includes:
#include <Component.h>
#include <Polymer.h>
#include <Spacer.h>
#include <LigandSet.h>

namespace Victor { namespace Biopool { 
    
    /**
    * @brief Is a container of Polymers (a component of components).
    *        Each chain in a PDB (a Protein object)
    *                   corresponds to a Polymers. Each Polymer corresponds to a Spacer and a LigandSet (if present).
    *  
    */
    
    class Protein : public Polymer {
    public:

        // CONSTRUCTORS/DESTRUCTOR:
        Protein();
        Protein(const Protein& orig);
        virtual ~Protein();

        // PREDICATES:

        virtual string getClassName() const {
            return "Protein";
        }
        Polymer& getPolymer(unsigned int n);
        Spacer* getSpacer(char c); //return a pointer to the Spacer with the correct chainID
        Spacer* getSpacer(unsigned int n);

        LigandSet* getLigandSet(char c); //return a pointer to the LigandSet with the correct chainID(NULL if not exist)
        LigandSet* getLigandSet(unsigned int n);

        const unsigned int sizeProtein() const;
        unsigned int getChainNum(char c);
        char getChainLetter(unsigned int i);
        vector <char> getAllChains();

        void save(Saver& s); // data saver                  

        // MODIFIERS:
        void addChain(char c);
        void insertComponent(Component* c);
        void setChainSelection();
        void removeComponent(Component* c);
        void deleteComponent(Component* c);

        void copy(const Protein& orig);
        void load(Loader& l); // data loader                 

        virtual Protein* clone();

        // OPERATORS:
        Protein& operator=(const Protein& orig);
        Polymer& operator[](unsigned int n);
        //const Polymer& operator[](unsigned int n) const;     

    protected:
        // HELPERS:
        void printComponents();
        // ATTRIBUTES 
    private:
        vector<char> chains;
        
    };
    // ---------------------------------------------------------------------------
    //                                    Protein
    // -----------------x-------------------x-------------------x-----------------

    //PREDICATES

    inline const unsigned int Protein::sizeProtein() const {
        return components.size();
    }

    inline vector<char>Protein::getAllChains() {
        return chains;
    }

    inline void Protein::save(Saver& s) {
        s.saveProtein(*this);
    }

    // MODIFIERS

    inline void Protein::addChain(char c) {
        chains.push_back(c);
    }

    inline void Protein::load(Loader& l) {
        l.loadProtein(*this);
    }

    // OPERATORS:
    inline Polymer&
    Protein::operator[](unsigned int n) {
        return getPolymer(n);
    }



    // HELPERS:

}}//namespace
#endif	/* PROTEIN2_H */

