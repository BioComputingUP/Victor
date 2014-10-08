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
#ifndef _SPACER_H_
#define _SPACER_H_


// Includes:
#include <Bond.h>
#include <Polymer.h>
#include <AminoAcid.h>
#include <Visitor.h>
#include <AminoAcidCode.h>
#include <map>
#include <set>


// Global constants, typedefs, etc. (to avoid):

namespace Biopool {

    /**@brief Implements a "Spacer" for a protein chain. Includes methods to obtain values from the atoms and its pdb information.
     * 
     *@Description ***Attention*** The current implementation allows for "1 to 1" spacers,
     *    ie. spacers composed of a single aminoacid chain. Class that develop the virtual methods from the parent class
     *   Polymer.
     * */
   

    class Spacer : public Polymer {
    public:
        //AK11
        friend class AtomEnergyCalculator;
        friend class SideChainPlacement;

        // CONSTRUCTORS/DESTRUCTOR:
        Spacer();
        Spacer(const Spacer& orig);
        virtual ~Spacer();

        // PREDICATES:

        int getStartOffset() {
            return startOffset;
        }

        int getAtomStartOffset() {
            return startAtomOffset;
        }
        int getIndexFromPdbNumber(int index); // in: pdb_num, out: array_num
        int getPdbNumberFromIndex(int index); // in: array_num, out: pdb_num

        int maxPdbNumber() {
            return startOffset + gaps.size() + sizeAmino();
        }
        bool isGap(int index);
        void printGaps();

        unsigned int sizeGaps() {
            return gaps.size();
        }

        virtual string getClassName() const {
            return "Spacer";
        }
        void save(Saver& s); // data saver
        AminoAcid& getAmino(unsigned int n);
        const AminoAcid& getAmino(unsigned int n) const;
        const unsigned int sizeAmino() const;

        Spacer& getSpacer(unsigned int n);
        const Spacer& getSpacer(unsigned int n) const;
        const unsigned int sizeSpacer() const;

        vector<pair<unsigned int, unsigned int> > getHelixData();
        vector<pair<unsigned int, unsigned int> > getStrandData();

        // MODIFIERS:

       
        
        void insertAminoAfter(string code, unsigned int n = 9999, double phi = -62,
                double psi = -41, double omega = 180,
                double NToCaLen = 1.458, double CaToCLen = 1.52, double CToOLen = 1.231,
                double atCaToCAng = 111.6, double CaToOAng = 120.80,
                double CtoNLen = 1.33, double atCToNAng = 116.4,
                double OToNAng = 123.2, double atNToCaAng = 121.9);
        
       
        void insertAminoAfterWithGaps(string code, unsigned int n = 9999, double phi = -62,
                double psi = -41, double omega = 180,
                int beginHole = 0, int endHole = 0,
                vgVector3 <double> ca1 = vgVector3<double>(), vgVector3 <double> ca2 = vgVector3<double>(),
                string target = "", Spacer* refSpacer = 0, Spacer* pOriginalSpacer = 0,
                double NToCaLen = 1.458, double CaToCLen = 1.52, double CToOLen = 1.231,
                double atCaToCAng = 111.6, double CaToOAng = 120.80,
                double CtoNLen = 1.33, double atCToNAng = 116.4,
                double OToNAng = 123.2, double atNToCaAng = 121.9);
        
       
        
        void insertAminoBefore(string code, unsigned int p = 0, double phi = -62,
                double psi = -41, double omega = 180,
                double NToCaLen = 1.458, double CaToCLen = 1.52, double CToOLen = 1.231,
                double atCaToCAng = 111.6, double CaToOAng = 120.80,
                double CtoNLen = 1.33, double atCToNAng = 116.4,
                double OToNAng = 123.2, double atNToCaAng = 121.9);
        void setStartOffset(int _offset);

        /**
         * @Description Set the atom offset
         * @param _offset
         */
        void setAtomStartOffset(int _offset) {
            startAtomOffset = _offset;
        }

        void addGap(int index);
        void removeGap(int index);
        void removeAllGaps();

        void insertSubSpacerAfter(Spacer* sp, unsigned int pos);

        void insertFirstSpacer(Spacer* sp);


        void insertComponent(Component* c);
        void removeComponent(Component* c);
        //-------------
        void removeComponentFromIndex(unsigned int i);
        //--------------
        void deleteComponent(Component* c);

        void mergeSpacer(Spacer* s);
        Spacer* splitSpacer(unsigned int index, unsigned int count);

        void copy(const Spacer& orig);

        void load(Loader& l); // data loader

        void setTrans(vgVector3<double> t);
        void addTrans(vgVector3<double> t);
        void setRot(vgMatrix3<double> r);
        void addRot(vgMatrix3<double> r);
        void sync();

        void setStateFromSecondary(string sec);
        void setStateFromTorsionAngles();
        void setDSSP(bool verbose);

        vector< set<char> > getDSSP() {
            return ss;
        };


        const Spacer& getInBond(unsigned int n) const;
        Spacer& getInBond(unsigned int n);
        const Spacer& getOutBond(unsigned int n) const;
        Spacer& getOutBond(unsigned int n);

        void makeFlat();
        void groupLikeStateCode();

        virtual const Atom& getOpenInBondRef(unsigned int n = 0) const;
        virtual Atom& getOpenInBondRef(unsigned int n = 0);
        virtual const Atom& getOpenOutBondRef(unsigned int n = 0) const;
        virtual Atom& getOpenOutBondRef(unsigned int n = 0);

        virtual Component* clone();

        virtual void acceptCalculator(EnergyVisitor* v);
        virtual void acceptOptimizer(OptimizationVisitor* v);

        // OPERATORS:
        Spacer& operator=(const Spacer& orig);

        // TESTERS AND DEVELOPERS
        void printSubSpacerList();
        void printComponents();
        vector< pair<int, int> > getHoles();

        static double BOND_ANGLE_AT_CPRIME_TO_N;
        static double BOND_LENGTH_CPRIME_TO_N;
        static double BOND_ANGLE_CA_TO_O;
        static double BOND_LENGTH_C_TO_O;
        static double BOND_LENGTH_CALPHA_TO_CPRIME;
        static double BOND_LENGTH_N_TO_CALPHA;
        static double BOND_ANGLE_AT_CALPHA_TO_CPRIME;
        static double BOND_ANGLE_O_TO_N;
        static double BOND_ANGLE_AT_N_TO_CALPHA;

        // 
        static bool NMRGetMinimumCADistanceVector(const string &strFileName, vector<double> *pNMR);

    protected:

        // HELPERS:
        void resetBoundaries();
        void modifySubSpacerList(Spacer*, int);
        void updateSubSpacerList();
        void getBackboneHbonds(); // backbone H bonds (for SS)

        // ATTRIBUTES

        int startOffset; //number (as reported in the pdb file) of first amminoacid loaded 
        int startAtomOffset; //number of first atom of the first amminoacid
        vector<int> gaps;

        pair<unsigned int, unsigned int> getSubSpacerListEntry(unsigned int);
        pair<unsigned int, unsigned int> getSubSpacerListEntry(unsigned int) const;

        bool **backboneHbonds;
        vector<set<char > > ss;

    private:
        vector<pair<unsigned int, unsigned int> > subSpacerList;
    };

    // ---------------------------------------------------------------------------
    //                                    Spacer
    // -----------------x-------------------x-------------------x-----------------

    // PREDICATES:

    /**
     *@Description Saves the information from the saver into the Spacer
     *@param s, reference to the saver
     */
    inline void Spacer::save(Saver& s) {
        s.saveSpacer(*this);
    }

    /**
     *@Description Allows to know if there is a spacer in the protein or not. 
     *@return quantity of spacers, possible values are 0 or 1
     */
    inline const unsigned int Spacer::sizeSpacer() const {
        return subSpacerList.size();
    }

    // MODIFIERS:

    /**
     *@Description Removes the corresponding component
     *@param c, pointer to the component to remove, usually the object that is calling the method.
     */
    inline void Spacer::removeComponent(Component* c) {
        Polymer::removeComponent(c);
        setModified();
    }

    /**
     *@Description Removes the corresponding component in the corresponding index
     *@param i, index to the component to remove, usually the index of the object that is calling the method.
     */
    inline void Spacer::removeComponentFromIndex(unsigned int i) {
        Polymer::removeComponentFromIndex(i);
        setModified();
    }

             /**
     *@Description Similar to removeComponent but also free the component's memory space
     *@param c, pointer to the Component to delete
     */inline void Spacer::deleteComponent(Component* c) {
        Polymer::deleteComponent(c);
        setModified();
    }

    /**
     *@Description Loads the Spacer from a Loader
     *@param l, reference of the Loader
     */
    inline void Spacer::load(Loader& l) {
        l.loadSpacer(*this);
        resetBoundaries();
    }

    /**
     *@Description Sets the translation vector
     *@param t, the translation vector to set
     */
    inline void Spacer::setTrans(vgVector3<double> t) {
        if (sizeAmino())
            getAmino(0).setTrans(t);
    }

    /**
     *@Description Adds a the translation vector
     *@param t, the translation vector to add
     */
    inline void Spacer::addTrans(vgVector3<double> t) {
        if (sizeAmino())
            getAmino(0).addTrans(t);
    }

    /**
     *@Description Sets the rotation matrix
     *@param r, the rotation matrix to set
     */
    inline void Spacer::setRot(vgMatrix3<double> r) {
        if (sizeAmino())
            getAmino(0).setRot(r);
    }

    /**
     *@Description Adds a rotation matrix
     *@param r, the rotation matrix to add
     */
    inline void Spacer::addRot(vgMatrix3<double> r) {
        if (sizeAmino())
            getAmino(0).addRot(r);
    }

    /**
     *@Description Returns the reference to the spacer that contains the inBond for the n Atom
     *@param n, index of the atom
     *@return reference to the spacer that contains the inBond information
     */
    inline const Spacer& Spacer::getInBond(unsigned int n) const {
        return dynamic_cast<const Spacer&> (Bond::getInBond(n));
    }

    /**
     *@Description Returns the reference to the spacer that contains the inBond for the n Atom
     *@param n, index of the atom
     *@return reference to the spacer that contains the inBond information
     */
    inline Spacer& Spacer::getInBond(unsigned int n) {
        return dynamic_cast<Spacer&> (Bond::getInBond(n));
    }

    /**
     *@Description Returns the reference to the spacer that contains the outBond for the n Atom
     *@param n, index of the atom
     *@return reference to the spacer that contains the outBond information
     */
    inline const Spacer& Spacer::getOutBond(unsigned int n) const {
        return dynamic_cast<const Spacer&> (Bond::getOutBond(n));
    }

    /**
     *@Description Returns the reference to the spacer that contains the outBond for the n Atom
     *@param n, index of the atom
     *@return reference to the spacer that contains the outBond information 
     */
    inline Spacer& Spacer::getOutBond(unsigned int n) {
        return dynamic_cast<Spacer&> (Bond::getOutBond(n));
    }

    /**
     *@Description Sets the Spacer as part of the Energy Visitor object
     *@param v, pointer to the Energy visitor object
     */
    inline void Spacer::acceptCalculator(EnergyVisitor* v) {
        v->PrepareSpacer(*this);
    }

    /**
     *@Description Sets the Spacer as part of the Optimization Visitor object
     *@param v, pointer to the Optimization Visitor
     */
    inline void Spacer::acceptOptimizer(OptimizationVisitor* v) {
        v->PrepareSpacer(*this);
    }





    // OPERATORS:

} // namespace
#endif //_SPACER_H_
