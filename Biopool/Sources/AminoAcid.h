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


#ifndef _AMINOACID_H_
#define _AMINOACID_H_



// Includes:
#include <Group.h>
#include <SideChain.h>
#include <AminoAcidCode.h>
#include <IntCoordConverter.h>
#include <Visitor.h>


// Global constants, typedefs, etc. (to avoid):

namespace Biopool {

    /**@brief It mplements a simple amino acid.
     * 
     *@Description Includes methods that allow to get and set angles, connections, bonds, side chain info, etc.
     * NB: Angles are in degrees.
     * */
    
    
    
    class AminoAcid : public Group {
    public:

        // CONSTRUCTORS/DESTRUCTOR:
        AminoAcid();
        AminoAcid(const AminoAcid& orig);
        virtual ~AminoAcid();

        // PREDICATES:

        virtual string getClassName() const {
            return "AminoAcid";
        }
        virtual unsigned int getCode() const;
        char getType1L();
        unsigned int size() const;
        unsigned int sizeBackbone() const;

        double getPhi(bool override = false);
        double getPsi(bool override = false);
        double getOmega(bool override = false);
        double getChi(unsigned int n);
        vector<double> getChi();
        unsigned int getMaxChi();

        StateCode getState(); // returns AA state (H, S, T, C)

        SideChain& getSideChain();
        const SideChain& getSideChain() const;

        bool isMember(const AtomCode& ac) const;

        void save(Saver& s); // data saver

        // MODIFIERS:
        virtual void connectIn(AminoAcid* c, unsigned int offset = 1);
        virtual void connectOut(AminoAcid* c, unsigned int offset = 1);
        virtual AminoAcid* unconnectIn();
        virtual AminoAcid* unconnectOut();

        void copy(const AminoAcid& orig);
        void setType1L(char _name);
        void setType(string _name);

        void setPhi(double a);
        void setPsi(double a);
        void setOmega(double a);
        void setChi(unsigned int n, double a);
        void setChi(vector<double> cv);

        void setState(StateCode sc); // sets AA state (H, S, T, C)

        void setDefault();
        void setBonds(double NToCaLen, double CaToCLen, double CToOLen, double atCaToCAng, double CaToOAng);

        void constructSideChain(SideChain& sc, vector<double> chi);
        void constructSideChain(SideChain& sc, double chi1, double chi2,
                double chi3, double chi4, double chi5);
        void setSideChain(SideChain& sc);
        void setStateFromTorsionAngles(); // sets states from its phi/psi angles
        void adjustLeadingN(); // adjusts the translation of the N atom
        void addTerminalOXT(); // adds an OXT atom to the C terminus
        void addMissingO(); // adds an O atom where missing
        void removeHAtomsfromLeadingNH3(); // removes H atoms in excess from NH3+
        void removeSideChain();

        void patchBetaPosition(unsigned int n = 0);
        // adds a CB atom to the sidechain, if necessary
        void patchAminoAcidCode(); // determines the 3-letter code from the sidechain
        bool setBondsFromPdbCode(bool connect, AminoAcid* prev = NULL,
                bool permissive = false); // returns false if failed to connect
        virtual void sync(); // synchronize coords with structure
        virtual void setModified();

        void load(Loader& l); // data loader

        const AminoAcid& getInBond(unsigned int n) const;
        AminoAcid& getInBond(unsigned int n);
        const AminoAcid& getOutBond(unsigned int n) const;
        AminoAcid& getOutBond(unsigned int n);

        virtual const Atom& getOpenInBondRef(unsigned int n = 0) const;
        virtual Atom& getOpenInBondRef(unsigned int n = 0);
        virtual const Atom& getOpenOutBondRef(unsigned int n = 0) const;
        virtual Atom& getOpenOutBondRef(unsigned int n = 0);

        virtual Component* clone();

        virtual void acceptCalculator(EnergyVisitor* v);
        virtual void acceptOptimizer(OptimizationVisitor* v);

        // OPERATORS:
        bool operator==(const AminoAcid& other) const;
        bool operator!=(const AminoAcid& other) const;
        AminoAcid& operator=(const AminoAcid& orig);
        Atom& operator[](unsigned int n);
        const Atom& operator[](unsigned int n) const;
        Atom& operator[](const AtomCode& ac);
        const Atom& operator[](const AtomCode& ac) const;

    protected:

        // HELPERS:
        virtual void resetBoundaries();

    private:

        // ATTRIBUTES:
        static double BOND_ANGLE_N_TO_CB;
        static double BOND_ANGLE_CB_TO_C;
        static double BOND_LENGTH_CA_TO_CB;

        AminoAcidCode type; // AminoAcid type
        StateCode state; //  -- '' -- state (H, S, T, C)
        double phi, psi, omega; // torsion angles
        SideChain sideChain;
        // vector atoms (inherited from Group) contains the backbone
        IntCoordConverter icc; // should be static, but compiler won't accept it?
    };

    // ---------------------------------------------------------------------------
    //                                    AminoAcid
    // -----------------x-------------------x-------------------x-----------------

    // PREDICATES:

    inline unsigned int
    AminoAcid::getCode() const {
        return aminoAcidThreeLetterTranslator(id);
    }

    inline char
    AminoAcid::getType1L() {
        return threeLetter2OneLetter(id);
    }

    inline unsigned int
    AminoAcid::size() const {
        return (Group::size() + sideChain.size());
    }

    inline unsigned int
    AminoAcid::sizeBackbone() const {
        return Group::size();
    }

    inline double
    AminoAcid::getPhi(bool override) {
        if ((phi > 990) || (override))
            if (sizeInBonds())
                phi = RAD2DEG * icc.getTorsionAngle(getInBond(0)[C], (*this)[N],
                    (*this)[CA], (*this)[C]);
        return phi;
    }

    inline double
    AminoAcid::getPsi(bool override) {
        if ((psi > 990) || (override))
            if (sizeOutBonds())
                psi = RAD2DEG * icc.getTorsionAngle((*this)[N], (*this)[CA],
                    (*this)[C], getOutBond(0)[N]);
        return psi;
    }

    inline double
    AminoAcid::getOmega(bool override) {
        if ((omega > 990) || (override))
            if (sizeOutBonds())
                omega = RAD2DEG * icc.getTorsionAngle((*this)[CA], (*this)[C],
                    getOutBond(0)[N], getOutBond(0)[CA]);
        return omega;
    }

    inline double
    AminoAcid::getChi(unsigned int n) {
        return sideChain.getChi(n);
    }

    inline vector<double>
    AminoAcid::getChi() {
        return sideChain.getChi();
    }

    inline unsigned int
    AminoAcid::getMaxChi() {
        return sideChain.getMaxChi();
    }

    inline StateCode
    AminoAcid::getState() {
        return state;
    }

    inline const SideChain&
    AminoAcid::getSideChain() const {
        return sideChain;
    }

    inline SideChain&
    AminoAcid::getSideChain() {
        return sideChain;
    }

    inline bool
    AminoAcid::isMember(const AtomCode& ac) const {
        return (pGetAtom(ac) != NULL);
    }

    inline void
    AminoAcid::save(Saver& s) {
        s.saveAminoAcid(*this);
    }


    // MODIFIERS:

    inline void
    AminoAcid::setType1L(char _name) {
        type = aminoAcidOneLetterTranslator(_name);
        id.setName(oneLetter2ThreeLetter(_name));
    }

    inline void
    AminoAcid::setType(string _name) {
        type = aminoAcidThreeLetterTranslator(_name);
        id.setName(_name);
    }

    inline void
    AminoAcid::setChi(unsigned int n, double a) {
        sideChain.setChi(n, a);
    }

    inline void
    AminoAcid::setChi(vector<double> cv) {
        sideChain.setChi(cv);
    }

    inline void
    AminoAcid::setState(StateCode sc) {
        state = sc;
    }

    inline void
    AminoAcid::removeSideChain() {
        SideChain sc;
        sideChain = sc;
        resetBoundaries();
    }

    inline void
    AminoAcid::load(Loader& l) {
        l.loadAminoAcid(*this);
        resetBoundaries();
    }

    inline const AminoAcid&
    AminoAcid::getInBond(unsigned int n) const {
        return dynamic_cast<const AminoAcid&> (Bond::getInBond(n));
    }

    inline AminoAcid&
    AminoAcid::getInBond(unsigned int n) {
        return dynamic_cast<AminoAcid&> (Bond::getInBond(n));
    }

    inline const AminoAcid&
    AminoAcid::getOutBond(unsigned int n) const {
        return dynamic_cast<const AminoAcid&> (Bond::getOutBond(n));
    }

    inline AminoAcid&
    AminoAcid::getOutBond(unsigned int n) {
        return dynamic_cast<AminoAcid&> (Bond::getOutBond(n));
    }

    inline void
    AminoAcid::patchAminoAcidCode() { // determines the 3-letter code from the sidechain
        sideChain.patchAminoAcidCode();
        setType(sideChain.getType());
    }

    inline void
    AminoAcid::setModified() {
        Group::setModified();
        sideChain.setModified();
    }

    inline void
    AminoAcid::acceptCalculator(EnergyVisitor* v) {
        v->PrepareAminoAcid(*this);
    }

    inline void
    AminoAcid::acceptOptimizer(OptimizationVisitor* v) {
        v->PrepareAminoAcid(*this);
    }


    // OPERATORS:

    inline bool
    AminoAcid::operator==(const AminoAcid& other) const {
        return (dynamic_cast<const Identity*> (this)) ==
                dynamic_cast<const Identity*> (&other);
    }

    inline bool
    AminoAcid::operator!=(const AminoAcid& other) const {
        return (dynamic_cast<const Identity*> (this)) !=
                dynamic_cast<const Identity*> (&other);
    }


    inline Atom&
    AminoAcid::operator[](unsigned int n) {
        PRECOND(n < this->size(), exception);
        return ((n < Group::size()) ? Group::getAtom(n) : sideChain[n - Group::size()]);
    }

    inline const Atom&
    AminoAcid::operator[](unsigned int n) const {
        PRECOND(n < this->size(), exception);
        return ((n < Group::size()) ? Group::getAtom(n) : sideChain[n - Group::size()]);
    }

    inline Atom&
    AminoAcid::operator[](const AtomCode& ac) {
        Atom* a = pGetAtom(ac);
        INVARIANT(a != NULL, exception);
        return *a;
    }

    inline const Atom&
    AminoAcid::operator[](const AtomCode& ac) const {
        Atom* a = pGetAtom(ac);
        INVARIANT(a != NULL, exception);
        return *a;
    }

    // HELPERS:

} // namespace
#endif //_AMINOACID_H_

