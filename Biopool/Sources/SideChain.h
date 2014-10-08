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


#ifndef _SIDECHAIN_H_
#define _SIDECHAIN_H_
#include <vector>
#include <Group.h>
#include <AminoAcidCode.h>
#include <IntCoordConverter.h>
#include <Visitor.h>

namespace Biopool {

    /**
     * @brief This class implements a side chain
     * 
     * */
    class SideChain : public Group {
    public:
        // CONSTRUCTORS/DESTRUCTOR:
        SideChain();
        SideChain(const SideChain& orig);
        virtual ~SideChain();

        // PREDICATES:

        virtual string getClassName() const {
            return "SideChain";
        }
        virtual unsigned int getCode() const;
        virtual string getType() const; // returns AA 3-L-code
        virtual string getExtendedType() const; // returns full code

        double getChi(unsigned int n);
        vector<double> getChi();
        unsigned int getMaxChi();

        AminoAcid* getBackboneRef();

        bool isMember(const AtomCode& ac) const;

        void save(Saver& s); // data saver

        // MODIFIERS:
        void setChi(unsigned int n, double c);
        void setChi(vector<double> cv);
        void setBackboneRef(AminoAcid* br);

        void copy(const SideChain& orig);
        void load(Loader& l); // data loader

        void patchAminoAcidCode(); // determines the 3-letter code from the sidechain
        bool setBondsFromPdbCode(bool connect, bool permissive = false);
        // returns false if failed to connect

        virtual const SideChain& getInBond(unsigned int n) const;
        virtual SideChain& getInBond(unsigned int n);
        virtual const SideChain& getOutBond(unsigned int n) const;
        virtual SideChain& getOutBond(unsigned int n);

        virtual void acceptCalculator(EnergyVisitor* v);
        virtual void acceptOptimizer(OptimizationVisitor* v);

        // OPERATORS:
        SideChain& operator=(const SideChain& orig);
        Atom& operator[](unsigned int n);
        const Atom& operator[](unsigned int n) const;
        Atom& operator[](const AtomCode& ac);
        const Atom& operator[](const AtomCode& ac) const;

    protected:

    private:
        // HELPERS:
        bool findCode(char c, unsigned int& pos, unsigned int& pos2);
        void findHCode(char c, char c2, Atom& at, bool connect);
        void setMaxChiFromCode();
        void setMaxChi(unsigned int _maxChi);
        double calculateChi(unsigned int n);
        void convertChi(unsigned int n, double a);

        bool hasB();
        bool hasG();
        bool hasD();
        bool hasE();
        bool hasZ();
        bool hasH();

        Atom& firstG();
        Atom& firstD();
        Atom& firstE();
        Atom& firstZ();
        Atom& firstH();

        // ATTRIBUTES:
        unsigned int maxChi; // number of chi torsion angles 
        vector<double> chi; // torsion angles
        AminoAcid* backboneRef; // reference to backbone (aminoacid) of this

        friend class AminoAcid;
    };

    // ---------------------------------------------------------------------------
    //                                    Group
    // -----------------x-------------------x-------------------x-----------------

    // PREDICATES:

    inline unsigned int SideChain::getCode() const {
        return aminoAcidThreeLetterTranslator(getType());
    }

    inline string SideChain::getType() const {
        string tmp = id;
        return tmp.substr(0, 3);
    }

    inline string SideChain::getExtendedType() const {
        return id;
    }

    inline double SideChain::getChi(unsigned int n) {
        PRECOND(n < chi.size(), exception);
        if (chi[n] > 990)
            chi[n] = RAD2DEG * calculateChi(n);
        return chi[n];
    }

    inline unsigned int SideChain::getMaxChi() {
        return chi.size();
    }

    inline AminoAcid* SideChain::getBackboneRef() {
        if (backboneRef == NULL)
            DEBUG_MSG("Sidechain has no assigned backbone.");
        return backboneRef;
    }

    inline bool SideChain::isMember(const AtomCode& ac) const {
        return (pGetAtom(ac) != NULL);
    }

    inline void SideChain::save(Saver& s) {
        s.saveGroup(*this);
    }

    inline void SideChain::acceptCalculator(EnergyVisitor* v) {
        v->PrepareSideChain(*this);
    }

    inline void SideChain::acceptOptimizer(OptimizationVisitor* v) {
        v->PrepareSideChain(*this);
    }


    // MODIFIERS:

    inline void SideChain::setChi(unsigned int n, double c) {
        PRECOND((n < chi.size()) && (c >= -180) && (c <= 180), exception);
        convertChi(n, DEG2RAD * c);
        chi[n] = c;
    }

    inline void SideChain::load(Loader& l) {
        l.loadSideChain(*this);
        resetBoundaries();
    }

    inline const SideChain& SideChain::getInBond(unsigned int n) const {
        return dynamic_cast<const SideChain&> (Bond::getInBond(n));
    }

    inline SideChain& SideChain::getInBond(unsigned int n) {
        return dynamic_cast<SideChain&> (Bond::getInBond(n));
    }

    inline const SideChain& SideChain::getOutBond(unsigned int n) const {
        return dynamic_cast<const SideChain&> (Bond::getOutBond(n));
    }

    inline SideChain& SideChain::getOutBond(unsigned int n) {
        return dynamic_cast<SideChain&> (Bond::getOutBond(n));
    }

    // OPERATORS:
    inline Atom& SideChain::operator[](unsigned int n) {
        return (Group::getAtom(n));
    }

    inline const Atom& SideChain::operator[](unsigned int n) const {
        return (Group::getAtom(n));
    }

    inline Atom& SideChain::operator[](const AtomCode& ac) {
        Atom* a = pGetAtom(ac);
        INVARIANT(a != NULL, exception);
        return *a;
    }

    inline const Atom& SideChain::operator[](const AtomCode& ac) const {
        Atom* a = pGetAtom(ac);
        INVARIANT(a != NULL, exception);
        return *a;
    }

    // HELPERS:

} // namespace
#endif //_SIDECHAIN_H_

