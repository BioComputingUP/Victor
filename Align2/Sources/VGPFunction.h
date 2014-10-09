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


#ifndef __VGPFunction_H__
#define __VGPFunction_H__

#include <GapFunction.h>
#include <PdbLoader.h>
#include <SolvExpos.h>
#include <Spacer.h>
#include <vector3.h>
#include <math.h>

namespace Victor { namespace Align2{

    /** @brief   Implement VGP (Variable Gap Penalty) function. 
     * 
     * @Description  Some
     *                  explanations can be found in:
     *                 Madhusudhan MS., Marti-Renom MA., Sanchez R., Sali A.
     *                  Variable gap penalty for protein sequence-structure alignment.
     *                  Department of Biopharmaceutical Sciences and Pharmaceutical
     *                  Chemistry, University of California at San Francisco, 94143, USA.
     *                  PMID: 16423846 [PubMed - indexed for MEDLINE]
     *
     * @This 
     **/
    class VGPFunction : public GapFunction {
    public:

        // CONSTRUCTORS:

        /// Default constructor.
        VGPFunction(string pdbFileName, string chainID);

        /// Constructor assigning o, e and weights.
        VGPFunction(string pdbFileName, string chainID, double o, double e, unsigned int extType,
                double wH, double wS, double wB, double wC, double wD);

        /// Copy constructor.
        VGPFunction(const VGPFunction &orig);

        /// Destructor.
        virtual ~VGPFunction();


        // OPERATORS:

        /// Assignment operator.
        VGPFunction& operator =(const VGPFunction &orig);


        // PREDICATES:

        /// Return open gap penalty for template position p.
        virtual double getOpenPenalty(int p);

        /// Return extension gap penalty for template position p.
        virtual double getExtensionPenalty(int p);


        // MODIFIERS:

        /// Copy orig object to this object ("deep copy").
        virtual void copy(const VGPFunction &orig);

        /// Construct a new "deep copy" of this object.
        virtual VGPFunction* newCopy();

        /// Set open gap penalty.
        virtual void setOpenPenalty(double pen);

        /// Set extension gap penalty.
        virtual void setExtensionPenalty(double pen);


        // HELPERS:

        /// Extract structural infos from PDB template file.
        void pExtractPdbInfo(string pdbFileName, string chainID);


    protected:


    private:

        // ATTRIBUTES:

        double o; ///< Open gap penalty.
        double e; ///< Extension gap penalty.
        unsigned int extType; ///< Extension gap type.
        unsigned int extCounter; ///< Extension gap counter.
        vector<double> hContent; ///< Template helical content.
        vector<double> sContent; ///< Template strand content.
        vector<double> solvAccess; ///< Template solvent accessibility.
        vector<double> bbStraight; ///< Template backbone straightness.
        vector<double> spaceProx; ///< Template space proximity.
        double wH; ///< Weight for helical content.
        double wS; ///< Weight for strand content.
        double wB; ///< Weight for solvent accessibility.
        double wC; ///< Weight for backbone straightness.
        double wD; ///< Weight for space proximity.

    };

    // -----------------------------------------------------------------------------
    //                                 VGPFunction
    // -----------------------------------------------------------------------------

    // MODIFIERS:

    inline void
    VGPFunction::setOpenPenalty(double pen) {
        o = pen;
    }

    inline void
    VGPFunction::setExtensionPenalty(double pen) {
        e = pen;
    }

}} // namespace

#endif
