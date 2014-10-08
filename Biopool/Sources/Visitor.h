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

#ifndef _VISITOR_H_
#define _VISITOR_H_

// Includes:
#include <Debug.h>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {

    class Group;
    class SideChain;
    class AminoAcid;
    class Spacer;

    /**
     * @brief Base class implementing the visitor pattern  
     * 
     * @Description   Uses the implicit copy operator. 
     * This is mostly used for implementing optimizers not contained in the Biopool module (e.g. EnergyCalculator, SideChainPlacement).
     **/
    class EnergyVisitor {
    public:

        // CONSTRUCTORS/DESTRUCTOR:

        EnergyVisitor() {
        };

        virtual ~EnergyVisitor() {
        };

        // MODIFIERS:
        virtual void PrepareGroup(Group& group) = 0;
        virtual void PrepareSideChain(SideChain& node) = 0;
        virtual void PrepareAminoAcid(AminoAcid& node) = 0;
        virtual void PrepareSpacer(Spacer& node) = 0;

    protected:

    private:

    };

    /**@brief Base class Optimizacion Patter. 
     * 
     * @Description  Uses the implicit copy operator. This is mostly used for implementing 
     * optimizers not contained in  the Biopool module (e.g. EnergyCalculator, SideChainPlacement).
     **/
    class OptimizationVisitor {
    public:

        // CONSTRUCTORS/DESTRUCTOR:

        OptimizationVisitor() {
        };
        // 

        virtual ~OptimizationVisitor() {
        };

        // MODIFIERS:
        virtual void PrepareGroup(Group& group) = 0;
        virtual void PrepareSideChain(SideChain& node) = 0;
        virtual void PrepareAminoAcid(AminoAcid& node) = 0;
        virtual void PrepareSpacer(Spacer& node) = 0;

    protected:

    private:

    };

} // namespace
#endif //_SAVER_H_
