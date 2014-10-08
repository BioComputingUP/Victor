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
/**
 *@Description A module to determine the solvent exposure/accessibility of residues in a
 *                              protein fragment.
 */
#ifndef __SolvExpos_H__
#define __SolvExpos_H__

#include <Spacer.h>

namespace Biopool {

    class SolvExpos {
    public:

        enum SolvExposEnum {
            CORE, EXPOSED
        };
        SolvExpos();
        ~SolvExpos();
        
        SolvExposEnum getSolvExpos(Spacer &chain, const unsigned int tgt,
                const unsigned int start, const unsigned int end);

        vector<SolvExposEnum>* getSolvExposVec(Spacer &chain,
                const unsigned int tgtS, const unsigned int tgtE,
                const unsigned int envS, const unsigned int envE);

        double getSolvAccess(Spacer &chain, unsigned int tgt,
                unsigned int start, unsigned int end);

        vector<double> getSolvAccessVec(Spacer &chain,
                unsigned int tgtS, unsigned int tgtE,
                unsigned int envS, unsigned int envE);
    private:
        
    };
} // namespace

#endif
