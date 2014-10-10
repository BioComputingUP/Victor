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

#ifndef _STRUCTURALALIGNMENT_H_
#define _STRUCTURALALIGNMENT_H_

// Includes:
#include <string>
#include <vector>
#include <Debug.h>
#include <AlignmentBase.h>
#include <Spacer.h>
#include <vector3.h>
#include <matrix3.h>
using namespace Victor::Biopool;

using namespace Victor;
using namespace std;

namespace Victor { namespace Align2{

/** @brief   Class for structural alignments.
     * 
     * @Description  
     * @This 
     **/
    struct EData {

        EData(unsigned int _o = 0, double _d = 0.0) : other(_o), dist(_d) {
        }
        unsigned int other;
        double dist;
    };

    struct FData {

        FData(unsigned int _mi = 0, unsigned int _ma = 0, double _d = 0.0) :
        min(_mi), max(_ma), dist(_d) {
        }
        unsigned int min;
        unsigned int max;
        double dist;
    };

    class StructuralAlignment : public AlignmentBase {
    public:

        // CONSTRUCTORS/DESTRUCTOR:
        StructuralAlignment();
        StructuralAlignment(const StructuralAlignment& orig);
        virtual ~StructuralAlignment();

        // PREDICATES:

        Spacer& getTarget() {
            return spTarget;
        }

        Spacer& getTemplate() {
            return spTemplate;
        }

        // MODIFIERS:

        void setTarget(Spacer& sp) {
            spTarget = sp;
        }

        void setTemplate(Spacer& sp) {
            spTemplate = sp;
        }

        void loadCE(istream& input, Spacer& spNew);
        void buildEquivalenceNetwork();
        void buildFragmentNetwork(double maxDist = 4.0);
        // maximum CA distance for "equivalence"
        void writeData(double maxDist = 4.0);
        virtual void copy(const StructuralAlignment& orig);

        // OPERATORS:

    protected:

    private:

        // HELPERS: 

        void pExecStructAli(vgMatrix3<double> rot, vgVector3<double> trans);

        // ATTRIBUTES:
        Spacer spTarget;
        Spacer spTemplate;

    public:
        vector<EData> equivData;
        vector<FData> fragData;

    };

}} // namespace

#endif //_STRUCTURALALIGNMENT_H_
