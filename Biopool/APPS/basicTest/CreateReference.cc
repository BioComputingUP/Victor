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
#include <Atom.h>
#include <Spacer.h>
#include <Group.h>
#include <SideChain.h>
#include <vector3.h>
#include <SeqLoader.h>
#include <PdbSaver.h>
#include <PdbLoader.h>
#include <XyzSaver.h>
#include <SeqSaver.h>
#include <IoTools.h>
#include <IntCoordConverter.h>
#include <IntSaver.h>
#include <RelSaver.h>

using namespace Victor;using namespace Victor::Biopool;

int main() {
    cout << "Create Reference \n   $Revision: 1.2 $:\n";
    cout << "-------------------------------------------------------\n";

    IntCoordConverter icc;
    Spacer sp;

    ifstream inFile("allaminoacids.pdb");

    if (!inFile)
        ERROR("File not found.", exception);

    PdbLoader pl(inFile);
    sp.load(pl);

    cout << "-------------------------------------------------------\n";
    PdbSaver ps(cout);
    sp.save(ps);
    cout << "-------------------------------------------------------\n";
    ofstream outFile22("test2.xyz");

    if (!outFile22)
        ERROR("File not found.", exception);

    XyzSaver xss2(outFile22);

    cout << "-------------------------------------------------------\n";
    ofstream outFile("reference.ref");

    if (!outFile)
        ERROR("File not found.", exception);

    RelSaver iss(outFile);

    for (unsigned int i = 0; i < sp.sizeAmino(); i++) {
        sp.getAmino(i)[N].bindOut(sp.getAmino(i)[CA]);
        sp.getAmino(i)[C].bindOut(sp.getAmino(i)[O]);
        sp.getAmino(i)[CA].bindOut(sp.getAmino(i)[C]);

        AminoAcid aa = sp.getAmino(i);

        if (i < sp.sizeAmino() - 1) {
            Atom at = sp.getAmino(i + 1)[N];

            at.setCode(OXT);
            at.unbindIn(at.getInBond(0));
            at.bindIn(aa[C]);
            aa.addAtom(at);
            for (unsigned int j = 0; j < aa[OXT].sizeOutBonds(); j++)
                aa[OXT].unbindOut(aa[OXT].getOutBond(j));
            sp.getAmino(i + 1)[N].bindIn(sp.getAmino(i)[C]);
        }

        vgVector3<double> t1(0, 1, 0);
        vgMatrix3<double> res(1);
        if (aa[N].getTrans().length() != 0) {
            alignVectors(t1, aa[N].getTrans(), res);
            aa[N].addRot(res);
            res = vgMatrix3<double>::createRotationMatrix(t1, 180 * DEG2RAD);
            aa.sync();
        };

        aa.save(xss2); // removing this will make the calculations go wrong???
        aa.save(iss);
    }
    cout << "-------------------------------------------------------\n";

    return 0;
}
