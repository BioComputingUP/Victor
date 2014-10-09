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

#include <Spacer.h>
#include <Group.h>
#include <SideChain.h>
#include <vector3.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <XyzSaver.h>
#include <SeqSaver.h>
#include <IoTools.h>
#include <IntCoordConverter.h>
#include <IntSaver.h>

using namespace Victor;using namespace Victor::Biopool;

int main(int nArgs, char* argv[]) {
    if (nArgs != 2) {
        cout << "Pdb Corrector $Revision: 1.2 $ -- adds missing oxygen atoms to "
                << "protein structure backbones" << endl;
        cout << "  Usage: \t\t PdbCorrector <filename> \n";
        return 1;
    };

    Spacer sp;
    ifstream inFile(argv[1]);
    if (!inFile)
        ERROR("File not found.", exception);

    PdbLoader il(inFile);
    sp.load(il);

    for (unsigned int i = 0; i < sp.sizeAmino(); i++)
        sp.getAmino(i).addMissingO();

    ofstream outFile2(argv[1]);

    if (!outFile2)
        ERROR("Couldn't write file.", exception);

    PdbSaver pss2(outFile2);
    sp.save(pss2);
    
    return 0;
}
