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
#include <CifLoader.h>
#include <CifSaver.h>
#include <XyzSaver.h>
#include <SeqSaver.h>
#include <IoTools.h>
#include <IntCoordConverter.h>
#include <IntSaver.h>

using namespace Victor;
using namespace Victor::Biopool;

/*
 * 
 */
int main(int argc, char** argv) {

    if (argc != 2) {
        cout << "Cif Corrector $Revision: 0.1 $ -- adds missing oxygen atoms to "
                << "protein structure backbones" << endl;
        cout << "  Usage: \t\t CifCorrector <filename> \n";
        return 1;
    };

    ifstream inFile(argv[1]);
    if (!inFile)
        ERROR("File not found.", exception);

    CifLoader il(inFile);
    Protein prot;
    prot.load(il);
    unsigned int zero = 0;
    Spacer sp = *(prot.getSpacer(zero));

    for (unsigned int i = 0; i < sp.sizeAmino(); i++)
        sp.getAmino(i).addMissingO();
    
    inFile.close();

    ofstream outFile(argv[1]);

    if (!outFile)
        ERROR("Couldn't write file.", exception);

    CifSaver pss(outFile);
    sp.save(pss);
    
    outFile.close();
    
    return 0;
}

