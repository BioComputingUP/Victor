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
#include <SeqLoader.h>
#include <PdbSaver.h>
#include <XyzSaver.h>
#include <SeqSaver.h>
#include <IoTools.h>
#include <IntCoordConverter.h>
#include <SeqConstructor.h>
#include <IntSaver.h>

using namespace Biopool;

int main()
{ 
  cout << "Start" << endl;

  IntCoordConverter icc;

  Spacer sp;

  ifstream refFile("reference.ref");

  if (!refFile)
    ERROR("File not found.", exception);

  SeqConstructor sc(refFile);

  cout << "-------------------------------------------------------\n";

  sp = sc.makeSpacer("GLY", 12);

  cout << "-------------------------------------------------------\n";
  cout << "Writing XYZ file...\n";

  ofstream outFile("builder.xyz");

  if (!outFile)
    ERROR("File not found.", exception);

  XyzSaver xs(outFile);
  sp.save(xs);
  cout << "-------------------------------------------------------\n";
  cout << "Writing PDB file...\n";

  ofstream outFile2("builder.pdb");

  if (!outFile2)
    ERROR("File not found.", exception);

  PdbSaver ps(outFile2);
  sp.save(ps);
  cout << "-------------------------------------------------------\n";
  SeqSaver is(cout);
  sp.save(is);
  cout << "-------------------------------------------------------\n";


  return 0;
}
