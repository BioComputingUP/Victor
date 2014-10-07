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
#include <RelSaver.h>
#include <SeqSaver.h>
#include <IoTools.h>
#include <IntCoordConverter.h>
#include <IntSaver.h>

using namespace Biopool;

int main()
{ 
  cout << "Start\n";
  
  IntCoordConverter icc;
  Spacer sp;

  ifstream refFile("reference.ref");
  ifstream inFile("pdbTest2.seq");

  if (!refFile)
    ERROR("File not found.", exception);
  if (!inFile)
    ERROR("File not found.", exception);

  SeqLoader sl(inFile, refFile);
  sp.load(sl);

  cout << "-------------------------------------------------------\n";
  RelSaver iss(cout);
  sp.save(iss);
  cout << "-------------------------------------------------------\n";
  SeqSaver sss(cout);
  sss.setWriteChi(false);
  sp.save(sss);
  cout << "-------------------------------------------------------\n";

  ofstream pdbFile("pdbRef.ent");

  if (!pdbFile)
    ERROR("File not found.", exception);
  
  PdbSaver pss(pdbFile);
  sp.save(pss);

  cout << "-------------------------------------------------------\n";

  return 0;
}
