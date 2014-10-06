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
#include <AminoAcid.h>
#include <Group.h>
#include <vector3.h>
#include <XyzLoader.h>
#include <XyzSaver.h>
#include <SeqSaver.h>
#include <IoTools.h>
#include <IntSaver.h>

using namespace Biopool;

int main()
{ 
  cout << "Start" << endl;

  AminoAcid aa;
  ifstream inFile("amino.xyz");
  if (!inFile)
    ERROR("File not found.", exception);

  XyzLoader il(inFile);
  aa.load(il);

  cout << "-------------------------------------------------------\n";
  XyzSaver is(cout, 0);
  aa.save(is);

  cout << "-------------------------------------------------------\n";
  cout << "SideChain->backbone = " << aa.getSideChain().getBackboneRef()->getType() << "\n";

  cout << "-------------------------------------------------------\n";
  SeqSaver iss(cout);
  aa.save(iss);
  cout << "-------------------------------------------------------\n";
  cout << "superior: " << "\t aa[2]= " << aa[2].getSuperior().getType()
       << " (" << aa[2].getSuperior().getCode() << ") \t aa[8]= " 
       << aa[8].getSuperior().getType() 
       << " (" << aa[8].getSuperior().getCode() << ")" 
       << endl;

  cout << "-------------------------------------------------------\n";
  IntSaver intS(cout);
  aa.save(intS);
  cout << "-------------------------------------------------------\n";
  return 0;
}
