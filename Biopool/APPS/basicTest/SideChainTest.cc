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
@Description */
#include <Spacer.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <SideChainConstructor.h>
#include <GetArg.h>

using namespace Victor;using namespace Victor::Biopool;

int main(int nArgs, char* argv[])
{ 
  cout << "Start" << endl;

  string inputFile, outputFile;
  getArg( "i", inputFile, nArgs, argv, "test.pdb");
  getArg( "o", outputFile, nArgs, argv, "sctest.pdb");

  IntCoordConverter icc;

  Spacer sp;

  ifstream inFile(inputFile.c_str());

  if (!inFile)
    ERROR("File not found.", exception);

  PdbLoader pl(inFile);
  sp.load(pl);

  ifstream refFile("allaminoacids.pdb");

  if (!refFile)
    ERROR("File not found.", exception);

  SideChainConstructor scc(refFile);
  vector<SideChain> vsc;
  
  cout << "-------------------------------------------------------\n";

  for (unsigned int i = 0; i < sp.size(); i++)
    sp.getAmino(i).setSideChain(scc.makeSideChain(sp.getAmino(i).getType()));
  
  cout << "-------------------------------------------------------\n";
  ofstream outFile(outputFile.c_str());

  if (!outFile)
    ERROR("File not found.", exception);

  PdbSaver ps(outFile);
  sp.save(ps);
  cout << "-------------------------------------------------------\n";


  return 0;
}
