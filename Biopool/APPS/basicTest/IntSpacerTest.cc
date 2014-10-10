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
#include "Spacer.h"
#include "Group.h"
#include "SideChain.h"
#include "vector3.h"
#include "XyzLoader.h"
#include "PdbLoader.h"
#include "PdbSaver.h"
#include "XyzSaver.h"
#include "SeqSaver.h"
#include "IoTools.h"
#include "IntCoordConverter.h"
#include "IntSaver.h"
#include "IntLoader.h"

using namespace Victor;using namespace Victor::Biopool;

int main()
{ 
  cout << "Start" << endl;

  Spacer sp;

  ifstream intFile_in("polyglycin.int");
  if (!intFile_in)
    ERROR("File not found.", exception);

  IntLoader IntL(intFile_in);
  
  cout << "-------------------------------------------------------\n";
  cout << "loading...\n";
  sp.load(IntL);
  
  cout << "-------------------------------------------------------\n";


  ofstream out5File("testPolyglycin.int");
  ofstream out4File("testPolyglycin.xyz");

  if (!out4File)
    ERROR("Couldn't write file.", exception);
  DEBUG_MSG("Next step : XyzSaver seqX(out4File);");
  if (!out5File)
    ERROR("Couldn't write file.", exception);

  DEBUG_MSG("Next step : IntSaver seqI(out5File);");
  XyzSaver seqX(out4File);
  DEBUG_MSG("Next step :  sp.save(seqX);");
  sp.save(seqX);
  cout << "XYZ data written to testPolyglycin.xyz.\n";
  cout << "-------------------------------------------------------\n";

  IntSaver seqI(out5File);
  DEBUG_MSG("Next step :  sp.save(seqI);");
  sp.save(seqI);
  cout << "Int data written to testPolyglycin.int. \n";
  cout << "-------------------------------------------------------\n";
 
  return 0;
}
