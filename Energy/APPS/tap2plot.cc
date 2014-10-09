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
// -*- C++ -*-----------------------------------------------------------------

#include <string>
#include <GetArg.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <SolvationPotential.h>
#include <EffectiveSolvationPotential.h>
#include <RapdfPotential.h>
#include <TorsionPotential.h>
#include <PhiPsi.h>
#include <PhiPsiOmegaChi1Chi2.h>
#include <AminoAcidCode.h>
using namespace Victor::Biopool;
using namespace Victor::Energy;
using namespace Victor;

void sShowHelp(){
  cout << "tap2plot\n\n"
       << "\n";
}
int main(int nArgs, char* argv[]){ 
  if (getArg( "h", nArgs, argv)){
      sShowHelp();
      return 1;
    };

  unsigned int arcstep;
  getArg( "a", arcstep, nArgs, argv, 10);

  PhiPsiOmegaChi1Chi2 tors(arcstep);

  Spacer *sp;
  ifstream inFile2("../data/allaminoacids_new.pdb");
  if (!inFile2)
    ERROR("File not found.", exception);
  PdbLoader pl(inFile2);
  pl.setNoHAtoms();
  pl.setNoVerbose();
  
  pl.setPermissive();
  Protein prot;
  prot.load(pl);
  unsigned int i=0;
  sp=prot.getSpacer(i);
  
 
  for (unsigned int offset = 0; offset < 20; offset++)    {
      if (sp->getAmino(offset).getMaxChi() < 1)
	continue;

      for (int chi = -180; chi <= 60; chi += 120){
	  string s = sp->getAmino(offset).getType() + "_";
	  if (chi == -180)
	    s += "trans.dat";
	  else if (chi == -60)
	    s += "gminus.dat";
	  else
	    s += "gplus.dat";
	  
	  ofstream out(s.c_str());
	  if (!out)
	    ERROR("Could not open file for writing.", exception);

	  sp->getAmino(offset).setOmega(179);
	  sp->getAmino(offset).setChi(0,chi);
	  	  
	  for (int i = -180; i < 180; i += arcstep)
	    for (int j = -180; j < 180; j += arcstep)  {
		sp->getAmino(offset).setPhi(i+1);
		sp->getAmino(offset).setPsi(j+1);
		sp->sync();
		
		out << ( tors.calculateEnergy(sp->getAmino(offset)) 
			  / (-log(tors.pReturnMaxPropensities(offset))) )
		     << "\n";
	      }
	  out.close();
	}
    }
}
