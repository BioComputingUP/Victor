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
 * @Description pdb2totalenergy - PDB Torsion potential Energy\n\n. This program calculates a 
 * pseudo-energy to evaluate the quality of a given protein structural model, as expressed in a single "
 * */
#include <string>
#include <GetArg.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <Potential.h>
#include <TorsionPotential.h>
#include <PhiPsi.h>
#include <PhiPsiPreAngle.h>
#include <PhiPsiOmega.h>
#include <PhiPsiOmegaPreAngle.h>
#include <PhiPsiOmegaChi1.h>
#include <PhiPsiOmegaChi1PreAngle.h>
#include <PhiPsiOmegaChi1Chi2.h>
#include <PhiPsiOmegaChi1Chi2PreAngle.h>
#include <Omega.h>
#include <Chi1Chi2.h>
using namespace Victor::Biopool;
using namespace Victor::Energy;
using namespace Victor;

void sShowHelp()
{
  cout << "pdb2totalenergy - PDB Torsion potential Energy\n\n"
       << "This program calculates a pseudo-energy to evaluate the quality "
       << "of a given protein structural model, as expressed in a single "
       << "(real) number. Contact the author for further details.\n"
       << "Developed by Silvio Tosatto <silvio.tosatto@unipd.it> and "
       << "Alessandro Albiero <alex@cribi.unipd.it> \n"
       << "The program may consider different angles in a dependent" 
       << "or indipendent way. Without options there are considered all"
       << "angle in a dependent way.\n\n"
       << "INPUT OPTIONS:\n"
       << "\t-i <filename> \t\t Input  PDB file." <<"\n"
       << "\t-I <filename> \t\t Input a filelist of PDB file.\n"
       << "\t[-f]          \t\t Force model to be loaded.\n"
       << "\t-c  <id>  \t ID of chain to load from PDB file\n"
       << "\t[--nmr] \t\t Calculate average over NMR ensemble ( only for NMR model ).\n"
       << "\t[-v] \t\t\t Verbose mode\n"
       << "\t[-p] \t\t\t Per Residue value\n"
       << "\t[-s  <int> ]\t Starting Residue. (DEFAULT = the second aminoacid)\n"
       << "\t                 \t with prengle DAFAULT is the third.\n"
       << "\t[-e    <int> ]\t Ending Residue. (DEFAULT = the last -1)\n\n"
       << "OPTION FOR DIPENDENTLY ANGLES:\n"
       << "\t--pp      \t to consider Phi and Phi Angles.\n"
       << "\t--ppo     \t to consider Phi, Phi, Omega Angles.\n"
       << "\t--ppoc    \t to consider Phi, Psi, Omega and Chi1 Angles.\n"
       << "\t--ppocc   \t to consider Phi, Psi, Omega, Chi1 and Chi2 Angles.\n\n"
       << "\t--pangle\t to consider PreAngle phi and psi (adding to previous angles ).\n\n"
       << "\t--complete\t to consider Phi, Psi, Omega, Chi1, Chi2, Pre-Phi and Pre-Psi angles.(DEFAULT).\n\n"
       << "\nOPTIONS FOR INDIPENDENTLY ANGLES:\n"
       << "\t--omega   \t to consider Omega indipendently.\n"
       << "\t--chi     \t to consider Chi angles indipendently.\n"
       << "\nRANGE FOR PHI/PSI ANGLES:\n"
       << "\t--setangle  <integer> \t the bin for phi and psi angles (Default = 20�).\n"
       << "\t--setangle2 <integer> \t the bin for prephi and prepsi angles (Default = 40� ).\n\n"
       << "KNOWLEDGE OPTIONS:\n"
       << "\t--know <filename> \t\t Knowledge filename (Dafault = TOP500 ).\n"
       << "                         nmr and X-ray Knowlwdge must be used WITHOUT pre-angle!!!!\n\n"
       << "RELATIVE ENERGY: \n"
       << "\t-R                 \t to consider a the relative value for the structure considered.\n"
       << "\t--allchains        \t to consider all chain of a pdb file.\n"
       << "\t                   \t Relative energy in this casa will be the average.\n"
       << "\t--cutoff <integer> \t The length of chain that will not considered.\n"
       << "\t                   \t Usable only with --allchain option. DEFAULT 20.\n"
       << "\t                   \t For more information on Relative energy type --relative\n";
}

void sShowSpecial()
{
  cout << "\n\nRELATIVE ENERGY:\n"
       << "Relative Energy is a sort of normalized energy using the best "
       << "energy value according to the the type of aminoacid.\n"
       << "The best energy value depends on the number and type of angles " 
       << "considered. More problematic angles are the pre-Angle\n(in this "
       << "case PrePhi and PrePsi) and in this version they are not considered.\n"
       << "The best value are obtained considering all angles that is phi, psi, omega,"
       << "chi1 and chi2.\n"
       << "There are differences between the sigle relative energy value and the per"
       << "residue relative energy value.\nIn the first case it is a log and the number"
       << "represent the exponent tho give to the best for having the model observed.\n"
       << "In the second case the number obtained for every residue is not a log and it "
       << "represent a multiplicative factor to give\nto the best aminoacid for obtaining "
       << "the observed.\n"
       << "The --allchains option is used for considering all chain in a pdb file. "
       << "In this case relative energy for this file\nis a sort of average value for the "
       << "entire cristall. The cutoff option is used for not considering the chains with a\n"
       << "length <= of the cutoff specified. It is used for not including in the relative"
       << " energy small peptide present in\nthe cristall but not part of the structure.\n\n\n";
}

void Header(string inputFile, bool verbose, bool nmr, int i)
{
  if (verbose)
    {
      cout << inputFile;
      if (nmr)
	cout <<" mod. "<< i;
      cout <<"\t";
    }
}

int 
main(int nArgs, char* argv[])
{
   if (nArgs == 1)
   cout << "SPECIFY AT LEAST INPUT PDB FILE\n";

  if (getArg( "h", nArgs, argv))
    {
      sShowHelp();
      return 1;
    };

  if (getArg( "-relative", nArgs, argv))
    {
      sShowSpecial();
      return 1;
    };

  //read options from command line
  string inputFile, listFile, know, chainID;
  bool verbose = getArg( "v", nArgs, argv);
  bool force = getArg( "f", nArgs, argv);
  getArg( "i", inputFile, nArgs, argv, "!");
  getArg( "I", listFile, nArgs, argv, "!");
  bool perres = getArg("p", nArgs, argv);    //per residue value
  
  //knowledge filename, default is TOP500
  getArg( "-know", know, nArgs, argv, "!");
  
  //which angle to consider
  bool pp = getArg("-pp", nArgs, argv);      
  bool ppo = getArg("-ppo", nArgs, argv);
  bool ppoc = getArg("-ppoc", nArgs, argv);
  bool ppocc = getArg("-ppocc", nArgs, argv);
  bool pangle = getArg("-pangle", nArgs, argv);
  bool complete = getArg("-complete", nArgs, argv);
  
  //angles chosen in an independent way, (omega and chi only)
  bool omega = getArg("-omega", nArgs, argv);
  bool chi = getArg("-chi", nArgs, argv);
  int ARCSTEP1, ARCSTEP2;
  unsigned int CUTOFF;
  unsigned int START = 2, END = 9999;
  getArg( "s", START, nArgs, argv, 2 );
  getArg( "e", END, nArgs, argv, 0 );
  getArg( "-setangle", ARCSTEP1, nArgs, argv, 20);
  getArg( "-setangle2", ARCSTEP2, nArgs, argv, 40);
  
  bool relative = getArg("R", nArgs, argv);
  getArg( "c", chainID, nArgs, argv, " ");
  bool nmr = getArg( "-nmr", nArgs, argv);
  bool allchains = getArg( "-allchains", nArgs, argv);
  getArg( "-cutoff", CUTOFF, nArgs, argv, 20);
  int MEMORYEND = END;
  
  if ( !pangle )
    START -= 1;
  if ( START < 1 )
    ERROR("ERROR. Not possible to start before the second residue.",exception);
  
  if ( ( allchains ) && ( chainID != " " ) )
    ERROR("ERROR. You have choose two option not compatible: --allchains and --chain.", exception);
  
  //no angle chosen, default: choose all
  if ( (!pp) && (!ppo) && (!ppoc) && (!ppocc) && (!complete) && (!pangle) &&
       (!omega) && (!chi))
    complete = true;
  
  if ( ((pp) && (ppo)) || ((pp) && (ppoc)) || ((pp) && (ppocc)) || ((pp) && (complete)) ||
       ((ppo) && (ppoc)) || ((ppo) && (ppocc)) || ((ppo) && (complete)) || 
       ((ppoc) && (ppocc)) || ((ppoc) && (complete)) || ((ppocc) && (complete)) )
    ERROR("ERROR in choosing options. You have chosen too many set of angles.", exception);

  if ( ( (complete) && ( (omega) || (chi) || (pangle) ) ) || 
       ( ( (ppocc) || (ppoc) ) && ( ( (omega) || (chi) ) ) ) ||
       ( (ppo) && (omega) ))
    ERROR("ERROR. You have chosen the same angles in dipendent and indipendent way.", exception);
  
   
  if ( ( ( 360%ARCSTEP1 ) != 0 ) || (360%ARCSTEP2) != 0 )
    ERROR("The range of angles (phi, psi, pre-phi and/or pre-psi) must be divisor of 360.", exception);
  
  if ( ( relative) && ( ( omega ) || ( chi ) ) )
    ERROR("Relative Energy may be used only without indipendent Angle.", exception);

  TorsionPotential* pot = NULL;
  TorsionPotential* omegaPot = NULL;//for omega indipendent
  TorsionPotential* chiPot = NULL;// for chi indipendent
  
  //Default knowledge
  if ( know == "!" )
    {  char *victor=getenv("VICTOR_ROOT");
	if (victor  == NULL)
		ERROR("Environment variable VICTOR_ROOT was not found.\n Use the command:\n export VICTOR_ROOT=......", exception);
      string env;
      env = getenv("VICTOR_ROOT");
      string tmpknow = "data/tor.par";
      string defaultknow = env + tmpknow;
      know = defaultknow;
      if (defaultknow.length() < 3 )
	ERROR("Environment variable VICTOR_ROOT was not found.", exception);
      
    }
  
  //torsion potential to calculate
  if ( pp )
    {
      if ( !pangle )
	pot =  new PhiPsi(ARCSTEP1, know);
      else if ( pangle )
	pot = new PhiPsiPreAngle(ARCSTEP1, ARCSTEP2, know);
    }
  else if ( ppo )
    {
      if ( !pangle )
	pot = new PhiPsiOmega(ARCSTEP1, know);
      else if ( pangle )
	pot = new PhiPsiOmegaPreAngle(ARCSTEP1, ARCSTEP2, know);
    }
  else if ( ppoc )
    {
      if ( !pangle )
	pot = new PhiPsiOmegaChi1(ARCSTEP1, know);
      else if ( pangle )
	pot = new PhiPsiOmegaChi1PreAngle(ARCSTEP1, ARCSTEP2, know);
    }
  else if ( ppocc )
    {
      if ( !pangle )
	pot = new PhiPsiOmegaChi1Chi2(ARCSTEP1, know);
      else if ( pangle )
	pot = new PhiPsiOmegaChi1Chi2PreAngle(ARCSTEP1, ARCSTEP2, know);
    }
  else if ( complete )
    pot = new PhiPsiOmegaChi1Chi2PreAngle(ARCSTEP1, ARCSTEP2, know);
  
  if ( omega )
    omegaPot = new Omega(know);
  if ( chi )
    chiPot = new Chi1Chi2(know);
  
  ifstream inFile(listFile.c_str());
  if ((!inFile) && (listFile != "!"))
    ERROR("File " + listFile + " not found.", exception);
  
  long double tres = 0.00; // the TorsionPotential propensities calculated
  long double total = 0.00; // the theoric best TorsionPotential
  int nmrmodel = 0; // the number of nmr model.
  
  
  while ((inFile) || (inputFile != "!"))
    {
      if (listFile != "!")
	{
	  inFile >> inputFile;
	  if (!inFile)
	    break;
	}
      
      //initialization of a vector containing the IDs of all chains
      vector<char> totalchain;
      ifstream pdbfile2(inputFile.c_str());
      
      //load components in standard PDB format from pdbfile2
      PdbLoader pl(pdbfile2);
      if ( allchains )
	{
	  //insert all ID of chains present in PDB file
	  totalchain = pl.getAllChains();
	  if ( totalchain.size() == 0 )
	    totalchain.push_back(chainID.c_str()[0]); 
	}
      else
	{
	  totalchain.push_back(chainID.c_str()[0]);
	}
      //maximum number of nmr models
      unsigned int max;
      if (!nmr)
	max = 1;
      else
	{
	  max = pl.getMaxModels();
	  if ( max == 0)
	    {
	      max = 1;
	      if ( verbose )
		cerr <<"Warning: the file "<<inputFile<<" probably is not an nmr structure.\n";
	    }
	}
      
      vector<unsigned int> chainleng; // vector with all length of chains requested
      
      //for each nmr model
      for (unsigned int i = 1; i <= max; i++)
	{
	  chainleng.clear();
	  
	  //for each chain
	  for ( unsigned int ch = 0; ch < totalchain.size(); ch++ )
	    {
	      PdbLoader pl(pdbfile2);
	      pl.setChain(totalchain[ch]);
	      pl.setNoHAtoms();
	      pl.setNoVerbose();
	      pl.setPermissive();
	      nmrmodel = i;

	      if (!force)
		pl.setModel(i); 
	      Protein prot;
              prot.load(pl);
              Spacer *sp;
              sp=prot.getSpacer(totalchain[ch]);
	      	
	      if (!pl.isValid())
		{
		  if (verbose)
		    {
		      cerr << "Warning: Invalid PDB file found "<<inputFile<<" Chain: "<<totalchain[ch];
		      if (nmr)
			cerr <<" model: "<<nmrmodel<<" \n";
		      else
			cerr <<" \n";
		    }
		  if ( ( i == max ) && ( ch == ( totalchain.size()-1 ) ))
		    inputFile = "!";
		  continue;
		}
	      //CUTOFF length of chain not to be considered
	      if ( sp->sizeAmino() <= CUTOFF )
		{
		  cerr << "Warning: "<<inputFile<<" Chain: "<<totalchain[ch]<<" under CUTOFF\n";
		  chainleng.push_back(0);
		  continue;
		}
	      else
		//insert chainlength for each chain 
		chainleng.push_back(sp->sizeAmino());
	      
	      if ( END == 0 )
		//number of all aminoacids in the spacer
		END = ( sp->sizeAmino()-1 );
	      else if ( END >= sp->sizeAmino() )
		ERROR("END out of protein max length permitted..", exception);
		
	      cout.setf(ios::fixed, ios::floatfield);
	      
	      //non-normalized pseudo energy
	      if ( !relative )
		{
		  if (perres)
		    {
		      //for each residue 
		      for ( unsigned int i = START; i < END; i++)
			{
			  int num = sp->getPdbNumberFromIndex(i);
			  unsigned int j = i+1;
			  
			  //no independent angle
			  if ( ( !omega ) && ( !chi ) )
			    tres = pot->calculateEnergy(*sp, i, j);
			  //omega independent angle
			  if ( ( omega ) && ( !chi ) )
			    tres = ( pot->calculateEnergy(*sp, i, j)) 
			      + (omegaPot->calculateEnergy(*sp, i, j));
			  //chi independent angle
			  if ( ( !omega) && ( chi ) )
			    tres = ( pot->calculateEnergy(*sp, i, j)) 
			      + (chiPot->calculateEnergy(*sp, i, j));
			  //both independent angles
			  if ( (omega) && (chi) )
			    tres = ( pot->calculateEnergy(*sp, i, j)) 
			      + (chiPot->calculateEnergy(*sp, i, j))
			      + (omegaPot->calculateEnergy(*sp, i, j));
			  Header( inputFile, verbose, nmr, nmrmodel );
			  cout <<num<<"\t"<<tres;
			  if ( verbose )
			    cout <<"\t"<<totalchain[ch];
			  cout <<endl;
			}
		    }
		  else //whole spacer
		    {
		      if ( ( !omega ) && ( !chi ) )
			tres = pot->calculateEnergy(*sp, START, END);
		      if ( ( omega ) && ( !chi ) )
			tres = ( pot->calculateEnergy(*sp, START, END)) + (omegaPot->calculateEnergy(*sp, START, END));
		      if ( ( !omega) && ( chi ) )
			tres = ( pot->calculateEnergy(*sp, START, END)) + (chiPot->calculateEnergy(*sp, START, END));
		      if ( (omega) && (chi) )
			tres = ( pot->calculateEnergy(*sp, START, END)) + 
			  (chiPot->calculateEnergy(*sp, START, END)) + (omegaPot->calculateEnergy(*sp, START, END));
		      Header( inputFile, verbose, nmr, nmrmodel );
		      cout << setw(9) << setprecision(4) << (tres) << "\t"; 
		      if (verbose)
			cout <<"\t"<<totalchain[ch]<<"\tsize: "<<sp->sizeAmino();
		      cout << endl;
		    }
		}
	      // calculate normalized energy
	      else if (relative)  
		{
		  // ------------------------------------
		  // relative energy, i.e. TAP score case
		  // ------------------------------------

		  if (perres)
		    { 
		      if (!pangle)
			{
			  //for each residue
			  for (unsigned int i = START; i < END; i++)
			    {
			      AminoAcidCode code = 
				aminoAcidThreeLetterTranslator(sp->getAmino(i).getType()); 
			      int num = sp->getPdbNumberFromIndex(i);
			      unsigned int j = i+1;
			      tres = pot->calculateEnergy(*sp, i, j);
			      Header( inputFile, verbose, nmr, nmrmodel );
// 			    
			      cout << num <<"\t"<<tres/(-log(pot->pReturnMaxPropensities(code)));
			      if ( verbose )
				cout<<"\t"<<totalchain[ch];
			      cout <<endl;
			    }
			}
		      else // if (pangle)
			{
			  for (unsigned int i = START; i < END; i++)
			    {
			      AminoAcidCode code = 
				aminoAcidThreeLetterTranslator(sp->getAmino(i).getType()); 
			      int num = sp->getPdbNumberFromIndex(i);
			      unsigned int j = i+1;
			      tres = (-exp(pot->calculateEnergy(*sp, i, j)));
			      Header( inputFile, verbose, nmr, nmrmodel);
			      cout <<"\n" << num <<"\t" <<tres/pot->
				pReturnMaxPropensitiesPreAngle
				(code,
				 pot->sGetPropBin2(sp->getAmino(i).getInBond(0).getPhi()),
				 pot->sGetPropBin2(sp->getAmino(i).getInBond(0).getPsi()));
			      if ( verbose )
				cout<<"\t"<<totalchain[ch];
			      cout <<endl;
			    }
			}
		    }
		  else //(!perres)
		    {
		      tres = tres + pot->calculateEnergy(*sp, START, END);
		      if (pangle)
			{
			  for ( unsigned int i = START; i < END; i++)
			    {
			      unsigned int code = sp->getAmino(i).getCode();

			      total = total + (-log(pot->pReturnMaxPropensitiesPreAngle
						    (static_cast<int>(code),
						     pot->sGetPropBin2(sp->getAmino(i).getInBond(0).getPhi()),
						     pot->sGetPropBin2(sp->getAmino(i).getInBond(0).getPsi())
						     )));
			    }
			}
		      else //(!pangle)
			{
			  for ( unsigned int i = START; i < END; i++)
			    {
			      unsigned int code = sp->getAmino(i).getCode(); 
			      total = total + (-log(pot->pReturnMaxPropensities(static_cast<int>(code))));
			    }
			}
		    }
		}
	      END = MEMORYEND;
	    }
	  if ( ( relative ) && (!perres)  )
	    {
	      Header( inputFile, verbose, nmr, nmrmodel);
	      cout << setw(9) << setprecision(4) <<tres/total;
	      if ( verbose )
		{
		  if ( !allchains )
		    cout <<"\t"<<totalchain[0] ; // <<"\tsize: "<<chainleng[0];
		  else
		    {
		      cout <<"\t all chains\t size:";
		      for ( unsigned int i = 0; i < chainleng.size(); i++ )
			cout <<" "<<totalchain[i] ; // <<"("<<chainleng[i]<<")";
		    }
		}
	      cout << endl;
	    }
	  total = 0.00;
	  tres = 0.00;
	}
      
      pdbfile2.close();
      // reset variable to trigger break condition:
      inputFile = "!";
    }
  //  delete pot;
  //  delete omegaPot;
  //  delete chiPot;
}
