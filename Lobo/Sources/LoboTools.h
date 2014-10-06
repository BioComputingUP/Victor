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
 *@Description: This head file contains a set of methods that allows to manage the options for the utils
  *      
*/


#ifndef _LOBO_TOOLS_H_
#define _LOBO_TOOLS_H_
#include <string>
#include <LoopModel.h>
#include <GetArg.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <iostream>
using namespace Biopool;

bool silent = false;
/**
*@Description prints a line
*/
inline void fillLine(){
   std::cout << "---------------------------------------------------\n";
}
/**
*@Description prints the common options of the apps
*/
inline void showCommonOptions(){
     char *victor=getenv("VICTOR_ROOT");
	if (victor  == NULL)
		ERROR("Environment variable VICTOR_ROOT was not found.\n Use the command:\n export VICTOR_ROOT=......", exception);
    string path="data/aa";
   std::cout << "\t[--sol1 <n-solutions>] \t Number of primary solutions"
       << " (default= " << LoopModel::MAX_ITER_SOL <<")\n"
       << "\t[--maxWrite <n-sol.>] \t Number of solutions to write out"
       << " (default = sol1)\n"
       << "\t[--scwrl <filename>] \t Sequence output file for SCWRL (-s option)\n"
       << "\t[--table <name>] \t Path and basename of LUT files (default = " 
       << getenv("VICTOR_ROOT") +  path << ") \n"
       << "\t[--withOxygen] \t\t Include Oxygen atoms in RMSD calculation\n"
       << "\t[--verbose] \t\t Verbose mode\n"
       << "\t[--silent] \t\t Silent mode\n";
}

/**
*@Description prints the special options of the apps
*/
inline void showSpecialOptions(){
  LoopModel lm;
   std::cout << "\t[--sol2 <n-solutions>] \t Number of secondary solutions"
       << " (default= 1)\n"
       << "\t[--endrmsWeigth] \t Choose weigthing of endrms"
       << " (default= " << lm.getENDRMS_WEIGHT() << ")\n"
       << "\t[--energyWeigth] \t Choose weigthing of energy"
       << " (default= " << LoopModel::ENERGY_WEIGTH << ")\n"
       << "\t[--secprefWeigth] \t Choose weigthing of secondary preference"
       << " (default= " << LoopModel::SECPREF_WEIGTH << ")\n"
       << "\t[--secprefTol] \t\t Choose tolerance for secondary preference"
       << " (default= " << LoopModel::SECPREF_TOL << ")\n"
       << "\t[--packingWeigth] \t Choose weigthing of packing density"
       << " (default= " << LoopModel::PACKING_WEIGTH << ")\n"
       << "\t[--weigthEP] \t\t Choose weigthing of EP compared to ED & EN"
       << " (default= " << LoopTableEntry::LAMBDA_EP << ")\n"
       << "\t[--maxSearch] \t\t Choose max fraction of tables to search "
       << "for best result (default= " << LoopTable::MAX_FACTOR << ")\n"
       << "\t[--vdwLimit] \t\t Choose threshold for VDW filter"
       << " (default= " << LoopModel::VDW_LIMIT << ")\n"
       << "\t[--energyLimit] \t Choose threshold for energy filter"
       << " (default= " << LoopModel::ENERGY_LIMIT << ")\n"
       << "\t[--simLimit] \t\t Choose threshold for clustering similarity"
       << " (default= " << LoopModel::SIM_LIMIT << ")\n"
       << "\t[--interpol] \t\t Use interpolated RAPDF energy\n"
       << "\t[--scatter <filename>] \t Write a scatter plot file.\n";
}
/**
*@Description prints the special options for Lobo
*/
inline void showSpecialOptionsLobo(){
   std::cout << "\t[--cluster] \t\t Cluster results.\n"
       << "\t[--norank] \t\t Do not rank or filter results.\n"
       << "\t[--norankRms] \t\t Do not rank or filter results according "
       << "to end AA RMS.\n"
       << "\t[--noFullModel] \t Do not write full model, write only"
       << " loop instead.\n"
       << "\t[--noWrite]\t\t Do not write output file.\n"
       << "\t[--refine] \t\t Refines the solutions (ie. reduced endRms)\n"
       << "\t[--optimize] \t\t Local optimization of the top solution\n"
       << "\t[--optall] \t\t Local optimization of all solutions pre-ranking\n"
       << "\t[--opt1] \t\t Number of outer iterations for optimization\n"
       << "\t[--opt2] \t\t Number of inner iterations for optimization\n"
       << "\t[--optnum] \t\t Number of solutions to optimize\n"
       << "\t[--test] \t\t Test mode, no optimization\n";
}
/**
*@Description shows the help information 
*/
static inline void sLoboHelp(){
   std::cout << "LOBO - LOop Building & Optimization\n"
       << "This program generates a model for a single loop or indel.\n"
       << " Options: \n"
       << "\t-i <filename> \t\t Input file\n"
       << "\t-c <chain> \t\t pdb chain\n"
       << "\t-s <n-start> \t\t Aminoacid where loop starts\n" 
       << "\t-e <n-end> \t\t Aminoacid where loop ends\n"
       << "\n"
       << "\t[-o <filename>] \t Output file (default = test.pdb)\n"
       << "\t[--seq <filename>] \t Sequence file (in one letter code)\n";
}
/**
*@Description   shows the help for the help option
*/
inline void showLoboHelp(){
  sLoboHelp();
  showCommonOptions();
   std::cout << "\n\t--help   \t\t Get detailed help.\n"
       << endl;
}
/**
*@Description shows all the possible options for the apps
*/
inline void showMoreLoboHelp(){
  sLoboHelp();
  showCommonOptions();
  showSpecialOptions();
  showSpecialOptionsLobo();
   std::cout << endl;
}
/**
*@Description shows  the options for Lobo Auto
*/
static inline void sLoboAutoHelp(){
   std::cout << "LOBO - LOop Building & Optimization\n"
       << "This program generates an automatic model for a single indel.\n"
       << " Options: \n"
       << "\t-i <filename> \t\t Input pdb file\n"
       << "\t-c <filename> \t\t pdb chain\n"
       << "\t-s <n-start> \t\t N-terminal aminoacid where indel starts\n" 
       << "\t-l <length> \t\t Length of indel\n"
       << "\t--del    \t\t Modelling of deletion requested (default)\n"
       << "\t--ins    \t\t Modelling of insertion requested\n"
       << "\t--seq <filename> \t Sequence file (in one letter code)\n"
       << "\n"
       << "\t[-o <filename>] \t Output file (default = test.pdb)\n";
}
/**
*@Description shows  help for for Lobo Auto
*/
inline void showLoboAutoHelp(){
  sLoboAutoHelp();
  showCommonOptions();
   std::cout << "\n\t--help   \t\t Get detailed help.\n"
       << endl;
}
/**
*@Description shows all  the possible options for Lobo Auto
*/
inline void showMoreLoboAutoHelp(){
  sLoboAutoHelp();
  showCommonOptions();
  showSpecialOptions();
  showSpecialOptionsLobo();
   std::cout << endl;
}
/** 
*@Description shows  options for Lobo full 
*/
static inline void sLoboFullHelp(){
   std::cout << "LOBO FULL - LOop Building & Optimization for full proteins\n"
       << "This program generates all models for every k-mer fragment in "
       << "a single protein over a sliding window.\n"
       << " Options: \n"
       << "\t-i <filename> \t\t Input file\n"
       << "\t[-c <chain>] \t\t pdb chain\n"
       << "\t[-o <filename>] \t Output file basename\n"
       << "\t[-s <offset>] \t\t Size of *real* sliding window (default= 5)\n"
       << "\n"
       << "\t[--norankRms] \t\t Do not rank or filter results according "
       << "to end AA RMS.\n"
       << "\t[--fullRank] \t\t Write full ranking.\n";
}
/**
*@Description shows  the common option and the full help for Lobo Full
*/
inline void showLoboFullHelp(){
  sLoboFullHelp();
  showCommonOptions();
   std::cout << "\n\t--help   \t\t Get detailed help.\n"
       << endl;
}
/**
*@Description shows  the common, special option and the full help for Lobo Full
*/
inline void showMoreLoboFullHelp(){
  sLoboFullHelp();
  showCommonOptions();
  showSpecialOptions();
   std::cout << endl;
}
/**
*@Description gets the sequence from a file
*@param  name of the file that contains the sequence(string), flag to make it verbose
*@return  the sequence(string)
*/
inline string sGetSequence(string sequenceFile, bool silent){
  if (!silent)
     std::cout << "seqFile = " << sequenceFile << "\n";

  ifstream seqFile(sequenceFile.c_str());
  if (!seqFile)
    ERROR("Sequence file not found.", exception);

  string tmp;
  seqFile >> tmp;
  return tmp;
}
/**
*@Description manages the arguments used to call the app
*@param  number of arguments(int),pointer to the list of arguments(char* [])
*@return  changes are made internally(void)
*/
inline void treatSpecialOptions(int nArgs, char* argv[]){
  getOption(LoopModel::ENERGY_WEIGTH, "-energyWeigth", nArgs, argv,true);
  getOption(LoopModel::SECPREF_WEIGTH, "-secprefWeigth", nArgs, argv,true);
  getOption(LoopModel::SECPREF_TOL, "-secprefTol", nArgs, argv,true);
  getOption(LoopModel::PACKING_WEIGTH, "-packingWeigth", nArgs, argv,true);
  getOption(LoopTableEntry::LAMBDA_EP, "-weigthEP", nArgs, argv,true);
  getOption(LoopTableEntry::LAMBDA_ED, "-weigthED", nArgs, argv,true);
  getOption(LoopTableEntry::LAMBDA_EN, "-weigthEN", nArgs, argv,true);
  getOption(LoopTable::MAX_FACTOR, "-maxSearch", nArgs, argv,true);
  getOption(LoopModel::VDW_LIMIT, "-vdwLimit", nArgs, argv,true);
  getOption(LoopModel::ENERGY_LIMIT, "-energyLimit", nArgs, argv,true);
  getOption(LoopModel::SIM_LIMIT, "-simLimit", nArgs, argv,true);
}
/**
*@Description prints the sequence and the b factors of the amino acids 
*@param  indexes to define the start and end of the sequence portion(unsigned int, unsigned int), reference to the spacer(Spacer& )
*@return   changes are made internally(void)
*/
inline void printSeqAndBfactors(unsigned int index1, unsigned int index2,Spacer& sp){
  fillLine();
   std::cout << "Sequence & B-factors:\n";
 
  for (unsigned int i = index1+1; i < index2+1; i++)    {
       std::cout << "   " << sp.getAmino(i).getType();
      for (unsigned int j = 0; j < sp.getAmino(i).sizeBackbone(); j++)
	 std::cout << "   " << setw(5) << setprecision(3) 
	     << sp.getAmino(i)[j].getBFac();
       std::cout << "\n";
    }
}
/**
*@Description Manages the problem with proline because of its different structure
*@param  indexes to define the start and end of the sequence portion(unsigned int, unsigned int), reference to the spacer(Spacer& )
*@return   changes are made internally(void)
*/
inline void treatProlineBug(unsigned int index1, unsigned int index2,Spacer& sp){
  // take care of PROline sidechain bug: 
  if (sp.getAmino(index1+1).getCode() == PRO)    {
      if (sp.getAmino(index1+1).getSideChain().isMember(CD))
	sp.getAmino(index1+1).getSideChain().removeAtom(
		     sp.getAmino(index1+1).getSideChain()[CD]);
    }
  if (sp.getAmino(index2+1).getCode() == PRO)    {
      if (sp.getAmino(index2+1).getSideChain().isMember(CD))
	sp.getAmino(index2+1).getSideChain().removeAtom(
		     sp.getAmino(index2+1).getSideChain()[CD]);
    }
}
/**
*@Description Prints the secondary structure and the loop
*@param  indexes to define the start and end of the sequence portion(unsigned int, unsigned int), reference to the spacer(Spacer& )
*@return   changes are made internally(void)
*/
inline void printSecStructureAndLoop(unsigned int index1, unsigned int index2,Spacer& sp){
  for (unsigned int i = 0; i < sp.sizeAmino(); i++)    {
      if (sp.getAmino(i).getState() == HELIX)
	 std::cout << "H";
      else if (sp.getAmino(i).getState() == STRAND)
	 std::cout << "E";
      else
	 std::cout << "-";
      if ((i+1) % 60 == 0)
	 std::cout << "\n";
    }
   std::cout << "\n";
  for (unsigned int i = 0; i < sp.sizeAmino(); i++)    {
      if ((i >= index1) && (i <= index2))
	 std::cout << "#";
      else
	 std::cout << ".";
      if ((i+1) % 60 == 0)
	 std::cout << "\n";
    }
   std::cout << "\n";
  fillLine();
}
/**
*@Description Prints all the results
*@param  indexes to define the start and end of the sequence portion(unsigned int, unsigned int), 
 * reference of the loop model, reference to the spacer(Spacer& ),reference to vector containing the posible solutions( vector<Spacer>&),
 * flag to consider Oxigen(bool), value for the maximum number of values printed(unsigned int)
*@return   changes are made internally(void)
*/
inline void printResults(unsigned int index1, unsigned int index2,
LoopModel& lm, Spacer& sp, vector<Spacer>& vsp, bool withOxygen, unsigned int maxWrite = 9999){
  fillLine();
   std::cout << " Results:  \t\t\t\t\t\t    1.35     121     180\n";

  maxWrite = (maxWrite < vsp.size() ? maxWrite : vsp.size());
  for (unsigned int i = 0; i < maxWrite; i++)    {
       std::cout << setw(3) << i << "   ";
      lm.calculateRms(sp, index1, index2, vsp[i], true, withOxygen);
    }
   std::cout << endl;
}
/**
*@Description Saves the Loops 
*@param   indexes to define the start and end of the sequence portion(unsigned int, unsigned int), 
 * reference of the loop model, reference to the spacer(Spacer& ),reference to vector containing the posible solutions( vector<Spacer>&),
 * name of the output file, value for the maximum number of values writed in the file(unsigned int), flag to consider nofullModel
*@return   changes are made internally(void)
*/
inline void saveLoopEnsemble(unsigned int index1, unsigned int index2,
LoopModel& lm, Spacer& sp, vector<Spacer>& vsp, string outputFile, 
unsigned int maxWrite, bool noFullModel){
  unsigned int maxW = vsp.size() < maxWrite ? vsp.size() : maxWrite;
  for (unsigned int i = 0; i < maxW; i++)    {
      string tmpStr = outputFile + ".";
      
      if (i < 10)
	tmpStr += "00";
      else if (i < 100)
	tmpStr += "0";
      
      tmpStr += itos(i);
      
      ofstream outFile(tmpStr.c_str());
      if (!outFile)
	ERROR("Could not create file.", exception);
      PdbSaver ps(outFile);
      ps.setWriteAtomOnly();
      
      if (!noFullModel)	{
	  lm.setStructure(sp, vsp[i], index1, index2);
	  sp.save(ps);
	}
      else	{
	  vsp[i].save(ps);
	}
    }
  
}

#endif //_LOBOTOOLS_H_
