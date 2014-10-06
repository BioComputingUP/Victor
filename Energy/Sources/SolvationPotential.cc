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
 *@Class               SolvationPotential 
 *@Project     Victor
 *@Description:
 *    This class implements a simple knowledge-based solvation potential.
  *   This causes a runtime error on 64 bit OS, use the path variable to avoid it if you are using a 64 bit SO
 *    string SolvationPotential::SOLV_PARAM_FILE = "data/solv.par";
*/

// Includes:
#include <float.h>
#include <SolvationPotential.h>

using namespace Biopool;

// Global constants, typedefs, etc. (to avoid):
unsigned int SolvationPotential::MAX_BINS = 30;
//NOTE:
//This causes a runtime error on 64 bit OS, use the path variable to avoid it if you are using a 64 bit SO
//string SolvationPotential::SOLV_PARAM_FILE = "data/solv.par";

// CONSTRUCTORS/DESTRUCTOR:


/**
 *@Description    basic constructor, allocate the information from the solv.par file, and set min and max propensities values
 */
SolvationPotential::SolvationPotential(unsigned int _resol) : sum(), binResolution(_resol){
  char *victor=getenv("VICTOR_ROOT");
	if (victor  == NULL)
		ERROR("Environment variable VICTOR_ROOT was not found.\n Use the command:\n export VICTOR_ROOT=......", exception);
  string path = "data/solv.par";//comment this if you are using a no 64bits OC
  string inputFile = getenv("VICTOR_ROOT");
  if (inputFile.length() < 3)
    ERROR("Environment variable VICTOR_ROOT was not found.", exception);
     //inputFile += SOLV_PARAM_FILE;//uncomment this if you are using a no 64bits OC
  inputFile += path;//comment this if you are using a no 64bits OC
  ifstream input(inputFile.c_str());
  if (!input)
    ERROR("Could not read data file.", exception);
  for (unsigned int i = 0; i < AminoAcid_CODE_SIZE; i++){
      vector<int >* tmp = new vector<int >;
      (*tmp).reserve(MAX_BINS/binResolution);
      for (unsigned int i = 0; i < (MAX_BINS/binResolution)+1; i++)
	(*tmp).push_back(0);
      sum.push_back(*tmp);
    }
  while (input){
      unsigned int num;
      string type;
      input >> num >> type;
      if (!input)
	break;
      AminoAcidCode code = aminoAcidThreeLetterTranslator(type);
      sum[code][(MAX_BINS/binResolution)] = num;
      for (unsigned int i = 0; i < (MAX_BINS/binResolution); i++){
	  sum[code][i] = 0;
 	  for (unsigned int j = 0; j < binResolution; j++){
 	      int tmp;
 	      input >> tmp; 
 	      sum[code][i] += tmp;
	    }
 	}
    }

  input.close();
  pConstructMaxPropensities();
  pConstructMinPropensities();
}

// PREDICATES:
/**
 *@Description calculates the total solvation potential for the residue in the spacer    
 *@param   reference of a Spacer(Spacer&)
 *@return    value of the total solvation potential(long double)
 */
long double SolvationPotential::calculateSolvation(Spacer& sp){
  long double solv = 0.0;
  for (unsigned int i = 0; i < sp.sizeAmino(); i++)
    solv += calculateSolvation(sp.getAmino(i), sp);
  return solv;
}

/**
 *@Description calculates the total energy for a portion of a residue in the spacer    
 *@param   reference of a Spacer(Spacer&), index for the beginning and ending position of the residue portion(unsigned int)
 *@return    value of the total solvation potential(long double)
 */
long double SolvationPotential::calculateEnergy(Spacer& sp, unsigned int index1, unsigned int index2){
  long double solv = 0.0;
  for (unsigned int i = index1; i < index2; i++)
    solv += calculateSolvation(sp.getAmino(i), sp);
  return solv;
}

/**
 *@Description calculates the total solvation potential for a portion of residues in the spacer, 
 * and an specific residue. the relation between them should be minor to SOLVATION_CUTOFF_DISTANCE and not identical    
 *@param  reference of a residue (AminoAcid&), reference of a Spacer(Spacer&), index for the beginning and ending position of the residues portion(unsigned int)
 *@return    value of the total solvation potential(long double)
 */
long double SolvationPotential::calculateSolvation(AminoAcid& aa, Spacer& sp, unsigned int start, unsigned int end){
  if ((aa.getCode() == GLY) || (!aa.getSideChain().isMember(CB)))
    return 0.0;
  start = (start >= 0) ? start : 0;
  end = (end < sp.sizeAmino()) ? end : sp.sizeAmino()-1;
  unsigned int count = 0;
  for (unsigned int j = start; j <= end; j++)
      {
	if (sp.getAmino(j).getSideChain().isMember(CB))
	  if ((aa.getSideChain()[CB].distance(
	     sp.getAmino(j).getSideChain()[CB]) <= SOLVATION_CUTOFF_DISTANCE) 
	      && (aa.getSideChain()[CB].distance(
	      sp.getAmino(j).getSideChain()[CB]) > 0.0)) // check not identical
	    count++;
      }
  // adapt to binResolution:
  count /= binResolution;
  if (count >= (MAX_BINS/binResolution))
    count = (MAX_BINS/binResolution)-1;
  long double a = static_cast<long double>(sum[aa.getCode()][count]+1) / 
    static_cast<long double>(sum[aa.getCode()][MAX_BINS/binResolution]);
  long double b = static_cast<long double>(sum[sum.size()-1][count]+1) / 
    static_cast<long double>(sum[sum.size()-1][MAX_BINS/binResolution]);
  return -propCoeff() * log ( a / b );
}
/**
 *@Description    solvation potential of the residue, considering a given amino acid type.
 *@param   reference of a residue (AminoAcid&), amino acid type (aminoacidCode) ,reference of a Spacer(Spacer&)
 *@return    value of the solvation potential(long double)
 */
long double SolvationPotential::calculateSolvation(AminoAcid& resid, AminoAcidCode type, 
Spacer& sp) {

  /* return a default value of zero if no CB is defined for resid or type.
     This is the only part of code in SolvationPotential that should be changed
     if a finer treatment of the no-CB case will be developed */ 
  if((resid.getCode() == GLY) || (!resid.getSideChain().isMember(CB)))
    return 0;
  if(type == GLY)
    return 0;

  /* count the number of CB atoms that are within a cutoff distance
     of the CB atom of resid. */
  unsigned int count = 0;
  for (unsigned int j = 0; j < sp.sizeAmino(); j++)
      {
	if (sp.getAmino(j).getSideChain().isMember(CB))
	  if ((resid.getSideChain()[CB].distance(
	     sp.getAmino(j).getSideChain()[CB]) <= SOLVATION_CUTOFF_DISTANCE) 
	      && (resid.getSideChain()[CB].distance(
	      sp.getAmino(j).getSideChain()[CB]) > 0.0)) // check not identical
	    count++;
      }
  /* use such a count as a bin index in the solvation-potential table.
     Legend: a:amino acid type; c:count; o:occurrences; t:total; f:fraction */
  // adapt to binResolution:
  count /= binResolution;
  if (count >= (MAX_BINS/binResolution))
    count = (MAX_BINS/binResolution)-1;
  long double aco = static_cast<long double>(sum[type][count]+1);
  long double ato = static_cast<long double>(sum[type][MAX_BINS/binResolution]); 
  long double cto = static_cast<long double>(sum[sum.size()-1][count]+1);
  long double  to = static_cast<long double>(sum[sum.size()-1][MAX_BINS/binResolution]);
  long double caf = aco/ato;
  long double ctf = cto/to;
  return -propCoeff() * log(caf/ctf); 
}


 
/**
 *@Description    obtains the propensity value
 *@param  amino acid type (AminoAcidCode), count (unsigned int)
 *@return  solvation propensity for the (type, count) slot in the 'sum' table(long double)    
 */
long double SolvationPotential::pGetPropensity(const AminoAcidCode type, unsigned int count) const {
  // legend: a:amino acid type; c:count; o:occurrences; t:total; f:fraction
  long double aco = static_cast<long double>(sum[type][count]+1);
  long double ato = static_cast<long double>(sum[type][MAX_BINS/binResolution]); 
  long double cto = static_cast<long double>(sum[sum.size()-1][count]+1);
  long double to = static_cast<long double>(sum[sum.size()-1][MAX_BINS/binResolution]);
  long double caf = aco/ato;
  long double ctf = cto/to;

  return (caf/ctf); // this is the propensity
}

/**
 *@Description    gets maximum solvation propensity for that type.
 *@param  amino acid type (AminoAcidCode) 
 *@return the maximum propensity value (long double)    
 */
long double SolvationPotential::pGetMaxPropensity(const AminoAcidCode type) const{
  long double max = -100000;
  for(unsigned int count = 0; count < (MAX_BINS/binResolution); count++){
    long double prop = pGetPropensity(type, count);
    if (prop > max)
      max = prop;
  }
  return max;
}
/**
 *@Description    calculates  minimum solvation propensity for that type.
 *@param    amino acid type (AminoAcidCode) 
 *@return the minimum propensity value (long double)    
 */
long double SolvationPotential::pGetMinPropensity(const AminoAcidCode type) const{
  long double min = LDBL_MAX; // database independent
  for(unsigned int count = 0; count < (MAX_BINS/binResolution); count++){
    long double prop = pGetPropensity(type, count);
    if (prop < min)
      min = prop;
  }

  return min;
} 


 
/**
 *@Description    stores the maximum propensity of each amino acid type. 
 *    the max. propensity of the i-th amino acid type in the AminoAcidCode
 *    enumeration is stored in the i-th position of amino_max_propensities
 *    (i=0,...,19)
 *@param   none 
 *@return    changes the object internally (void)  
 */
void Biopool::SolvationPotential::pConstructMaxPropensities() {
  amino_max_propensities.clear();
  amino_max_propensities.reserve(AminoAcid_CODE_SIZE-1);
  
  for(AminoAcidCode i = ALA; i <= TYR; i++) {
    long double mpi = pGetMaxPropensity(i);
    amino_max_propensities.push_back(mpi);
  }
}



/**
 *@Description    stores the minimum propensity of each amino acid type. 
 * the min. propensity of the i-th amino acid type in the AminoAcidCode
*    enumeration is stored in the i-th position of amino_max_propensities
*    (i=0,...,19).
*
*    this function relies on the fact that the codes of the 20 natural amino
*    acid types are 20 consecutive integers.
 *@param   none 
 *@return     changes the object internally (void)  
 */
void Biopool::SolvationPotential::pConstructMinPropensities() {
  amino_min_propensities.clear();
  amino_min_propensities.reserve(AminoAcid_CODE_SIZE-1);
  
  for(AminoAcidCode i = ALA; i <= TYR; i++) {
    long double mpi = pGetMinPropensity(i);
    amino_min_propensities.push_back(mpi);
  }
}

