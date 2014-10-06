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
* @Class               Chi1Chi2
* @Project       Victor
*/

#ifndef _CHI1CHI2_H_
#define _CHI1CHI2_H_


// Includes:
#include <TorsionPotential.h>
#include <vector>
#include <Spacer.h>
#include <AminoAcidCode.h>
#include <IoTools.h>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
  /** @brief class manages the angle qualities and the energy 
 * 
* @Description This class implements a simple torsion potential based on the  statistical preference of aminoacid types for Chi1 and Chi2 angles.
* @This 
 * */
class Chi1Chi2 : public TorsionPotential
{
public: 

  // CONSTRUCTORS/DESTRUCTOR:
  Chi1Chi2 ( string knownledge = "data/tor.par"); // OTHER KNOWNLEDGE
  virtual ~Chi1Chi2() { pResetData(); } 

  // PREDICATES:
  virtual long double calculateEnergy(Spacer& sp);
  virtual long double calculateEnergy(AminoAcid& aa);
  virtual long double calculateEnergy(AminoAcid& aa, Spacer& sp)
  { return calculateEnergy(aa); }
  virtual long double calculateEnergy(Spacer& sp, unsigned int index1, unsigned int index2);
  virtual long double calculateEnergy(AminoAcid& diheds, AminoAcidCode code) 
    {ERROR("ERROR. NOT IMPLEMENTED FOR THIS CLASS.", exception)}
  virtual double pReturnMaxPropensities(int amino);
  virtual double pReturnMaxPropensitiesPreAngle(int amino, int prephi, int prepsi)
  { ERROR("ERROR. NOT IMPLEMENTED FOR THIS CLASS.", exception)}
  virtual int sGetPropBin2(double p)
  { ERROR("ERROR. NOT IMPLEMENTED FOR THIS CLASS.",exception)}

  // MODIFIERS:
  
  // OPERATORS:
  
protected:
  
  // HELPERS:
  virtual void pConstructData();
  virtual void pResetData();
  virtual double pGetMaxPropensities(int amino);
  //not implemented in this class. No pre-angle considered.
  virtual double pGetMaxPropensities(int amino, int prephi, int prepsi)
  { ERROR("ERROR. NOT IMPLEMENTED FOR THIS CLASS.", exception)}
  virtual void sAddProp(int code, int x, int y);
  virtual int sGetPropChiBin(double p);
  virtual void pConstructMaxPropensities();
  
private:
  
  // ATTRIBUTES:
  string TOR_PARAM_FILE;   // File with prop torsion angles
  int amino_count[AminoAcid_CODE_SIZE]; 
  // total number of entries for all amino acids
  vector<vector<vector<int>* >* > propensities;// the propensities table.
  vector<vector<int>* > all_propensities;// the sum of propropensities table.
  int CHI_RANGE;
  double total;
  vector<double> amino_max_propensities;//vector with max amino propensities
                                        // according to knowledge.

};

// ---------------------------------------------------------------------------
//                            Chi1Chi2
// -----------------x-------------------x-------------------x-----------------
} // namespace
#endif //_CHI1CHI2_H_


