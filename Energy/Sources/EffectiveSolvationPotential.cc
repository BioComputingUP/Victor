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
 *  
*@Class:              EffectiveSolvationPotential
 *@Project Name:       Victor
*@Description:
*    This class implements a knowledge-based solvation with polar/hydrophobic 
*    information potential.
*/
// Includes:
#include <EffectiveSolvationPotential.h>
#include <AminoAcidCode.h>

using namespace Biopool;

// Global constants, typedefs, etc. (to avoid):

// CONSTRUCTORS/DESTRUCTOR:
/**
 *@Description    basic constructor, sets the efective solvation potential for the aa
 */
EffectiveSolvationPotential::EffectiveSolvationPotential() { 

    solvCoeff.resize(AminoAcid_CODE_SIZE);

  // "hydrophobic":
  solvCoeff[ALA] = 71.55;
  solvCoeff[CYS] = 89.80;
  solvCoeff[GLY] = 69.09;
  solvCoeff[HIS] = 70.12;
  solvCoeff[ILE] = 86.90;
  solvCoeff[LEU] = 85.48;
  solvCoeff[MET] = 81.92;
  solvCoeff[PHE] = 86.82;
  solvCoeff[SER] = 61.94;
  solvCoeff[THR] = 65.91;
  solvCoeff[TRP] = 84.95;
  solvCoeff[TYR] = 80.38;
  solvCoeff[VAL] = 85.36;
  
  // "polar":
  solvCoeff[ARG] = 59.57;
  solvCoeff[ASN] = 58.77;
  solvCoeff[ASP] = 55.87;
  solvCoeff[GLN] = 57.46;
  solvCoeff[GLU] = 50.44;
  solvCoeff[LYS] = 44.30;
  solvCoeff[PRO] = 57.34;

  solvCoeff[XXX] = 60.0;

}


// PREDICATES:
/**
 *@Description Calculates Energy from the amino acids in the spacer
*@param spacer reference (spacer&)
 *@return energy value (long double)
*/
long double EffectiveSolvationPotential::calculateEnergy(Spacer& sp){
  long double solv = 0.0;

  for (unsigned int i = 0; i < sp.sizeAmino(); i++)
    solv += calculateEnergy(sp.getAmino(i), sp);
  
  return solv;
}

/**
 *@Description Calculates Energy
*@param spacer reference(spacer&),positions of the start and eend of the amino acids portion of the spacer(unsigned int,unsigned int)
 *@return energy value(long double)
*/
long double EffectiveSolvationPotential::calculateEnergy(Spacer& sp, unsigned int index1, 
unsigned int index2){
  long double solv = 0.0;

  for (unsigned int i = index1; i < index2; i++)
    solv += calculateEnergy(sp.getAmino(i), sp);
  
  return solv;
}

/**
 *@Description Calculates Fraction of the amino acids Buried in the spacer
*@param   central amino acid position to considers,(unsigned int ),spacer reference(spacer&)
 *@return the fraction of the amino acids buried (double)
*/
double EffectiveSolvationPotential::pCalcFracBuried(unsigned int index, Spacer& sp){
  double barea = 0.0;
  double fneb = 0.0;
  double scheck = 0.0;
  unsigned int count = 0;

  for (unsigned int j = 0; j < sp.sizeAmino(); j++){
      //only does this for the amino acids that are not in the group of 3, considering the index as the center of this 3 aas
      if ((j > index+1) || (j < index-1)){
	fneb += 4.0 / sp.getAmino(index)[CA].distance(sp.getAmino(j)[CA]);

	for (unsigned int k = 0; k < sp.sizeAmino(); k++)
	  if ((k > index+1) || (k < index-1)){
	      scheck = sp.getAmino(index)[CA].distance(sp.getAmino(j)[CA])
		/ (sp.getAmino(index)[CA].distance(sp.getAmino(k)[CA])
		+ sp.getAmino(j)[CA].distance(sp.getAmino(k)[CA]));

	      if (scheck > 0.7){
		  count++;

 
		    fneb *= (1-scheck);
		}
	    }
	barea += fneb;
      }
  }
  return (barea > 10.0 ? 10.0 : barea);
}

/**
 *@Description Calculates the energy for the amino acids in the Spacer
*@param  amino acid reference(not used AminoAcid&),spacer reference(Spacer&),
 *@return energy value (long double)
*/
long double EffectiveSolvationPotential::calculateEnergy(AminoAcid& aa, Spacer& sp){
  long double energy = 0.0;

  for (unsigned int i = 1; i < sp.sizeAmino()-1; i++){
      double fracBuried = pCalcFracBuried(i, sp);
      
      double sign = (isPolar(sp.getAmino(i)) ? -1.0 : 1.0);

      energy += sign * fracBuried 
	* (solvCoeff[sp.getAmino(i).getCode()]-60.0)/29.8;
    }
  return energy;
}

