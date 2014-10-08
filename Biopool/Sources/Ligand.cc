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

//Includes:
#include <Ligand.h>
#include <Debug.h>
#include <IntCoordConverter.h>
#include <limits.h>
#include <float.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;

// CONSTRUCTORS/DESTRUCTOR:

Ligand::Ligand() : Group(0, 0) {
}

Ligand::Ligand(const Ligand& orig) {
    this->copy(orig);
}

Ligand::~Ligand() {
    PRINT_NAME;
}

// PREDICATES:

/**
 *@Description 
 *@param 
 */
bool Ligand::isMetalCompound() //not used
{
    return ((size() == 1) && (
            (getType() == "1CU") ||
            (getType() == "2MO") || (getType() == "3CO") ||
            (getType() == "4MO") || (getType() == "6MO") ||

            (getType() == "AG") || (getType() == "AL") ||
            (getType() == "ALF") || (getType() == "ATH") ||
            (getType() == "AU") || (getType() == "AUC") ||

            (getType() == "BA") || (getType() == "BR") ||

            (getType() == "CA") || (getType() == "CD") ||
            (getType() == "CD1") || (getType() == "CE") ||
            (getType() == "CL") || (getType() == "CO") ||
            (getType() == "CO5") || (getType() == "CR") ||
            (getType() == "CS") || (getType() == "CU") ||
            (getType() == "CU1") || (getType() == "CUA") ||
            (getType() == "CUZ") ||

            (getType() == "EU") || (getType() == "EU3") ||

            (getType() == "F") || (getType() == "FE") ||
            (getType() == "FE2") || (getType() == "GA") ||

            (getType() == "HG") || (getType() == "HGC") ||

            (getType() == "IN") || (getType() == "IOD") ||
            (getType() == "IR") || (getType() == "IUM") ||

            (getType() == "K") || (getType() == "LA") ||
            (getType() == "LCO") || (getType() == "LCP") ||
            (getType() == "LI") || (getType() == "LU") ||

            (getType() == "MG") || (getType() == "MN") ||
            (getType() == "MN3") || (getType() == "MN5") ||
            (getType() == "MO1") || (getType() == "MO2") ||
            (getType() == "MO3") || (getType() == "MO4") ||
            (getType() == "MO5") || (getType() == "MO6") ||
            (getType() == "MOO") || (getType() == "MOS") ||
            (getType() == "MW1") || (getType() == "MW2") ||
            (getType() == "NA") || (getType() == "NA5") ||
            (getType() == "NA6") || (getType() == "NAW") ||
            (getType() == "NCO") || (getType() == "ND4") ||
            (getType() == "NI") || (getType() == "NI1") ||
            (getType() == "NI2") || (getType() == "O4M") ||

            (getType() == "OC1") || (getType() == "OC2") ||
            (getType() == "OC3") || (getType() == "OC4") ||
            (getType() == "OC5") || (getType() == "OC6") ||
            (getType() == "OC7") || (getType() == "OCL") ||
            (getType() == "OCM") || (getType() == "OCN") ||
            (getType() == "OCO") || (getType() == "OF1") ||
            (getType() == "OF3") || (getType() == "OS") ||

            (getType() == "PB") || (getType() == "PBM") ||
            (getType() == "PT") || (getType() == "PTN") ||

            (getType() == "RB") || (getType() == "SB") ||
            (getType() == "SE4") || (getType() == "SM") ||
            (getType() == "SO3") || (getType() == "SO4") ||
            (getType() == "SR") || (getType() == "SUL") ||

            (getType() == "TB") || (getType() == "TCN") ||
            (getType() == "TL") || (getType() == "VO4") ||
            (getType() == "WO4") || (getType() == "Y1") ||
            (getType() == "YB") || (getType() == "YT3") ||

            (getType() == "ZN") || (getType() == "ZN1") ||
            (getType() == "ZN2") || (getType() == "ZN3") ||
            (getType() == "ZNO") || (getType() == "1CU")

            ) ? true : false);
}

/**
 *@Description 
 */
bool Ligand::isCommonMetal() {
    return ((size() == 1) && (
            (getType() == "1CU") || (getType() == "3CO") ||

            (getType() == "AG") || (getType() == "AL") ||
            (getType() == "AU") ||

            (getType() == "CA") || (getType() == "CD") ||
            (getType() == "CD1") || (getType() == "CO") ||
            (getType() == "CO5") || (getType() == "CR") ||
            (getType() == "CU") ||
            (getType() == "CU1") || (getType() == "CUA") ||
            (getType() == "CUZ") ||

            (getType() == "F") || (getType() == "FE") ||
            (getType() == "FE2") || (getType() == "HG") ||

            (getType() == "IOD") || (getType() == "K") ||
            (getType() == "LI") ||

            (getType() == "MG") || (getType() == "MN") ||
            (getType() == "MN3") || (getType() == "MN5") ||
            (getType() == "MO1") || (getType() == "MO2") ||
            (getType() == "MO3") || (getType() == "MO4") ||
            (getType() == "MO5") || (getType() == "MO6") ||
            (getType() == "MW1") || (getType() == "MW2") ||
            (getType() == "NA") || (getType() == "NA5") ||
            (getType() == "NA6") || (getType() == "NAW") ||
            (getType() == "ND4") || (getType() == "NI") ||
            (getType() == "NI1") || (getType() == "NI2") ||

            (getType() == "OC1") || (getType() == "OC2") ||
            (getType() == "OC3") || (getType() == "OC4") ||
            (getType() == "OC5") || (getType() == "OC6") ||
            (getType() == "OC7") || (getType() == "OCL") ||
            (getType() == "OCM") || (getType() == "OCN") ||
            (getType() == "OCO") || (getType() == "OF1") ||
            (getType() == "OF3") || (getType() == "PB") ||
            (getType() == "PT") || (getType() == "SB") ||
            (getType() == "SE4") ||

            (getType() == "ZN") || (getType() == "ZN1") ||
            (getType() == "ZN2") || (getType() == "ZN3") ||
            (getType() == "ZNO") || (getType() == "1CU")

            ) ? true : false);
}

/**
 *@Description 
 *@param 
 */
bool Ligand::isWater() {
    return (getType() == "H2O");
}

/**
 *@Description 
 *@param 
 */
bool Ligand::isCofactor() {
    return ( !(Ligand::isWater()) && !(Ligand::isSimpleMetalIon()));
}

/**
 *@Description 
 *@param 
 */
bool Ligand::isSimpleMetalIon() {
    return (
            (getType() == "3CO") || (getType() == "3NI") ||
            (getType() == "4MO") || (getType() == "6MO") ||
            (getType() == "AG") || (getType() == "AL") ||
            (getType() == "AU") || (getType() == "AU3") ||
            (getType() == "BA") || (getType() == "CA") ||
            (getType() == "CD") || (getType() == "CO") ||
            (getType() == "CR") || (getType() == "CS") ||
            (getType() == "CU") || (getType() == "CU1") ||
            (getType() == "F") || (getType() == "FE") ||
            (getType() == "FE2") || (getType() == "GA") ||
            (getType() == "HG") || (getType() == "IN") ||
            (getType() == "IR") || (getType() == "IR3") ||
            (getType() == "LA") || (getType() == "LI") ||

            (getType() == "K") || (getType() == "Y1") ||
            (getType() == "YT3") || (getType() == "W") ||

            (getType() == "MG") || (getType() == "MN") ||
            (getType() == "MN3") || (getType() == "MO") ||
            (getType() == "NA") || (getType() == "ND4") ||
            (getType() == "NI") || (getType() == "OS") ||
            (getType() == "OS4") || (getType() == "PB") ||
            (getType() == "PD") || (getType() == "PT") ||
            (getType() == "PT4") || (getType() == "RB") ||
            (getType() == "RH") || (getType() == "RH3") ||
            (getType() == "RU") || (getType() == "SB") ||
            (getType() == "SR") || (getType() == "TL") ||
            (getType() == "V") || (getType() == "ZN")
            );
}
// MODIFIERS:

/**
 *@Description 
 *@param 
 */
void Ligand::copy(const Ligand& orig) {
    PRINT_NAME;
    Group::copy(orig);
}

// OPERATORS:

/**
 *@Description 
 *@param 
 */
Ligand& Ligand::operator=(const Ligand& orig) {
    PRINT_NAME;
    if (&orig != this)
        copy(orig);
    return *this;
}

// HELPERS:

