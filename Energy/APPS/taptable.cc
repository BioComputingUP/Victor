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
 @Description A function that performes statistics over tap energy of every amino.
*      We print out a table with MIN, MAX, AVG and SD.
* NOTE: This is only for each beta region in the Ramachandran's map.
 */
#include <string>
#include <GetArg.h>
#include <AminoAcidCode.h>
#include <StatTools.h>
#include <PhiPsi.h>

using namespace std;
using namespace Biopool;

 
void printAlpha( PhiPsi &pp, long phiAlphaBegin, long phiAlphaEnd, long psiAlphaBegin, 	long psiAlphaEnd ){		
	for ( int i=0; i<AminoAcid_CODE_SIZE-1; ++i ) 	{
		cout <<aminoAcidThreeLetterTranslator( static_cast<AminoAcidCode>(i) ) <<"\t";
		vector<double> data;
		for ( long phi=phiAlphaBegin; phi<=phiAlphaEnd; phi += 1 )
			for ( long psi=psiAlphaBegin; psi<=psiAlphaEnd; psi += 1 ) {
				double tmp=pp.getEnergyFromPhiPsi( static_cast<AminoAcidCode>(i), phi, psi );
				data.push_back( tmp );
			}
		cout <<"\t\t" <<minimum(data) <<"\t\t" 
			<<maximum(data) <<"\t\t" 
			<<average(data) <<"\t\t"
			<<standardDeviation(data) <<endl;
	}
}

 
void printBeta( PhiPsi &pp, 	long phiBetaBegin,long phiBetaEnd, long psiBetaBegin, long psiBetaEnd ){
	for ( int i=0; i<AminoAcid_CODE_SIZE-1; ++i ) 	{
		cout <<aminoAcidThreeLetterTranslator( static_cast<AminoAcidCode>(i) ) <<"\t";
		vector<double> data;
		for ( long phi=phiBetaBegin; phi<=phiBetaEnd; phi += 1 )
			for ( long psi=psiBetaBegin; psi<=psiBetaEnd; psi += 1 ) {
				double tmp=pp.getEnergyFromPhiPsi( static_cast<AminoAcidCode>(i), phi, psi );
				data.push_back( tmp );
			}
		cout <<"\t\t" <<minimum(data) 
			<<"\t\t" <<maximum(data) <<"\t\t" 
			<<average(data) <<"\t\t"
			<<standardDeviation(data) <<endl;
	}
}


 
bool isInside( 	long value, long begin, long end ){
	if ( value>=begin && value<=end )
		return true;
	return false;
}

 

void printCoil( PhiPsi &pp, long phiAlphaBegin, long phiAlphaEnd, long psiAlphaBegin, 
	long psiAlphaEnd, long phiBetaBegin, long phiBetaEnd, long psiBetaBegin, long psiBetaEnd ){		
	for ( int i=0; i<AminoAcid_CODE_SIZE-1; ++i )	{
		cout <<aminoAcidThreeLetterTranslator( static_cast<AminoAcidCode>(i) ) <<"\t";
		vector<double> data;
		for ( long phi=-180; phi<=180; phi += 1 )
			for ( long psi=-180; psi<=180; psi += 1 ){
				if ( !isInside(phi,phiAlphaBegin,phiAlphaEnd) && 
					 !isInside(phi,phiBetaBegin,phiBetaEnd) && 
					 !isInside(psi,psiAlphaBegin,psiAlphaEnd) && 
					 !isInside(psi,psiBetaBegin,psiBetaEnd) )	{
					double tmp=pp.getEnergyFromPhiPsi( static_cast<AminoAcidCode>(i), phi, psi );
					data.push_back( tmp );
				}
			}

		cout <<"\t\t" <<minimum(data) <<"\t\t" 
			<<maximum(data) <<"\t\t" 
			<<average(data) <<"\t\t"
			<<standardDeviation(data) <<endl;
	}
}


// ----------------------------------------------------------------------------
// MAIN PROGRAM
// ----------------------------------------------------------------------------
int main(int nArgs, char* argv[])
{
	
	// ----------------
	// 0. Treat options
	// ----------------
	 char *victor=getenv("VICTOR_ROOT");
	if (victor  == NULL)
		ERROR("Environment variable VICTOR_ROOT was not found.\n Use the command:\n export VICTOR_ROOT=......", exception);
	string vct = getenv("VICTOR_ROOT");
	if ( vct.empty() )
		ERROR("Environment variable VICTOR_ROOT is not set.", exception);
	
	long phiAlphaBegin, phiAlphaEnd, psiAlphaBegin, psiAlphaEnd, phiBetaBegin, phiBetaEnd, psiBetaBegin, psiBetaEnd;
	
	getArg( "-phiAlphaBegin", phiAlphaBegin, nArgs, argv, -160); 
	getArg( "-phiAlphaEnd", phiAlphaEnd, nArgs, argv, -60); 
	getArg( "-psiAlphaBegin", psiAlphaBegin, nArgs, argv, -65); 
	getArg( "-psiAlphaEnd", psiAlphaEnd, nArgs, argv, -40); 
	getArg( "-phiBetaBegin", phiBetaBegin, nArgs, argv, -160); 
	getArg( "-phiBetaEnd", phiBetaEnd, nArgs, argv, -60); 
	getArg( "-psiBetaBegin", psiBetaBegin, nArgs, argv, 90); 
	getArg( "-psiBetaEnd", psiBetaEnd, nArgs, argv, 175);
		
	PhiPsi pp;
	cout <<endl<<"ALPHA"<<endl;
	printAlpha( pp, phiAlphaBegin, phiAlphaEnd, psiAlphaBegin, psiAlphaEnd ); 
	cout <<endl<<"BETA"<<endl;
	printBeta( pp, phiBetaBegin, phiBetaEnd, psiBetaBegin, psiBetaEnd );
	cout <<endl<<"COIL"<<endl;	
	printCoil( pp, phiAlphaBegin, phiAlphaEnd, psiAlphaBegin, psiAlphaEnd, phiBetaBegin, phiBetaEnd, psiBetaBegin, psiBetaEnd );
	
	return 0;
}
