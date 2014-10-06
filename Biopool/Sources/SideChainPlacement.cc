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
 *@Class:             SideChainPlacement
 *@Base Class(es):    EnergyCalculator
 *@Description: 
*    is Singleton(unique), which calculate the energy of Atom based chemical 
*    Structures. 
*    There' s only one existing SideChainPlacement Object, but it will be 
*    changed by changing the forcefield via User-Input.   
*    "Singleton"-Design Pattern is used (Gamma et al. Entwurfsmuster,
*    p.140 - 148, Addison-Wesley, 1996)
*/

// Includes:
#include <float.h>
#include <SideChainPlacement.h>
#include <AminoAcidCode.h>
#include <Group.h>
#include <AminoAcid.h>
#include <PdbSaver.h>

using namespace std;
using namespace Biopool;
// Global constants, typedefs, etc. (to avoid):

#define INFINITY DBL_MAX-1
//static Member:
SideChainPlacement* SideChainPlacement::mSideChainPlacer=NULL;
//unsigned int PurgePriorityQueue<QueueElem>::maxSize = 1000000;
template <class T> typename PurgePriorityQueue<T>::size_type PurgePriorityQueue<T>::maxSize=1000000;

const double ENTROPY_THRESHOLD = 100;


bool QueueElem::operator<(const QueueElem& fCompElem) const{
  //"wrong" operator for "right" topping
  return EstimatedEnergyCostMinimum > fCompElem.EstimatedEnergyCostMinimum;
}

QueueElem::QueueElem(const QueueElem& fOrig)
  :RotamerList(fOrig.RotamerList),ActualEnergyCosts(fOrig.ActualEnergyCosts),
   EstimatedEnergyCostMinimum(fOrig.EstimatedEnergyCostMinimum),
   NodeDepth(fOrig.NodeDepth){
}

QueueElem&  QueueElem::operator=(const QueueElem& fOrig){
  PRINT_NAME;
  
  if (&fOrig != this)
  {
    ActualEnergyCosts=fOrig.ActualEnergyCosts;
    NodeDepth=fOrig.NodeDepth;
    EstimatedEnergyCostMinimum=fOrig.EstimatedEnergyCostMinimum;
    RotamerList=fOrig.RotamerList;
  }
  return *this;
}

void SkipSpace(string& fSetIt,string& fGetIt){
  //skips first occurence of spacecharacters in a line and extracts the next 
  //word assigning it to another string:
  fSetIt='\0';
  bool MyIteratorExhausted=false; 

  //ATTENTION, if iterator exhausted, you don' t know where it is:
  string::iterator MyPos=fGetIt.begin();
  if(fGetIt.begin()!=fGetIt.end())    {
      while ((*MyPos!=' ')&&!(MyIteratorExhausted)) {
	if (MyPos != fGetIt.end()) MyPos++;
	else MyIteratorExhausted = true;
      }

      while((*MyPos==' ')&&!(MyIteratorExhausted))
	if (MyPos != fGetIt.end()) MyPos++;
	else MyIteratorExhausted=true;

      fSetIt.assign(MyPos,fGetIt.end());
      MyPos=fSetIt.begin();
      MyIteratorExhausted=false;
      while((*MyPos!=' ')&&!(MyIteratorExhausted)) {
	if(MyPos!=fSetIt.end()) MyPos++;
	else MyIteratorExhausted=true;
      }
    }
  else
    fSetIt.assign(MyPos,fGetIt.end());

  return;
}

int SideChainPlacement::GetAllBetaInTemplate(Spacer& fSpacer,int fSpacerPos){
  switch (fSpacer.getAmino(fSpacerPos).getSideChain().getCode())
    {
    case GLY:
      return 0;

    case VAL:    case ILE:    case THR:
      return 2;
    case SER:    case LEU:    case PHE:
    case TYR:    case TRP:    case CYS:
    case MET:    case ASP:    case GLU:
    case ASN:    case GLN:    case LYS:
    case ARG:    case HIS:    case PRO:
      return 3;
    case ALA:
      return 4;
    default:
      return 0;
    }
}

void
SideChainPlacement::RMSCalculate(DeviationMeasure& dist, 
SideChain& fSideChain, int fSpacerPos){
  //function could be implemented more effective , each case could be a method
  //in class Biopool::SideChain:

  SideChain& aa = mSpacer->getAmino(fSpacerPos).getSideChain();
  AminoAcidCode aaCode = static_cast<AminoAcidCode>(
			   mSpacer->getAmino(fSpacerPos).getCode());

  // ****** temporary! include CB in RMS calculation
  int offset = 0;
    if (fSideChain.isMember(CB))      {
	dist.addCounter(aaCode, fSideChain[CB].distance(aa[CB]));
	offset++;
      }
  // ****** temporary! include CB in RMS calculation

  double MyMinFirst=0.0;

  switch(fSideChain.getCode())
    {
    case CYS:
      dist.addCounter(aaCode, fSideChain[SG].distance(aa[SG]));
      dist.addSum(aaCode, 1+offset);
      return;

    case ASP:
      dist.addCounter(aaCode, fSideChain[CG].distance(aa[CG]));
      MyMinFirst += fabs(fSideChain[OD1].distance(aa[OD1]));
      MyMinFirst += fabs(fSideChain[OD2].distance(aa[OD2]));
      if(MyMinFirst > fabs(fSideChain[OD1].distance(aa[OD2]))
	 + fabs(fSideChain[OD2].distance(aa[OD1])))	{
	  dist.addCounter(aaCode, fSideChain[OD1].distance(aa[OD2]));
	  dist.addCounter(aaCode, fSideChain[OD2].distance(aa[OD1]));
	}
      else	{
	  dist.addCounter(aaCode, fSideChain[OD1].distance(aa[OD1]));
	  dist.addCounter(aaCode, fSideChain[OD2].distance(aa[OD2]));
	}
      dist.addSum(aaCode, 3+offset);
      return;

    case ASN:
      dist.addCounter(aaCode, fSideChain[CG].distance(aa[CG]));
      dist.addCounter(aaCode, fSideChain[OD1].distance(aa[OD1]));
      dist.addCounter(aaCode, fSideChain[ND2].distance(aa[ND2]));
      dist.addSum(aaCode, 3+offset);
      return;

    case GLU:
      dist.addCounter(aaCode, fSideChain[CG].distance(aa[CG]));
      dist.addCounter(aaCode, fSideChain[CD].distance(aa[CD]));
      MyMinFirst += fabs(fSideChain[OE1].distance(aa[OE1]));
      MyMinFirst += fabs(fSideChain[OE2].distance(aa[OE2]));
      if(MyMinFirst > fabs(fSideChain[OE1].distance(aa[OE2]))
	 + fabs(fSideChain[OE2].distance(aa[OE1])))	{
	  dist.addCounter(aaCode, fSideChain[OE1].distance(aa[OE2]));
	  dist.addCounter(aaCode, fSideChain[OE2].distance(aa[OE1]));
	}
      else	{
	  dist.addCounter(aaCode, fSideChain[OE1].distance(aa[OE1]));
	  dist.addCounter(aaCode, fSideChain[OE2].distance(aa[OE2]));
	}
      dist.addSum(aaCode, 4+offset);
      return;

    case GLN:
      dist.addCounter(aaCode, fSideChain[CG].distance(aa[CG]));
      dist.addCounter(aaCode, fSideChain[CD].distance(aa[CD]));
      dist.addCounter(aaCode, fSideChain[OE1].distance(aa[OE1]));
      dist.addCounter(aaCode, fSideChain[NE2].distance(aa[NE2]));
      dist.addSum(aaCode, 4+offset);
      return;

    case PHE:
      dist.addCounter(aaCode, fSideChain[CG].distance(aa[CG]));
      MyMinFirst += fabs(fSideChain[CD1].distance(aa[CD1]));
      MyMinFirst += fabs(fSideChain[CD2].distance(aa[CD2]));
      MyMinFirst += fabs(fSideChain[CE1].distance(aa[CE1]));
      MyMinFirst += fabs(fSideChain[CE2].distance(aa[CE2]));
      if(MyMinFirst > fabs(fSideChain[CD1].distance(aa[CD2]))
	 + fabs(fSideChain[CE2].distance(aa[CE1]))
	 + fabs(fSideChain[CD2].distance(aa[CD1]))
	 + fabs(fSideChain[CE1].distance(aa[CE2])))	{
	  dist.addCounter(aaCode, fSideChain[CD1].distance(aa[CD2]));
	  dist.addCounter(aaCode, fSideChain[CE2].distance(aa[CE1]));
	  dist.addCounter(aaCode, fSideChain[CD2].distance(aa[CD1]));
	  dist.addCounter(aaCode, fSideChain[CE1].distance(aa[CE2]));
	}
      else	{
	  dist.addCounter(aaCode, fSideChain[CD1].distance(aa[CD1]));
	  dist.addCounter(aaCode, fSideChain[CE2].distance(aa[CE2]));
	  dist.addCounter(aaCode, fSideChain[CD2].distance(aa[CD2]));
	  dist.addCounter(aaCode, fSideChain[CE1].distance(aa[CE1]));
	}
      dist.addCounter(aaCode, fSideChain[CZ].distance(aa[CZ]));
      dist.addSum(aaCode, 6+offset);
      return;
      
    case HIS:
      dist.addCounter(aaCode, fSideChain[CG].distance(aa[CG]));
      dist.addCounter(aaCode, fSideChain[ND1].distance(aa[ND1]));
      dist.addCounter(aaCode, fSideChain[CD2].distance(aa[CD2]));
      dist.addCounter(aaCode, fSideChain[CE1].distance(aa[CE1]));
      dist.addCounter(aaCode, fSideChain[NE2].distance(aa[NE2]));
      dist.addSum(aaCode, 5+offset);
      return;

    case ILE:
      dist.addCounter(aaCode, fSideChain[CG1].distance(aa[CG1]));
      dist.addCounter(aaCode, fSideChain[CG2].distance(aa[CG2]));
      dist.addCounter(aaCode, fSideChain[CD1].distance(aa[CD1]));
      dist.addSum(aaCode, 3+offset);
      return;

    case LYS:
      dist.addCounter(aaCode, fSideChain[CG].distance(aa[CG]));
      dist.addCounter(aaCode, fSideChain[CD].distance(aa[CD]));
      dist.addCounter(aaCode, fSideChain[CE].distance(aa[CE]));
      dist.addCounter(aaCode, fSideChain[NZ].distance(aa[NZ]));
      dist.addSum(aaCode, 4+offset);
      return;

    case LEU:
      dist.addCounter(aaCode, fSideChain[CG].distance(aa[CG]));
      MyMinFirst += fabs(fSideChain[CD1].distance(aa[CD1]));
      MyMinFirst += fabs(fSideChain[CD2].distance(aa[CD2]));
      if (MyMinFirst > fabs(fSideChain[CD1].distance(aa[CD2]))
	 + fabs(fSideChain[CD2].distance(aa[CD1])))	{
	  dist.addCounter(aaCode, fSideChain[CD1].distance(aa[CD2]));
	  dist.addCounter(aaCode, fSideChain[CD2].distance(aa[CD1]));
	}
      else	{
	  dist.addCounter(aaCode, fSideChain[CD1].distance(aa[CD1]));
	  dist.addCounter(aaCode, fSideChain[CD2].distance(aa[CD2]));
	}
      dist.addSum(aaCode, 3+offset);
      return;

    case MET:
      dist.addCounter(aaCode, fSideChain[CG].distance(aa[CG]));
      dist.addCounter(aaCode, fSideChain[SD].distance(aa[SD]));
      dist.addCounter(aaCode, fSideChain[CE].distance(aa[CE]));
      dist.addSum(aaCode, 3+offset);
      return;

    case PRO:
      dist.addSum(aaCode, offset);
      return;

    case ARG:
      dist.addCounter(aaCode, fSideChain[CG].distance(aa[CG]));
      dist.addCounter(aaCode, fSideChain[CD].distance(aa[CD]));
      dist.addCounter(aaCode, fSideChain[NE].distance(aa[NE]));
      dist.addCounter(aaCode, fSideChain[CZ].distance(aa[CZ]));
      MyMinFirst += fabs(fSideChain[NH1].distance(aa[NH1]));
      MyMinFirst += fabs(fSideChain[NH2].distance(aa[NH2]));
      if(MyMinFirst > fabs(fSideChain[NH1].distance(aa[NH2]))
	 + fabs(fSideChain[NH2].distance(aa[NH1])))	{
	  dist.addCounter(aaCode, fSideChain[NH1].distance(aa[NH2]));
	  dist.addCounter(aaCode, fSideChain[NH2].distance(aa[NH1]));
	}
      else	{
	  dist.addCounter(aaCode, fSideChain[NH1].distance(aa[NH1]));
	  dist.addCounter(aaCode, fSideChain[NH2].distance(aa[NH2]));
	}
      dist.addSum(aaCode, 6+offset);
      return;

    case SER:
      dist.addCounter(aaCode, fSideChain[OG].distance(aa[OG]));
      dist.addSum(aaCode, 1+offset);
      return;

    case THR:
      dist.addCounter(aaCode, fSideChain[OG1].distance(aa[OG1]));
      dist.addCounter(aaCode, fSideChain[CG2].distance(aa[CG2]));
      dist.addSum(aaCode, 2+offset);
      return;

    case VAL:
      MyMinFirst += fabs(fSideChain[CG1].distance(aa[CG1]));
      MyMinFirst += fabs(fSideChain[CG2].distance(aa[CG2]));
      if(MyMinFirst > fabs(fSideChain[CG1].distance(aa[CG2]))
	 + fabs(fSideChain[CG2].distance(aa[CG1])))	{
	  dist.addCounter(aaCode, fSideChain[CG1].distance(aa[CG2]));
	  dist.addCounter(aaCode, fSideChain[CG2].distance(aa[CG1]));
	}
      else
	{
	  dist.addCounter(aaCode, fSideChain[CG1].distance(aa[CG1]));
	  dist.addCounter(aaCode, fSideChain[CG2].distance(aa[CG2]));
	}
      dist.addSum(aaCode, 2+offset);
      return;

    case TRP:
      dist.addCounter(aaCode, fSideChain[CG].distance(aa[CG]));
      dist.addCounter(aaCode, fSideChain[CD1].distance(aa[CD1]));
      dist.addCounter(aaCode, fSideChain[CD2].distance(aa[CD2]));
      dist.addCounter(aaCode, fSideChain[NE1].distance(aa[NE1]));
      dist.addCounter(aaCode, fSideChain[CE2].distance(aa[CE2]));
      dist.addCounter(aaCode, fSideChain[CE3].distance(aa[CE3]));
      dist.addCounter(aaCode, fSideChain[CZ2].distance(aa[CZ2]));
      dist.addCounter(aaCode, fSideChain[CZ3].distance(aa[CZ3]));
      dist.addCounter(aaCode, fSideChain[CH2].distance(aa[CH2]));
      dist.addSum(aaCode, 9+offset);
      return;

    case TYR:
      dist.addCounter(aaCode, fSideChain[CG].distance(aa[CG]));
      MyMinFirst += fabs(fSideChain[CD1].distance(aa[CD1]));
      MyMinFirst += fabs(fSideChain[CD2].distance(aa[CD2]));
      MyMinFirst += fabs(fSideChain[CE1].distance(aa[CE1]));
      MyMinFirst += fabs(fSideChain[CE2].distance(aa[CE2]));
      if(MyMinFirst > fabs(fSideChain[CD1].distance(aa[CD2]))
	 + fabs(fSideChain[CE2].distance(aa[CE1])) 
	 + fabs(fSideChain[CD2].distance(aa[CD1]))
	 + fabs(fSideChain[CE1].distance(aa[CE2])) )	{
	  dist.addCounter(aaCode, fSideChain[CD1].distance(aa[CD2]));
	  dist.addCounter(aaCode, fSideChain[CE2].distance(aa[CE1]));
	  dist.addCounter(aaCode, fSideChain[CD2].distance(aa[CD1]));
	  dist.addCounter(aaCode, fSideChain[CE1].distance(aa[CE2]));
	}
      else	{
	  dist.addCounter(aaCode, fSideChain[CD1].distance(aa[CD1]));
	  dist.addCounter(aaCode, fSideChain[CE2].distance(aa[CE2]));
	  dist.addCounter(aaCode, fSideChain[CD2].distance(aa[CD2]));
	  dist.addCounter(aaCode, fSideChain[CE1].distance(aa[CE1]));
	}
      dist.addCounter(aaCode, fSideChain[CZ].distance(aa[CZ]));
      dist.addCounter(aaCode, fSideChain[OH].distance(aa[OH]));
      dist.addSum(aaCode, 7+offset);
      return;

    case ALA:
      dist.addSum(aaCode, offset);

    default:
      return;

    }
  return;
}

//the hidden public constructor:
SideChainPlacement* 
SideChainPlacement::SideChainPlacer()
{
  if(mSideChainPlacer==NULL)    {
      mSideChainPlacer=new SideChainPlacement;
    }
  return mSideChainPlacer; 
}

double SideChainPlacement::DeadEndElimination(double fEnergyCut){
  mNodePositionInSpacer.clear();
  mPositionEntropies.clear();
  mEntropyCombinedLeachLemon.clear();

  // Maximale Knotentiefe
  mMaxNodeDepth=mSpacer->size();

  unsigned int MyPosInSpacer=0;
  int MyCopiedPosInSpacer=0;
  do  {
    if(mAllPossibleRotamerSettings[MyPosInSpacer]==-1)      {
	//no DEE needed:
	mMaxNodeDepth--;
      }
    else      {
	// Schleife ueber alle Rotamer settings und checke, ob 
	// MyFirst durch die Plazierung von MySecond (oder umgekehrt)
	// herausgeworfen werden kann.
	mNodePositionInSpacer.push_back(MyCopiedPosInSpacer);
	for(int MyFirst=mAllPossibleRotamerSettings[MyPosInSpacer];
	    MyFirst<GetLastPos(MyCopiedPosInSpacer);
	    MyFirst++)	  {
	    for(int MySecond=MyFirst+1;
		MySecond<GetLastPos(MyCopiedPosInSpacer);
		MySecond++)	      {
		if(!(mTemplateInteractions[MyFirst].DeadEndingRotamer))  {
		    //this statement-block contains one DEE-step:
		    //double MyMinSum=0.0;
		    //  double MyMaxSum=0.0;
		    //check DEE--Inequality!
		    // this is the DEE Energy:
		    double MyDeadEndingFirst
		      =mTemplateInteractions[MyFirst].InteractionEnergy
		      -mTemplateInteractions[MySecond].InteractionEnergy
		      -fEnergyCut;
		    double MyDeadEndingSecond
		      =mTemplateInteractions[MySecond].InteractionEnergy
		      -mTemplateInteractions[MyFirst].InteractionEnergy
		      -fEnergyCut;
		    //This is Goldsteins minimization
		    
		    for(unsigned int j=0;j<MyPosInSpacer;j++)    {
			double MyIrItJsMin=0.0;
			double MyIrItJsMax=0.0;
			if(mAllPossibleRotamerSettings[j]!=-1)	 {
			    MyIrItJsMin= INFINITY;
			    MyIrItJsMax=-INFINITY;
			    for(int S=mAllPossibleRotamerSettings[j];
				S<GetLastPos(j);S++)			  {
				if(mResidueAllInteractions[MyFirst][S]-
				   mResidueAllInteractions[MySecond][S]
				   <MyIrItJsMin)
				  MyIrItJsMin
				    =mResidueAllInteractions[MyFirst][S]-
				    mResidueAllInteractions[MySecond][S];
				if(mResidueAllInteractions[MyFirst][S]-
				   mResidueAllInteractions[MySecond][S]
				   >MyIrItJsMax)
				  MyIrItJsMax
				    =mResidueAllInteractions[MyFirst][S]-
				    mResidueAllInteractions[MySecond][S];
			      }
			  }
			MyDeadEndingFirst+=MyIrItJsMin;
			MyDeadEndingSecond=MyDeadEndingSecond-MyIrItJsMax;
		      }
		    for(unsigned int j=MyPosInSpacer+1;j<mSpacer->size();j++)    {
			double MyIrItJsMin=0.0;
			double MyIrItJsMax=-0.0;
			if(mAllPossibleRotamerSettings[j]!=-1)  {
			    MyIrItJsMin= INFINITY;
			    MyIrItJsMax=-INFINITY;
			    for(int S=mAllPossibleRotamerSettings[j];
				S<GetLastPos(j);S++)     {
				if(mResidueAllInteractions[MyFirst][S]-
				   mResidueAllInteractions[MySecond][S]
				   <MyIrItJsMin)
				  MyIrItJsMin
				    =mResidueAllInteractions[MyFirst][S]-
				    mResidueAllInteractions[MySecond][S];
				if(mResidueAllInteractions[MyFirst][S]-
				   mResidueAllInteractions[MySecond][S]
				   >MyIrItJsMax)
				  MyIrItJsMax
				    =mResidueAllInteractions[MyFirst][S]-
				    mResidueAllInteractions[MySecond][S];
			      }
			  }
			MyDeadEndingFirst+=MyIrItJsMin;
			MyDeadEndingSecond=MyDeadEndingSecond-MyIrItJsMax;
		      }
		    //MySecond:sum over all Maxima E_ij from mResidueMaxima
		    
		    if(MyDeadEndingFirst>0) {
			mTemplateInteractions[MyFirst].DeadEndingRotamer=true;
		      }
		    if(MyDeadEndingSecond>0)    {
			mTemplateInteractions[MySecond].DeadEndingRotamer=true;
		      }
		  }
		//after this parenthesis one DEE--step is executed
	      }
	  }
      }
    MyPosInSpacer++;
    MyCopiedPosInSpacer++;
  }
  while(MyPosInSpacer<mSpacer->size());
  MyPosInSpacer=0;
  MyCopiedPosInSpacer=0;
  
  do  {
    double MyNorm=0.0;
    //calculating norm
    for(int MyFirst=mAllPossibleRotamerSettings[MyPosInSpacer];
	MyFirst<GetLastPos(MyCopiedPosInSpacer);
	MyFirst++)    {
	if(!mTemplateInteractions[MyFirst].DeadEndingRotamer)
	  MyNorm+=mTemplateInteractions[MyFirst].RelFrequency;
	else
	  mTemplateInteractions[MyFirst].RelFrequency=0.0;
      }
    //norming frequencies
    for(int MyFirst=mAllPossibleRotamerSettings[MyPosInSpacer];
	MyFirst<GetLastPos(MyCopiedPosInSpacer);
	MyFirst++)    {
	mTemplateInteractions[MyFirst].RelFrequency=
	  mTemplateInteractions[MyFirst].RelFrequency/MyNorm;
      }
    //calculate Entropies
    double MyEntropy=0.0;
    for(int MyFirst=mAllPossibleRotamerSettings[MyPosInSpacer];
	MyFirst<GetLastPos(MyCopiedPosInSpacer);
	MyFirst++)   {
	if(!mTemplateInteractions[MyFirst].DeadEndingRotamer)
	  {
	    MyEntropy=MyEntropy
	      +(-mTemplateInteractions[MyFirst].RelFrequency
		*log(mTemplateInteractions[MyFirst].RelFrequency));
	  }
      }
    mPositionEntropies.push_back(MyEntropy);
    MyPosInSpacer++;
    MyCopiedPosInSpacer++;
  }
  while(MyPosInSpacer<mSpacer->size());

  //claculate entropies reduced per LeachLemonIndicators
  for(unsigned int i=0;i<mPositionEntropies.size();i++)  {
    if(mLeachLemonIndicators[i]==0)
      if(mPositionEntropies[i]!=0)
	mEntropyCombinedLeachLemon.push_back(INFINITY);
      else
	mEntropyCombinedLeachLemon.push_back(0);
    else
      mEntropyCombinedLeachLemon.push_back
	(mPositionEntropies[i]/mLeachLemonIndicators[i]);
  }

  
  return 0.0;
}
  

bool SideChainPlacement::DEEPreCalculateAllEnergies(int fCalculationStage){
  //this function performs all Energy Calculations necessary for 
  //DEE--A* Algorithm:
  vector<int> MyAllPossibleRotamerSettings;
  mTemplateInteractions.clear();
  mLeachLemonIndicators.clear();
  mPositionInteractionCounters.clear();
  mTemplateInteractionMinima.clear();

  //Kopie des Spacers mit Testrechnung:
  Spacer MyPrivateSpacer = *mSpacer;

  AtomEnergyCalculator* AECPointer = AtomEnergyCalculator::AtomEnergy();
  double MyCompleteEnergy = 0;

  mNumberOfAllRotamers=0;
  int MyLast=0; 
  //skip this lines if not the first PrecalculatioStage
  if (fCalculationStage == 0) {
    for (unsigned int i=0; i<mSpacer->size(); i++) {
      mFinalResult.push_back(0);
    }
    
    //first set for each position number of available rotamers
    for (unsigned int i=0; i<MyPrivateSpacer.size(); i++) {
      int aaCode_i = MyPrivateSpacer.getAmino(i).getCode();
      if ((aaCode_i != CYS) && (aaCode_i != PRO)) {
	MyAllPossibleRotamerSettings.push_back
	  (allRotamersWithFrequencies[aaCode_i].getNumberOfRotamers());
      } 
	else {
	//Searching Cystine SG--SG Bonds
	if(MyPrivateSpacer.getAmino(i).getCode()==CYS)	  {
		  // Erstes Cystein gefunden
		  bool MySSBond=false;
		  for(unsigned int j=0;j<i;j++)
		    {
		      if(MyPrivateSpacer.getAmino(j).getCode()==CYS)			{
			  // Zweites Cystein gefunden
			  if((MyPrivateSpacer.getAmino(i).
			      getSideChain()[SG].distance
			      (MyPrivateSpacer.getAmino(j).getSideChain()[SG])<2.3)
			     &&(MyPrivateSpacer.getAmino(i).
				getSideChain()[SG].distance
				(MyPrivateSpacer.getAmino(j).
				 getSideChain()[SG])>1.8))   {
			      MyAllPossibleRotamerSettings.push_back(0);
			      //must quit loop
			      j=MyPrivateSpacer.size()+1;
			      MySSBond=true;
			    }
			}
		    }
		  if(!MySSBond)    {
		      for(unsigned int j=i+1;j<MyPrivateSpacer.size();j++){
			  if(MyPrivateSpacer.getAmino(j).getCode()==CYS)
			    {
			      // Zweites Cystein gefunden
			      if((MyPrivateSpacer.getAmino(i).
				  getSideChain()[SG].distance
				  (MyPrivateSpacer.getAmino(j).getSideChain()[SG])<2.3)
				 &&(MyPrivateSpacer.getAmino(i).
				    getSideChain()[SG].distance
				    (MyPrivateSpacer.getAmino(j).
				     getSideChain()[SG])>1.8)){
				  MyAllPossibleRotamerSettings.push_back(0);
				  //must quit loop
				  j=MyPrivateSpacer.size()+1;
				  MySSBond=true;
				}
			    }
			}
		    }
		  //no SSBond detected, treat Cystine as anormal case
		  if(!MySSBond)		    {
		      MyAllPossibleRotamerSettings.push_back
		(allRotamersWithFrequencies[MyPrivateSpacer.getAmino(i).getCode()].getNumberOfRotamers());
		    }
		}
	      else {
	      // PRO detected
	      //set Proline to -1
	      MyAllPossibleRotamerSettings.push_back(0);
	    }
	    }
	  
	  //Testausgabe:
	  // cout<<i+1<<".Aminosaure "<<MyAllPossibleRotamerSettings[i]<<" Rotamere\n";
	}
      
      for(unsigned int i=0;i<MyPrivateSpacer.size();i++)	{
	  if(MyAllPossibleRotamerSettings[i]==0)	    {
	      //cerr<<"Mist1 "<<i<<endl;
	      mAllPossibleRotamerSettings.push_back(-1);
	    }
	  else	    {
	      mAllPossibleRotamerSettings.push_back
		(MyLast/*+MyAllPossibleRotamerSettings[i]*/);
	      MyLast+=MyAllPossibleRotamerSettings[i];
	    }
	  //Testausgabe:
	  /*
	  cout<<"Erstes Rotamer an Position "<<i+1
	      <<"steht im Feld aller Rotamere an der"
	      <<mAllPossibleRotamerSettings[i]<<" .ten Stelle\n";
	  */
	}
    }
  else    {
      for(unsigned int i=0;i<MyPrivateSpacer.size();i++){
	  if(mAllPossibleRotamerSettings.size()==0){
		  cerr<<"Error: please set CalculationStage to zero!";
		  exit(1);  
		}
	  if(mAllPossibleRotamerSettings[i]!=-1)    {
	      
	      MyLast+=allRotamersWithFrequencies[MyPrivateSpacer.getAmino(i).getCode()].getNumberOfRotamers();
	      MyAllPossibleRotamerSettings.push_back
	      (allRotamersWithFrequencies[MyPrivateSpacer.getAmino(i).getCode()].getNumberOfRotamers());
	    }
	  else	    {
	      MyAllPossibleRotamerSettings.push_back(0);
	    }
	  if(mFinalResult[i])    {
	      chi_helper(MyPrivateSpacer.getAmino(i), mFinalResult[i]-1);
	    }	
	}
    }
  AECPointer->mCalculationAtoms.clear();
  MyPrivateSpacer.acceptCalculator(AECPointer);
  MyCompleteEnergy=AECPointer->CalculateEnergy
    (AECPointer->mCalculationAtoms,"AMBER",0);
  //if setchi() works uncorrect, else deactivate TurningBias!!!
  double MyTurningBias=0;
  cout<<"Energie der gedrehten Kopie: "<<MyCompleteEnergy<<endl;
  if(MyCompleteEnergy>100000)    {
      MyTurningBias=MyCompleteEnergy;
    }
  //Anzahl insgesamt waehlbarer Rotamere fuer das Protein:
  mNumberOfAllRotamers=MyLast;
  //mTemplateInteractions.reserve(mNumberOfAllRotamers);
  //to insert: for each Rotamer Calculate Interactions
  //AECP means AtomEnergyCalculatorPointer
  //first Calculate mTemplateEnergyBias;
  mCalculationAtoms.clear();
  for(unsigned int k=0;k<MyPrivateSpacer.size();k++)    {
      if((mFinalResult[k]==0)&&(mAllPossibleRotamerSettings[k]!=-1))	{
	  ExpandInteractingTemplate(MyPrivateSpacer,k,-1);
	}
      else	{
	  ExpandInteractingTemplate(MyPrivateSpacer,k,1024);
	}
    }
  
  mTemplateBias=AECPointer->CalculateEnergy(mCalculationAtoms,"AMBER",0);
  if(mNumberOfAllRotamers<1)    {
      cout<<"no more calculations needed!"<<endl;

      return false;
    }
  AtomEnergyCalculator* AECP=AtomEnergyCalculator::AtomEnergy();

  //initialize InteractionMatrices for IrJs -- Min/Max 
  
  mResidueAllInteractions=new double*[mNumberOfAllRotamers+1];
  if(!mResidueAllInteractions)
    return NULL;
  for(int i=0;i<mNumberOfAllRotamers+1;i++)    {
      mResidueAllInteractions[i]=new double[mNumberOfAllRotamers+1];
      if(!mResidueAllInteractions[i])
	return NULL;
    }

  //building minima and maxima interacting matrices for cost estimation needed
  mIrJsEnergyMinima=new double*[mNumberOfAllRotamers+1];
  if(!mIrJsEnergyMinima)
    return NULL;
  for(int i=0;i<mNumberOfAllRotamers+1;i++)    {
      mIrJsEnergyMinima[i]=new double[MyPrivateSpacer.size()+1];
      if(!mIrJsEnergyMinima[i])
	return NULL;
    }
  for(int i=0;i<mNumberOfAllRotamers+1;i++)
    for(unsigned int k=0;k<MyPrivateSpacer.size();k++)      {
	//WARNING: this initialization seems to be wrong, but it isn't:
	//it ensures, that ((no Interaction)== (0.0)) and in 
	//Calculating step later the minimum will initialized with 
	//INFINITY 
	mIrJsEnergyMinima[i][k]=0.0;
      }
  
  //building maxima analogue minima:
  mIrJsEnergyMaxima=new double*[mNumberOfAllRotamers+1];
  if(!mIrJsEnergyMaxima)
    return NULL;
  for(int i=0;i<mNumberOfAllRotamers+1;i++)    {
      mIrJsEnergyMaxima[i]=new double[MyPrivateSpacer.size()+1];
      if(!mIrJsEnergyMaxima[i])
	return NULL;
    }

  for(int i=0;i<mNumberOfAllRotamers+1;i++)
    for(int unsigned k=0;k<MyPrivateSpacer.size();k++)     {
	//setting 0.0 if you want to use CBCutoffs
	mIrJsEnergyMaxima[i][k]=0.0;
      }
  
  for(int i=0;i<mNumberOfAllRotamers+1;i++)
    for(int k=0;k<mNumberOfAllRotamers+1;k++)     {
	mResidueAllInteractions[i][k]=0.0;
      }
  
  //build the map using CB:
  // i: Position
  for(unsigned int i=0;i<MyPrivateSpacer.size();i++)    {
      int MyInteractionCounter=0;
      // R: Rotamer
      for(int R=0;R<MyAllPossibleRotamerSettings[i];R++){
	mCalculationAtoms.clear();
	//BindInRotamerI: to be implemented
	int MyFirst=mAllPossibleRotamerSettings[i]+R;

	vector<double> MyAllCheese
	  =MyPrivateSpacer.getAmino(i).getSideChain().getChi();

	//skip first component of Rotamer in library which is not an angle
	//but relative frequency
	chi_helper(MyPrivateSpacer.getAmino(i), R);

		MyAllCheese=MyPrivateSpacer.getAmino(i).getSideChain().getChi();

	mCalculationAtoms.clear();
	int MyBetaBlocker=ExpandInteractingTemplate(MyPrivateSpacer,i,R);
	unsigned int MYSPTEMP
	  =MyPrivateSpacer.getAmino(i).getSideChain().size()
	  -MyBetaBlocker;
	//ExpandInteractingTemplate(*mSpacer,i,R);

	mIrJsInteractionAtoms.clear();
	ExpandInteractingRotamers(MyPrivateSpacer,i,R);
	double IREnergy=AECP->CalculateEnergy
	  (mIrJsInteractionAtoms,"AMBER",0); 
	mRotOwnEnergy.push_back(IREnergy);
	      	      
	for(unsigned int j = 0; j < MyPrivateSpacer.size(); j++)	  {
	    double MyIrJsMax = -INFINITY;
	    double MyIrJsMin = INFINITY;
	    if (j != i)	      {
		//loop will be skipped if Residue j ==ALA/GLY
		bool MyPosJAddedToTemplate=false;
		
		for(int S=0;S<MyAllPossibleRotamerSettings[j];S++)		  {
		    int MySecond=mAllPossibleRotamerSettings[j]+S;
		    //this if for building map via CB:bug!
		    //if(TestCBCutoff(i,j))
		    
		    if((MyPrivateSpacer.getAmino(i).getSideChain()[CB]).distance
		       (MyPrivateSpacer.getAmino(j).getSideChain()[CB])<
		       AECP->mCACutOffs[MyPrivateSpacer.getAmino(i).getCode()]
		       [MyPrivateSpacer.getAmino(j).getCode()])      {
			chi_helper(MyPrivateSpacer.getAmino(j), S);
			
			if (!MyPosJAddedToTemplate)  {
			    ExpandInteractingTemplate(MyPrivateSpacer,j,-1);
			    //ExpandInteractingTemplate(*mSpacer,j,-1);
			    MyPosJAddedToTemplate=true;
			    //cout<<"Groesse des Template nach Einfuegen von"
				//<<" Saeure "<<j+1<<" ist: "
			    //<<mCalculationAtoms.size()<<endl;
			    //doing this only for first Romtamer in Pos i&&j:
			    if((S==0)&&(R==0))
			      MyInteractionCounter++;
			  }
			
			if(i<j)  {
			    mIrJsInteractionAtoms.clear();
			    ExpandInteractingRotamers(MyPrivateSpacer,i,R);
			    unsigned int MYSP=mIrJsInteractionAtoms.size();
			    ExpandInteractingRotamers(MyPrivateSpacer,j,S);
			    //cout<<"Laenge der interagierenden Rotamere: "
			    //<<mIrJsInteractionAtoms.size()<<endl;

			    
			    double IrJsInteractionEnergy=
			      AECP->CalculateEnergy
			      (mIrJsInteractionAtoms,"AMBER",MYSP);
			   
			    mResidueAllInteractions[MyFirst][MySecond]
			      =IrJsInteractionEnergy;
			    //hier evt Spiegeln:
			    mResidueAllInteractions[MySecond][MyFirst]
			      =IrJsInteractionEnergy;
			    
			    if(IrJsInteractionEnergy<MyIrJsMin)      {
				mIrJsEnergyMinima[MyFirst][j]=IrJsInteractionEnergy;
				MyIrJsMin=IrJsInteractionEnergy;
			      }
			    if(IrJsInteractionEnergy>MyIrJsMax)     {
				mIrJsEnergyMaxima[MyFirst][j]=IrJsInteractionEnergy;
				MyIrJsMax=IrJsInteractionEnergy;
			      }
			  }
		      }
		  }
				
		//only TemplateExpansion:
		//how to expand if CYS/PRO or other fixed residues???
		if(mAllPossibleRotamerSettings[j]==-1)  {
		    if(MyPrivateSpacer.getAmino(j).getCode()!=GLY)     {
			//cout<<"The ALA or fixed CASE"<<i+1<<" "<<j+1<<endl;
			if((MyPrivateSpacer.getAmino(i).getSideChain()[CB]).distance
			   (MyPrivateSpacer.getAmino(j).getSideChain()[CB])<
			   AECP->
			   mCACutOffs[MyPrivateSpacer.getAmino(i).getCode()]
			   [MyPrivateSpacer.getAmino(j).getCode()])  {
			    ExpandInteractingTemplate(MyPrivateSpacer,j,1024);
			    //ExpandInteractingTemplate(*mSpacer,j,1024);
			    if(R==0)
			      MyInteractionCounter++;
			  }
		      }
		    else     {
			if((MyPrivateSpacer.getAmino(i).getSideChain()[CB]).distance
			   (MyPrivateSpacer.getAmino(j)[CA])<
			   AECP->
			   mCACutOffs[MyPrivateSpacer.getAmino(i).getCode()]
			   [MyPrivateSpacer.getAmino(j).getCode()])  {
			    ExpandInteractingTemplate(MyPrivateSpacer,j,1024);
			    //ExpandInteractingTemplate(*mSpacer,j,1024);
			    if(R==0)
			      MyInteractionCounter++;
			  }
		      }
		  }
	      }
	  }

	double iRTemplateInteractionEnergy=
	  AECP->CalculateEnergy
	  (mCalculationAtoms,"AMBER",MYSPTEMP);
	//Prepare E_i_R,template for DEE--Algorithm:
	RotamerBackboneInteractionElem MyNext;
	//deactivate MyTurningBias if setchi() works correct:
	MyNext.InteractionEnergy=
	  iRTemplateInteractionEnergy;
	MyNext.RelFrequency=allRotamersWithFrequencies[MyPrivateSpacer.getAmino(i).getCode()].getRotamer(R).getFrequency();
	MyNext.DeadEndingRotamer=false;
	mTemplateInteractions.push_back(MyNext);
	}
      mPositionInteractionCounters.push_back(MyInteractionCounter);	  
    }
  //Testausgaben fuer alle iR jS rotamerische WEWIs
  //Calculate Minima for nonmirrorable Interactions

  for(unsigned int i=0;i<MyPrivateSpacer.size();i++)    {
      for(int R=0;R<MyAllPossibleRotamerSettings[i];R++)	{
	  int MyFirst=mAllPossibleRotamerSettings[i]+R;
	  for(unsigned int j=0;j<i;j++)	    {
	      double MyIrJsMax=-INFINITY;
	      double MyIrJsMin=INFINITY;
	
	      for(int S=0;S<MyAllPossibleRotamerSettings[j];S++){
		  int MySecond=mAllPossibleRotamerSettings[j]+S;
		  //this if for building map via CB:bug!
		  //if(TestCBCutoff(i,j))
		  
		  if((MyPrivateSpacer.getAmino(i).getSideChain()[CB]).distance
		     (MyPrivateSpacer.getAmino(j).getSideChain()[CB])<
		     AECP->
		     mCACutOffs[MyPrivateSpacer.getAmino(i).getCode()]
		     [MyPrivateSpacer.getAmino(j).getCode()])  {
		      if(mResidueAllInteractions[MyFirst][MySecond]<MyIrJsMin){
			  mIrJsEnergyMinima[MyFirst][j]
			    =mResidueAllInteractions[MyFirst][MySecond];
			  MyIrJsMin
			    =mResidueAllInteractions[MyFirst][MySecond];
			}
		      if(mResidueAllInteractions[MyFirst][MySecond]>MyIrJsMax)	{
			  mIrJsEnergyMaxima[MyFirst][j]
			    =mResidueAllInteractions[MyFirst][MySecond];
			  MyIrJsMax=
			    mResidueAllInteractions[MyFirst][MySecond];  
			}
		    }
		}
	    }
	}
    }
  
  
  //Setting mTemplateInteractionMinima used for EstimatedEnergyCostMinima in
  //AStarSearch
  
  for(unsigned int i=0;i<mSpacer->size();i++)    {
      mTemplateInteractionMinima.push_back(0.0);
    }
  for(unsigned int i=0;i<mSpacer->size();i++)    {
      if(mAllPossibleRotamerSettings[i]!=-1)	{
	  double MyMinimum=INFINITY;
	  int R;
	  for(R=mAllPossibleRotamerSettings[i];R<GetLastPos(i);R++)	    {
	      if(mTemplateInteractions[R].InteractionEnergy<MyMinimum)
		MyMinimum=mTemplateInteractions[R].InteractionEnergy;
	    }
	  // R==mAllPossibleRotamerSettings[i]) => "Dumpfbacke"
	  mTemplateInteractionMinima[i]=MyMinimum;
	}
      else	{
	  // cerr<<"Mist"<<i<<endl;
	}
    }

  //Setting mIJRotamericPositionInteractionMinima used for Estimation  
  
  mIJRotamericPositionInteractionMinima=new double*[mSpacer->size()+1];
  if(!mIJRotamericPositionInteractionMinima)
    return NULL;
  for(unsigned int i=0;i<mSpacer->size()+1;i++)    {
      mIJRotamericPositionInteractionMinima[i]=
	new double[mSpacer->size()+1];
      if(!mIJRotamericPositionInteractionMinima[i])
	return NULL;
    }
  
  //PreInitialization
  for(unsigned int i=0;i<mSpacer->size();i++)
    for(unsigned int j=0;j<mSpacer->size();j++)     {
	mIJRotamericPositionInteractionMinima[i][j]=0.0;
      }
  
  
  for(unsigned int i=0;i<mSpacer->size();i++)    {
      for(unsigned int j=i+1;j<mSpacer->size();j++)
	{
	  if((mAllPossibleRotamerSettings[i]!=-1)
	     &&(mAllPossibleRotamerSettings[j]!=-1))	    {
	      double MyMinimum=INFINITY;
	      for(int R=mAllPossibleRotamerSettings[i];R<GetLastPos(i);R++){
		  if(mIrJsEnergyMinima[R][j]<MyMinimum)		    {
		      MyMinimum=mIrJsEnergyMinima[R][j];
		    }
		}
	      mIJRotamericPositionInteractionMinima[i][j]=MyMinimum;
	      mIJRotamericPositionInteractionMinima[j][i]=MyMinimum;
	    }
	}
    }

  //calculate v_i_r like Leach and Lemon
  vector<double> MyLeachLemonHelpers;
  for(int i=0;i<mNumberOfAllRotamers;i++)    {
      double MySum=0;
      MySum+=mTemplateInteractions[i].InteractionEnergy;
      for(unsigned int j=0;j<mSpacer->size();j++)	{
	  MySum=MySum+mIrJsEnergyMinima[i][j];
	}
      MyLeachLemonHelpers.push_back(MySum);
    }

  for(unsigned int i=0;i<mSpacer->size();i++)    {
      if(mAllPossibleRotamerSettings[i]!=-1)	{
	  double MyFirstMin=INFINITY;
	  double MySecondMin=INFINITY;
	  for(int R=mAllPossibleRotamerSettings[i];R<GetLastPos(i);R++)	    {
	      //cout<<"Position i"<<R+1<<" wird gemint!"<<endl;
	      
	      if(MyLeachLemonHelpers[R]<MyFirstMin)		{
		  MySecondMin=MyFirstMin;
		  MyFirstMin
		    =MyLeachLemonHelpers[R];
		}
	      else		{
		  if(MyLeachLemonHelpers[R]    <MySecondMin)  {
		      MySecondMin
			=MyLeachLemonHelpers[R];
		    }
		}
	    }
	  double MyDifference=(MySecondMin-MyFirstMin);
	  mLeachLemonIndicators.push_back(MyDifference);
	}
      else	{
	  mLeachLemonIndicators.push_back(0.0);
	}
    }
 
 return true;
}

int SideChainPlacement::ExpandInteractingTemplate(Spacer& fTemplateSource,
 int fInSequence, int fRotamericChoice){
  //eventuell die Schleifenenden optimieren:unsigned int MyEnd;

  int MyBetaBlocker=0;
  if(fRotamericChoice!=-1)
    {
      //Put complete Acid in  the Template:(first or fixed or GLY residue)
      for(unsigned int i=fTemplateSource.getAmino(fInSequence).sizeBackbone();
	  i<fTemplateSource.getAmino(fInSequence).size();i++)
	{
	  if(!((fTemplateSource.getAmino(fInSequence)[i].getCode()==CB) ||
	       (fTemplateSource.getAmino(fInSequence)[i].
		getCode()==HB)||
	       (fTemplateSource.getAmino(fInSequence)[i].
		getCode()==HB1)||
	       (fTemplateSource.getAmino(fInSequence)[i].
		getCode()==HB2)||
	       (fTemplateSource.getAmino(fInSequence)[i].
		getCode()==HB3)||
	       (fTemplateSource.getAmino(fInSequence)[i].
		getCode()==HA)||
	       (fTemplateSource.getAmino(fInSequence)[i].
		getCode()==HA1)||
	       (fTemplateSource.getAmino(fInSequence)[i].
		getCode()==HA2)))
	    mCalculationAtoms.push_back
	      (&(fTemplateSource.getAmino(fInSequence)[i]));
	}
    
      for(unsigned int i=fTemplateSource.getAmino(fInSequence).sizeBackbone();
	  i<fTemplateSource.getAmino(fInSequence).size();i++)
	{
	  if((fTemplateSource.getAmino(fInSequence)[i].getCode()==CB) || 
	     (fTemplateSource.getAmino(fInSequence)[i].
	      getCode()==HB)||
	     (fTemplateSource.getAmino(fInSequence)[i].
	      getCode()==HB1)||
	     (fTemplateSource.getAmino(fInSequence)[i].
	      getCode()==HB2)||
	     (fTemplateSource.getAmino(fInSequence)[i].
	      getCode()==HB3)||
	     (fTemplateSource.getAmino(fInSequence)[i].
	      getCode()==HA)||
	     (fTemplateSource.getAmino(fInSequence)[i].
	      getCode()==HA1)||
	     (fTemplateSource.getAmino(fInSequence)[i].
	      getCode()==HA2))
	    {
	      mCalculationAtoms.push_back
		(&(fTemplateSource.getAmino(fInSequence)[i]));
	      MyBetaBlocker++;
	    }
	}
    }
  
  //put BackBone in Template 
  for(unsigned int i=0;
      i<fTemplateSource.getAmino(fInSequence).sizeBackbone();i++)
    {
      mCalculationAtoms.push_back(&(fTemplateSource.getAmino(fInSequence)[i]));
    }

  if(fRotamericChoice==-1)
    {
      for(unsigned int i=fTemplateSource.getAmino(fInSequence).sizeBackbone();
	  i<fTemplateSource.getAmino(fInSequence).size();i++)
	{
	  if((fTemplateSource.getAmino(fInSequence)[i].getCode()==CB)
	     ||
	     (fTemplateSource.getAmino(fInSequence)[i].
	      getCode()==HB)||
	     (fTemplateSource.getAmino(fInSequence)[i].
	      getCode()==HB1)||
	     (fTemplateSource.getAmino(fInSequence)[i].
	      getCode()==HB2)||
	     (fTemplateSource.getAmino(fInSequence)[i].
	      getCode()==HB3)||
	     (fTemplateSource.getAmino(fInSequence)[i].
	      getCode()==HA)||
	     (fTemplateSource.getAmino(fInSequence)[i].
	      getCode()==HA1)||
	     (fTemplateSource.getAmino(fInSequence)[i].
	      getCode()==HA2))
	    {
	      mCalculationAtoms.push_back
		(&(fTemplateSource.getAmino(fInSequence)[i]));
	    }
	}
    }
  return MyBetaBlocker;
}


void SideChainPlacement::ExpandInteractingRotamers(Spacer& fRotamerSource,
 int fInSequence, int fRotamericChoice){
  //Vorsicht: Nur SideChainAtoms, da BackboneWeWi in ExpandInterActimgTemplate
  //eventuell die Schleifenenden optimieren:unsigned int MyEnd;
  if(fRotamericChoice>-1)    {
      for(unsigned int i=0;
	  i<fRotamerSource.getAmino(fInSequence).getSideChain().size();i++)	{
	  if((fRotamerSource.getAmino(fInSequence).getSideChain()[i].
	      getCode()!=CB)&&
	     (fRotamerSource.getAmino(fInSequence).getSideChain()[i].
	      getCode()!=HB)&&
	     (fRotamerSource.getAmino(fInSequence).getSideChain()[i].
	      getCode()!=HB1)&&
	     (fRotamerSource.getAmino(fInSequence).getSideChain()[i].
	      getCode()!=HB2)&&
	     (fRotamerSource.getAmino(fInSequence).getSideChain()[i].
	      getCode()!=HB3)&&
	     (fRotamerSource.getAmino(fInSequence).getSideChain()[i].
	      getCode()!=HA)&&
	     (fRotamerSource.getAmino(fInSequence).getSideChain()[i].
	      getCode()!=HA1)&&
	     (fRotamerSource.getAmino(fInSequence).getSideChain()[i].
	      getCode()!=HA2))
	    mIrJsInteractionAtoms.push_back
	      (&(fRotamerSource.getAmino(fInSequence).getSideChain()[i]));
	}
      return;
    }
  cout<<"FEHLER in Rotamerauswahl!\n";
  return;
}

int SideChainPlacement::GetLastPos(int fFirstPos){
  int MyLastPos=-1;

  if(mAllPossibleRotamerSettings[fFirstPos]==-1)     {
      return 0;
    } 

  if(static_cast<unsigned int>(fFirstPos)>mAllPossibleRotamerSettings.size())    {
      cout<<"ERROR: in function GetLastPos\n";
      return MyLastPos;
    }
  //Special case for last Position in Spacer is needed:
  if(static_cast<unsigned int>(fFirstPos)==mAllPossibleRotamerSettings.size()-1)    {
      MyLastPos=mAllPossibleRotamerSettings[fFirstPos]
	     +allRotamersWithFrequencies[mSpacer->getAmino(fFirstPos).getCode()].getNumberOfRotamers();
      return MyLastPos;
    }

  do{
      fFirstPos++;
  } while((static_cast<unsigned int>(fFirstPos)
	 <mAllPossibleRotamerSettings.size()-1)
	&&(mAllPossibleRotamerSettings[fFirstPos]==-1));

  if(mAllPossibleRotamerSettings[fFirstPos]==-1)    {
      MyLastPos=0;
      for(unsigned int k=0;k<mSpacer->size();k++){
	  //this if is necessary for hierarchical DEE
	  if(mAllPossibleRotamerSettings[k]!=-1)
	    MyLastPos+=allRotamersWithFrequencies[mSpacer->getAmino(k).getCode()].getNumberOfRotamers();
	}
      return MyLastPos;
    }

  // normaler Fall
  MyLastPos=mAllPossibleRotamerSettings[fFirstPos];
  return MyLastPos;
}

void SideChainPlacement::PrepareSpacer(Spacer& fSpacer){
  this->mSpacer=&fSpacer;
  for(unsigned int i=0;i<fSpacer.components.size();i++)    {
      fSpacer[i].acceptOptimizer(this);  
    }
  return;
}

void SideChainPlacement::PrepareGroup(Group& fGroup){
  vector<Atom> MyGroupAtoms=fGroup.giveAtoms();
  mCalculationAtoms.reserve(mCalculationAtoms.size()+MyGroupAtoms.size());
  for(unsigned int i=0;i<MyGroupAtoms.size();i++)    {
      mCalculationAtoms.push_back(&MyGroupAtoms[i]);
    }
}

void SideChainPlacement::PrepareAminoAcid(AminoAcid& fAminoAcid){
  for(unsigned int i=0;i<fAminoAcid.sizeBackbone();i++)    {
      mAllBackBoneAtoms.push_back(&fAminoAcid[i]);
    }
  PrepareSideChain(fAminoAcid.getSideChain());
}

void SideChainPlacement::PrepareSideChain(SideChain& fSideChain){
  //GLY all HA atoms treated as Backbone:
  if(fSideChain.size()!=0)    {
      if(fSideChain[0].getSuperior().getCode()==GLY)	{
	  for(unsigned int i=0;i<fSideChain.size();i++)	    {
	      mAllBackBoneAtoms.push_back(&fSideChain[i]);
	    }
	  return;
	}
    }
  //normal case:
  for(unsigned int i=0;i<fSideChain.size();i++)    {
      mAllSideChainAtoms.push_back(&fSideChain[i]);
    }
}

void SideChainPlacement::BuildStartNode(){
  //starts building A*---SearchTree:
  QueueElem MyFirstQueueElem;
  MyFirstQueueElem.ActualEnergyCosts=0.0+mTemplateBias;
  MyFirstQueueElem.NodeDepth=0;
  for(unsigned int i=0;i<mMaxNodeDepth;i++)    {
      //Caution 0 is the first Rotamer!
      MyFirstQueueElem.RotamerList.push_back(-1);
    }

  MyFirstQueueElem.EstimatedEnergyCostMinimum=0.0+mTemplateBias;
  //first add for each SideChainPosition its MinimumTemplateInteraction
  for(unsigned int i=0;i<mSpacer->size();i++)    {
      if(mAllPossibleRotamerSettings[i]!=-1)	{
	  MyFirstQueueElem.EstimatedEnergyCostMinimum+=
	    mTemplateInteractionMinima[i];
	}
    }

  //first add for each SideChainPosition its MinimumSideChainSideChain
  //Interaction 
  // i and k are Spacerpositions  
  for(unsigned int i=0;i<mNodePositionInSpacer.size();i++)     {
      for(unsigned int k=i+1;k<mNodePositionInSpacer.size();k++)	{
	  MyFirstQueueElem.EstimatedEnergyCostMinimum+=
	    mIJRotamericPositionInteractionMinima
	   [mNodePositionInSpacer[i]][mNodePositionInSpacer[k]];
	}
    }

  //mAStarPriorityQueue.push_back(MyFirstQueueElem);
  //for more iterations, the priority Queue must resetted:
  while (!(mBStarPriorityQueue.empty()))    {
      mBStarPriorityQueue.pop();
    }
  mBStarPriorityQueue.push(MyFirstQueueElem);
  return;
}

unsigned int SideChainPlacement::ChooseExpandingPos(const QueueElem& fExpandNode,
				       unsigned int fMaxSearchDepth){
  //"Brute force"--method:
  unsigned int i=0;
  bool MyChoosed=false;

  while((i<mMaxNodeDepth)&&(!MyChoosed))    {
      if((fExpandNode.RotamerList[i])!=-1)	{
	  i++;
	}
      else	{
	  // i wird ausgewaehlt.
	  MyChoosed=true;
	}
    }
  return i;
}

unsigned int SideChainPlacement::ChooseExpandingPosViaEntropy(const QueueElem& fExpandNode, unsigned int fSearchDepth){
  //"Brute force"--method:
  unsigned int i=0;
  //this Minimum assuming equal distribution of 28 ARGININE--residues, it has
  //to be changed if more than 28 Rotamers of one sidechain exist
  //MyMin=-log(1/N), where N is max Number of available Rotamers :
  double MyMin=INFINITY;
  int MyMinPos=0;

  while(i<mMaxNodeDepth)    {

      if((fExpandNode.RotamerList[i])!=-1)	{
	}
      else	{
	  // i wird ausgewaehlt
	  if(MyMin>mPositionEntropies[mNodePositionInSpacer[i]])    {
	      MyMin=mPositionEntropies[mNodePositionInSpacer[i]];
	      MyMinPos=i;
	    }
	}
      i++;
    }
  return MyMinPos;
}

unsigned int SideChainPlacement::ChooseExpandingPosViaEntropyCombinedLeachLemon
(const QueueElem& fExpandNode, unsigned int fSearchDepth){
  //"Brute force"--method:
  unsigned int i=0;
  //this Minimum assuming equal distribution of 28 ARGININE--residues, it has
  //to be changed if more than 28 Rotamers of one sidechain exist
  //MyMin=-log(1/N), where N is max Number of available Rotamers :
  double MyMin=INFINITY;
  int MyMinPos=0;

  while(i<mMaxNodeDepth)    {
      if((fExpandNode.RotamerList[i])!=-1)	{
	}
      else	{
	  // i wird ausgewaehlt
	  if(MyMin>mEntropyCombinedLeachLemon[mNodePositionInSpacer[i]])  {
	      MyMin=mEntropyCombinedLeachLemon[mNodePositionInSpacer[i]];
	      MyMinPos=i;
	    }
	}
      i++;
    }
  return MyMinPos;
}

unsigned int SideChainPlacement::ChooseExpandingPosViaLeachLemon
(const QueueElem& fExpandNode, unsigned int fSearchDepth){
  //"Brute force"--method:
  unsigned int i=0;
  //this Minimum assuming equal distribution of 28 ARGININE--residues, it has
  //to be changed if more than 28 Rotamers of one sidechain exist
  //MyMin=-log(1/N), where N is max Number of available Rotamers :
  double MyMax=-INFINITY;
  int MyMaxPos=0;

  while(i<mMaxNodeDepth){
      if((fExpandNode.RotamerList[i])!=-1){
	}
      else{
	  if(MyMax<mLeachLemonIndicators[mNodePositionInSpacer[i]]) {
	    MyMax=mLeachLemonIndicators[mNodePositionInSpacer[i]];
	    MyMaxPos=i;
	    }
	}
      i++;
    }
  return MyMaxPos;
}

/*This Method implements one Step in Searching the AStarTree*/
void SideChainPlacement::ExpandAStarTree(bool fIterative, unsigned int fNumberOfIterations, 
				    double fCutOff, int fMaxResultSize){
  //fCutOff in kcal/mole
  unsigned int MyMaxSearchDepth=mMaxNodeDepth;
  if(fIterative==true)    {
      if(fNumberOfIterations<mMaxNodeDepth)
	MyMaxSearchDepth=fNumberOfIterations;
      else
	MyMaxSearchDepth=mMaxNodeDepth;
    }

  if(MyMaxSearchDepth>0)    {
      bool MyContinue=true;
      int MyThrowAways=0;
      mResultQueue.setMaxSize(65535);
      int MyActualResultSize=0;
      while(!(mBStarPriorityQueue.empty())&&(MyActualResultSize<fMaxResultSize)
	    &&(MyContinue))	{
	  const QueueElem MyTopElem=mBStarPriorityQueue.top();
	  mBStarPriorityQueue.pop();
	  
	  int MyNextRotamericPos = ChooseExpandingPosViaEntropy
	    (MyTopElem,MyMaxSearchDepth);
	  
      int MyNewQueueSize = allRotamersWithFrequencies[mSpacer->getAmino(mNodePositionInSpacer[MyNextRotamericPos]).getCode()].getNumberOfRotamers();
      
      for(int i=0;i<MyNewQueueSize;i++)	{
	  int MyNewRotamer
	    =mAllPossibleRotamerSettings[mNodePositionInSpacer
					[MyNextRotamericPos]]+i;
	  if(!mTemplateInteractions[MyNewRotamer].DeadEndingRotamer)	    {	      
	      QueueElem MyNewQueueElem;
	      //inserting new nodes in searchtree:
	      MyNewQueueElem.RotamerList.clear();
	      for(unsigned int k=0;k<MyTopElem.RotamerList.size();k++)	{
		  MyNewQueueElem.RotamerList.push_back
		    (MyTopElem.RotamerList[k]);
		}
	      MyNewQueueElem.RotamerList[MyNextRotamericPos]=i;
	      
	      MyNewQueueElem.ActualEnergyCosts=
		MyTopElem.ActualEnergyCosts
		+mTemplateInteractions[MyNewRotamer].InteractionEnergy
		+mRotOwnEnergy[MyNewRotamer];
	      //add new TemplateInteraction:
	      
	      //add all Interactions with Rotamers from MyTopElem:
	      for(int k=0;k<MyNextRotamericPos;k++)		{
		  int MyOldRotamer;
		  if(MyNewQueueElem.RotamerList[k]!=-1)		    {
		      MyOldRotamer=
			mAllPossibleRotamerSettings[mNodePositionInSpacer
						   [k]]+MyNewQueueElem.RotamerList[k];
		      MyNewQueueElem.ActualEnergyCosts+=mResidueAllInteractions
			[MyNewRotamer][MyOldRotamer];
		    }
		}
	      for(int k=MyNextRotamericPos+1;
		  k<static_cast<int>(mMaxNodeDepth);k++){
		  int MyOldRotamer;
		  if(MyNewQueueElem.RotamerList[k]!=-1)    {
		      MyOldRotamer=
			mAllPossibleRotamerSettings[mNodePositionInSpacer
						   [k]]+MyNewQueueElem.RotamerList[k];
		      MyNewQueueElem.ActualEnergyCosts+=mResidueAllInteractions
			[MyNewRotamer][MyOldRotamer];
		    }
		}
	      MyNewQueueElem.NodeDepth=MyTopElem.NodeDepth+1;
	      
	      MyNewQueueElem.EstimatedEnergyCostMinimum=
	      MyTopElem.EstimatedEnergyCostMinimum
	      -mTemplateInteractionMinima
		[mNodePositionInSpacer[MyNextRotamericPos]]
	      +mTemplateInteractions[MyNewRotamer].InteractionEnergy
		+mRotOwnEnergy[MyNewRotamer];

		      //next two loops are very sophisticated estimations
	      for(int k=0;k<MyNextRotamericPos;k++)		{
		  if(MyNewQueueElem.RotamerList[k]==-1)		    {
		      MyNewQueueElem.EstimatedEnergyCostMinimum=
		      MyNewQueueElem.EstimatedEnergyCostMinimum
		      -mIJRotamericPositionInteractionMinima
		      [mNodePositionInSpacer[k]]
		      [mNodePositionInSpacer[MyNextRotamericPos]]
		      +mIrJsEnergyMinima[MyNewRotamer]
		      [mNodePositionInSpacer[k]];
		    }
		  else    {
		      
		      MyNewQueueElem.EstimatedEnergyCostMinimum=
		      MyNewQueueElem.EstimatedEnergyCostMinimum
		      -mIrJsEnergyMinima[MyNewQueueElem.RotamerList[k]
		      		  +mAllPossibleRotamerSettings[mNodePositionInSpacer[k]]]
		      [mNodePositionInSpacer[MyNextRotamericPos]]
		      +mResidueAllInteractions
		      [MyNewQueueElem.RotamerList[k]
		      +mAllPossibleRotamerSettings[mNodePositionInSpacer[k]]]
		      [MyNewRotamer];
		      
		    }
		}

	     
	      for(int k=MyNextRotamericPos+1;
		  k<static_cast<int>(mMaxNodeDepth);k++)		{
		  if(MyNewQueueElem.RotamerList[k]==-1)		    {
		      
		      MyNewQueueElem.EstimatedEnergyCostMinimum=
		      MyNewQueueElem.EstimatedEnergyCostMinimum
		      -mIJRotamericPositionInteractionMinima
		      [mNodePositionInSpacer[k]]
		      [mNodePositionInSpacer[MyNextRotamericPos]]
		      +mIrJsEnergyMinima[MyNewRotamer]
		      [mNodePositionInSpacer[k]];
		      
		    }
		  else    {
		      //MyNewQueueElem.EstimatedEnergyCostMinimum+=
		      //mResidueAllInteractions
		      //[MyNewQueueElem.RotamerList[k]
		      //+mAllPossibleRotamerSettings[mNodePositionInSpacer[k]]]
		      //[MyNewRotamer];
		      MyNewQueueElem.EstimatedEnergyCostMinimum=
		      MyNewQueueElem.EstimatedEnergyCostMinimum
		      -mIrJsEnergyMinima[MyNewQueueElem.RotamerList[k]
		      		  +mAllPossibleRotamerSettings[mNodePositionInSpacer[k]]]
		      [mNodePositionInSpacer[MyNextRotamericPos]]
		      +mResidueAllInteractions
		      [MyNewQueueElem.RotamerList[k]
		      +mAllPossibleRotamerSettings
		      [mNodePositionInSpacer[k]]]			
		      [MyNewRotamer];
		    }
		}
	      if(MyNewQueueElem.NodeDepth<MyMaxSearchDepth){
		  
		    mBStarPriorityQueue.push(MyNewQueueElem);
		  
		}
	      else{
		  if(MyNewQueueElem.EstimatedEnergyCostMinimum<fCutOff* mSpacer->size()) {
		      //MyNewQueueElem.EstimatedEnergyCostMinimum=
		      //MyNewQueueElem.ActualEnergyCosts;
		      mResultQueue.push(MyNewQueueElem);
		      MyActualResultSize++;
		     
		    }
		  else {
		      MyThrowAways++;
		      // cerr<<"Wirfs weg "<<endl;
		    }
		  
		}
	    }
	  
	} 
	}
    }

//here comes the great delete, big fun!:
 if (mNumberOfAllRotamers!=0)   {
     for (int i=0;i<mNumberOfAllRotamers+1;i++)       {
	 delete[] mResidueAllInteractions[i];
       }
     delete[] mResidueAllInteractions;

     for (int i=0;i<mNumberOfAllRotamers+1;i++)       {
	 delete[] mIrJsEnergyMinima[i];
       }
     delete[] mIrJsEnergyMinima;
     
     delete[] mIrJsEnergyMaxima;
    
     delete[] mIJRotamericPositionInteractionMinima;
   }

 if((MyMaxSearchDepth==0)&&(fIterative))   {
     cerr<<"Ziel wurde bereits erreicht!"<<endl;
     exit(1);
   }

 return;
}

void SideChainPlacement::AStarSearch(){
  BuildStartNode();

  ExpandAStarTree(1,15,8200,2048); 
}

void SideChainPlacement::ValidateResults(){
  Spacer MyPrivateSpacer = *mSpacer;
  int i = 0;

  //push first result in mFinalResult:
  if(!mResultQueue.empty())    {
      const QueueElem MyResult = mResultQueue.top();
      for(unsigned int R = 0; R < MyResult.RotamerList.size(); R++)
	  mFinalResult[mNodePositionInSpacer[R]] = MyResult.RotamerList[R]+1;

      //and prepare a possible next calculation stage:
      int MyReduceLevel = 0;
      for(unsigned int R = 0; R < mSpacer->size(); R++)
	if((mFinalResult[R] != 0) && (mAllPossibleRotamerSettings[R] != -1))	  {
	      mAllPossibleRotamerSettings[R] = -1;
	      MyReduceLevel = MyReduceLevel
                + allRotamersWithFrequencies[mSpacer->getAmino(R).getCode()].getNumberOfRotamers();
	  }
	else	  {
	    if(mAllPossibleRotamerSettings[R] != -1)
	      mAllPossibleRotamerSettings[R] =
		mAllPossibleRotamerSettings[R] - MyReduceLevel;
	  }
    }

  while( !(mResultQueue.empty()) && (i < 200) )    {
      const QueueElem MyResult = mResultQueue.top();
      mResultQueue.pop();
      i++;
              
      for (unsigned int k= 0; k < mFinalResult.size(); k++)	{
	  //Copy from spacer to make turns on it:

	  if (mFinalResult[k] != 0)	    {
	      chi_helper(MyPrivateSpacer.getAmino(k), mFinalResult[k]-1);
	    }//each Rotamer is turned (i.e. the complete Spacer)
	  //calculate RMS
	}
      
      if(i == 1)
	calculateResults(MyPrivateSpacer);
    }
  while(!(mResultQueue.empty()))
    mResultQueue.pop();
}

void SideChainPlacement::calculateChi1Correct(DeviationMeasure& currChi1, 
DeviationMeasure& currChi12, AminoAcid& aa, int aminoOffset){
  if (aa.getSideChain().getMaxChi() < 1)
    return;

  currChi1.addSum(static_cast<AminoAcidCode>(aa.getCode()));
  currChi12.addSum(static_cast<AminoAcidCode>(aa.getCode()));

  if ((fabs(aa.getChi(0) - mSpacer->getAmino(aminoOffset).getChi(0)) <= 40.0)
      || (fabs(aa.getChi(0) - mSpacer->getAmino(aminoOffset).getChi(0)) >= 320.0))    {
      currChi1.addCounter(static_cast<AminoAcidCode>(aa.getCode()));
      if ( (aa.getSideChain().getMaxChi() < 2) 
	   || ( fabs(aa.getChi(1) 
		- mSpacer->getAmino(aminoOffset).getChi(1)) <= 40.0)
	   || ( fabs(aa.getChi(1) 
		- mSpacer->getAmino(aminoOffset).getChi(1)) >= 320.0) )
	currChi12.addCounter(static_cast<AminoAcidCode>(aa.getCode()));
    }
}


void SideChainPlacement::calculateResults(Spacer MyPrivateSpacer){
  AtomEnergyCalculator* AECPointer = 
    AtomEnergyCalculator::AtomEnergy();
  AECPointer->mCalculationAtoms.clear();
  MyPrivateSpacer.acceptCalculator(AECPointer);
  double MyCompleteEnergy = AECPointer->CalculateEnergy
    (AECPointer->mCalculationAtoms, "AMBER", 0);  
  cout<<"Energie des Ergebnisses: " << MyCompleteEnergy 
      << endl;
  
  DeviationMeasure rms(true);
  DeviationMeasure chi1;
  DeviationMeasure chi12;

    ofstream MySavedFile("Erg.pdb");
    PdbSaver TestSaver(MySavedFile);
    TestSaver.saveSpacer(MyPrivateSpacer);
  for(unsigned int j = 0; j < MyPrivateSpacer.size(); j++)    {
      RMSCalculate(rms, MyPrivateSpacer.getAmino(j).getSideChain(), j);
      calculateChi1Correct(chi1, chi12, MyPrivateSpacer.getAmino(j), j);
    }
  
  cout << "----------------------------------------------\n";
  cout << "RMS: \n" << rms;
  cout << "----------------------------------------------\n";
  cout << "Chi1: \n" << chi1;
  cout << "----------------------------------------------\n";
  cout << "Chi12: \n" << chi12;
  cout << "----------------------------------------------\n";
}


void SideChainPlacement::findBestFit(){
  Spacer localSpacer = *mSpacer;
  Spacer resultSpacer = *mSpacer;

  for (unsigned int i = 0; i < localSpacer.sizeAmino(); i++)    {
      if ( (localSpacer.getAmino(i).getCode() == GLY) || 
	   (localSpacer.getAmino(i).getCode() == ALA) || 
	   (localSpacer.getAmino(i).getCode() == CYS) || 
	   (localSpacer.getAmino(i).getCode() == PRO) )
	continue;
      DeviationMeasure partialRms(true);
      partialRms.addCounter(ALA, 10000.0);
      partialRms.addSum(ALA);
      
      for (unsigned int r = 0; r < allRotamersWithFrequencies[ 
	     localSpacer.getAmino(i).getCode()].getNumberOfRotamers(); r++)	{
	  DeviationMeasure currentRms(true);
	  chi_helper(localSpacer.getAmino(i), r);
	  RMSCalculate(currentRms, localSpacer.getAmino(i).getSideChain(), i);
	  if (currentRms.getTotalAverage() < partialRms.getTotalAverage())    {
	      partialRms = currentRms;
	      chi_helper(resultSpacer.getAmino(i), r);
	    }
	}
    }
  
  cout << "----------------------------------------------------\n";
  cout << "  Best fit:\n";
  calculateResults(resultSpacer);
  cout << "----------------------------------------------------\n";
}


void SideChainPlacement::calculateEntropyFit(){
  Spacer localSpacer = *mSpacer;
  AtomEnergyCalculator* AECP = AtomEnergyCalculator::AtomEnergy();
  AECP->mCalculationAtoms.clear();
  localSpacer.acceptCalculator(AECP);

  for (unsigned int i = 0; i < localSpacer.sizeAmino(); i++)    {
      if ( (localSpacer.getAmino(i).getCode() == GLY) || 
	   (localSpacer.getAmino(i).getCode() == ALA) || 
	   (localSpacer.getAmino(i).getCode() == CYS) || 
	   (localSpacer.getAmino(i).getCode() == PRO) )
	continue;
      double origEnergy = AECP->CalculateEnergy(
			AECP->mCalculationAtoms, "AMBER", 0);  
      for (unsigned r = 1; r < allRotamersWithFrequencies[
	   localSpacer.getAmino(i).getCode()].getNumberOfRotamers();
	   r++)	{
	  chi_helper(localSpacer.getAmino(i), r);

	  double newEnergy = AECP->CalculateEnergy(
			     AECP->mCalculationAtoms, "AMBER", 0);  
	  if ( fabs(origEnergy - newEnergy) < ENTROPY_THRESHOLD)
	    break;
	}
    }
  cout << "----------------------------------------------------\n";
  cout << "  Greedy entropy fit:\n";
  calculateResults(localSpacer);
  cout << "----------------------------------------------------\n";
}


void SideChainPlacement::optimizeResults(Spacer localSpacer){
  AtomEnergyCalculator* AECP = AtomEnergyCalculator::AtomEnergy();
  
  for (unsigned int i = 0; i < localSpacer.sizeAmino(); i++)    {
      AECP->mCalculationAtoms.clear();
      localSpacer.acceptCalculator(AECP);
      if (localSpacer.getAmino(i).getSideChain().getMaxChi() > 0)	{
	  pOptimizeChi(localSpacer, AECP, i, 0);
	  if (localSpacer.getAmino(i).getSideChain().getMaxChi() > 1)
	    pOptimizeChi(localSpacer, AECP, i, 1);
	    if (localSpacer.getAmino(i).getSideChain().getMaxChi() > 2)
	      pOptimizeChi(localSpacer, AECP, i, 2);
	  }
    }
    

  cout << "--------------------------------------\n";
  cout << "  Optimized result: \n";
  calculateResults(localSpacer);
  cout << "--------------------------------------\n";
}


void SideChainPlacement::pOptimizeChi(Spacer& localSpacer, 
AtomEnergyCalculator* AECP, int i, int currChi){
  PRECOND( localSpacer.getAmino(i).getSideChain().getMaxChi() >= currChi,
	   exception);
  double origEnergy = AECP->CalculateEnergy(
			     AECP->mCalculationAtoms, "AMBER", 0);  
  double opt = 0.0;
  double oldChi = localSpacer.getAmino(i).getChi(currChi);
  
  for (double offset = -5.0; offset <= 5.0; offset += 10.0)    {
      localSpacer.getAmino(i).setChi(currChi, oldChi + offset);
      double partialEnergy = AECP->CalculateEnergy(
 			      AECP->mCalculationAtoms, "AMBER", 0);  
      if (partialEnergy < origEnergy)	{
	  origEnergy = partialEnergy;
	  opt = offset;
	}
    }
  
  localSpacer.getAmino(i).setChi(currChi, oldChi + opt);
}

// This function reads a definition file for rotamers which contains
// probabilities for each conformation.
// An example for a file follows:
//
//  # this denotes a comment till the end of the line
//  # The header of a definition sequence for an Aminoacid starts with
//  # section {Aminoacid name}
//  section ARG  # ARG stands for Arginine
//  
//  # now follows the definition of single conformations and their
//  # relative probability:
//  # {probability} {chi1} {chi2} ...
//  1.12   58.7  -179.8    69.1  -175.4  
//  1.37   63.7   177.8   178.4    86.2  
//  1.88   62.4  -179.1   178.6  -179.1 
//
//  # Note that the number of given chi angles is checked against the
//  # number that the function SideChain::getMaxChi() returns.
//  # Any violation results in an error message
void SideChainPlacement::LoadRotamerLibrary(istream &in) {
  // allocate enough memory...
  allRotamersWithFrequencies.resize(20);

  // now start to read the file.
  while (in.good()) {
    // skip comments and whitespace
    eatComment(in);

    // now (hopefully) read in the Aminoacid section statement
    string token;
    in >> token;

    if (token != "section") {
      ERROR("Error in rotamer library. Expected \"section\", read "
	    + token, exception);
    }

    // We've got a section statment here. Read the aminoacid code...
    in >> token;
    // ... and translate it into an integer:
    int aminoAcidCode = aminoAcidThreeLetterTranslator(token);

    if (aminoAcidCode == XXX) {
      ERROR("LoadRotamerLibrary: Error in rotamer library. \
Unknown aminoacid code:" + token, exception);
    }

    // OK, we found the start of a new section and got the corresponding
    // Aminoacid code.
    // now try to read all rows of rotamer definitions, for one AminoAcid.
    allRotamersWithFrequencies[aminoAcidCode].read(in);
  }
}


void SideChainPlacement::chi_helper(AminoAcid &aa, int i) {
  // the side chain of this amino acid is set to the
  // rotamer from the library

  const Rotamer &rotamer = allRotamersWithFrequencies[aa.getCode()][i];
  aa.getSideChain().setChi(rotamer.getChiAngles());
  aa.getSideChain().sync();
}
