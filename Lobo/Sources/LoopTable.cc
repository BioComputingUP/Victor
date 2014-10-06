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
 *@Class:           LoopTable    
 * 
 *@Description:Defines a table of possible amino chain end points and end directions 
 *  after k amino acids have been concatenated. 
 *      
*/

// Includes:
#include <LoopTable.h>
#include <LoopTableEntry.h>
#include <IntCoordConverter.h>
#include <VectorTransformation.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;

unsigned int LoopTable::MAX_BINS = 128;
unsigned int LoopTable::MAX_FACTOR = 5; // ie. applied as / MAX_FACTOR 
double LoopTable::BOND_ANGLE_AT_CPRIME_TO_N = 112.7;
 
double LoopTable::BOND_LENGTH_N_TO_CALPHA = 1.46;
double LoopTable::BOND_ANGLE_AT_N_TO_CALPHA = 121.0;

const double SIM_WEIGTH = 0.0001;
/**
*@Description Calculates the square value
*@param  value to find the square(float)
*@return resulting value (static float)
*/
static float sqr(float x){
  return x*x;
}
/**
*@Description Calculates the distance between two points
*@param  vectors containing the coords for each point( vgVector3<float>, vgVector3<float>)
*@return  corresponding value(static double)
*/
static double sDistance( vgVector3<float> xv,  vgVector3<float> yv){
  return  sqrt(sqr(xv.x-yv.x) + sqr(xv.y-yv.y) + sqr(xv.z-yv.z));
}


// CONSTRUCTORS/DESTRUCTOR:
/**
*@Description basic constructor
*/
LoopTable::LoopTable() : rama(NULL), nAminoAcid(1), lowerLimit(1000.0), 
  upperLimit(-1000.0), stepLimit(0),  nextIndex(0), stepWidth(1){ 
  vgVector3<float> tmp(1000, 1000, 1000);
  for (unsigned int i = 0; i < 6; i++)    {
      max[i] = -tmp;
      min[i] = tmp;
    }

  entry.resize(MAX_BINS);
 
  // compute a constant expression needed if nAminoAcid = 1:
  vgVector3<float> tmpEndPoint(1.403962, 2.015868, 0);
  vgVector3<float> tmp2(0,BOND_LENGTH_N_TO_CALPHA,0);
  psiNormal = -(tmpEndPoint - tmp2);
  psiNormal.normalize();

}

/**
*@Description constructor base from the copy of another object
*@param reference to the original object (const LoopTable&)
*/
LoopTable::LoopTable(const LoopTable& orig){
 this->copy(orig);
}
/**
*@Description basic destructor
*/
LoopTable::~LoopTable(){
  PRINT_NAME;
}


// PREDICATES:
/**
*@Description Prints the loop table distribution
*@param  none
*/
void LoopTable::showDistribution(){
  for (unsigned int i = 0; i < entry.size(); i++)    {
      cout << " " << entry[i].size() << "\t";
      if ((i+1) % 5 == 0)
	cout << "\n";
    }
  cout << "\n";
}
 
/**
*@Description Find a matching entry with "minimal" deviation from destination 
*    (EP, ED, EN) and return it.
*@param  reference to the destination loop table entry(const LoopTableEntry&), index of the current selection (unsigned int)
*@return  the new loop table entry(LoopTableEntry)
*/
LoopTableEntry LoopTable::getClosest(const LoopTableEntry& dest, 
unsigned int currentSelection){ 
  PRECOND( (entry.size() > 0), exception); 

  priority_queue<solutionQueueElem> solutionQueue;
  unsigned int index = pGetBin(dest.endPoint);
  unsigned int max = size() / MAX_FACTOR;

  pAddToSolutionQueue(index, solutionQueue, dest);

  unsigned int offset = 0;
  while ((solutionQueue.size() < max) || (offset > MAX_BINS/2))    {
      offset++;
      if (index >= offset) 
	pAddToSolutionQueue(index-offset, solutionQueue, dest);
      if (index+offset < MAX_BINS)
	pAddToSolutionQueue(index+offset, solutionQueue, dest);
    }

  if (solutionQueue.size() > currentSelection)
    for (unsigned int i = 0; i < currentSelection; i++)
      solutionQueue.pop();

  return entry[solutionQueue.top().index1][solutionQueue.top().index2];
}
/**
*@Description Gets the closest N 
*@param  reference to the table( LoopTableEntry&), value of N(unsigned int ), number of amino acids (unsigned int )
*@return  vector containing the loop entry tables(vector<LoopTableEntry> )
*/
vector<LoopTableEntry> LoopTable::getNClosest(const LoopTableEntry& dest, unsigned int num, 
unsigned int nAmino){
  PRECOND( (entry.size() > 0), exception); 

  vector<LoopTableEntry> result;
  priority_queue<solutionQueueElem> solutionQueue;
  unsigned int index = pGetBin(dest.endPoint);
  unsigned int max = size() / MAX_FACTOR;

  pAddToSolutionQueue(index, solutionQueue, dest);

  unsigned int offset = 0;
  while ((solutionQueue.size() < max) || (offset > MAX_BINS/2))    {
      offset++;
      if (index >= offset) 
	pAddToSolutionQueue(index-offset, solutionQueue, dest);
      if (index+offset < MAX_BINS)
	pAddToSolutionQueue(index+offset, solutionQueue, dest);
    }

  if (solutionQueue.size() == 0)
    ERROR("No valid entries found.", exception);
  
  double maxDev = 20 * solutionQueue.top().dev;
 
  while ((solutionQueue.size() > 0) && (result.size() < num))    {
      if (solutionQueue.top().dev > maxDev)
	break;

     
      if (nAmino <= 5)      {
      vgVector3<float> offsetEP = dest.endPoint - 
        entry[solutionQueue.top().index1][solutionQueue.top().index2].endPoint;
      entry[solutionQueue.top().index1][solutionQueue.top().index2].midPoint 
	  += (offsetEP/2);
      entry[solutionQueue.top().index1][solutionQueue.top().index2].endPoint 
	  = dest.endPoint;
      }
      
      result.push_back(entry[solutionQueue.top().index1]
		       [solutionQueue.top().index2]);
      solutionQueue.pop();
    }
  
  return result;

}


// MODIFIERS:
/**
*@Description copies the original object into a new one
*@param  reference to the original object(const LoopTable& )
*@return  changes are made internally(void)
*/
void LoopTable::copy(const LoopTable& orig){
  PRINT_NAME; 
  rama = orig.rama;
  nAminoAcid = orig.nAminoAcid;

  entry.clear();
  entry.resize(orig.entry.size());
  for (unsigned int i = 0; i < orig.entry.size(); i++)
    for (unsigned int j = 0; j < orig.entry[i].size(); j++)
      entry[i].push_back(orig.entry[i][j]);

  lowerLimit = orig.lowerLimit;
  upperLimit = orig.upperLimit;
  stepLimit = orig.stepLimit;

  min = orig.min;
  max = orig.max;
  psiNormal = orig.psiNormal;
  nextIndex = orig.nextIndex;
  stepWidth = orig.stepWidth;
}

 
/**
*@Description Defaults the protein table to the simple case where nAminoacid = 1.
*@param  none
*@return   changes are made internally(void)
*/
void LoopTable::setToSingleAminoAcid(){
 nAminoAcid = 1;
 LoopTableEntry tmpEntry;
 tmpEntry.setToSingleAminoAcid();

 entry[0].clear();
 entry[0].push_back(tmpEntry);
}

 
/**
*@Description The function concatenates two protein tables (src1, src2) and stores 
*    the result in this. Adds the src to the vector entry.
*@param  reference to the loop tables(LoopTable& ,LoopTable& )src sizes(unsigned long,unsigned long)
*@return   changes are made internally(void)
*/
void LoopTable::concatenate(LoopTable& src1, LoopTable& src2, 
		       unsigned long nSrc1, unsigned long nSrc2){
  src1.initOccurrence(nSrc1);

  // nSrc1 statistical cases of src1:
  for (unsigned long i = 0; i < nSrc1; i++)
    {
      LoopTableEntry selectedSrc1 = src1.selectOccurrence();
      src2.initOccurrence(nSrc2);
        
      // nSrc2 statistical cases of src2:
      for (unsigned long j = 0; j < nSrc2; j++)	{
	  LoopTableEntry selectedSrc2 = src2.selectOccurrence();
	  LoopTableEntry result = selectedSrc1.concatenate(selectedSrc2, 
		  (src1.nAminoAcid+1));  // concatenate the selected entries 
	  store(result); // write result into entry table
	}
    }
  
  nAminoAcid = src1.nAminoAcid + src2.nAminoAcid;

  adjustTable();
}

/**
*@Description Adjusts the hash table representation, distributing its contents among
*    all bins, according to a min -> max distribution with uniform sampling.
*@param  none
*@return   changes are made internally(void)
*/
void LoopTable::adjustTable(){ // adjusts the internal hash table 

  vector<vector<LoopTableEntry> > tmpEntry;
  tmpEntry.resize(MAX_BINS);

  // set boundaries:

  stepLimit = (upperLimit-lowerLimit) / (MAX_BINS-1);
  
  for (unsigned int i = 0; i < entry.size(); i++)
      for (unsigned int j = 0; j < entry[i].size(); j++)
	pInsertElem(entry[i][j], tmpEntry);

  entry = tmpEntry;
}

/**
*@Description Reads the table contents from a file.
*@param  reference to the file name(const string&)
*@return   changes are made internally(void)
*/
void LoopTable::read(const string& filename){ 
  ifstream is(filename.c_str(), ios::in | ios::binary );
  if (!is)    {
      cout << "Table data file " << filename 
	   << " not found. Please construct it first.\n";
      ERROR("Table data file " + filename + " not found.", exception);
    }

  // first read the global table parameters
  is.read((char*)&nAminoAcid, sizeof(unsigned int));

  unsigned int nBins = 0;
  unsigned int nEntry = 0;

  is.read((char*)&nBins, sizeof(unsigned int));

  entry.clear();
  entry.resize(nBins);

  if (nBins != MAX_BINS)
    ERROR("LoopTable::read() : nBins out of scope.", exception);

  //load array into heap (ie. faster):
  unsigned int* tmpBin = new unsigned int[nBins];
  is.read((char*)tmpBin, sizeof(unsigned int) * nBins);

  for (unsigned int i = 0; i < nBins; i++)    {
      nEntry += tmpBin[i];
      entry[i].resize(tmpBin[i]);
    }
  
  is.read((char*)&lowerLimit, sizeof(double));
  is.read((char*)&upperLimit, sizeof(double));
  stepLimit = (upperLimit-lowerLimit) / (nBins-1);

  // then read the max and min values:
  for (unsigned int i = 0; i < 6; i++)
    for (unsigned int j = 0; j < 3; j++)
      is.read((char*)&max[i][j], sizeof(double));
      
  for (unsigned int i = 0; i < 6; i++)
    for (unsigned int j = 0; j < 3; j++)
      is.read((char*)&min[i][j], sizeof(double));

  // load the complete table into the heap because it's faster
  CompressedLoopTableEntry* compTable = new CompressedLoopTableEntry[nEntry]; 
  is.read((char*)compTable, sizeof(CompressedLoopTableEntry) * nEntry);

  // decode the table from memory itself:
  unsigned long compIndex = 0;

  for (unsigned int index = 0; index < nBins; index++)    {
      for (unsigned  i = 0; i < tmpBin[index]; i++)	{
	  for (unsigned int j = 0; j < 6; j++)
	    for (unsigned int k = 0; k < 3; k++)
	      (entry[index][i])[j][k] = decode(compTable[compIndex][j][k], 
					       min[j][k], max[j][k]);
	  compIndex++;
	}
    }
  
  delete tmpBin;
  delete compTable;
  is.close();
}
 
/**
*@Description Clusters entries below cutoff, leaving only a single representative.
*@param  cutoff value(double)
*@return   changes are made internally(void)
*/
void LoopTable::cluster(double cutoff){

  if (entry.size() == 0)
    return;
  
  vector<vector<LoopTableEntry> > tmpEntry; // initialize tmpEntry:
  tmpEntry.resize(MAX_BINS);

  for (unsigned int i = 0; i < MAX_BINS; i++)
    if (entry[i].size() > 0)      {
	tmpEntry[i].push_back(entry[i][0]);
	
	for (unsigned int ii = 1; ii < entry[i].size(); ii++)	  {
	    bool contained = false;
	    for (unsigned int j = 0; j < tmpEntry[i].size(); j++)  {
		if (sDistance(tmpEntry[i][j].endPoint, 
			      entry[i][ii].endPoint)
		    + sDistance(tmpEntry[i][j].endDirection, 
				entry[i][ii].endDirection)
		    + sDistance(tmpEntry[i][j].endNormal, 
				entry[i][ii].endNormal)
		    <= cutoff)   {
		    contained = true;
		    break;
		  }  
	      }
	    
	    if (!contained)
	      tmpEntry[i].push_back(entry[i][ii]);
	  }
      }

entry = tmpEntry;
}
 
/**
*@Description Writes the table contents to file.
*@param  reference to the file name(const string& )
*@return   changes are made internally(void)
*/
void LoopTable::write(const string& filename){
  ofstream os(filename.c_str(), ios::out | ios::binary);
  
  PRECOND(os, exception);
  
  // first write the global table parameters:
  os.write((char*)&nAminoAcid, sizeof(unsigned int));

  unsigned int nBin = entry.size();

  os.write((char*)&nBin, sizeof(unsigned int));
  for (unsigned int i = 0; i < entry.size(); i++) {
      nBin = entry[i].size();
      os.write((char*)&nBin, sizeof(unsigned int));
    }

  os.write((char*)&lowerLimit, sizeof(double));
  os.write((char*)&upperLimit, sizeof(double));

  for (unsigned int i = 0; i < 6; i++)
    for (unsigned int j = 0; j < 3; j++)
      os.write((char*)&max[i][j], sizeof(double));
  
  for (unsigned int i = 0; i < 6; i++)
    for (unsigned int j = 0; j < 3; j++)
      os.write((char*)&min[i][j], sizeof(double));
  
  // next write the table itself, compressing it on the fly:

  INVARIANT(os, exception);
  for (unsigned int index = 0; index < entry.size(); index++)
    for (unsigned int i = 0; i < entry[index].size(); i++)
      for (unsigned int j = 0; j < 6; j++)
	for (unsigned int k = 0; k < 3; k++)	  {
	    unsigned short tmp = code((entry[index][i])[j][k], 
				      min[j][k], max[j][k]);
	    os.write( (char*)&tmp, sizeof(unsigned short));
	  }

  os.close();
}
 
/**
*@Description Writes the table contents to file in ASCII format.
//    (Useful for viewing with GNU Plot)
*@param  reference to the file name (const string& ), size(unsigned long), entry weigth( unsigned int)
 *    ,entry dimension( unsigned int)
*@return   changes are made internally(void)
*/
void LoopTable::writeASCII(const string& filename, unsigned long num, 
			 unsigned int wEntry, unsigned int wDim){
  ofstream os(filename.c_str(), ios::out);
  PRECOND(os, exception);
  
  if (num == 0)
    num = size();

  initOccurrence(num);
  
  unsigned int w;
  unsigned int v;
  
  if (wDim == 1)     { 
      w = 0; 
      v = 2; 
    }
  else if (wDim == 2)    { 
      w = 1; 
      v = 2;
    }
  else    { 
      w = 0; 
      v = 1; 
    }
  
  for (unsigned long i = 0; i < num; i++)    {
      LoopTableEntry tmpEntry = selectOccurrence();
      os << tmpEntry[wEntry][w] << " ";
      os << tmpEntry[wEntry][v] << "\n";
    }

  os.close();
}



// OPERATORS:
/**
*@Description Assigns one object into another
*@param  reference to the original object(const LoopTable&)
*@return    the reference to the the one(LoopTable&)
*/
LoopTable& LoopTable::operator=(const LoopTable& orig){
  if (&orig != this)
    copy(orig);
  return *this;
}


// HELPERS:
/**
*@Description calculates the esptimate similarity cutoff
*@param reference to the loop table entry (const LoopTableEntry& ),value to consider (unsigned int)
*@return  the corresponding value(double)
*/
double LoopTable::pEstimateSimilarityCutoff(const LoopTableEntry& le, 
unsigned int num){
  double base = 0.01;

  double distanceFactor = 1.0;
  vgVector3<float> tmp(0,0,0);

  double dist = sDistance(tmp, le.endPoint);

  if (dist > (num * 1.25))
    distanceFactor = 0.65;
  else if (dist < (num / 1.25))
    distanceFactor = 1.5;

  return base * distanceFactor * sqr(num);
}
 
/**
*@Description  Adds the src to the vector entry.  stores the results in the table
*@param  reference to the loop table entry(const LoopTableEntry&)
*@return   changes are made internally(void)
*/
void LoopTable::store(const LoopTableEntry& src){
  entry[0].push_back(src);
  
  double tmp = sqrt(sqr(src.endPoint.x) 
		    + sqr(src.endPoint.y) 
		    + sqr(src.endPoint.z));
  if (tmp < lowerLimit)
    lowerLimit = tmp;
  if (tmp > upperLimit)
    upperLimit = tmp;
 
  //update max & min where necessary:
  for (unsigned int i = 0; i < 6; i++)
    for (unsigned int j = 0; j < 3; j++)
      {
	if (src[i][j] > max[i][j])
	  max[i][j] = src[i][j];
	if (src[i][j] < min[i][j])
	  min[i][j] = src[i][j];
      }
}

 
/**
*@Description The function initializes the SuS selection mechanism.
*@param  (const unsigned long)
*@return   changes are made internally(void)
*/
void LoopTable::initOccurrence(const unsigned long nSelect){
  if (nAminoAcid == 1)
    return;
  
  nextIndex = rand() % size();

  stepWidth = 0;
}

 
/**
*@Description The function selects the next (random) table entry to be employed, 
*    using a SuS selection mechanism. A special case occurs if 
*    nAminoAcid = 1, see below.
*@param  none
*@return new table entry (LoopTableEntry)
*/
LoopTableEntry LoopTable::selectOccurrence(){
  PRECOND(nAminoAcid >= 1, exception);
  
  // if the chain is composed of a single amino acid the table entry 
  // needs to be rotated along a  random psi & phi  torsion angles to 
  // sample all solutions.
  if (nAminoAcid == 1)    { 
      INVARIANT(rama != NULL, exception);
      LoopTableEntry tmpEntry;
      tmpEntry.setToSingleAminoAcid();

      vgMatrix3<float> rotationMatrix = 
	  vgMatrix3<float>::createRotationMatrix( -psiNormal.normalize(), 
		  DEG2RAD * rama->getRandomPsi(1) );

      // NB: axis has to be -psiNormal in order to produce the correct angles,
      // +psiNormal would produce -angle, which is wrong.
      tmpEntry.endDirection = rotationMatrix * tmpEntry.endDirection;
      tmpEntry.endNormal = rotationMatrix * tmpEntry.endNormal;
      
      // rotate tmpEntry by randomPhi around (0,1,0): (update EP, ED, N)
      vgVector3<float> axis(0, 1,0);
      tmpEntry.rotate(axis, DEG2RAD * rama->getRandomPhi());

      return tmpEntry;
    }
  
  // now return the selected entry:
  INVARIANT(nextIndex < size(), exception);

  LoopTableEntry tmpEntry;

      tmpEntry = (*this)[nextIndex];
      nextIndex = rand() % size();


  return tmpEntry;
}

 
/**
*@Description The function codes an ieee32 value as the 1 byte wide 
*    offset between min and max.
*@param  value to consider ( const unsigned ),min and max value( const doubl,const doubl))
*@return the corresponding value (unsigned short) 
*/
unsigned short LoopTable::code(const double value, const double min, 
		const double max){
  PRECOND(((value >= min) && (value <= max)), exception);

  if ((max - min) == 0)
    return 0;
  double tmpResult = (value - min) * 65500 / (max - min);
  
  // round the resulting offset:
  if ( (tmpResult - static_cast<int>(tmpResult)) >= 0.5)
    return static_cast<unsigned short>(tmpResult) + 1;
  else
    return static_cast<unsigned short>(tmpResult);
}
 
/**
*@Description The three functions below decode a ProteinTableEntry, vector 
*		or ieee32 respectively into its uncompressed form. 
*@param value to consider ( const unsigned ),min and max value( const doubl,const doubl))
*@return the corresponding value (double)
*/
double LoopTable::decode(const unsigned short value, const double min, 
		  const double max){
  PRECOND((min <= max), exception);
  if ((max - min) == 0)
    return min;

  return ( min + (value * ((max - min) / 65500)) );
}


// TESTERS:
/**
*@Description print the table statistics (min, max):
*@param  index entry(unsigned int)
*/
void LoopTable::printTable(unsigned int num){
int tmpNEntry = size();

 
 cout << "nAminoAcid = " << static_cast<int>(nAminoAcid) 
     << "\n\n" << "Min:\n" << "\t EP: "<< setw(5) << setprecision(4) 
     << min.endPoint.x << "\t ED: " << setw(5) << setprecision(4) 
     << min.endDirection.x << "\t N: " << setw(5) << setprecision(4)
     << min.endNormal.x << "\t MP: " << setw(5) << setprecision(4)  
     << min.midPoint.x << "\t MD: " << setw(5) << setprecision(4)  
     << min.midDirection.x << "\t MN: " << setw(5) << setprecision(4) 
     << min.midNormal.x << "\n" << "\t EP: " << setw(5) << setprecision(4) 
     << min.endPoint.y << "\t ED: " << setw(5) << setprecision(4)  
     << min.endDirection.y << "\t N: " << setw(5) << setprecision(4)
     << min.endNormal.y << "\t MP: " << setw(5) << setprecision(4) 
     << min.midPoint.y << "\t MD: " << setw(5) << setprecision(4) 
     << min.midDirection.y << "\t MN: " << setw(5) << setprecision(4) 
     << min.midNormal.y << "\n" << "\t EP: " << setw(5) << setprecision(4)  
     << min.endPoint.z << "\t ED: " << setw(5) << setprecision(4) 
     << min.endDirection.z << "\t N: " << setw(5) << setprecision(4) 
     << min.endNormal.z << "\t MP: " << setw(5) << setprecision(4)  
     << min.midPoint.z << "\t MD: " << setw(5) << setprecision(4) 
     << min.midDirection.z << "\t MN: " << setw(5) << setprecision(4) 
     << min.midNormal.z << "\n"  << "Max:\n" << "\t EP: "<< setw(5) 
     << setprecision(4) << max.endPoint.x << "\t ED: " << setw(5) 
     << setprecision(4) << max.endDirection.x << "\t N: " 
     << setw(5) << setprecision(4)  << max.endNormal.x << "\t MP: " 
     << setw(5) << setprecision(4)  << max.midPoint.x << "\t MD: " 
     << setw(5) << setprecision(4)  << max.midDirection.x 
     << "\t MN: " << setw(5) << setprecision(4) << max.midNormal.x 
     << "\n" << "\t EP: " << setw(5) << setprecision(4)  << max.endPoint.y 
     << "\t ED: " << setw(5) << setprecision(4) << max.endDirection.y 
     << "\t N: " << setw(5) << setprecision(4)  << max.endNormal.y 
     << "\t MP: " << setw(5) << setprecision(4)  << max.midPoint.y 
     << "\t MD: " << setw(5) << setprecision(4)  << max.midDirection.y 
     << "\t MN: " << setw(5) << setprecision(4) << max.midNormal.y 
     << "\n" << "\t EP: " << setw(5) << setprecision(4) << max.endPoint.z
     << "\t ED: " << setw(5) << setprecision(4) << max.endDirection.z 
     << "\t N: " << setw(5) << setprecision(4)  << max.endNormal.z 
     << "\t MP: " << setw(5) << setprecision(4)  << max.midPoint.z 
     << "\t MD: " << setw(5) << setprecision(4)  << max.midDirection.z 
     << "\t MN: " << setw(5) << setprecision(4) << max.midNormal.z 
     << "\n";

cout << "----------------------------\n";

if (num > 0)
  tmpNEntry = num;

// print the first tmpNEntry table entries:
for (int i = 0; i < tmpNEntry; i++)
  (*this)[i].printTable(i);

 cout << "----------------------------\n";
}

