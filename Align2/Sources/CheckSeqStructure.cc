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
// --*- C++ -*------x---------------------------------------------------------
//
// Description:     This class is used to control the correspondence between
//                  a Fasta sequence file and a Pdb structure file. 
//                  It assign a progressive number to the aa of the sequence
//                  according to the aa of the structures.
// 
// -----------------x-------------------x-------------------x-----------------

#include <CheckSeqStructure.h>
#include <iostream>


namespace Biopool
{

// CONSTRUCTORS:

CheckSeqStructure::CheckSeqStructure(const Spacer& sp, const string& seq) : 
  sp(sp), seq(seq)
{
  setConfidence(15);
  pConstructArray();
}


CheckSeqStructure::CheckSeqStructure(const CheckSeqStructure& orig)
{
  copy(orig);
}


CheckSeqStructure::~CheckSeqStructure()
{ }


// OPERATORS:

CheckSeqStructure&
CheckSeqStructure::operator = (const CheckSeqStructure& orig)
{
  if (&orig != this)
    copy(orig);
  return *this;
}


// PREDICATES:

void
CheckSeqStructure::SetStartOffset(string& seq, Spacer& sp)
{
  if ( ( sp.sizeAmino() < conf) || ( seq.length() < conf ))
    ERROR("Sequence or structure too short!!!", exception);
  
  for ( unsigned int j = 0; j < (sp.sizeAmino()- conf); j++ )
    for ( unsigned int i = 0; i < (seq.length() - conf); i++ )
      {
	unsigned int tmp = i;
	unsigned int amino = j;
	
	while ( tmp <  (i + conf))
	  if ( seq[tmp] == sp.getAmino(amino).getType1L() )
	    {
	      tmp++;
	      amino++;
	      if ( tmp == (i + conf) )
		{
		  startoffset[0] = i;
		  startoffset[1] = j;
		  return;
		}
	    }
	  else 
	    break;
      }
}


void
CheckSeqStructure::SetEndOffset( string& seq, Spacer& sp )
{
  for ( unsigned int j = sp.sizeAmino()-1; j > conf; j-- )
    for (unsigned int i = seq.length()-1; i > conf; i--)
      {
	unsigned int tmp = i;
	unsigned int amino = j; 
	
	while ( tmp > (i - conf))
	  if ( (seq[tmp]) == (sp.getAmino(amino).getType1L()) )
	    {
	      tmp--;
	      amino--;
	      if ( (i-tmp) == conf )
		{
		  endoffset[0] = i;
		  endoffset[1] = j;
		  return;
		}
	      
	    }
	  else 
	    break;
      }
}


void
CheckSeqStructure::pConstructArray()
{
  // inizializzazione vettore
  for ( unsigned int i = 0; i < seq.length(); i++)
    pdbseqnumber.push_back(0);
  
  //riempimento vettore
  
  SetStartOffset(seq, sp); 
  SetEndOffset(seq, sp);
  
  /* caso1: la sequenza e' totalmente contenuta nella struttura */
  /* endoffset - startoffset == lunghezza sequenza */
  if ( (sp.getPdbNumberFromIndex(endoffset[1])) - (sp.getPdbNumberFromIndex(startoffset[1])) == static_cast<int>(seq.length()-1) )
    {
      {
	unsigned int j = 0;
	
	for ( unsigned int i = 0; i < seq.length(); i++ )
	  {
	    if (sp.isGap(sp.getPdbNumberFromIndex(startoffset[1]) + j ))
	      {
		pdbseqnumber[i] = -5;
		j++;
	      }
	    else
	      {
		pdbseqnumber[i] = sp.getPdbNumberFromIndex(startoffset[1]) + i;
		j++;
	      }
	  }
	
	//CONTROLLO
	unsigned int z = 0;
	unsigned int zz = 0;
	while ( z < pdbseqnumber.size() )
	  {
	    if ( pdbseqnumber[z] < 0 )
	      { 
		z++; 
	      }
	    else
	      { 
		if ( (seq[z]) == (sp.getAmino(startoffset[1] + zz).getType1L()) )
		  {
		    zz++;
		    z++;
		  }
		else
		  {
		    for ( unsigned int i = 0; i < pdbseqnumber.size(); i++)
		      pdbseqnumber[i] = 0;
		    pSpecialCase();
		    break;
		  }
	      }
	  }
	
      }
      
    }
  
  /* caso 2: l'inizio( o la fine) della sequenza manca nella struttura e la fine( o l'inizio) della sequenza e' nella struttura.
     viene rifatto il ciclo per controllare eventuali mancanze di amminoacidi interni della struttura rispetto alla sequenza. */
  /* endoffset - startoffset < lunghezza sequenza */
  
  if ( ( sp.getPdbNumberFromIndex(endoffset[1])) - (sp.getPdbNumberFromIndex(startoffset[1]) ) < static_cast<int>( seq.length() ) )
    {
      int count = 0;
      int count2 = startoffset[1];
      if (( endoffset[0] < static_cast<int>( (seq.length()-1)) ) && (startoffset[0] == 0))
	{	
	  for (int i = 0; i <= endoffset[0]; i++)
	    {
	      
	      if (sp.isGap(sp.getPdbNumberFromIndex(startoffset[1]) + count))
		{
		  
		  pdbseqnumber[i] = -5;
		  count++;
		  count2++;
		}
	      else
		{
		  pdbseqnumber[i] = sp.getPdbNumberFromIndex(startoffset[1])+i;
		  count2++;
		  count++;
		}
	    }
	  
	  for ( unsigned int i =  endoffset[0]+1; i < seq.length(); i++)
	    pdbseqnumber[i] = (-1);
	  
	}
      
      
      else if (( startoffset[0] > 0 ) && ( endoffset[0] == static_cast<int>(seq.length()-1) ))
	{
	  int count = 0;
	  int count2= startoffset[1];
	  {
	    for ( int i = startoffset[0]; i < static_cast<int>(seq.length()); i++)
	      {
		
		if ( sp.isGap(sp.getPdbNumberFromIndex(startoffset[1]) + count))
		  {
		    pdbseqnumber[i] = -5;
		    count++;
		    count2++;
		  }
		else
		  {
		    pdbseqnumber[i] = sp.getPdbNumberFromIndex(startoffset[1])+i - startoffset[0];
		    count2++;
		    count++;
		  }
		
	      }
	    for ( int i = 0; i <  startoffset[0]; i++)
	      pdbseqnumber[i] = (-1);
	    
	  }
	}
      
      else if (( startoffset[0] > 0 ) && ( endoffset[0] < static_cast<int>((seq.length()-1)) ))
	{
	  {
	    int count = 0;
	    int count2= startoffset[1];
	    for ( int i = startoffset[0]; i <= endoffset[0]; i++ )
	      {
		if ( sp.isGap(sp.getPdbNumberFromIndex(startoffset[1]) + count))
		  {
		    pdbseqnumber[i] = -5;
		    count2++;
		    count++;
		  }
		else
		  {
		    pdbseqnumber[i] =sp.getPdbNumberFromIndex(startoffset[1])+i - startoffset[0];
		    count2++;
		    count++;
		  }
		
	      }
	    
	    for ( int i = 0; i < startoffset[0]; i++)
	      pdbseqnumber[i] = (-1);
	    for ( unsigned int i =  endoffset[0]+1; i < seq.length(); i++)
	      pdbseqnumber[i] = (-1);
	    
	  }
	  
	  //CONTROLLO:
	  unsigned int z = 0;
	  unsigned int zz = 0;
	  while ( z < pdbseqnumber.size() )
	    {
	      if ( pdbseqnumber[z] < 0 )
		{ 
		  z++; 
		}
	      else
		{ 
		  if ( seq[z] == sp.getAmino(startoffset[1] + zz).getType1L() )
		    {
		      zz++;
		      z++;
		    }
		  else
		    {
		      for ( unsigned int i = 0; i < pdbseqnumber.size(); i++)
			pdbseqnumber[i] = 0;
		      pSpecialCase();
		      break;
		    }
		}
	    }
	}
    }
  
  /* caso 3 la struttura e' molto spezzettata o viene passata la sequenza di un dominio contro una struttura */
  
  if ( ( sp.getPdbNumberFromIndex(endoffset[1])) - (sp.getPdbNumberFromIndex(startoffset[1]) ) > static_cast<int>(seq.length()) )
    pSpecialCase();
}


void   	 
CheckSeqStructure::pSpecialCase() 
{
  unsigned int j = 0;
  unsigned int i = 0;
  string seq2 = seq;
  unsigned int s = startoffset[0];
  unsigned int count = startoffset[0];

  while ( i < seq2.length() )
    {
      if (sp.isGap(sp.getPdbNumberFromIndex(startoffset[1]) + j ))
	{
	  pdbseqnumber[s] = -5;
	  j++;
	  s++;
	  i++;
	}
      else
	{
	  if ( seq2[count + i] == sp.getAmino(startoffset[1] + j).getType1L() )
	    {
	      pdbseqnumber[s] = sp.getPdbNumberFromIndex(startoffset[1]) + i;
	      j++;
	      i++;
	      s++;
	      
	    }
	  else
	    {
	      string seq3;
	      unsigned int tmp = count + i;
	      if ( ( tmp < 10 ) && ( tmp >= 3 ) )
		setConfidence(tmp);
	      else if ( tmp >= 10 )
		setConfidence(10);
	      else if ( tmp < 3 )
		return;
	      while ( tmp < seq2.length() )
		{
		  seq3 = seq3 + seq2[tmp];
		  tmp++;
		}
	      seq2.erase();
	      seq2 = seq3;
	      SetStartOffset(seq2, sp); 
	      j = 0;
	      i = 0;
	      s = s + startoffset[0];
	      count = startoffset[0];
	      seq3.erase();
	      
	    }
	}
    }
  for ( unsigned int ii = 0; ii < pdbseqnumber.size(); ii++ )
    {
      if (pdbseqnumber[ii] == 0)
	{
	  pdbseqnumber[ii] = -5;
	}
    }
       
}


void 
CheckSeqStructure::copy(const CheckSeqStructure& orig)
{
  sp = orig.sp;  
  seq = orig.seq; 
  endoffset[2] = orig.endoffset[2];
  startoffset[2] = orig.startoffset[2];
  for ( unsigned int i = 0; i < seq.length(); i++)
    pdbseqnumber.push_back(orig.pdbseqnumber[i]);
  conf = orig.conf;  
}

} // namespace
