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
 

#ifndef __Profile_H__
#define __Profile_H__

#include <Alignment.h>
#include <AminoAcidCode.h>
#include <Debug.h>
#include <IoTools.h>
#include <string>
#include <vector>

namespace Biopool
{
/** @brief  Calculate a frequency profile or PSSM.
 * 
* @Description  
* @This 
 **/
class Profile
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	Profile();

	/// Copy constructor.
	Profile(const Profile &orig);

	/// Destructor.
	virtual ~Profile();


// OPERATORS:

	/// Assignment operator.
	Profile& operator = (const Profile &orig);


// PREDICATES:

	/// Return the frequency of the aminoacid amino for position i.
	virtual double getAminoFrequencyFromCode(AminoAcidCode amino,
		unsigned int i);

	/// Return the frequency of the aminoacid amino for position i.
	virtual double getAminoFrequency(char amino, unsigned int i);

	/// Return the frequency of the most frequent aminoacid for position i.
	virtual double getFreqMaxAminoFrequency(unsigned int i);

	/// Return the most frequent aminoacid for position i.
	virtual AminoAcidCode getAminoMaxFrequencyCode(unsigned int i);

	/// Return the most frequent aminoacid for position i.
	virtual char getAminoMaxFrequency(unsigned int i);

	/// Return the number of gaps for position i.
	virtual unsigned int getNumGap(unsigned int i);

	/// Return the number of sequences in the profile.
	virtual unsigned int getNumSequences();

	/// Return the lenght of sequences in the profile.
	virtual unsigned int getSequenceLength();

	/// Return the master sequence.
	virtual const string getSeq();

	/// Return the consensus of the profile.
	virtual string getConsensus();


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const Profile &orig);

	/// Construct a new "deep copy" of this object.
	virtual Profile* newCopy();

	/// Set the frequency of the aminoacid amino for position i.
	virtual void setFrequency(double freq, AminoAcidCode amino, int i);

	/// Set the frequency of the aminoacid amino for position i.
	virtual void setFrequency(double freq, char amino, int i);

	/// Set the number of gaps for position i.
	virtual void setNumGap(int numGap, int j);

	/// Set the number of sequences in the profile.
	virtual void setNumSequences(int i);

	/// Set the master sequence.
	virtual void setSeq(string master);

	/// Set the profile with or without gaps in the master sequence.
	virtual void setProfile(Alignment &ali);

	/// Set the profile with or without gaps in the master sequence.
	virtual void setProfile(Alignment &ali, istream &is);

	/// Set wether to include/exclude gaps in the master sequence.
	virtual void setAllowGaps(bool g);

	/// Reverse profile.
	virtual void reverse();


// HELPERS:

	/// Calculate the raw (ie. unnormalized) aminoacids frequencies for position i.
	virtual void pCalculateRawFrequency(vector<double> &freq, double &gapFreq,
		Alignment &ali, unsigned int i);

	/// Construct data from alignment.
	virtual void pConstructData(Alignment &ali);

	/// Reset all data.
	virtual void pResetData();


// ATTRIBUTES:

	vector< vector<double> > profAliFrequency;    ///< Aminoacids frequencies.
	vector<double> gapFreq;                       ///< Gaps frequencies.
	string seq;                                   ///< Master sequence.
	unsigned int seqLen;                          ///< Lenght of sequences.
	unsigned int numSeq;                          ///< Number of sequences.
	bool gap;                                     ///< If true, consider gaps in the master sequence.


protected:


private:

};

// -----------------------------------------------------------------------------
//                                   Profile
// -----------------------------------------------------------------------------

// PREDICATES:

inline double
Profile::getAminoFrequencyFromCode(AminoAcidCode amino, unsigned int i)
{
	return profAliFrequency[i][amino];
}


inline double
Profile::getAminoFrequency(char amino, unsigned int i)
{
	return getAminoFrequencyFromCode(aminoAcidOneLetterTranslator(amino), i);
}


inline double
Profile::getFreqMaxAminoFrequency(unsigned int i)
{
	return profAliFrequency[i][getAminoMaxFrequencyCode(i)];
}


inline AminoAcidCode
Profile::getAminoMaxFrequencyCode(unsigned int i)
{
	AminoAcidCode amino = XXX;
	double max = 0;

	for (AminoAcidCode j = ALA; j <= TYR; j++)
	{
		double tmp = profAliFrequency[i][j];
		if (tmp > max)
		{
			amino = j;
			max = tmp;
		}
	}

	return amino;
}


inline char
Profile::getAminoMaxFrequency(unsigned int i)
{
	return aminoAcidOneLetterTranslator(getAminoMaxFrequencyCode(i));
}


inline unsigned int
Profile::getNumGap(unsigned int i)
{
	return static_cast<int>(gapFreq[i]);
}


inline unsigned int
Profile::getNumSequences()
{
	return numSeq;
}


inline unsigned int
Profile::getSequenceLength()
{
	return profAliFrequency.size();
}


inline const string
Profile::getSeq()
{
	return seq;
}


inline string
Profile::getConsensus()
{
	string consensus;

	for (unsigned int i = 0; i < getSequenceLength(); i++)
	{
		char amino = getAminoMaxFrequency(i);
		consensus += amino;
	}

	return consensus;
}


// MODIFIERS:

inline void
Profile::setFrequency(double freq, AminoAcidCode amino, int i)
{
	profAliFrequency[i][amino] = freq;
}


inline void
Profile::setFrequency(double freq, char amino, int i)
{
	setFrequency(freq, aminoAcidOneLetterTranslator(amino), i);
}


inline void
Profile::setNumGap(int numGap, int j)
{
	gapFreq[j] = numGap;
}


inline void
Profile::setNumSequences(int i)
{
	numSeq = i;
}


inline void
Profile::setSeq(string master)
{
	seq = master;
}


inline void
Profile::setAllowGaps(bool g)
{
	gap = g;
}

} // namespace

#endif
