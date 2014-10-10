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


#ifndef _LOOPTABLE_H_
#define _LOOPTABLE_H_

// Includes:
#include <string>
#include <vector3.h>
#include <matrix3.h>
#include <Debug.h>
#include <RamachandranData.h>
#include <LoopTableEntry.h>
#include <queue>

namespace Victor { namespace Lobo {

    // Global constants, typedefs, etc. (to avoid):

    /**@brief  Struct that contains a queue element.
     * 
     * */
    struct solutionQueueElem {
        double dev;
        unsigned int index1, index2;

        bool operator<(const solutionQueueElem& _sqe) const {
            return (_sqe.dev < dev);
        }
    };

    /**
     * @brief  Defines a table of possible amino chain end points and end directions 
     *  after k amino acids have been concatenated. 
     * 
     * */
    class LoopTable {
    public:

        // CONSTRUCTORS/DESTRUCTOR:
        LoopTable();
        LoopTable(const LoopTable& orig);
        virtual ~LoopTable();

        // PREDICATES:
        LoopTableEntry getMin();
        LoopTableEntry getMax();
        RamachandranData* getRama();
        virtual unsigned int size();
        unsigned int getLength();
        unsigned int getMaxBins();
        void showDistribution();
        virtual LoopTableEntry getClosest(const LoopTableEntry& le,
                unsigned int currentSelection = 1);
        virtual vector<LoopTableEntry> getNClosest(const LoopTableEntry& le,
                unsigned int num, unsigned int nAmino);

        // MODIFIERS:
        void setRama(RamachandranData* r);
        void setLength(unsigned int l);

        void copy(const LoopTable& orig);

        // initializes the table for nAminoAcid = 1
        void setToSingleAminoAcid();  

         // concatenates two tables
        void concatenate(LoopTable&, LoopTable&, unsigned long, unsigned long);
       

        virtual void adjustTable();
        virtual void read(const string&);
        virtual void cluster(double cutoff);
        virtual void write(const string&);
        void writeASCII(const string&, unsigned long num = 0,
                unsigned int wEntry = 0, unsigned int wDim = 0);

        // OPERATORS:
        LoopTable& operator=(const LoopTable& orig);
        virtual LoopTableEntry& operator[](unsigned int n);
        virtual const LoopTableEntry& operator[](unsigned int n) const;


        // TESTERS:
        void printTable(unsigned int);

        static unsigned int MAX_FACTOR;

    protected:

        // HELPERS: 
        double pEstimateSimilarityCutoff(const LoopTableEntry& le,
                unsigned int num);

        virtual void store(const LoopTableEntry&);
        void pInsertElem(LoopTableEntry& elem,
                vector<vector<LoopTableEntry> >& table);
        unsigned int pGetBin(const vgVector3<float>& e);
        void pAddToSolutionQueue(unsigned int offset,
                priority_queue<solutionQueueElem>& solutionQueue,
                const LoopTableEntry& dest);

        //methods used to concatenate two tables:
        void initOccurrence(const unsigned long);
        LoopTableEntry selectOccurrence();

        //methods used to compress and uncompress a table:
        unsigned short code(const double, const double, const double);
        double decode(const unsigned short, const double, const double);


    private:

        // ATTRIBUTES:
        RamachandranData* rama;

        // number of amino acids represented by this table, ie. chain segment length
        unsigned int nAminoAcid;
        

        // table of (sorted) bins containing the chain information
        vector<vector<LoopTableEntry> > entry;
        

        double lowerLimit; // information about the lower and 
        double upperLimit; // upper bin limits
        double stepLimit; // bin step size

        LoopTableEntry max, min; // min and max values for each coordinate

        // constant expression required for selectOccurrence:
        vgVector3<float> psiNormal;

        unsigned long nextIndex; // data required for SuS selection of 
        unsigned long stepWidth; // next occurrence

        static unsigned int MAX_BINS;
        static double BOND_ANGLE_AT_CPRIME_TO_N;
        static double BOND_LENGTH_N_TO_CALPHA;
        static double BOND_ANGLE_AT_N_TO_CALPHA;

    };

    void printTable(LoopTable&, int num = -1);

    // ---------------------------------------------------------------------------
    //                                  LoopTable
    // -----------------x-------------------x-------------------x-----------------


    // PREDICATES:

    /**
     *@Description returns the minimum value for each coordinate
     *@param  none
     *@return  corresponding loop entry table( LoopTableEntry)
     */
    inline LoopTableEntry LoopTable::getMin() {
        return min;
    }

    /**
     *@Description returns the maximum value for each coordinate
     *@param  none
     *@return  corresponding loop entry table( LoopTableEntry)
     */
    inline LoopTableEntry LoopTable::getMax() {
        return max;
    }

    /**
     *@Description return the   Ramachandran Data
     *@param  none
     *@return  pointer to the data( RamachandranData* )
     */
    inline RamachandranData* LoopTable::getRama() {
        return rama;
    }

    /**
     *@Description return the entries total size
     *@param  none
     *@return  corresponding value ( unsigned int)
     */
    inline unsigned int LoopTable::size() {
        unsigned int count = 0;
        for (unsigned int i = 0; i < entry.size(); i++)
            count += entry[i].size();
        return count;
    }

    /**
     *@Description return the amino acids number 
     *@param  none
     *@return  corresponding value ( unsigned int)
     */
    inline unsigned int LoopTable::getLength() {
        return nAminoAcid;
    }

    /**
     *@Description return the maximum bindings
     *@param  none
     *@return  corresponding value ( unsigned int)
     */
    inline unsigned int LoopTable::getMaxBins() {
        return MAX_BINS;
    }

    // MODIFIERS:

    /**
     *@Description Sets the Ramachandran Data
     *@param  pointer to the RamachandranData(RamachandranData*)
     *@return  changes are made internally(void)
     */
    inline void LoopTable::setRama(RamachandranData* r) {
        rama = r;
    }

    /**
     *@Description Sets the length
     *@param  length value(unsigned int)
     *@return   changes are made internally(void)
     */
    inline void LoopTable::setLength(unsigned int l) {
        nAminoAcid = l;
    }


    // OPERATORS:
    /**
     *@Description returns the corresponding table entry for the given index
     *@param  index(unsigned int )
     *@return  reference to the loop entry table( LoopTableEntry&)
     */
    inline LoopTableEntry& LoopTable::operator[](unsigned int n) {

        unsigned int count = MAX_BINS + 1;
        for (unsigned int i = 0; i < entry.size(); i++) {
            if (n >= entry[i].size())
                n -= entry[i].size();
            else {
                count = i;
                break;
            }
        }

        if ((count >= MAX_BINS + 1) || (entry.size() == 0))
            ERROR("LoopTable::operator[] : Argument out of scope.", exception);

        return entry[count][n];
    }
    /**
     *@Description returns the corresponding table entry for the given index
     *@param  index(unsigned int )
     *@return  reference to the loop entry table( LoopTableEntry&)
     */
    inline const LoopTableEntry& LoopTable::operator[](unsigned int n) const {
        unsigned int count = MAX_BINS + 1;
        for (unsigned int i = 0; i < entry.size(); i++) {
            if (n >= entry[i].size())
                n -= entry[i].size();
            else {
                count = i;
                break;
            }
        }
        if ((count >= MAX_BINS + 1) || (entry.size() == 0))
            ERROR("LoopTable::operator[] : Argument out of scope.", exception);
        return entry[count][n];
    }


    // HELPERS:

    /**
     *@Description Inserts a single element in table.
     *@param  referecnce to the Loop Entry table(LoopTableEntry&), reference to the loop table entry table ( vector<vector<LoopTableEntry>&) 
     *@return   changes are made internally(void)
     */
    inline void LoopTable::pInsertElem(LoopTableEntry& elem, vector<vector<LoopTableEntry> >& table) {
        unsigned int numBin = pGetBin(elem.endPoint);
        if (numBin >= MAX_BINS)
            ERROR("LoopTable::pInsertElem : element out of scope.", exception);

        table[numBin].push_back(elem);
    }

    /**
     *@Description returns the corresponding binding
     *@param  reference for the coordenate values (const vgVector3<float>&)
     *@return  corresponding binding value( unsigned int)
     */
    inline unsigned int LoopTable::pGetBin(const vgVector3<float>& e) {
        double dist = sqrt(sqr(e.x) + sqr(e.y) + sqr(e.z));
        if (dist < lowerLimit)
            return 0;
        if (dist >= upperLimit)
            return MAX_BINS - 1;

        return static_cast<unsigned int> ((dist - lowerLimit) / stepLimit);
    }

    /**
     *@Description Adds the data into the solution queue
     *@param  offset(unsigned int),reference to the solutions 
     * queue(priority_queue<solutionQueueElem>&), reference to the loop entry table(const LoopTableEntry&)
     *@return   changes are made internally(void)
     */
    inline void LoopTable::pAddToSolutionQueue(unsigned int offset,
            priority_queue<solutionQueueElem>& solutionQueue, const LoopTableEntry& dest) {
        solutionQueueElem sqe;
        for (unsigned int i = 0; i < entry[offset].size(); i++) {
            sqe.dev = entry[offset][i].calculateDeviation(dest, nAminoAcid);
            sqe.index1 = offset;
            sqe.index2 = i;
            solutionQueue.push(sqe);
        }
    }

}} // namespace

#endif //_LOOPTABLE_H_


