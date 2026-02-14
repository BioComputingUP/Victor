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


#ifndef _SeqNodeGraph_h_
#define _SeqNodeGraph_h_

// Includes:

#include <string>
#include <vector>
#include <iostream>
#include <Alignment.h>


using std::string;
using std::vector;



// Global constants, typedefs, etc. (to avoid):
namespace Victor { namespace Phylo {

    /**@brief  Class for create node of a graph of sequences
     * 
     *@Description  Use for create ClustalW, for create a graph of token
     * */
    class SeqNodeGraph{
    public:

        // CONSTRUCTORS/DESTRUCTOR:
    	SeqNodeGraph();
    	SeqNodeGraph(unsigned int indexS,unsigned int indexF,unsigned int totNumS,vector <string> tokenS,unsigned int numSeqInS2);

        ~SeqNodeGraph();


    	// PREDICATES:
        int getIndexStart();
        int getIndexFinish();
        unsigned int getTotNumSeq();
        string getTokenSeq(unsigned int index);
        double getTaxEdgeInPosition(unsigned int i);
        char getCharOfTokenSeq(unsigned int seqNum,unsigned int pos);
        unsigned int getTokenSize();
        double getAverageTax();
        unsigned int returnBestEdge();
        int returnBestEdgeAfterIndex(unsigned int index);

	    // OPERATORS:
        void setIndexStart(unsigned int i);
        void setIndexFinish(unsigned int i);
        void setTotNumSeq(int i);
        void setTaxEdgeInPosition(unsigned int i, double tax);
        void setTokenSize(int size);
        void printTaxEdge();
        void calculateAverageTax();
        static void setNode(SeqNodeGraph* node, vector <SeqNodeGraph*> vNode);
        static void setNodeWithWeigth(SeqNodeGraph* node, vector <SeqNodeGraph*> vNode,vector <double> vWeigth1,vector <double> vWeigth2);




    protected:
        // ATTRIBUTES:
       unsigned int indexStart;
       unsigned int indexFinish;
       unsigned int totNumSeq;
       unsigned int numTokenInSeq2;
       unsigned int numSeq;
       vector<double> taxEdge;
       vector<string> tokenSeq;
       unsigned int tokenSize;
       double averageTax;


    private:



    };


}}
#endif

