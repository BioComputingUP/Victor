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


#ifndef _PhyloSupport_h_
#define _PhyloSupport_h_

// Includes:

#include <string>
#include <vector>
#include <iostream>
#include <Alignment.h>

using namespace Victor::Align2;
using std::string;
using std::vector;


// Global constants, typedefs, etc. (to avoid):
namespace Victor { namespace Phylo {

    /**@brief  Methods to manages the global statistic data.
     * 
     *@Description  Use for create Tree during phylogeny and controll the parameter in NJ, UPGMA and CLUSTALW
     * */
    class PhyloSupport {
    public:

        // CONSTRUCTORS/DESTRUCTOR:
    	PhyloSupport();


        ~PhyloSupport();


        //STATIC:

        //Calc multi Align
        static vector<Alignment> calcAlignmentV(Alignment *aliSec, vector<vector<double> > &distance, bool ktuples=false, bool verbose=false);
        //Calc distance from 2 seq of DNA or protein
        static double distanceCalcTwoSeq(string seq1,string seq2);
        //Calc distance from 2 seq of DNA or protein use k-tuples method
        static double distanceCalcTwoSeqktuples(string seq1,string seq2);
        //for split string
        static std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
        static std::vector<std::string> split(const std::string &s, char delim);
        //for calculate divergent usefull for NJ
        static double calcDivR(vector<double> vDist);
        //for print matrix of double
        static void printMatrix(vector<vector<double> > &distance);
        //aling 2 sequence
        static vector<string> AlingSvsS(string seq1,string seq2,bool verbose=false);
        //aling 2 or more sequence use method of graph
        static vector<string> AlingMultiSvsMultiS(vector <string> seq1,vector <string> seq2,vector <double> vWeigth1,vector <double> vWeigth2,bool verbose=false,int tokenSize=1);
        //insert a gap in certain position in a string
        static string insertGapPosition(string seq, int position);
        //cast a int to string
        static string intToString( int num );
        //to merge 2 vector of string, order of merge is v=v1[0],v1[1]....v2[0]....v2[n-1]
        static vector<string> mergeStringVector(vector <string> v1,vector <string> v2 );
        //to merge 2 vector of double, order of merge is v=v1[0],v1[1]....v2[0]....v2[n-1]
        static vector<double> mergeDoubleVector(vector <double> v1,vector <double> v2 );

        //static parameter for NJ UPGMA and ClustalW:

        //gap for align
        static double openGapPenalty;
        static double extensionGapPenalty;
		//Cseq for align
        static double cSeq;
		//for align
        static double downs, downa, ups, upa;
		static unsigned int weightingScheme;
		//ClustalW
		static int tokenSize;
		//Substitution matrix
		static string matrix;
		//verbose output
		static bool verbose;



    protected:


    private:


    };


}}
#endif

