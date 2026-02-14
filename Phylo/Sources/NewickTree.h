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


#ifndef _NewickTree_h_
#define _NewickTree_h_

// Includes:

#include <string>
#include <vector>
#include <iostream>
#include <Alignment.h>


using std::string;
using std::vector;

//rappresent the node of a phylogenic tree
struct iNode
{
	string name;
	int nseq;
	string seq;//pure seq
	vector<string> allignSeq;
	bool isRigth;
	bool isLeft;
	iNode *left;
	iNode *right;
	iNode *parent;
	int numberOfChildLeaf;
	bool isLeaf;
	double branchLength;
	double divergenceR;//for NJ
	double weigth;//for ClustalW
	vector<double> weigthV;//for ClustalW
	vector<string> nameV;//for ClustalW
	bool ClustalW;//for ClustalW

};

// Global constants, typedefs, etc. (to avoid):
namespace Victor { namespace Phylo {

    /**@brief  Methods to manages the global statistic data.
     * 
     *@Description  Use for create Tree during phylogeny
     * */
    class NewickTree {
    public:

        // CONSTRUCTORS/DESTRUCTOR:
    	NewickTree();
		NewickTree(int position,string name,string pureSeq,double divR);
		NewickTree(NewickTree *rTree, NewickTree *lTree);
		NewickTree(int position,string name,string pureSeq);

        ~NewickTree();


    	// PREDICATES:
    	iNode *getRoot();
    	iNode *getRightChild();
    	iNode*getLeftChild();
    	iNode*getParent();

    	string getName();
        string getPureSeq();

        double getValueOFChildBranchLength();
        double getRdiv();
        double getWeigth();

        bool isLeaf();

        int getNumOfChildren();

        unsigned int getNumberOfLeaf();

        NewickTree getLeafInPosition(unsigned int index);




	    // OPERATORS:
        void setBranchLength(double val);
	    void neighborJoining(Align2::Alignment ali,bool ktuples=false, bool verbose=false);
	    void upgma(Align2::Alignment ali, bool ktuples=false, bool verbose=false);
	    void setRdiv(double val);
	    void setParent(NewickTree pTree);
	    void setWeigth(double w);
	    void setPosition(unsigned int pos);
	    void setAlingSeq(vector<string> vSeq);


	    // Print in a string the tree using Newick format
        string printNewickTree();

    protected:
        // HELPERS:
        string printNewickTreeNode(iNode node);
        double calcWeigth(iNode* seq);

    // ATTRIBUTES:
        iNode *root;
        vector<NewickTree> leafList;


    private:
        // HELPERS:
        void destroy_tree(iNode *leaf);


    };


}}
#endif

