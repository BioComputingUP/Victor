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


// Includes:
#include <SubMatrix.h>
#include <Alignment.h>
#include <AlignmentData.h>
#include <SequenceData.h>
#include <NWAlign.h>
#include <AGPFunction.h>
#include <ScoringS2S.h>

#include <NewickTree.h>
#include <PhyloSupport.h>

#include <sstream>
#include <string>
#include <vector>

// Global constants, typedefs, etc. (to avoid):
using namespace Victor::Align2;
using namespace Victor;
using namespace std;

namespace Victor { namespace Phylo{

	// CONSTRUCTORS:


std::vector<double> distance;

/**
 *@Description for create void node
 */
	NewickTree::NewickTree() {
		//root->weigthV;
		//root->nameV;
		//root->allignSeq;
		root=new iNode;
		root->name="NULL";
		root->left=NULL;
		root->right=NULL;
		root->parent=NULL;
		root->seq="noPURESeq";
		root->branchLength=-1;
		root->divergenceR=0;
		root->nseq=-1;
		root->ClustalW=false;
		root->isRigth=false;
		root->isLeft=false;

		root->isLeaf=true;
		root->numberOfChildLeaf=0;
	}

	/**
	 *@Description for create node of tree
	 *@param int position numeber of sequenze in FASTA file
	 *@param string name name of sequence in FASTA file
	 *@param string pureSeq the  Sequence with no modification like in FASTA file
	 *@param double divR divergence usefull inNJ method
	 */
	NewickTree::NewickTree(int position,string name,string pureSeq,double divR){
		//not need
		//root->allignSeq;
		//root->nameV;
		//root->weigthV;
		//void
		root=new iNode;
		root->left=NULL;
		root->right=NULL;
		root->parent=NULL;
		root->branchLength=-1;
		root->ClustalW=false;
		root->isRigth=false;
		root->isLeft=false;
		//important
		root->name=name;
		root->divergenceR=divR;
		root->seq=pureSeq;
		root->nseq=position;
		root->isLeaf=true;
		root->numberOfChildLeaf=0;
	}

	/**
	 *@Description for create node of tree
	 *@param int position numeber of sequenze in FASTA file
	 *@param string name name of sequence in FASTA file
	 *@param string pureSeq the  Sequence with no modification like in FASTA file
	 */
	NewickTree::NewickTree(int position,string name,string pureSeq){
		//void
		//root->weigthV;
		//root->nameV;
		//root->allignSeq;
		root=new iNode;
		root->left=NULL;
		root->right=NULL;
		root->parent=NULL;
		root->branchLength=0;
		root->ClustalW=false;
		root->isRigth=false;
		root->isLeft=false;
		//important
		root->name=name;
		root->divergenceR=0;
		root->seq=pureSeq;
		root->nseq=position;
		root->isLeaf=true;
		root->numberOfChildLeaf=0;
	}



	/**
	 *@Description Create new Tree. Used for union two different Tree.
	 *@param NewickTree *rTree the right child
	 *@param NewickTree *lTree the left child
	 */
	NewickTree::NewickTree( NewickTree *rTree, NewickTree *lTree) {
		//not need
		//root->weigthV;
		//root->nameV;
		//root->allignSeq;
		//void
		root=new iNode;
		root->branchLength=-1;
		root->ClustalW=false;
		root->isRigth=false;
		root->isLeft=false;
		//important
		root->isLeaf=false;//internal
		int  tmp=0;
		if(rTree->isLeaf())
			tmp++;
		if(lTree->isLeaf())
			tmp++;
		root->numberOfChildLeaf=tmp+rTree->getNumOfChildren()+lTree->getNumOfChildren();
		root->name="";
		root->seq="";

		rTree->getRoot()->isRigth=true;
		root->right=rTree->getRoot();

		lTree->getRoot()->isLeft=true;
		root->left=lTree->getRoot();

		root->parent=NULL;

		root->divergenceR=0;
	}


	/**
	 *@Description Basic destructor
	 */
	NewickTree::~NewickTree() {
	}


	// PREDICATES:
	iNode *NewickTree::getRoot(){
		return root;
	}

	iNode *NewickTree::getRightChild(){
		return root->right;
    }

	iNode*NewickTree::getLeftChild(){
		return root->left;
    }
	iNode*NewickTree::getParent(){
		return root->parent;
    }

    string NewickTree::getName(){
    	return root->name;
    }

    string NewickTree::getPureSeq(){
    	return root->seq;
    }


    double NewickTree::getRdiv(){
    	return root->divergenceR;
    }


    bool NewickTree::isLeaf(){
    	return root->isLeaf;
    }

    int NewickTree::getNumOfChildren(){
    	return root->numberOfChildLeaf;
    }

    double NewickTree::getWeigth(){
 	   return root->weigth;
    }

    NewickTree NewickTree::getLeafInPosition(unsigned int index){
    	return leafList[index];
    }

    unsigned int NewickTree::getNumberOfLeaf(){
    	return leafList.size();
    }

    // OPERATORS:
	/**
	 *@Description for create phylogenic Tree use NJ method
	 *@param Align2::Alignment ali alignament of starting sequence
	 *@param bool ktuples decide to use ktuples method for calculate identity of 2 sequence
	 *@param bool verbose for more information on video during computing.
	 */
    void NewickTree::neighborJoining(Align2::Alignment ali,bool ktuples, bool verbose){

		int Nseq=ali.size()+1;//number of seq
		vector<vector<double> > distance(Nseq, std::vector<double>(Nseq, 0));

		if(verbose)
			cout<<"---------Creating Distance Matrix--------------/"<<endl;
		vector<Alignment>vet=PhyloSupport::calcAlignmentV(&ali,distance,ktuples,verbose);



		if (verbose){
			cout<<"---------Print Distance Matrix--------------/"<<endl;
			PhyloSupport::printMatrix(distance);
		}

		vector<NewickTree> trees(Nseq);

		//Create all node
		for(unsigned int i=0;i<trees.size();i++){
			if(i==trees.size()-1){
				trees[i] = NewickTree(i,PhyloSupport::split(ali.getTargetName(), ' ')[0],Alignment::getPureSequence(ali.getTarget()),PhyloSupport::calcDivR(distance[distance.size()-1]));
			}
			else
				trees[i] = NewickTree(i,PhyloSupport::split(ali.getTemplateName(i), ' ')[0],Alignment::getPureSequence(ali.getTemplate(i)),PhyloSupport::calcDivR(distance[i]));

		}


		vector<NewickTree> seqList(trees.size());//never modify, used after create guide tree
		int count=0;//count number of seq insert in seqList


		//memorize the nearest seq
		unsigned int minj=0;
		unsigned int mini=1;
		//vector for new distance
		vector<double> newMediaD ( trees.size());

		//for creat new matrix;
		unsigned int newi=0;
		unsigned int newj=0;


		while(trees.size()!=3){// end 3
			if(verbose)
				cout<<"----------------While num= "<<(Nseq-trees.size())<<"---------------"<<endl;

			//calculate a new distanceMatrix using for each pair of seq
			//Example: M(ij)=d(ij) -([r(i) + r(j)]/N-2) or in the case of the pair A,B:
			//Example: M(AB)=d(AB) -([(r(A) + r(B)]/N-2)
			// r=divergenceR, d(AB)=distanceA-B, M(AB)=newDistnace, N=total node in this moment
			//matrix Mij for calculate minimun distance
			vector<vector<double> > Mij(trees.size(), std::vector<double>(trees.size(), 0));

			for(unsigned int i=0; i<trees.size(); i++){
				for(unsigned int j=0; j<trees.size(); j++){
					if(i!=j){
						Mij[i][j]=distance[i][j] - (  (trees[i].getRdiv()+trees[j].getRdiv())/(trees.size()-2)  );
						if(verbose)
							cout<<"MIJ "<<i<<" "<<j<<" = "<<distance[i][j]<<" - ( "<<trees[i].getRdiv()<<" + "<<trees[j].getRdiv()<<" ) = "<<Mij[i][j]<<endl;
					}
					else
						Mij[i][j]=0;//diagonal useless
				}
			}

			if(verbose){
				cout<<"----------------------------MIJ MATRIX--------------------------------"<<endl;
				PhyloSupport::printMatrix(Mij);
				cout<<endl<<"----------------------------find Min RDiv--------------------------------"<<endl;
			}

			//reset
			mini=0;
			minj=1;

			//find min i e min j
			for(unsigned int i=0; i<trees.size(); i++){
				for(unsigned int j=0; j<trees.size(); j++){
					if( Mij[mini][minj]  > Mij[i][j] && i!=j)
					{
						minj=j;
						mini=i;
					}
				}
			}

			if(verbose)
				cout<<endl<<"----------------------------min i j find! "<<mini<<" "<<minj<<" --------------------------------"<<endl;


			//distance of new node U respect A B the children
			//S(AU) =(d(AB)/2) + ([r(A)-r(B)]/2(n-2) ))
			//S(BU) =(d(AB)/2) + ([r(B)-r(A)]/2(n-2) ))
			trees[mini].setBranchLength((distance[mini][minj]/2) +  ((trees[mini].getRdiv()-trees[minj].getRdiv() )/(2*(trees.size()-2))  ) );
			trees[minj].setBranchLength((distance[mini][minj]/2) +  ((trees[minj].getRdiv()-trees[mini].getRdiv() )/(2*(trees.size()-2))  ) );

			vector<vector<double> > tmpdistance(trees.size()-1, std::vector<double>(trees.size()-1, 0));

			//reset
			newi=0;
			newj=0;
			//copy old matrix in new matrix, but not mini (A) and minj (B) and not the new node (U)
			for(unsigned int i=0; i<distance.size();i++){
				for(unsigned int j=0; j<distance.size();j++){
					if(i!=mini && j!=minj && i!=minj && j!=mini){
						if(i<mini && i<minj){
							newi=i;
						}
						else if(i>mini && i<minj){
							newi=i-1;
						}
						else if(i>mini && i>minj){
							newi=i-2;
						}
						if(j<mini && j<minj){
							newj=j;
						}
						else if(j>mini && j<minj){
							newj=j-1;
						}
						else if(j>mini && j>minj){
							newj=j-2;
						}
						if(verbose)
							cout<<endl<<"mini e minj "<<mini<<" "<<minj<<" i e j "<<i<<" "<<j<<" newi e new j "<<newi<<" "<<newj<<" "<<"distance value  "<<distance[i][j]<<endl;
						tmpdistance[newi][newj]=distance[i][j];
					}
				}//end for j
			}//end for i

			//insert the new distance
			/*	example:
				d(CU) = (d(AC) + d(BC) - d(AB) )/2
				d(DU) = (d(AD) + d(BD) - d(AB) )/2
				d(EU) = (d(AE) + d(BE) - d(AB) )/2
				d(FU) = (d(AF) + d(BF) - d(AB) )/2
			*/
			newi=0;
			for(unsigned int i=0; i<distance.size();i++){
				if(i!=mini && i!=minj){
					if(i<mini && i<minj){
						newi=i;
					}
					else if(i>mini && i<minj){
						newi=i-1;
					}
					else if(i>mini && i>minj){
						newi=i-2;
					}
				tmpdistance[newi][tmpdistance.size()-1]=(distance[mini][i]+distance[minj][i]-distance[mini][minj])/2;
				if(verbose)
					cout<<endl<<"("<<distance[mini][i]<<" + "<<distance[minj][i]<<" - "<<distance[mini][minj]<<")/2 = "<<tmpdistance[newi][tmpdistance.size()-1]<<endl;
				tmpdistance[tmpdistance.size()-1][newi]=tmpdistance[newi][tmpdistance.size()-1];
				}
			}

			//ccreate ne distance matrix dimension-1
			distance=tmpdistance;
			if(verbose)
				PhyloSupport::printMatrix(distance);

			//insert new tree
			trees.push_back(NewickTree(&trees[mini],&trees[minj]));

			//insert parent
			trees[mini].setParent(trees[trees.size()-1]);
			trees[minj].setParent(trees[trees.size()-1]);

			//save in array if leaf
			if(trees[mini].isLeaf()){
				seqList[count]=trees[mini];
				trees[mini].setPosition(count);
				count++;
			}
			if(trees[minj].isLeaf()){
				seqList[count]=trees[minj];
				trees[minj].setPosition(count);
				count++;
			}

			//erase trees
			trees.erase(trees.begin()+mini);
			//minus 1, because  length was decreased
			trees.erase(trees.begin()+minj-1);
			if(verbose)
				cout<<"the distance of tree vector after erase"<<trees.size()<<endl;

			//SET NEW r
			for(unsigned int i=0; i<trees.size();i++)
				trees[i].setRdiv(PhyloSupport::calcDivR(distance[i]));


		}//end while
		NewickTree *tmpT;

		//The root is placed by a "mid-point" method (15) at a position
		//where the means of the branch lengths on either side of the root
		//are equal
		//now use upgma for decide root
		if(distance[0][1]<distance[1][2] || distance[0][2]<distance[1][2]){
			if(distance[0][1]<distance[0][2]){//min distance 0,1
				trees[0].setBranchLength(distance[0][1]/2-trees[0].getValueOFChildBranchLength());
				trees[1].setBranchLength(distance[0][1]/2-trees[1].getValueOFChildBranchLength());
				tmpT=new NewickTree(&trees[0],&trees[1]);
				//parent
				trees[0].setParent(*tmpT);
				trees[1].setParent(*tmpT);

				trees[2].setBranchLength((distance[1][2]/2+distance[0][2]/2)/2-trees[2].getValueOFChildBranchLength());
				tmpT->setBranchLength((distance[1][2]/2+distance[0][2]/2)/2-tmpT->getValueOFChildBranchLength());
				tmpT=new NewickTree(tmpT,&trees[2]);

				//parent
				trees[2].setParent(*tmpT);

			}
			else{//min distance 0,2
				trees[0].setBranchLength(distance[0][2]/2-trees[0].getValueOFChildBranchLength());
				trees[2].setBranchLength(distance[0][2]/2-trees[2].getValueOFChildBranchLength());
				tmpT=new NewickTree(&trees[0],&trees[2]);
				//parent
				trees[0].setParent(*tmpT);
				trees[2].setParent(*tmpT);

				trees[1].setBranchLength((distance[1][2]/2+distance[1][0]/2)/2-trees[1].getValueOFChildBranchLength());
				tmpT->setBranchLength((distance[1][2]/2+distance[1][0]/2)/2-tmpT->getValueOFChildBranchLength());
				tmpT=new NewickTree(tmpT,&trees[1]);
				//parent
				trees[1].setParent(*tmpT);

			}

		}
		else{//min distance 1,2
			trees[1].setBranchLength(distance[1][2]/2-trees[1].getValueOFChildBranchLength());
			trees[2].setBranchLength(distance[1][2]/2-trees[2].getValueOFChildBranchLength());
			tmpT=new NewickTree(&trees[1],&trees[2]);
			//parent
			trees[1].setParent(*tmpT);
			trees[2].setParent(*tmpT);

			trees[0].setBranchLength((distance[0][2]/2+distance[0][1]/2)/2-trees[0].getValueOFChildBranchLength());
			tmpT->setBranchLength((distance[0][1]/2+distance[0][2]/2)/2-tmpT->getValueOFChildBranchLength());
			tmpT=new NewickTree(tmpT,&trees[0]);
			//parent
			trees[0].setParent(*tmpT);

			}

		if(trees[0].isLeaf()){
			seqList[count]=trees[0];
			trees[0].setPosition(count);
			count++;
		}
		if(trees[1].isLeaf()){
			seqList[count]=trees[1];
			trees[1].setPosition(count);
			count++;
		}
		if(trees[2].isLeaf()){
			seqList[count]=trees[2];
			trees[2].setPosition(count);
			count++;
		}

		//calculate weigth
		//The weights are dependent upon the distance from the root of the tree
		//but sequences which have a common branch with other sequences share the weight
		//derived from the shared branch
		if(verbose)
			cout<<endl<<"Use Tree Guide for calculate Weigth"<<endl;
		for(unsigned int i=0; i<seqList.size();i++)
		{
			seqList[i].setWeigth(calcWeigth(seqList[i].getRoot()));


			if(verbose)
				cout<<seqList[i].getName()<<" "<<seqList[i].getWeigth()<<endl;
		}


		leafList=seqList;


		root=tmpT->getRoot();


	}//end method


    // OPERATORS:
	/**
	 *@Description for create phylogenic Tree use UPGMA method
	 *@param Align2::Alignment ali alignament of starting sequence
	 *@param bool ktuples decide to use ktuples method for calculate identity of 2 sequence
	 *@param bool verbose for more information on video during computing.
	 */
    void NewickTree::upgma(Align2::Alignment ali,bool ktuples, bool verbose){

    		int Nseq=ali.size()+1;//number of seq
    		vector<vector<double> > distance(Nseq, std::vector<double>(Nseq, 0));

    		vector<Alignment>vet=PhyloSupport::calcAlignmentV(&ali,distance,ktuples,verbose);
       		vector<NewickTree> trees(Nseq);


    		if(verbose){
    			cout<<"---------Print Distance Matrix--------------/"<<endl;
    			PhyloSupport::printMatrix(distance);
    		}


    		//Create all node
    		for(unsigned int i=0;i<trees.size();i++){
    			if(i==trees.size()-1)
    				trees[i] = NewickTree(i,PhyloSupport::split(ali.getTargetName(), ' ')[0],ali.getTarget());

    			else
    				trees[i] = NewickTree(i,PhyloSupport::split(ali.getTemplateName(i), ' ')[0],ali.getTemplate(i));

    		}

    		//memorize the nearest seq
    		unsigned int minj=0;
    		unsigned int mini=1;
    		//vector for new distance
    		vector<double> newMediaD ( trees.size());

    		//for creat new matrix;
    		unsigned int newi=0;
    		unsigned int newj=0;


    		while(trees.size()!=1){

    			mini=0;
    			minj=1;
    			//search nearest seq
    			for(unsigned int i=0; i<trees.size(); i++){
    				for(unsigned int j=0; j<trees.size(); j++){
    					if( distance[mini][minj]  > distance[i][j] && i!=j)
    					{
    						minj=j;
    						mini=i;
    					}
    				}
    			}

    			if(verbose)
    				cout<<endl<<"----------------------------min i j find! "<<mini<<" "<<minj<<" --------------------------------"<<endl;


    			//distance of new node U respect A B the children
    			trees[mini].setBranchLength((distance[mini][minj]/2) - trees[mini].getValueOFChildBranchLength());
    			trees[minj].setBranchLength((distance[mini][minj]/2) - trees[mini].getValueOFChildBranchLength());


    			vector<vector<double> > tmpdistance(trees.size()-1, std::vector<double>(trees.size()-1, 0));
    			newi=0;
    			newj=0;

    			//copy old matrix in new matrix, but not mini (A) and minj (B) and not the new node (U)
    			for(unsigned int i=0; i<distance.size();i++){
    				for(unsigned int j=0; j<distance.size();j++){
    					if(i!=mini && j!=minj && i!=minj && j!=mini){
    						if(i<mini && i<minj){
    							newi=i;
    						}
    						else if(i>mini && i<minj){
    							newi=i-1;
    						}
    						else if(i>mini && i>minj){
    							newi=i-2;
    						}
    						if(j<mini && j<minj){
    							newj=j;
    						}
    						else if(j>mini && j<minj){
    							newj=j-1;
    						}
    						else if(j>mini && j>minj){
    							newj=j-2;
    						}
    						if(verbose)
    							cout<<"mini e minj "<<mini<<" "<<minj<<" i e j "<<i<<" "<<j<<" newi e new j "<<newi<<" "<<newj<<" "<<"value "<<distance[i][j]<<endl;
    						tmpdistance[newi][newj]=distance[i][j];
    					}
    				}//end for j
    			}//end for i

    			//insert the new distance
    			/*	example:
    				d(CU) = (d(AC) + d(BC))/2
    				d(DU) = (d(AD) + d(BD))/2
    				d(EU) = (d(AE) + d(BE))/2
    				d(FU) = (d(AF) + d(BF))/2
    			*/
    			newi=0;
    			for(unsigned int i=0; i<distance.size();i++){
    				if(i!=mini && i!=minj){
    					if(i<mini && i<minj){
    						newi=i;
    					}
    					else if(i>mini && i<minj){
    						newi=i-1;
    					}
    					else if(i>mini && i>minj){
    						newi=i-2;
    					}
    				tmpdistance[newi][tmpdistance.size()-1]=(distance[mini][i]+distance[minj][i])/2;
    				if(verbose)
    					cout<<distance[mini][i]<<" + "<<distance[minj][i]<<"/2 = "<<tmpdistance[newi][tmpdistance.size()-1]<<endl;
    				tmpdistance[tmpdistance.size()-1][newi]=tmpdistance[newi][tmpdistance.size()-1];
    				}
    			}


    			distance=tmpdistance;
    			if(verbose){
    				cout<<"The new distance matrix------------------------"<<endl;
    				PhyloSupport::printMatrix(distance);
    				cout<<"the distance of tree vector "<<trees.size()<<endl;
    				cout<<endl<<"---------------Create distance of other node from new node!--------------------------"<<endl;
    			}
    			//insert new tree
    			trees.push_back(NewickTree(&trees[mini],&trees[minj]));

    			//erase trees
    			trees.erase(trees.begin()+mini);
    			//minus 1, because  length was decreased
    			trees.erase(trees.begin()+minj-1);
    			if(verbose)
    				cout<<"the distance of tree vector "<<trees.size()<<endl;


    		}//end while
    		root=trees[0].getRoot();

    		cout<<"------------------end UPGMA-------------------------------"<<endl;
    	}//end method


    // OPERATORS:
	/**
	 *@Description Print in a string the tree using Newick format
	 */
    string NewickTree::printNewickTree()
    {
    	return printNewickTreeNode(*root)+";";
    }

	/**
	 *@Description helper method for recursion look NewickTree::printNewickTree()
	 */
    string NewickTree::printNewickTreeNode(iNode node)
    {
    	string printTree="";
    	std::string varAsString="";
    	std::ostringstream sstream;
    	if(node.right != NULL)
    		printTree+= "(" +printNewickTreeNode(*node.left);
    	if(node.left != NULL)
    		printTree+= "," + printNewickTreeNode(*node.right)+")";
    	if(node.branchLength<=0)
		{
    		return printTree+node.name;
		}
    	else
    		sstream << node.branchLength;
    		varAsString = sstream.str();
    		return printTree+node.name+":"+varAsString;

    }

	/**
	 *@Description calculate the wiegth of seq usign guide tree
	 *	very important for ClustalW
	 */
    double NewickTree::calcWeigth(iNode* seq){
    	double w=seq->branchLength;
    	if(seq->parent!=NULL){
    		w+=( calcWeigth(seq->parent)/(seq->parent->numberOfChildLeaf) );
    	}
    	return w;
    }


    /**
     *@Description Return max values of this node children branchLength.
     * Return zero if left child and right are NULL
     */

    double NewickTree::getValueOFChildBranchLength(){
    	double value1=0;
    	double value2=0;
    	if(root->left!=NULL)
    	{
    		value1=root->left->branchLength;
    	}
    	if(root->right!=NULL)
    	{
    		value2=root->right->branchLength;
    	}
    	if(value1>value2)
    		return value1;
    	return value2;
    }


	/**
	* Set the value of node if is -1 its escape during print
	* @param val
	*
	*/
	void NewickTree::setBranchLength(double val){
		if(val>0)
			root->branchLength=val;
		else
			root->branchLength=0;
	}

	/**
	* Set the value of divergence will use in NJ if is -1 its escape during print
	* @param val
	*
	*/
    void NewickTree::setRdiv(double val){
    	root->divergenceR=val;
    }


   void NewickTree::setParent(NewickTree pTree) {
		root->parent=pTree.getRoot();
	}

   void NewickTree::setWeigth(double w){
	   root->weigth=w;
   }


   void NewickTree::setPosition(unsigned int pos){
	   root->nseq=pos;
   }

   void NewickTree::setAlingSeq(vector<string> vSeq){
	   root->allignSeq=vSeq;
   }


}} // namespace
