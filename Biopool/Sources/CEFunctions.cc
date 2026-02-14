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
/* 
 * Author: Alberto Tilocca
 *
 */
#include <Spacer.h>
#include <Group.h>
#include <SideChain.h>
#include <vector3.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <XyzSaver.h>
#include <SeqSaver.h>
#include <IoTools.h>
#include <IntCoordConverter.h>
#include <IntSaver.h>
#include <GetArg.h>
#include "XyzLoader.h"
#include <cmath>
#include <Kabsch.cc>
#include "CEFunctions.h"


using namespace Victor;using namespace Victor::Biopool;
/**
 * Constructor
 * @param none
 */
CEFunctions::CEFunctions() {

}

/**
 * Destructor
 * @param none
 */
CEFunctions::~CEFunctions() {

}
/**
 * calculate the distances of CA atoms of a spacer.
 * @param Spacer* vector<vector<double> >&
 */
void CEFunctions::innerDistance(Spacer*sp ,vector<vector<double> >&dist ){
    
    int len=sp->sizeAmino();
    for (int row=0;row<len;row++)
       for (int col=0;col<len;col++)
        dist[row][col] = sp->getAmino(row)[CA].distance(sp->getAmino(col)[CA]);


}
/**
 * calculate the similarity matrix using a full set of inter-residue distaces, where all possibile distances except those for neighboring residues are evaluated.
 * @param Spacer* 
 * @param Spacer* 
 * @param vector<vector<double> >&
 * @param vector<vector<double> >&
 * @param vector<vector<double> >&
 *  */
void CEFunctions::createSMatrix(Spacer *sp1,Spacer *sp2, int combination,int m,vector<vector<double> > &score,vector<vector<double> > &distA ,vector<vector<double> > &distB)
    {
    
    int dimSpacer1=sp1->sizeAmino();
    int dimSpacer2=sp2->sizeAmino(); 
    
    for (int i=0;i<dimSpacer1;i++)
        for (int j=0;j<dimSpacer2;j++)
            score[i][j]=-1;
    
    //apply formula 7
        for (int i=0;i<dimSpacer1-(m+1);i++){
           for (int j=0;j<dimSpacer2-(m+1);j++){
             for (int k=0; k<m;k++)
                 for(int l=k+2; l<m;l++)
                   score[i][j]+= fabs(distA[i+k][i+l]-distB[j+k][j+l]);
                        
     
             score[i][j]= score[i][j]/combination;
                
            }}

}
/**
 * search into the similarty matrix the best alignment of AFP according to the CE criteria, including gaps in one of the pair if needed.
 * create a buffer with the best 30 paths.
 * 
 * @param Spacer* 
 * @param Spacer* 
 * @param vector<vector<afp> >&
 * @param vector<vector<double> >&
 * @param vector<vector<double> >&
 * @param int
 * @param int
 * @param int
 * @param vector<vector<double> >
 * @param vector<vector<double> >
 * @param vector<int> 
 *  */
  void CEFunctions::CESearch(Spacer *sp1,Spacer *sp2, vector<vector<afp> >&tempPathBuffer,vector<double>  &scoreBuffer,vector<vector<double> >&score,int smaller,int combination,int m,vector<vector<double> > &distA,vector<vector<double> > &distB,vector<int> &lenBuffer){
  
    int D0=3;         //threshold
    int D1=4;         //threshold
    double bestScore= 1e6;
    int MAXPATH=30;
    int gapMax=30;
    double bestScoreGap =1e6;
    int bestPositionGap=-1;
    vector<afp> bestPath(smaller);
    vector<afp> pathCorr(smaller);
    int bufferIndex=0;
    int bufferSize=0;
    
    int dimSpacer1=sp1->sizeAmino();
    int dimSpacer2=sp2->sizeAmino();    
  
    
    for ( int i = 0; i < smaller; i++ ) {
		bestPath[i].first = -1;
		bestPath[i].second = -1;
	}
    for ( int i = 0; i < MAXPATH; i++ ) {
        for(int j=0;j<smaller;j++){
		// initialize the paths
		scoreBuffer[i] = 1e6;
		lenBuffer[i] = 0;
		
        }
	}
for ( int i = 0; i < smaller; i++ ) {
    for(int j=0;j<smaller;j++){
		tempPathBuffer[i][j].first = 0;
                tempPathBuffer[i][j].second = 0;
    }
	}
    
    vector<int> winCache(smaller);
    for (int i=0;i<smaller;i++)
        winCache[i] = (i+1)*i*m/2 + (i+1)*combination;
    
    vector<int> index(smaller);
    //in this matrix there will be the score releted to the gaps
   vector<vector<double> > partialScoreGap(smaller);
   
   for (int i = 0; i< smaller; i++)
      partialScoreGap[i].resize(gapMax*2+1);
   
   //inizialize
	for ( int i = 0; i < smaller; i++ ) 
		for ( int j = 0; j < gapMax*2+1; j++ )
			partialScoreGap[i][j] = 1e6;
     
       int bestLength=0;
      
     //here start the searching
   
        for(int i=0;i<dimSpacer1;i++){
            if(i>dimSpacer1-m*(bestLength-1)){
      
                break;
            }
            for(int j=0;j<dimSpacer2;j++){
               
                if(score[i][j]>=D0 )
                    continue;
                 
                if(score[i][j]==-1 )
                    continue;
               
                if(j> dimSpacer2 - m*(bestLength-1))
                    break;
                 
                //reset corrent path
                for (int z=0;z<smaller;z++){
                    pathCorr[z].first=-1;
                    pathCorr[z].first=-1;
                }
                pathCorr[0].first = i;
                pathCorr[0].second=j;
                int pathLengthCorr=1;
                index[pathLengthCorr-1]=0;
                
                double allScoreCorr=0.0;
                
                //check all the possible path starting from i , j
                bool ok=false;
                
                while(!ok){
                
                     bestScoreGap =1e6;
                     bestPositionGap=-1;
                     //according to the fact that gaps can stay just in one of the two afp use the even and odds method
                    
                     for(int g=0;g<(gapMax*2)+1;g++){
                         //gli afp successivi rispettivamente
                         int succ1=pathCorr[pathLengthCorr-1].first+m;
                         int succ2=pathCorr[pathLengthCorr-1].second+m;
                         if ((g+1)%2 ==0){
                             succ1+=(g+1)/2 ;//gap in the first
                         }
                         else{
                             succ2+=(g+1)/2; //gap in the second
                         }
                         
                         // check if going out of matrix, can be not necessary
					if ( succ1> dimSpacer1-m-1 || succ2 > dimSpacer2-m-1 ){
						continue;
					}
					//se non è buono con i gap.
					if ( score[succ1][succ2] > D0 )
						continue;
					// sicurezza nel caso non abbia modificato tutta la mtrice
					if ( score[succ1][succ2] == -1.0 )
						continue;
                                         
                          double scoreCorr=0.0;
                                    
                     //calculate formula 6
                       for(int s=0;s<pathLengthCorr;s++){
                           scoreCorr += fabs(distA[pathCorr[s].first][succ1]-distB[pathCorr[s].second][succ2]);
                           scoreCorr += fabs(distA[pathCorr[s].first + (m-1)][succ1+(m-1)]-distB[pathCorr[s].second + (m-1)][succ2+(m-1)]);
                           
                           for(int k=0;k<m-2;k++)
                               scoreCorr+=fabs(distA[pathCorr[s].first + k][succ1 + (m-1)-k] - distB[pathCorr[s].second +k][succ2 +(m-1)-k]);
                           
                       }
                          // formula 10
                           scoreCorr /=(double)m * (double)pathLengthCorr;
                           
                           if (scoreCorr >= D1){
                            continue;
                           }
                           //keep the best score with gaps
                           
                           if( scoreCorr < bestScoreGap ){
                               pathCorr[pathLengthCorr].first =succ1;
                               pathCorr[pathLengthCorr].second=succ2;
                               bestScoreGap = scoreCorr;
                               bestPositionGap = g;
                               partialScoreGap[pathLengthCorr -1][g] = scoreCorr;
                           }                        
                     }//end of gapping
                     allScoreCorr =0.0;
                     double score1=0.0;
                     double score2=0.0;
                     
                 
                     int addGap;
                     int gap1,gap2;
                     if (bestPositionGap != -1){
                         addGap=(bestPositionGap +1)/2;
                         if ((bestPositionGap +1)%2==0){ //gap nel primo
                             gap1=pathCorr[pathLengthCorr-1].first + m + addGap;
                             gap2=pathCorr[pathLengthCorr-1].second + m;
                         }
                         else{ //gap nel secondo
                         gap1=pathCorr[pathLengthCorr-1].first + m ;
                         gap2=pathCorr[pathLengthCorr-1].second + m + addGap;
                         }
                         
                         score1= (partialScoreGap[pathLengthCorr-1][bestPositionGap]
                                 *m*pathLengthCorr + score[gap1][gap2]*combination)/
                                 (m*pathLengthCorr+combination);
                         score2= ((pathLengthCorr>1 ? (partialScoreGap[pathLengthCorr-2][index[pathLengthCorr-1]])
                                 :score[i][j])*winCache[pathLengthCorr-1]+ score1 
                                 * (winCache[pathLengthCorr]-winCache[pathLengthCorr-1]))
                                 /winCache[pathLengthCorr];
                         
                         allScoreCorr=score2;
                         
                         if(allScoreCorr >D1){
                             ok=true;
                             bestPositionGap = -1;
                             break;
                         
                         }
                         else {
                             partialScoreGap[pathLengthCorr-1][bestPositionGap] = allScoreCorr;
                             index[pathLengthCorr]=bestPositionGap;
                             pathLengthCorr++;
                         }
                     }
                     else{
                         //no good path with gaps
                         ok=true;
                         pathLengthCorr--;
                         break;
                                              }
                    
                     //if the gapped path is longer or same size but higher score keep the the gapped one.
                     if(pathLengthCorr > bestLength || (pathLengthCorr == bestLength && allScoreCorr<bestScore)){
                     bestLength=pathLengthCorr;
                     bestScore=allScoreCorr;
                     //corrent path become best path
                     
                     for (int i=0;i<smaller;i++)
                     {
                     bestPath[i].first = pathCorr[i].first;
                     bestPath[i].second = pathCorr[i].second;
                     }
                    
                     }
                     
                }//while
            //here I already have the best path
                
            //keep the 20 best path for the last optimization.
                
                if (bestLength > lenBuffer[bufferIndex] ||  
                    ( bestLength == lenBuffer[bufferIndex] &&
					   bestScore < scoreBuffer[bufferIndex] )){
                    bufferIndex = ( bufferIndex == MAXPATH-1 ) ? 0 : bufferIndex+1;
	       	bufferSize = ( bufferSize < MAXPATH ) ? bufferSize+1 : MAXPATH;
                
                vector<afp> pathCopy(smaller);
                for (int i=0;i<smaller;i++)
                {
                    pathCopy[i].first=bestPath[i].first;
                    pathCopy[i].second=bestPath[i].second;
                }
                
                 if ( bufferIndex == 0 && bufferSize == MAXPATH ) {
                for(int i=0;i<smaller;i++){
                        tempPathBuffer[MAXPATH-1][i].first=-1;
                        tempPathBuffer[MAXPATH-1][i].second=-1;
                            }                    
                for(int i=0;i<smaller;i++){
                         
                        tempPathBuffer[MAXPATH-1][i].first=pathCopy[i].first;
                        tempPathBuffer[MAXPATH-1][i].second=pathCopy[i].second;
                        }
                scoreBuffer[MAXPATH-1] = bestScore;
                 lenBuffer[MAXPATH-1] = bestLength;
                 
                 }
                 else{
                 for(int i=0;i<smaller;i++){
                         
                        tempPathBuffer[bufferIndex-1][i].first=-1;
                        tempPathBuffer[bufferIndex-1][i].second=-1;
                        
                 }    
                     
                 for(int i=0;i<smaller;i++){
                         
                        tempPathBuffer[bufferIndex-1][i].first=pathCopy[i].first;
                        tempPathBuffer[bufferIndex-1][i].second=pathCopy[i].second;
                        
                 }
                 scoreBuffer[bufferIndex-1] = bestScore;
                 lenBuffer[bufferIndex-1] = bestLength;
                 
                 }
                 }
              
               
                
             
                 for (int i=0;i<smaller;i++){
                                            pathCorr[i].first=-1;
                                            pathCorr[i].second=-1;
                                        }
                bestLength=0; 
                
            }
        
        }
  }
/**
 * for each of the 30 best path first superimpose the path then calulate the RMSD, this is done using the Kabsch algorithm, then store the path 
 * with the lower RMSD and finally write the pdb file.
 * 
 * @param Spacer* 
 * @param Spacer* 
 * @param vector<vector<afp> >&
 * @param int
 * @param vector<int> 
 * @param int
 * @param string
 * @param int
 *  */
int CEFunctions::finalSteps(Spacer *sp1,Spacer *sp2,vector<vector<afp> >&tempPathBuffer,int smaller,vector<int> &lenBuffer,int m,string outputFile,string outputFile2 ){
    
    int MAXPATH=30;
    Eigen::Matrix3Xd Matrix1;
    Eigen::Matrix3Xd Matrix2;
    Matrix1.resize(3,smaller); 
    Matrix2.resize(3,smaller); 
    double minRMSD=1e06;
    int indexFinalPath;
    Spacer *bestSp1;
    Spacer *bestSp2;
    
    for (int p=0;p<MAXPATH;p++) {
       Spacer* copia1 = new Spacer(*sp1);
       Spacer* copia2 = new Spacer(*sp2);            
       int l=0;
       int l2=0;                
       int cnt=0;
       int cnt2=0;
       int iter=0;
       int iter2=0;
       bool end1;
       bool end2;
                   
       for (int i=0;i<lenBuffer[p];i++){                   
            //deleting the componets  
            
            //chain 1
            end1 =false;               
            while(end1==false)
              {
                if (iter>=(int)(sp1->sizeAmino()-1))
                     end1=true;
                      if (l<tempPathBuffer[p][i].first){                
                           copia1->removeComponentFromIndex(cnt*m);
                           iter++;
                           l++;
                       }
                       else{                             
                            l=l+m+cnt;
                            cnt++;
                            iter++;
                            end1=true;
                       }
                   }
                      
               //chain 2
            end2 =false;
            while(end2==false)
              {
                if (iter2>=(int)(sp2->sizeAmino()-1))
                     end2=true;
                      if (l2<tempPathBuffer[p][i].second){
                           copia2->removeComponentFromIndex(cnt2*m);
                           iter2++;
                           l2++;
                       }
                       else{     
                            l2=l2+m+cnt2;
                            cnt2++;
                            iter2++;
                            end2=true;
                       }
                   }
              }                   
        double RMSD=0.0;
        vgVector3<double>CaCoordsSp1;
        vgVector3<double>CaCoordsSp2;

        for (int c=0; c< smaller;c++){
             CaCoordsSp1 = copia1->getAmino(c)[CA].getCoords();
             CaCoordsSp2 =  copia2->getAmino(c)[CA].getCoords();
             Matrix1(0,c)=CaCoordsSp1.x;
             Matrix1(1,c)=CaCoordsSp1.y;
             Matrix1(2,c)=CaCoordsSp1.z;
             Matrix2(0,c)=CaCoordsSp2.x;
             Matrix2(1,c)=CaCoordsSp2.y;
             Matrix2(2,c)=CaCoordsSp2.z;                
              }    

            vgMatrix3<double> vgMat;
            vgVector3<double>vgVet;
            //apply kabasch for the first time
            Eigen::Affine3d A= Find3DAffineTransform(Matrix2,Matrix1);
            //coordinates of rototaion matrix in victor format
            vgMat.x.x=A.rotation().coeff(0,0);
            vgMat.x.y=A.rotation().coeff(0,1);
            vgMat.x.z=A.rotation().coeff(0,2);
            vgMat.y.x=A.rotation().coeff(1,0);
            vgMat.y.y=A.rotation().coeff(1,1);
            vgMat.y.z=A.rotation().coeff(1,2);
            vgMat.z.x=A.rotation().coeff(2,0);
            vgMat.z.y=A.rotation().coeff(2,1);
            vgMat.z.z=A.rotation().coeff(2,2);
            
            //coordinates of rototaion matrix in victor format
            vgVet.x=A.translation().coeff(0);
            vgVet.y=A.translation().coeff(1);
            vgVet.z=A.translation().coeff(2);
             
            //apply rotation
             copia2->getAmino(0)[CA].addRot(vgMat);
             copia2->sync();

            //new coordinates after ONLY rotation
            for (int c=0; c< smaller;c++)
              {
                CaCoordsSp2 =  copia2->getAmino(c)[CA].getCoords();
                Matrix2(0,c)=CaCoordsSp2.x;
                Matrix2(1,c)=CaCoordsSp2.y;
                Matrix2(2,c)=CaCoordsSp2.z;
              }
          
             //again Kabash after applied rotation
            A= Find3DAffineTransform(Matrix2,Matrix1);
            
            //Victor format
            vgVet.x=A.translation().coeff(0);
            vgVet.y=A.translation().coeff(1);
            vgVet.z=A.translation().coeff(2);
            
            //apply translation
             copia2->getAmino(0)[CA].addTrans(vgVet);
             copia2->sync();
           //calculate RMSD
            for (int i=0; i<smaller; i++)
                RMSD += pow(copia1->getAmino(i)[CA].distance(copia2->getAmino(i)[CA]),2);          
                     RMSD= sqrt((RMSD/smaller));
          
            if (RMSD<minRMSD){
                minRMSD=RMSD;
                indexFinalPath=p;
            //save the best spacer
                 bestSp1 = new Spacer(*copia1);
                 bestSp2 = new Spacer(*copia2);
                 
            }
            
            delete copia1;
            delete copia2;
           }
            //write files
            ofstream outFile(outputFile.c_str());             
            ofstream outFile2(outputFile2.c_str());  
            
            outFile <<"";
            outFile.close();
            outFile.open(outputFile.c_str(),ios::app);
            PdbSaver ps(outFile);       // Pdb format
            PdbSaver ps2(outFile2);
            outFile <<"MODEL    1"<<endl;
            bestSp1->save(ps);
            outFile <<"ENDMDL"<<endl;
            outFile <<"MODEL    2"<<endl;
            bestSp2->save(ps);
            outFile <<"ENDMDL"<<endl;
            bestSp1->save(ps2);
            outFile.close();
            cout<<"The pdb file has been created."<<endl;
            cout << "RMSD: "<<minRMSD<<endl;
            return indexFinalPath;
}
