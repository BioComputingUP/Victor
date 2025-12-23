#include <PdbLoader.h>
#include <Protein.h>
#include <iostream>
#include <Eigen/Geometry> 
#include <StructuralSuperimposition.h>

using namespace Victor::Biopool;
using namespace Victor; 
using namespace std;


 
int main( int argc, char* argv[] ) 
{
	double d0= 3.5; //algorithm threshold
	int L=4; //algorithm seeds size
	string inputFile1, inputFile2, outputFile, score;

	// Check the number of parameters
    if (argc < 5)
    {
        // Tell the user how to run the program
        cout<<endl<<endl;
        cout<<"STRUCTURAL SUPERIMPOSITION utility:"<<endl;
        cout<<"This program evaluates the optimal rigid superimposition of two sequence equivalent proteic structures."<<endl;
        cout<<"It uses an heuristic iterative algorithm (Siew, Elofsson, Rychlewski, Fischer) to find the largest subset "<<endl;
        cout<<"of model (A) residues that best superimpose with their corresponding residues in the experimental structure (B)"<<endl;
        cout<<"From the resulting subset it's possible to select which superimposition quality score to view."<<endl;
        cout<<"The avaiable score are: RMSD, GDT_TS, TM-Score, MaxSub."<<endl<<endl<<endl;
        cout<<"Usage: "<<endl<<endl;
        cout<<"\t"<<argv[0]<<"\t<in_file1>\t<in_file2>\t<out_file>\t<seed_size>\t<score>\t[OPTIONS]"<<endl<<endl;
        cout<<"Files required:"<<endl<<endl;
        cout<<"\tin_file1\t PDB file of the first protein structure (the \"native\" structure)"<<endl;
        cout<<"\tin_file2\t PDB file of the second protein structure (the \"experimental\" structure)"<<endl;
        cout<<"\tout_file\t PDB file of the experimental structure optimally superimposed"<<endl<<endl;
        cout<<"Seed: "<<endl<<endl;
        cout<<"\t<integer>\t Select the initial seed matches size (must be >=4, default is 4)."<<endl<<endl;
        cout<<"Scores: "<<endl<<endl;
        cout<<"\t[-r <double>]\tUse RMSD score with custom d0 threshold (default is 3.5 A)"<<endl;
        cout<<"\t[-t]\t\tUse TM-Score"<<endl;
        cout<<"\t[-m <double>]\tUse MaxSub score with custom d0 threshold (default is 3.5 A)"<<endl;
        cout<<"\t[-g]\t\tUse GDT-TS score"<<endl;
        cout<<"\t[-a]\t\tAll"<<endl<<endl;
        
        return 1;
    }
    else
    {
    	//argument parsing 
    	inputFile1 = argv[1];
    	inputFile2 = argv[2];
    	outputFile = argv[3];
    	L= strtol(argv[4],NULL,10);
    	score = argv[5];

    	if(L<4)
    	{
    		cout<<"\nError on parsing seed_size: input must be >= 4."<<endl;
    		return 1;
    	}

    	//options parsing
	    if(argc > 6)
	    {
		    d0 = strtod (argv[6], NULL);
		    if(d0 <= 0.0)
		    {
		    	cout<<"\nError on parsing the d0 threshold: input is not a number and can't be <= 0"<<endl;
		    	return 1;
		    }
	    }


    	if(score=="-r" || score=="-t" || score=="-m" || score=="-g" || score=="-a")
	    {
	    	
	    	ifstream inFile1(inputFile1.c_str());
			ifstream inFile2(inputFile2.c_str());
			
			Protein prot1, prot2;   
			Spacer *sp1, *sp2;

	
		    PdbLoader pl1(inFile1);
		    PdbLoader pl2(inFile2);

		    pl1.setNoHAtoms();
		    pl2.setNoHAtoms();

			cout<<"Loading protein A ("<<inputFile1<<")"<<endl;
			prot1.load(pl1);       
			cout<<endl<<endl<<"Loading protein B ("<<inputFile2<<")"<<endl;
			prot2.load(pl2);       
			cout<<endl<<endl;
			

			vector<char> chain1 = prot1.getAllChains();
			vector<char> chain2 = prot2.getAllChains();
			
			sp1=prot1.getSpacer(chain1[0]);
			sp2=prot2.getSpacer(chain2[0]);
			
			StructuralSuperimposition s1, s2= StructuralSuperimposition();
			
			Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> A=StructuralSuperimposition::convertSpacerToEigen(sp1);

			Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> B=StructuralSuperimposition::convertSpacerToEigen(sp2);

			if(L>A.rows() || L>B.rows())
    		{
    			cout<<"\nError: L is bigger than the protein size ("<<A.rows()<<")"<<endl;
    			return 1;
    		}
			
			StructuralSuperimposition ssi= StructuralSuperimposition();

	    	if(score=="-r")
	    	{
	    		double rmsdscore= ssi.getRmsdScore(A,B,L,d0);
	    		cout<<"RMSD score: "<<rmsdscore<<endl;
	    		ssi.printPDBFile(outputFile.c_str(),*sp1,*sp2, A);
	    	}

	    	if(score=="-t")
	    	{
				double tmscore= ssi.getTMScore(A, B, L);
				cout<<"TM-score: "<<tmscore<<endl;
				ssi.printPDBFile(outputFile.c_str(),*sp1,*sp2, A);
	    	}
	    	if(score=="-m")
	    	{
	    		double msscore= ssi.getMaxSubScore(A, B, L, d0);
	    		cout<<"MaxSub score: "<<msscore<<endl;
	    		ssi.printPDBFile(outputFile.c_str(),*sp1,*sp2, A);
	    	}

	    	if(score=="-g")
	    	{
	    		double gdttsscore= ssi.getGDTTSScore(A,B,L);
				cout<<"GDT-TS score: "<<gdttsscore<<endl;
				ssi.printPDBFile(outputFile.c_str(),*sp1,*sp2, A);
	    	}

	    	if(score=="-a")
	    	{
	    		double rmsdscore= ssi.getRmsdScore(A,B,L,d0);
	    		cout<<"\nRMSD score: "<<rmsdscore<<endl<<endl;

	    		ssi.printPDBFile((outputFile+".RMSD.pdb").c_str(),*sp1,*sp2, A);

	    		double tmscore= ssi.getTMScore(A, B, L);
	    		cout<<"TM-score: "<<tmscore<<endl<<endl;
	    		ssi.printPDBFile( (outputFile+".TMS.pdb").c_str(),*sp1,*sp2, A);

	    		double msscore= ssi.getMaxSubScore(A, B, L, d0);
	    		cout<<"MaxSub score: "<<msscore<<endl<<endl;

	    		ssi.printPDBFile( (outputFile+".MSS.pdb").c_str(),*sp1,*sp2, A);

	    		double gdttsscore= ssi.getGDTTSScore(A,B,L);
				cout<<"GDT-TS score: "<<gdttsscore<<endl<<endl;
				ssi.printPDBFile( (outputFile+".GDTTS.pdb").c_str(),*sp1,*sp2, A);
	    	}
	    }
		else
		{
			cout<<"\nInvalid score selection, the valid ones are: -r (RMSD), -t (TM-Score), -m (MaxSub), -g (GDT-TS)."<<endl;
			return 1;
		}
	}	
	cout<<endl<<endl;
}
