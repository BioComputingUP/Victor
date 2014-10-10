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
 @Description This program calculates a pseudo-energy to evaluate the quality  
 *of a given protein structural model, as expressed in a single  
 *(real) number.
 */
#include <string>
#include <GetArg.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <Potential.h>
#include <TorsionPotential.h>
#include <PhiPsi.h>
#include <PhiPsiOmega.h>
#include <PhiPsiOmegaChi1.h>
#include <PhiPsiOmegaChi1Chi2.h>
//#include <Omega.h>
//#include <Chi1Chi2.h>
#include <StatTools.h>
using namespace Victor::Biopool;
using namespace Victor::Energy;
using namespace Victor;

void sShowHelp()
{
	cout << "pdb2tap - PDB Torsion angle propensity\n\n"
	<< "This program calculates a pseudo-energy to evaluate the quality "
	<< "of a given protein structural model, as expressed in a single "
	<< "(real) number. Contact the author for further details.\n"
	<< "Developed by Silvio Tosatto <silvio.tosatto@unipd.it> \n"
	<< "\nOPTIONS:\n"
	<< "\t-i <filename> \t\t Input  PDB file." <<"\n"
	<< "\t-I <filename> \t\t Input a filelist of PDB file.\n"
	<< "\n"
	<< "\t[-c <id>]  \t\t  ID of chain to load from PDB file\n"
	<< "\t[--nmr] \t\t Calculate average over NMR ensemble.\n"
	<< "\t[-s <int> ] \t\t Start residue (DEFAULT = second AA)\n"
	<< "\t[-e <int> ] \t\t End residue (DEFAULT = last -1)\n"
	<< "\n\t--acc   \t\t Display information on estimated accuracy of stucture\n"
	<< "\n\t--help  \t\t Extended help for the method\n"
	<< endl;
}

void sShowSpecial()
{
	sShowHelp();
	
	cout << "\nADDITIONAL OPTIONS:\n"
		<< "\t[--casp]    \t\t CASP-like output mode\n"
		<< "\t[-v]    \t\t Verbose mode\n"
		<< "\t[-f]    \t\t Force utilisation of first model disregarding its number\n"
		<< "\t[-P] \t\t\t Per Residue output file (NB: quick hack)\n"
		//       << "\t[-p] \t\t\t Per Residue value\n"
		<< "\t[-r]        \t\t Rank tap values of file list models\n"
		<< "\t     \t\t\t Should be used with -f flag if models are CONCOORD output\n"
		<< "\t[--know <filename>] \t Knowledge filename (Default = tor.par).\n"
		<< "\t[--allchains]     \t to consider all chain of a pdb file.\n"
		<< "\t         \t\t Relative energy in this casa will be the average.\n"
		<< "\t[--cutoff <integer>] \t Length of chain that will not considered.\n"
		<< "\t         \t\t Usable only with --allchain option. DEFAULT 20.\n"
		<< "\n"
		<< "\t[--pp]      \t\t consider Phi and Psi Angles\n"
		<< "\t[--ppo]     \t\t consider Phi, Psi, Omega Angles\n"
		<< "\t[--ppoc]    \t\t consider Phi, Psi, Omega, and Chi1 Angles\n"
		<< "\t[--ppocc]   \t\t consider Phi, Psi, Omega, Chi1, and Chi2 Angles (DEFAULT) \n"
		<< "\t[--bin <integer>] \t the bin for Phi and Psi angles (Default = 10 degrees).\n"
		<< "\t[--maxOnly] \t\t Use only maximum energy normalization instead of min/max.\n"
		<< "\t[--raw] \t\t Do _NOT_ use energy normalization.\n"
		<< endl;
}


double sGetAcc(double tap, unsigned int c)
{
	double acc3[27] = {0.766781097, 0.766893126,	0.767100263,	0.767212396,	0.76770041,	
		0.768864469,	0.770813432,	0.775362319,	0.78161348,	0.790236173,	0.803274908,	
		0.821307468,	0.843368559,	0.869625281,	0.897667069,	0.925998758,	0.950525287,	
		0.968298807,	0.980840259,	0.989392085,	0.994704325,	0.995485327,	0.982608696,	
		0.962962963,	0.916666667,	0.75};
	double acc2[27] = {0.420495216,	0.420556651,	0.420710319,	0.420771817,	0.421079579,	
		0.421684982,	0.422955397,	0.425539781,	0.429018789,	0.434090395,	0.442266298,	
		0.453197697,	0.470222004,	0.493869798,	0.526907705,	0.572448768,	0.636696799,	
		0.706842436,	0.771468144,	0.82619339,	0.87466902,	0.900677201,	0.886956522,	
		0.888888889,	0.833333333,	0.5};
	double acc1[27] = {0.085603681,	0.085616188,	0.085647471,	0.085659991,	
		0.085722645,	0.085860806,	0.08611948,	0.086660751,	0.087384432,	0.088432808,	
		0.090098401,	0.092421733,	0.096010486,	0.101105163,	0.108374384,	
		0.120057959,	0.139139995,	0.167922159,	0.210064635,	0.263157895,	
		0.337157988,	0.426636569,	0.365217391,	0.444444444,	0.333333333,	0.25};
	
	double tmp = (tap - 0.64) * 100;
	if (tmp < 0)
		tmp = 0;
	if (tmp > 29)
		tmp = 29;
	unsigned int data = static_cast<unsigned int>(tmp);
	
	if (c == 3)
		return 100 * acc3[data];
	else if (c == 2)
		return 100 * acc2[data];
	else
		return 100 * acc1[data];
	
	return 0.0;
}

double sGetCov(double tap, unsigned int c)
{
	double cov3[27] = {1,	1,	0.999904744,	0.999904744,	0.999809488,	0.999714231,	0.99923795,	
		0.998856925,	0.998571156,	0.997618594,	0.995332444,	0.992093732,	0.980662983,	
		0.959420842,	0.919984759,	0.852257573,	0.741188798,	0.587730996,	0.404743761,	
		0.23099638,	0.107353782,	0.042008002,	0.010763955,	0.002476662,	
		0.001047819,	0.000285769};
	double cov2[27] = {1,	1,	1,	1,	1,	0.999826298,	0.999826298,	0.999652597,	0.999478895,	
		0.999305194,	0.999305194,	0.998262984,	0.997047073,	0.993573042,	0.984714261,	
		0.960743443,	0.905332639,	0.782351919,	0.580510683,	0.351745701,	0.172138266,	
		0.069306931,	0.017717561,	0.004168838,	0.001737016,	0.000347403};
	double cov1[27] = {1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	0.999146758,	0.994880546,	
		0.989761092,	0.971843003,	0.912969283,	0.776450512,	0.550341297,	0.325938567,	
		0.161262799,	0.035836177,	0.010238908,	0.003412969,	0.000853242};
	
	double tmp = (tap - 0.64) * 100;
	if (tmp < 0)
		tmp = 0;
	if (tmp > 29)
		tmp = 29;
	unsigned int data = static_cast<unsigned int>(tmp);
	
	if (c == 3)
		return 100 * cov3[data];
	else if (c == 2)
		return 100 * cov2[data];
	else
		return 100 * cov1[data];
	
	return 0.0;
}

void sShowAccuracy(double tap)
{
	cout << "\nBased on an analysis of almost 14,000 X-ray structures from the PDB, the \n"
	<< "following confidence estimates could be derived for the TAP scores of X-ray\n"
	<< "structures at various experimental quality levels:\n\n"
	<< "Experimental quality cutoffs:\n\n" 
	<< "               \t Resolution \t R-free \t Fraction of X-ray structures\n"
	<< "Medium (ME)    \t <= 2.5 A  \t <= 0.3   \t     76.7% \n"
	<< "High (HI)      \t <= 2.2 A  \t <= 0.25  \t     42.2% \n"
	<< "Very high (VH) \t <= 1.7 A  \t <= 0.22  \t      8.7% \n"
	<< "\n\n"
	<< "Confidence estimates based on TAP score cutoffs:\n\n"
	<< "Accuracy (acc) measures the probability (in %) that the current TAP score\n"
	<< "corresponds to a structure of the given experimental quality class.\n\n"
	<< "Coverage (cov) measures the fraction (in %) of structures for a given \n"
	<< "experimental quality class that have (at least) the current TAP score.\n\n"
	<< endl;
}

inline void writeAccuracy(double tap)
{
	double acc1, acc2, acc3;
	double cov1, cov2, cov3;
	
	acc1 = sGetAcc(tap, 1);
	acc2 = sGetAcc(tap, 2);
	acc3 = sGetAcc(tap, 3);
	cov1 = sGetCov(tap, 1);
	cov2 = sGetCov(tap, 2);
	cov3 = sGetCov(tap, 3);
	
	cout << "ME: acc= " << setw(5) << setprecision(1) << acc3 << "  cov= " 
		<< setw(5) << setprecision(1) << cov3 << " \t";
	cout << "HI: acc= " << setw(5) << setprecision(1) << acc2 << "  cov= " 
		<< setw(5) << setprecision(1) << cov2 << " \t";
	cout << "VH: acc= " << setw(5) << setprecision(1) << acc1 << "  cov= " 
		<< setw(5) << setprecision(1) << cov1 << " \t";
}


inline void Header(string inputFile, bool verbose, bool nmr, int i)
{
	if (verbose)
		cout << inputFile << "\t";
	
	if (nmr)
		cout << "Model# "<< i <<"\t";
}

static void qsort(vector<double> a, vector<string> b, int lo, int hi)
{
	int h, l;
	double p,t;
	string r;
	
	if (lo < hi)
    {
		l = lo;
		h = hi;
		p = a[hi];
		// t = b[hi];
		
		do 
		{
			while ((l < h) && (a[l] <= p))
				l = l+1;
			while ((h > l) && (a[h] >= p))
				h = h-1;
			if (l < h)
			{
				t = a[l];
				a[l] = a[h];
				a[h] = t;
				
				r = b[l];
				b[l] = b[h];
				b[h] = r;
			}
		} while (l < h);
		
		t = a[l];
		a[l] = a[hi];
		a[hi] = t;
		
		r = b[l];
		b[l] = b[hi];
		b[hi] = r;
		
		qsort(a, b, lo, l-1);
		qsort(a, b, l+1, hi);
    } // if
}

int 
main(int nArgs, char* argv[])
{
	if (nArgs == 1)
    {
		cerr << "Missing argument specification. Use -h to get help.\n";
		return 1;
    }
	
	if (getArg( "h", nArgs, argv))
    {
		sShowHelp();
		return 1;
    };
	
	if (getArg( "-help", nArgs, argv))
    {
		sShowSpecial();
		return 1;
    };
	
	//read options from command line
	string inputFile, listFile, know, chainID, perresFile;
	bool verbose = getArg( "v", nArgs, argv);
	bool caspVerbose = getArg( "-casp", nArgs, argv);
	if (verbose)
		caspVerbose = true;
	
	getArg( "i", inputFile, nArgs, argv, "!"); 
	getArg( "I", listFile, nArgs, argv, "!");
	bool perres = getArg("p", nArgs, argv);    // per residue value
	bool showAccuracy = getArg("-acc", nArgs, argv);  
	bool maxOnly = getArg("-maxOnly", nArgs, argv);  
	bool rawOnly = getArg("-raw", nArgs, argv);  
	
	getArg( "P", perresFile, nArgs, argv, "!"); // per-residue hack
	
	if ((inputFile == "!") && (listFile == "!") )
        {
                    cerr << "Missing input file. Use -h to get help.\n";
                    return 1;
        }
	
	bool force = getArg("f", nArgs, argv);
	bool rank = getArg("r", nArgs, argv);
	//knowledge filename, default is TOP500
	getArg( "-know", know, nArgs, argv, "!");
	
	//which angle to consider
	bool pp = getArg("-pp", nArgs, argv);      
	bool ppo = getArg("-ppo", nArgs, argv);
	bool ppoc = getArg("-ppoc", nArgs, argv);
	bool ppocc = getArg("-ppocc", nArgs, argv);
	
	int ARCSTEP1;
	unsigned int CUTOFF;
	unsigned int START = 0, END = 9999;
	unsigned int START1 = 0, END1 = 9999;
	
	//set ranges for phi, psi angles
	getArg( "-bin", ARCSTEP1, nArgs, argv, 10);
	
	getArg( "c", chainID, nArgs, argv, " ");
	bool nmr = getArg( "-nmr", nArgs, argv);
	bool allchains = getArg( "-allchains", nArgs, argv);
	getArg( "-cutoff", CUTOFF, nArgs, argv, 20);
	int MEMORYEND = END;
	
	getArg( "s", START1, nArgs, argv, 9999 );
	getArg( "e", END1, nArgs, argv, 9999 );
	
	double tap = 9999.9;
	
	vector<double> tapNmr; // vector containing the per-model TAP scores, for statistical purposes
	
	if ( START1 < 1 )
		ERROR("Impossible to start before the second residue.",exception);
	
	if ( (allchains) && (chainID != " ") )
		ERROR("ERROR. You have choose two incompatible options: --allchains and --chain.", 
			  exception);
	
	//incompatible sets of angles
	if ( ((pp) && (ppo)) || ((pp) && (ppoc)) || ((pp) && (ppocc)) ||
		 ((ppo) && (ppoc)) || ((ppo) && (ppocc)) ||  
		 ((ppoc) && (ppocc)) )
		ERROR("ERROR in choosing options. You have chosen too many set of angles.", exception);
	
	//verify range of angles phi, psi
	if ( ( ( 360%ARCSTEP1 ) != 0 ) )
		ERROR("The range of angles (phi, psi,) must be divisor of 360.", exception);
	
	//energy potential, assignment will be done according to the chosen angles
	TorsionPotential* pot = NULL; 
	
	//Set default knowledge
	if ( know == "!" )
        {           char *victor=getenv("VICTOR_ROOT");
                    if (victor  == NULL)
                        ERROR("Environment variable VICTOR_ROOT was not found.\n Use the command:\n export VICTOR_ROOT=......", exception);
                    string env;
                    env = getenv("VICTOR_ROOT");
                    string tmpknow = "data/tor.par";
                    string defaultknow = env + tmpknow;
                    know = defaultknow;
                    if (defaultknow.length() < 3 )
                            ERROR("Environment variable VICTOR_ROOT was not found.", exception);

        }
	
	//torsion potential to calculate
	if ( pp )
		pot =  new PhiPsi(ARCSTEP1, know);
	else if ( ppo )
		pot = new PhiPsiOmega(ARCSTEP1, know);
	else if ( ppoc )
		pot = new PhiPsiOmegaChi1(ARCSTEP1, know);
	else // default 
		pot = new PhiPsiOmegaChi1Chi2(ARCSTEP1, know);
    
	ifstream inFile(listFile.c_str());
	if ((!inFile) && (listFile != "!"))
		ERROR("File " + listFile + " not found.", exception);
	
	long double tres = 0.00; // the TorsionPotential propensities calculated
	long double total = 0.00; // the theoric best TorsionPotential
	long double totalMin = 0.00; // the theoric worst TorsionPotential
	int nmrmodel = 0; // the number of nmr model.
	
	vector<double> tapRank;
	vector<double> resTap; // TAP per-residue hack
	vector<string> models;
	
	while ((inFile) || (inputFile != "!"))
	{
		if (listFile != "!")
		{
			inFile >> inputFile;
			if (!inFile)
				break;
		}
		
		//initialization of a vector containing the IDs of all chains
		vector<char> totalchain;
		ifstream pdbfile2(inputFile.c_str());
		if ( !pdbfile2 )
			ERROR( "File not found", exception );
		
		//load components in standard PDB format from pdbfile2
		PdbLoader pl(pdbfile2);
		if ( allchains )
		{
			//insert all ID of chains present in PDB file
			totalchain = pl.getAllChains();
			if ( totalchain.size() == 0 )
				totalchain.push_back(chainID.c_str()[0]); 
		}
		else
		{
                    if (chainID!=" "){
			totalchain.push_back(chainID.c_str()[0]);
                    }
                    else{
                        totalchain.push_back(pl.getAllChains()[0]);
                    }
		}
                
                
                
		//maximum number of nmr models
		unsigned int max;
		if (!nmr)
			max = 1;
		else
		{
			max = pl.getMaxModels();
			
			if ( max == 0)
			{
				max = 1;
				if (verbose)
					cerr <<"Warning: the file " << inputFile 
						<< " is probably not an nmr structure.\n";
			}
		}
		
		vector<unsigned int> chainleng; // vector with all lengths of chains requested
		
		//for each nmr model
		for (unsigned int i = 1; i <= max; i++)
		{
			chainleng.clear();
			
			//for each chain
			for ( unsigned int ch = 0; ch < totalchain.size(); ch++ )
			{
				PdbLoader pl(pdbfile2);
				pl.setChain(totalchain[ch]);
				pl.setNoHAtoms();
				pl.setNoVerbose();
				pl.setPermissive();
				if (!force)
				{		
					nmrmodel = i;	     
					pl.setModel(i);
				} 
				Spacer *sp;
				
				Protein prot;
                                prot.load(pl);
                                sp=prot.getSpacer(totalchain[ch]);
				//get internal array number from PDB aminoacid number
				if (START1 < 9999)
					START = sp->getIndexFromPdbNumber(START1);
				else
					START = 1;
				
				if (END1 < 9999)
					END = sp->getIndexFromPdbNumber(END1);
				else
					END = sp->sizeAmino()-1;
				
				if (!pl.isValid())
				{
					if (verbose)
					{
						cerr << "Warning: Invalid PDB file found "<<inputFile<<" Chain: "<<totalchain[ch];
						if (nmr)
							cerr <<" model: "<<nmrmodel<<" \n";
						else
							cerr <<" \n";
					}
					if ( ( i == max ) && ( ch == ( totalchain.size()-1 ) ))
						inputFile = "!";
					continue;
				}
				//CUTOFF length of chain not to be considered
				if ( sp->sizeAmino() <= CUTOFF )
				{
					cerr << "Warning: "<<inputFile<<" Chain: "<<totalchain[ch]<<" under CUTOFF\n";
					chainleng.push_back(0);
					continue;
				}
				else
					//insert chainlength for each chain 
					chainleng.push_back(sp->sizeAmino());
				
				if ( END == 0 )
					
					//number of all aminoacids in the spacer
					END = ( sp->sizeAmino()-1 );
				
				else if ( END >= sp->sizeAmino() )
					ERROR("END out of protein max length permitted..", exception);
				
				cout.setf(ios::fixed, ios::floatfield);
				
				// ------------------------------------
				// relative energy, i.e. TAP score case
				// ------------------------------------
				
				if (perres)
				{
					
					ERROR("This feature is currently unsupported.", exception);
					// NB: Have to fix TAP calculation for this case first. -- ST, 03/2007
					
					//for each residue
					for (unsigned int i = START; i < END; i++)
					{
						AminoAcidCode code = 
						aminoAcidThreeLetterTranslator(sp->getAmino(i).getType()); 
						int num = sp->getPdbNumberFromIndex(i);
						unsigned int j = i+1;
						tres = pot->calculateEnergy(*sp, i, j);
						Header( inputFile, caspVerbose, nmr, nmrmodel );
						cout << num <<"\t"<<tres/(-log(pot->pReturnMaxPropensities(code)));
						if ( verbose )
							cout<<"\t"<<totalchain[ch];
						cout <<endl;
					}
				}
				else //(!perres)
				{
					tres = tres + pot->calculateEnergy(*sp, START, END);
					
					for ( unsigned int i = START; i < END; i++)
					{
						AminoAcidCode code = 
						aminoAcidThreeLetterTranslator(sp->getAmino(i).getType());
						
						double tmpMin = (-log(pot->pReturnMinPropensities(code)));
						double tmpMax = (-log(pot->pReturnMaxPropensities(code)));
						
						totalMin += tmpMin;
						total += tmpMax;
						unsigned int j=i+1;
						double zero=pot->calculateEnergy(*sp,i, j) - tmpMin;
						if ( zero==0.0 )
						{
							cout <<"Zero detected" <<endl;
							resTap.push_back(zero);
						}
						else
						{
							// this is a quick hack to print out (hidden) per-residue values:
							double tmp = ( zero ) / ( tmpMax - tmpMin );
							resTap.push_back(tmp);
						}
					}
					
				}
				
				END = MEMORYEND;
			}
			
			//print out results in (!perres) case
			if (!perres)
			{
				tapRank.push_back(tres/total);
				models.push_back(inputFile);
				
				if (!rank)
				{
					Header( inputFile, caspVerbose, nmr, nmrmodel);
					
					if (rawOnly)
						tap = tres;
					else if (maxOnly)
						tap = tres / total;
					else
						tap = (tres - totalMin) / (total - totalMin);
					
					cout << setw(9) << setprecision(4) << tap << "\t";
					
					tapNmr.push_back(tap);
					
					if (showAccuracy)
						writeAccuracy(tap);
					
					if ( verbose )
					{
						if (maxOnly)
							cout << "max\t";
						
						if ( !allchains )
							cout <<totalchain[0]<<"\tsize: "<<chainleng[0] << "\n";
						else
						{
							cout <<"all chains\t size:";
							for ( unsigned int i = 0; i < chainleng.size(); i++ )
								cout <<" "<<totalchain[i]<<"("<<chainleng[i]<<")";
						}
					}
					else
						cout << "\n";
				}
			}
			totalMin = 0.00;
			total = 0.00;
			tres = 0.00;
			
		}
		
		pdbfile2.close();
		// reset variable to trigger break condition:
		inputFile = "!";
	}
	
	if (rank)
	{
		vector<double> ranks;
		vector<string> mods;
		
		for (unsigned int i = 0; i < tapRank.size(); i++)
		{
			ranks.push_back(tapRank[i]);
			mods.push_back(models[i]);
		}
		
		qsort(ranks, mods, 0, tapRank.size()-1);
		
		for (unsigned int pos = 0; pos < mods.size(); pos++)
		{
			cout <<  mods[pos] << ":" << "\t";
			cout << setw(9) << setprecision(4) << ranks[pos] << "\n" ;
		}
	}
	
	if (nmr)
    {
		double avgTap = average(tapNmr);
		double sdTap = standardDeviation(tapNmr, avgTap);
		double minTap = minimum(tapNmr);
		double maxTap = maximum(tapNmr);
		cout << "Average TAP = \t" << setw(9) << setprecision(4) << avgTap << "\tSD ="
			<< setw(9) << setprecision(4) << sdTap << "\tmin ="
			<< setw(9) << setprecision(4) << minTap << "\tmax ="
			<< setw(9) << setprecision(4) << maxTap << "\n";
    }
	else if (perresFile != "!") // write (hidden) prer-residue data only for X-ray structures
    {
		ofstream outFile(perresFile.c_str());
		if (!outFile)
			ERROR("File " + perresFile + " could not be created.", exception);
		
		for (unsigned int i = 0; i < resTap.size(); i++)
			outFile << setw(4) << (i+1) << "\t" << setw(9) << setprecision(4) << resTap[i] << "\n";
    }
	
	if (showAccuracy)
		sShowAccuracy(tap);
}
