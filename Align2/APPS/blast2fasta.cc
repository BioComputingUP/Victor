// --*- C++ -*------x-----------------------------------------------------------
//
// $Id: blast2fasta.cc,v 1.7 2008-08-03 21:26:33 biocomp Exp $
//
// Author:          Enrico Negri, Silvio Tosatto, Eckart Bindewald
//
// Project name:    Victor/Align
//
// Date:            01/2008
//
// Description:     Format converter from BLAST M6 to FASTA.
//
// -----------------x-----------------------------------------------------------

#include <Alignment.h>
#include <AlignmentBase.h>
#include <GetArg.h>
#include <IoTools.h>
#include <string>


using namespace Biopool;


/// Show command line options and help text.
void
sShowHelp()
{
	cout << "\nBLAST TO FASTA"
	     << "\nFormat converter from BLAST M6 to FASTA.\n"
	     << "\nOptions:"
	     << "\n * [--in <name>]     \t Name of input BLAST mode 6 file"
	     << "\n * [--out <name>]    \t Name of output FASTA file"
	     << "\n   [--eVal <double>] \t E-value significance threshold (default = 0.01)"
	     << "\n   [-t <double>]     \t % threshold for acceptance, relative to max score (default = 0.50 -> 50%)"
	     << "\n   [--max <int>]     \t Maximum number of sequences for alignment (default = all)"
	     << "\n" << endl;
}


int
main(int argc, char **argv)
{
	string inputFileName, outputFileName;
	double eVal, threshold;
	unsigned int max;


	// --------------------------------------------------
	// 0. Treat options
	// --------------------------------------------------

	if (getArg("h", argc, argv))
	{
		sShowHelp();
		return 1;
	}

	getArg("-in", inputFileName, argc, argv, "!");
	getArg("-out", outputFileName, argc, argv, "!");
	getArg("-eVal", eVal, argc, argv, 0.01);
	getArg("t", threshold, argc, argv, 0.50);
	getArg("-max", max, argc, argv, 9999);


	// --------------------------------------------------
	// 1. Read alignment
	// --------------------------------------------------

	string examplesPath = getenv("VICTOR_ROOT");
	if (examplesPath.length() < 3)
		cout << "Warning: environment variable VICTOR_ROOT is not set." << endl;
	//examplesPath += "examples/";
	examplesPath="";

	Alignment ali;
	if (inputFileName != "!")
	{
		inputFileName = examplesPath + inputFileName;
		ifstream inputFile(inputFileName.c_str());
		if (!inputFile)
			ERROR("Error opening input BLAST mode 6 file.", exception);
		ali.loadBlastMode6(inputFile);
	}
	else
		ERROR("blast2fasta needs input BLAST mode 6 file.", exception);


	// --------------------------------------------------
	// 2. Edit alignment
	// --------------------------------------------------

	if (ali.size() == 0)
		ERROR("No sequences in input BLAST mode 6 file.", exception);

	double maxScore = ali.getScore(0);

	for (unsigned int i = 0; i < ali.size(); i++)
		if (ali.getEvalue(i) > eVal)
		{
			ali.cutTemplate(i);
			break;
		}

	for (unsigned int i = 0; i < ali.size(); i++)
		if (ali.getScore(i) < (maxScore * threshold))
		{
			ali.cutTemplate(i);
			break;
		}

	if ((max < 9999) && (ali.size() > max))
		ali.cutTemplate(max);


	// --------------------------------------------------
	// 3. Output alignment statistics
	// --------------------------------------------------

	fillLine(cout);
	cout << "\n   PDB Code   Bit Score   E-value   Conservation (%)\n" << endl;
	for (unsigned int i = 0; i < ali.size(); i++)
	{
		double conservation = 0;
		for (unsigned int j = 0; j < ali.getLength(); j++)
			if (ali.isConserved(j, i))
				conservation++;

		cout << setw(11) << ali.getTemplateName(i)
		     << setw(12) << ali.getScore(i)
		     << setw(10) << ali.getEvalue(i)
		     << setw(19) << setprecision(3)
		     << (ali.calculatePairwiseIdentity(ali.getTarget(), ali.getTemplate(i)) * 100)
		     << endl;
	}
	cout << endl;
	fillLine(cout);


	// --------------------------------------------------
	// 4. Write alignment to disk
	// --------------------------------------------------

	if (outputFileName != "!")
	{
		outputFileName = examplesPath + outputFileName;
		ofstream outputFile(outputFileName.c_str());
		if (!outputFile)
			ERROR("Error creating output FASTA file.", exception);
		ali.saveFasta(outputFile);
	}
	else
		ERROR("blast2fasta needs output FASTA file.", exception);

	return 0;
}
