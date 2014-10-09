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

// --*- C++ -*------x-----------------------------------------------------------
//
//
// Description:     This program calculates sub-optimal alignments for two
//                  sequences. Alignments can be either global, local or
//                  free-shift; both with or without structural informations.
//
// -----------------x-----------------------------------------------------------

#include <ScoringS2S.h>
#include <ScoringP2S.h>
#include <ScoringP2P.h>
#include <HenikoffProfile.h>
#include <PSICProfile.h>
#include <SeqDivergenceProfile.h>
#include <CrossProduct.h>
#include <LogAverage.h>
#include <DotPFreq.h>
#include <DotPOdds.h>
#include <EDistance.h>
#include <Pearson.h>
#include <JensenShannon.h>
#include <AtchleyDistance.h>
#include <AtchleyCorrelation.h>
#include <NWAlign.h>
#include <SWAlign.h>
#include <FSAlign.h>
#include <SubMatrix.h>
#include <AGPFunction.h>
#include <VGPFunction.h>
#include <VGPFunction2.h>
#include <Sec.h>
#include <Ss2.h>
#include <Prof.h>
#include <Alignment.h>
#include <AlignmentBase.h>
#include <SequenceData.h>
#include <SecSequenceData.h>
#include <GetArg.h>
#include <iostream>
#include <ctime>


using namespace Victor::Align2;
using namespace Victor::Biopool;
using namespace Victor;



/// Show command line options and help text.

void
sShowHelp() {
    cout << "\nSUBOPTIMAL ALIGNMENT GENERATOR"
            << "\nThis program calculates sub-optimal alignments for two sequences. Alignments can be either global,"
            << "\nlocal or free-shift; both with or without structural informations. It is possible to select between"
            << "\nsequence-to-sequence, profile-to-sequence or profile-to-profile alignments. Profiles (optional)"
            << "\nare extracted from BLAST M4 output files, which can be further subjected to a number of sequence"
            << "\nsimilarity filters to optimize similarity.\n"
            << "\nOptions:"
            << "\n"
            << "\n * [--in <name>]     \t Name of input FASTA file"
            << "\n   [--pro1 <name>]   \t Name of target profile (psiBLAST M4 format) file"
            << "\n   [--pro2 <name>]   \t Name of template profile (psiBLAST M4 format) file"
            << "\n   [--out <name>]    \t Name of output FASTA file (default = to screen)"
            << "\n   [--fasta]         \t Use FASTA format to load profiles"
            << "\n   [-d <double>]     \t Min. master vs all seq. id. filter threshold (suggested = 0.00)"
            << "\n   [-D <double>]     \t Min. all vs all seq. id. filter threshold (suggested = 0.00)"
            << "\n   [-u <double>]     \t Max. master vs all seq. id. filter threshold (suggested = 1.00)"
            << "\n   [-U <double>]     \t Max. all vs all seq. id. filter threshold (suggested = 1.00)"
            << "\n   [--ws <0|1|2|3>]  \t Weighting scheme for profiles (default = 0, i.e. no weighting scheme)"
            << "\n                     \t --ws=0: No weighting scheme (default)."
            << "\n                     \t --ws=1: Calculate a frequency profile or PSSM using Henikoff weighting scheme."
            << "\n                     \t --ws=2: Calculate a frequency profile or PSSM using PSIC weighting scheme."
            << "\n                     \t --ws=3: Calculate a frequency profile or PSSM using SeqDivergence weighting scheme."
            << "\n   [--sf <0|...|8>]  \t Scoring function for profile-to-profile alignments (default = 1, i.e. LogAverage)."
            << "\n                     \t --sf=0: CrossProduct."
            << "\n                     \t --sf=1: LogAverage."
            << "\n                     \t --sf=2: DotPFreq."
            << "\n                     \t --sf=3: DotPOdds."
            << "\n                     \t --sf=4: EDistance."
            << "\n                     \t --sf=5: Pearson."
            << "\n                     \t --sf=6: JensenShannon."
            << "\n                     \t --sf=7: AtchleyDistance."
            << "\n                     \t --sf=8: AtchleyCorrelation."
            << "\n"
            << "\n   [--global]        \t Needleman-Wunsch global alignment (default)"
            << "\n   [--local]         \t Smith-Waterman local alignment"
            << "\n   [--freeshift]     \t Free-shift alignment"
            << "\n   [-n <int>]        \t Number of suboptimal alignments (default = 1)"
            << "\n   [-p <double>]     \t Penalty multiplier for suboptimal alignments (default = 1.00)"
            << "\n   [-a <double>]     \t Penalty subtractor for suboptimal alignments (default = 1.00)"
            << "\n"
            << "\n   [-m <name>]       \t Name of substitution matrix file (default = blosum62.dat)"
            << "\n   [-M <name>]       \t Name of structural substitution matrix file (default = secid.dat)"
            << "\n"
            << "\n   [--gf <0|1|2>]    \t Gap function (default = 0, i.e. AGP)"
            << "\n                     \t --gf=0: AGP (Affine Gap Penalty) function (default)."
            << "\n                     \t --gf=1: VGP (Variable Gap Penalty) function."
            << "\n                     \t --gf=2: VGP2 (Variable Gap Penalty) function."
            << "\n   [-o <double>]     \t Open gap penalty (default = 12.00)"
            << "\n   [-e <double>]     \t Extension gap penalty (default = 3.00)"
            << "\n   [--eType <0|1|2>] \t Extension gap type (default = 0, i.e. constant)"
            << "\n                     \t --eType=0: constant (default)."
            << "\n                     \t --eType=1: One sort of fraction value."
            << "\n                     \t --eType=2: We use a power function."
            << "\n   [--pdb <name>]    \t Name of template PDB file"
            << "\n   [--c <id>    ]    \t Chain identifier to read(default is first chain)"
            << "\n   [--wH <double>]   \t Weight for helical content (default = 1.00)"
            << "\n   [--wS <double>]   \t Weight for strand content (default = 1.00)"
            << "\n   [--wB <double>]   \t Weight for solvent accessibility (default = 1.00)"
            << "\n   [--wC <double>]   \t Weight for backbone straightness (default = 1.00)"
            << "\n   [--wD <double>]   \t Weight for space proximity (default = 1.00)"
            << "\n"
            << "\n   [--str <0|1|2|3>] \t Structural informations type (default = 0, i.e. no structural informations)"
            << "\n                     \t --str=0 No structural information (default)."
            << "\n                     \t --str=1 Calculate structural scores with info derived from secondary structure."
            << "\n                     \t --str=2 Calculate structural scores with info derived from PSI-PRED."
            << "\n                     \t --str=3 Calculate structural scores with info derived from PHD."
            << "\n   [--sec <name>]    \t Name of secondary structure FASTA file"
            << "\n   [--psi1 <name>]   \t Name of SS2 file for target sequence"
            << "\n   [--psi2 <name>]   \t Name of SS2 file for template sequence"
            << "\n   [--prof1 <name>]  \t Name of PROF file for target sequence"
            << "\n   [--prof2 <name>]  \t Name of PROF file for template sequence"
            << "\n   [--cSeq <double>] \t Coefficient for sequence alignment (default = 0.80)"
            << "\n   [--cStr <double>] \t Coefficient for structural alignment (default = 0.20)"
            << "\n"
            << "\n   [--verbose]       \t Verbose mode"
            << "\n" << endl;
}

int
main(int argc, char **argv) {
    string inputFileName, pro1FileName, pro2FileName, outputFileName, matrixFileName, matrixStrFileName;
    string secFileName, psi1FileName, psi2FileName, prof1FileName, prof2FileName, pdbFileName, chainID;
    string seq1Name, seq2Name, seq1, seq2, sec1, sec2;
    double downs, downa, ups, upa;
    double suboptPenaltyMul, suboptPenaltyAdd;
    double openGapPenalty, extensionGapPenalty;
    double weightHelix, weightStrand, weightBuried, weightStraight, weightSpace;
    double cSeq, cStr;
    unsigned int weightingScheme, scoringFunction, suboptNum, gapFunction, extensionType, structure;
    bool fasta, global, local, freeshift, verbose;
    struct tm* newtime;
    time_t t;

    // --------------------------------------------------
    // 0. Treat options
    // --------------------------------------------------

    if (getArg("h", argc, argv)) {
        sShowHelp();
        return 1;
    }

    getArg("-in", inputFileName, argc, argv, "!");
    getArg("-pro1", pro1FileName, argc, argv, "!");
    getArg("-pro2", pro2FileName, argc, argv, "!");
    getArg("-out", outputFileName, argc, argv, "!");
    fasta = getArg("-fasta", argc, argv);
    getArg("d", downs, argc, argv, 999.9);
    getArg("D", downa, argc, argv, 999.9);
    getArg("u", ups, argc, argv, 999.9);
    getArg("U", upa, argc, argv, 999.9);
    getArg("-ws", weightingScheme, argc, argv, 0);
    getArg("-sf", scoringFunction, argc, argv, 1);

    global = getArg("-global", argc, argv);
    local = getArg("-local", argc, argv);
    freeshift = getArg("-freeshift", argc, argv);
    if (!local && !freeshift)
        global = true;
    getArg("n", suboptNum, argc, argv, 1);
    getArg("p", suboptPenaltyMul, argc, argv, 1.00);
    getArg("a", suboptPenaltyAdd, argc, argv, 1.00);

    getArg("m", matrixFileName, argc, argv, "blosum62.dat");
    getArg("M", matrixStrFileName, argc, argv, "secid.dat");

    getArg("-gf", gapFunction, argc, argv, 0);
    getArg("o", openGapPenalty, argc, argv, 12.00);
    getArg("e", extensionGapPenalty, argc, argv, 3.00);
    getArg("-eType", extensionType, argc, argv, 0);
    getArg("-pdb", pdbFileName, argc, argv, "!");
    getArg("-c", chainID, argc, argv, " ");
    getArg("-wH", weightHelix, argc, argv, 1.00);
    getArg("-wS", weightStrand, argc, argv, 1.00);
    getArg("-wB", weightBuried, argc, argv, 1.00);
    getArg("-wC", weightStraight, argc, argv, 1.00);
    getArg("-wD", weightSpace, argc, argv, 1.00);

    getArg("-str", structure, argc, argv, 0);
    getArg("-sec", secFileName, argc, argv, "!");
    getArg("-psi1", psi1FileName, argc, argv, "!");
    getArg("-psi2", psi2FileName, argc, argv, "!");
    getArg("-prof1", prof1FileName, argc, argv, "!");
    getArg("-prof2", prof2FileName, argc, argv, "!");
    getArg("-cSeq", cSeq, argc, argv, 0.80);
    getArg("-cStr", cStr, argc, argv, 0.20);

    verbose = getArg("-verbose", argc, argv);


    // --------------------------------------------------
    // 1. Load data
    // --------------------------------------------------

    string path = getenv("VICTOR_ROOT");
    if (path.length() < 3)
        cout << "Warning: environment variable VICTOR_ROOT is not set." << endl;

    string examplesPath;

    string dataPath = path + "data/";

    if (inputFileName != "!") {
        inputFileName = examplesPath + inputFileName;
        ifstream inputFile(inputFileName.c_str());
        if (!inputFile)
            ERROR("Error opening input FASTA file.", exception);
        Alignment ali;
        ali.loadFasta(inputFile);
        if (ali.size() < 1)
            ERROR("Input FASTA file must contain two sequences.", exception);
        seq1Name = ali.getTargetName();
        seq2Name = ali.getTemplateName();
        seq1 = Alignment::getPureSequence(ali.getTarget());
        seq2 = Alignment::getPureSequence(ali.getTemplate());

    } else
        ERROR("subali needs input FASTA file.", exception);


    matrixFileName = dataPath + matrixFileName;
    ifstream matrixFile(matrixFileName.c_str());
    if (!matrixFile)
        ERROR("Error opening substitution matrix file.", exception);

    matrixStrFileName = dataPath + matrixStrFileName;
    ifstream matrixStrFile(matrixStrFileName.c_str());
    if (!matrixStrFile)
        ERROR("Error opening structural substitution matrix file.", exception);


    if (pdbFileName != "!")
        pdbFileName = examplesPath + pdbFileName;


    if (secFileName != "!") {
        secFileName = examplesPath + secFileName;
        ifstream secFile(secFileName.c_str());
        if (!secFile)
            ERROR("Error opening secondary structure FASTA file.", exception);
        Alignment aliSec;
        aliSec.loadFasta(secFile);
        if (aliSec.size() < 1)
            ERROR("Secondary structure FASTA file must contain two sequences.", exception);
        sec1 = Alignment::getPureSequence(aliSec.getTarget());
        sec2 = Alignment::getPureSequence(aliSec.getTemplate());

    }


    Ss2Input *psipred1 = 0;
    if (psi1FileName != "!") {
        psi1FileName = examplesPath + psi1FileName;
        ifstream psi1File(psi1FileName.c_str());
        if (!psi1File)
            ERROR("Error opening SS2 file for target sequence.", exception);
        psipred1 = new Ss2Input(psi1File);
    }

    Ss2Input *psipred2 = 0;
    if (psi2FileName != "!") {
        psi2FileName = examplesPath + psi2FileName;
        ifstream psi2File(psi2FileName.c_str());
        if (!psi2File)
            ERROR("Error opening SS2 file for template sequence.", exception);
        psipred2 = new Ss2Input(psi2File);
    }


    ProfInput *phd1 = 0;
    if (prof1FileName != "!") {
        prof1FileName = examplesPath + prof1FileName;
        ifstream prof1File(prof1FileName.c_str());
        if (!prof1File)
            ERROR("Error opening PROF file for target sequence.", exception);
        phd1 = new ProfInput(prof1File);
    }

    ProfInput *phd2 = 0;
    if (prof2FileName != "!") {
        prof2FileName = examplesPath + prof2FileName;
        ifstream prof2File(prof2FileName.c_str());
        if (!prof2File)
            ERROR("Error opening PROF file for template sequence.", exception);
        phd2 = new ProfInput(prof2File);
    }


    // --------------------------------------------------
    // 2. Output test data
    // --------------------------------------------------

    if (verbose) {
        fillLine(cout);
        if (secFileName != "!")
            cout << "\nTarget sequence:\n" << seq1
                << "\n\nTarget secondary structure:\n" << sec1
                << "\n\nTemplate sequence:\n" << seq2
                << "\n\nTemplate secondary structure:\n" << sec2
                << "\n" << endl;
        else
            cout << "\nTarget sequence:\n" << seq1
                << "\n\nTemplate sequence:\n" << seq2
                << "\n" << endl;
    }
    fillLine(cout);


    // --------------------------------------------------
    // 3. Select alignment mode
    // --------------------------------------------------

    SubMatrix sub(matrixFile);
    SubMatrix subStr(matrixStrFile);

    AlignmentData *ad;
    Structure *str;
    ScoringFunction *sf;
    ScoringScheme *ss;
    GapFunction *gf;

    switch (structure) {
        case 1:
            ad = new SecSequenceData(4, seq1, seq2, sec1, sec2, seq1Name, seq2Name);
            str = new Sec(&subStr, ad, cStr);

            break;
        case 2:
            ad = new SecSequenceData(4, seq1, seq2, sec1, sec2, seq1Name, seq2Name);
            str = new Ss2(&subStr, ad, psipred1, psipred2, cStr);
            break;
        case 3:
            ad = new SecSequenceData(4, seq1, seq2, sec1, sec2, seq1Name, seq2Name);
            str = new Prof(&subStr, phd1, phd2, cStr);
            break;
        default:
            ad = new SequenceData(2, seq1, seq2, seq1Name, seq2Name);
            str = 0;
            cSeq = 1.00;
            break;
    }


    if (pro1FileName != "!") {
        // Construct first profile

        pro1FileName = examplesPath + pro1FileName;
        ifstream pro1File(pro1FileName.c_str());
        if (!pro1File)
            ERROR("Error opening target profile (BLAST M6 format) file.", exception);

        Alignment ali1;
        if (fasta)
            ali1.loadFasta(pro1File);
        else
            ali1.loadPsiBlastMode4(pro1File);
        if (downs < 999.9)
            ali1.RemoveLowerSimple(downs);
        else
            if (downa < 999.9)
            ali1.RemoveLowerAll(downa);
        if (ups < 999.9)
            ali1.RemoveUpperSimple(ups);
        else
            if (upa < 999.9)
            ali1.RemoveUpperAll(upa);


        Profile *pro1;
        switch (weightingScheme) {
            case 1:
                pro1 = new HenikoffProfile();
                cout << "switch weightingScheme: Henikoff\n";
                break;
            case 2:
                pro1 = new PSICProfile();
                cout << "switch weightingScheme: PSICProfile\n";
                break;
            case 3:
                pro1 = new SeqDivergenceProfile();
                break;
            default:
                pro1 = new Profile();
                cout << "switch weightingScheme: no scheme\n";
                break;
        }
        pro1->setProfile(ali1);
        time(&t);
        newtime = localtime(&t);
        cout << "weightingScheme setted on prof1 " << newtime->tm_hour << "/" << newtime->tm_min << endl;

        if (pro2FileName != "!") {
            // --------------------------------------------------
            // 3.1. Profile-to-Profile case
            // --------------------------------------------------

            // Construct second profile

            pro2FileName = examplesPath + pro2FileName;
            ifstream pro2File(pro2FileName.c_str());
            if (!pro2File)
                ERROR("Error opening template profile (BLAST M6 format) file.", exception);

            Alignment ali2;
            if (fasta)
                ali2.loadFasta(pro2File);
            else
                ali2.loadPsiBlastMode4(pro2File);
            if (downs < 999.9)
                ali2.RemoveLowerSimple(downs);
            else
                if (downa < 999.9)
                ali2.RemoveLowerAll(downa);
            if (ups < 999.9)
                ali2.RemoveUpperSimple(ups);
            else
                if (upa < 999.9)
                ali2.RemoveUpperAll(upa);

            Profile *pro2;
            switch (weightingScheme) {
                case 1:
                    pro2 = new HenikoffProfile();
                    cout << "switch weightingScheme: Henikoff\n";
                    break;
                case 2:
                    pro2 = new PSICProfile();
                    cout << "switch weightingScheme: PSICProfile\n";
                    break;
                case 3:
                    pro2 = new SeqDivergenceProfile();
                    cout << "switch weightingScheme: SeqDivergenceProfile\n";
                    break;
                default:
                    pro2 = new Profile();
                    cout << "switch weightingScheme: no scheme\n";
                    break;
            }
            pro2->setProfile(ali2);
            time(&t);
            newtime = localtime(&t);
            cout << "weightingScheme setted on prof2 " << newtime->tm_hour << "/" << newtime->tm_min << endl;
            switch (scoringFunction) {
                case 1:
                    sf = new LogAverage(&sub, pro1, pro2);
                    cout << "switch scoring function: LogAverage\n";
                    break;
                case 2:
                    sf = new DotPFreq(pro1, pro2);
                    cout << "switch scoring function: DotPFreq\n";
                    break;
                case 3:
                    sf = new DotPOdds(pro1, pro2);
                    cout << "switch scoring function: DotPOdds\n";
                    break;
                case 4:
                    sf = new EDistance(pro1, pro2);
                    cout << "switch scoring function: EDistance\n";
                    break;
                case 5:
                    sf = new Pearson(pro1, pro2);
                    cout << "switch scoring function: Pearson\n";
                    break;
                case 6:
                    sf = new JensenShannon(pro1, pro2);
                    cout << "switch scoring function: JensenShannon\n";
                    break;
                case 7:
                    sf = new AtchleyDistance(pro1, pro2);
                    cout << "switch scoring function: AtchleyDistance\n";
                    break;
                case 8:
                    sf = new AtchleyCorrelation(pro1, pro2);
                    cout << "switch scoring function: AtchleyCorrelation\n";
                    break;
                default:
                    sf = new CrossProduct(&sub, pro1, pro2);
                    cout << "switch scoring function: CrossProduct\n";
                    break;
            }
            ss = new ScoringP2P(&sub, ad, str, pro1, pro2, sf, cSeq);
        } else
            ss = new ScoringP2S(&sub, ad, str, pro1, cSeq);
    } else
        ss = new ScoringS2S(&sub, ad, str, cSeq);


    switch (gapFunction) {
        case 1:
            gf = new VGPFunction(pdbFileName, chainID, openGapPenalty, extensionGapPenalty,
                    extensionType, weightHelix, weightStrand, weightBuried, weightStraight,
                    weightSpace);
            cout << "switch gapfunction: VGPFunction\n";
            break;
        case 2:
            gf = new VGPFunction2(secFileName, openGapPenalty, extensionGapPenalty,
                    extensionType, weightHelix, weightStrand);
            cout << "switch gapfunction: VGPFunction(no pdb)\n";
            break;
        default:
            gf = new AGPFunction(openGapPenalty, extensionGapPenalty);
            cout << "switch gapfunction: AGPFunction\n";
            break;
    }

    // --------------------------------------------------
    // 4. Calculate alignments
    // --------------------------------------------------

    Align *a;

    if (global) {
        cout << "\nSuboptimal Needleman-Wunsch alignments:\n" << endl;
        a = new NWAlign(ad, gf, ss);
    } else
        if (local) {
        cout << "\nSuboptimal Smith-Waterman alignments:\n" << endl;
        a = new SWAlign(ad, gf, ss);
    } else {
        cout << "\nSuboptimal free-shift alignments:\n" << endl;
        try {
            a = new FSAlign(ad, gf, ss);
        } catch (const char* a) {
            cout << "FSAlign error!\n";
        }
    }
    time(&t);
    newtime = localtime(&t);
    cout << "object FSAlign created " << newtime->tm_hour << "/" << newtime->tm_min << endl;

    // --------------------------------------------------
    // 5. Output alignments
    // --------------------------------------------------
    cout << "Preparing alignments...\n";
    a->setPenalties(suboptPenaltyMul, suboptPenaltyAdd);
    vector<Alignment> a2 = a->generateMultiMatch(suboptNum);
    if (a2.size() == 0)
        ERROR("No output alignments generated.", exception);
    a2[0].cutTemplate(1);
    Alignment a3 = a2[0];
    for (unsigned int i = 1; i < a2.size(); i++) {
        a2[i].cutTemplate(1);
        a3.addAlignment(a2[i]);
    }

    if (outputFileName != "!") {
        outputFileName = examplesPath + outputFileName;
        ofstream outputFile(outputFileName.c_str());
        if (!outputFile)
            ERROR("Error opening output FASTA file.", exception);
        cout << "Saving output to FASTA file: " << outputFileName << endl;
        a3.saveFasta(outputFile);
    } else
        a3.saveFasta(cout);
    cout << endl;
    fillLine(cout);

    time(&t);
    newtime = localtime(&t);
    cout << "Done, subali step is finished " << newtime->tm_hour << "/" << newtime->tm_min << "\n\n";
    return 0;
}
