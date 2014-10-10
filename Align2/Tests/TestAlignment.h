/*
 * TestAlignment.cpp
 *
 *  Created on: Oct 6th, 2014
 *      Author: Manuel Giollo
 */

#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>

#include <AlignmentBase.h>
#include <Alignment.h>
using namespace std;
using namespace Victor;

class TestAlignment  : public CppUnit::TestFixture {
private:
    Alignment *testAlignment;
    
public:

    TestAlignment() : testAlignment(NULL) {
    }

    virtual ~TestAlignment() {
        delete testAlignment;
    }

    static CppUnit::Test *suite() {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestAlignmentBase");

        suiteOfTests->addTest(new CppUnit::TestCaller<TestAlignment>("Test1 - Loads alignment evaluates size.",
                &TestAlignment::testAlignment_A));
        suiteOfTests->addTest(new CppUnit::TestCaller<TestAlignment>("Test2 - Loads alignment using psiblast output as input.",
                &TestAlignment::testAlignment_B));
        suiteOfTests->addTest(new CppUnit::TestCaller<TestAlignment>("Test3 - loading alignments from psiblast output.",
                &TestAlignment::testAlignment_C));

        return suiteOfTests;
    }

    /// Setup method

    void setUp() {
    }

    /// Teardown method

    void tearDown() {
    }

protected:

    void testAlignment_A() {
        //Using the same sequence as target and template
        string path = getenv("VICTOR_ROOT");
        string dataPath = path + "Align2/Tests/data/";
        string inputFileName="test.fasta";
        inputFileName = dataPath + inputFileName;
        ifstream inputFile(inputFileName.c_str());
        if (!inputFile)
            ERROR("Error opening input FASTA file.", exception);
        Alignment ali;
        ali.loadFasta(inputFile);
        CPPUNIT_ASSERT((ali.getTargetName()== ali.getTemplateName())&&(ali.getTarget()==ali.getTemplate()));
    }

    void testAlignment_B() {
        //Using psiblast output as input
        string path = getenv("VICTOR_ROOT");
        string dataPath = path + "Align2/Tests/data/";
        string inputFileName="target.out";
        inputFileName = dataPath + inputFileName;
        ifstream inputFile(inputFileName.c_str());
        if (!inputFile)
            ERROR("Error opening input FASTA file.", exception);
        Alignment ali;
        ali.loadPsiBlastMode4(inputFile);
        //from the input file
        //ref|ZP_06622438.1|  phosphopyruvate hydratase [Turicibacter sp. P...    912   0.0  
        //ref|YP_001664731.1|  phosphopyruvate hydratase [Thermoanaerobacte...    907   0.0  
        CPPUNIT_ASSERT((ali.getScore(0)==912)&&(907==ali.getScore(1))&& (ali.getEvalue(0)==ali.getEvalue(1))); 
    }

    void testAlignment_C() {
        //Using psiblast output as input to load the alignments
        string path = getenv("VICTOR_ROOT");
        string dataPath = path + "Align2/Tests/data/";
        string inputFileName="target.out";
        inputFileName = dataPath + inputFileName;
        ifstream inputFile(inputFileName.c_str());
        if (!inputFile)
            ERROR("Error opening input FASTA file.", exception);
        Alignment ali;
        ali.loadPsiBlastMode4(inputFile);
        //from the input file
        //Query_1       1    IVKIIGREIIDSRGNPTVEAEVHLEGGFVGMAAAPSGASTGSREALELRDGDKSRFLGKG  60
        //ZP_06622438   4    IVDVYAREVLDSRGNPTVEVEVTTESGSFGRALVPSGASTGIYEAVELRDGDKSRYLGKG  63
        //Query_1       61   VTKAVAAVNGPIAQALI--G--K-DAKDQAGIDKIMIDLDGTENKSKFGANAILAVSLAN  115
        //ZP_06622438   64   VLNAVKNVNDIIAPELV--G--M-DVTDQCGIDRLMIALDGTKNKGKLGANAILGVSMAV  118
        //Query_1       116  AKAAAAAKGMPLYEHIAELNGTPGKYSMPVPMMNIINGGEHADNNVDIQEFMIQPVGAKT  175
        //ZP_06622438   119  AHAAADFVGLPLYRYLGGFNSKE----LPTPMMNIINGGEHADNNIDFQEFMIMPVGAPT  174
        CPPUNIT_ASSERT((ali.getTarget()=="IVKIIGREIIDSRGNPTVEAEVHLEGGFVGMAAAPSGASTGSREALELRDGDKSRFLGKGVTKAVAAVNGPIAQALI--G--K-DAKDQAGIDKIMIDLDGTENKSKFGANAILAVSLANAKAAAAAKGMPLYEHIAELNGTPGKYSMPVPMMNIINGGEHADNNVDIQEFMIQPVGAKTVKEAIRMGSEVFHHLAKVLKAKGM--NTA-VGDEGGYAP-NLGS-NAEALAVIAEAVKAAGYELG----------K----------DITLAMDCAASEF--Y--K---D--G--K-----Y---V---L-AG----E----------GNKAFTSEEFTHFLEELTKQY-P-IVSIEDGLDESDWDGFAYQTKVLGDKIQLVGDDLFVTNTKILKEGIEKGIANSILIKFNQIGSLTETLAAIKMAKDAGYTAVISHRSGETEDATIADLAVGTAAGQIKTGSMSRSDRVAKYNQLIRIEEALG-----EKAPYNGRKEIKG")
        && (ali.getTemplate()          =="IVDVYAREVLDSRGNPTVEVEVTTESGSFGRALVPSGASTGIYEAVELRDGDKSRYLGKGVLNAVKNVNDIIAPELV--G--M-DVTDQCGIDRLMIALDGTKNKGKLGANAILGVSMAVAHAAADFVGLPLYRYLGGFNSKE----LPTPMMNIINGGEHADNNIDFQEFMIMPVGAPTFKEAIRMGAEVFHALKSVLHGMGL--NTA-VGDEGGFAP-NLES-NEAAIKVILEAIEKAGYVPG----------K----------DVMIAMDVASSEF--Y--K---D--G--K-----Y---V---L-AG----E----------GGKVFTSEELCDFYAELCEKY-P-IISIEDGLDQDDWAGWDYLTKKIGDKVQLVGDDFFVTNTERLAEGIEKNVANSILIKVNQIGTLTETFEAIEMAKKAGYTAVVSHRSGETEDATIADIAVATNAGQIKTGSMSRTDRIAKYNQLLRIEDELG-----QQAVYNGVKSF--"));
        

    }

};
