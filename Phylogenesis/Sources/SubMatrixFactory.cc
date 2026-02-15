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

#include <SubMatrixFactory.h>
#include <SubMatrix.h>
#include <fstream>

namespace Victor {
namespace Phylogenesis {

SubMatrixFactory::SubMatrixFactory(bool positive) {

	string path = getenv("VICTOR_ROOT");
	path += + "data/";
	string matrixFileName;

	matrixFileName = path + "blosum30.dat";
	ifstream b30(matrixFileName.c_str());
	blosum30 = new SubMatrix(b30);

	matrixFileName = path + "blosum45.dat";
	ifstream b45(matrixFileName.c_str());
	blosum45 = new SubMatrix(b45);

	matrixFileName = path + "blosum62.dat";
	ifstream b62(matrixFileName.c_str());
	blosum62 = new SubMatrix(b62);

	matrixFileName = path + "blosum80.dat";
	ifstream b80(matrixFileName.c_str());
	blosum80 = new SubMatrix(b80);

	matrixFileName = path + "pam20.dat";
	ifstream p20(matrixFileName.c_str());
	pam20 = new SubMatrix(p20);

	matrixFileName = path + "pam60.dat";
	ifstream p60(matrixFileName.c_str());
	pam60 = new SubMatrix(p60);

	matrixFileName = path + "pam120.dat";
	ifstream p120(matrixFileName.c_str());
	pam120 = new SubMatrix(p120);

	matrixFileName = path + "pam350.dat";
	ifstream p350(matrixFileName.c_str());
	pam350 = new SubMatrix(p350);

	if(positive){
		rescore(blosum30);
		rescore(blosum45);
		rescore(blosum62, 2.0);
		rescore(blosum80);
		rescore(pam20);
		rescore(pam60);
		rescore(pam120);
		rescore(pam350);
	}
}

void SubMatrixFactory::rescore(SubMatrix* matrix, double scalar){
	string residues = matrix->getResidues();
	double res = residues.size()-2; //avoids '-'
	for(int i=0;i<=res;i++){
		char res1 = residues[i];
		char res1_bis = res1;
		if(res1>=65 && res1<=90){
			res1_bis = res1 + 32;
		}
		if(res1>=97 && res1<=122){
			res1_bis = res1 - 32;
		}
		for(int j=0;j<=i;j++){
			char res2 = residues[j];
			char res2_bis = res2;
			if(res2>=65 && res2<=90){
				res2_bis = res2 + 32;
			}
			if(res2>=97 && res2<=122){
				res2_bis = res2 - 32;
			}
    		double rescoredValue = matrix->score[res1][res2] - matrix->getMinScore();
    		matrix->score[res1][res2] = matrix->score[res2][res1] =
			matrix->score[res1][res2_bis] = matrix->score[res2_bis][res1] =
			matrix->score[res1_bis][res2] = matrix->score[res2][res1_bis] =
			matrix->score[res1_bis][res2_bis] = matrix->score[res2_bis][res1_bis] = rescoredValue * scalar;
		}
    }
}

SubMatrixFactory::~SubMatrixFactory() {
	delete blosum30;
	delete blosum45;
	delete blosum62;
	delete blosum80;
	delete pam20;
	delete pam60;
	delete pam120;
	delete pam350;
}

SubMatrix* SubMatrixFactory::getSubMatrix(SubMatrixName name){
	switch(name){
		case(BLOSUM30):
			return blosum30;
		case(BLOSUM45):
			return blosum45;
		case(BLOSUM62):
			return blosum62;
		case(BLOSUM80):
			return blosum80;
		case(PAM20):
			return pam20;
		case(PAM60):
			return pam60;
		case(PAM120):
			return pam120;
		case(PAM350):
			return pam350;
	}
	return NULL;
}

}}
