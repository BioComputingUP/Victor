/* 	This file is part of Victor

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
 * @file TmScore.cc
 * @author Riccardo Zanella
 * @date Lug 2015
 * @version 0.1
 */

//Includes:
#include <TmScore.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <Debug.h>
#include <stdio.h>
#include <errno.h>

using namespace Victor::Mobi;
using namespace Victor::Biopool;
using namespace std;

const string TMTMP_OUT = "TMScore.pdb";

/**
 * @brief  Default Constructor
 * @param _binary (string) full path to binary TMScore file, must have execution permission
 * @param _tmp (string) full path to temp dir, must have write permission
 * @param _tmp bool verbose
 */
TmScore::TmScore(string _binary, string output, bool verbose) :
		binary(_binary), outputTmScore(output + TMTMP_OUT), verbose(verbose) {
}

//TmScore::TmScore(string _binary, string output, bool verbose) :
//		binary(_binary), tmp(output.substr(output.length() - 1, 1) == "/" ?
//						output : output + "/"), verbose(verbose) {
//}

TmScore::~TmScore() {
}

/**
 *
 * @param _binary (string) full path to binary TMScore file, must have execution permission
 * @param _tmp string, path of model file, must have read permission
 * @param _tmp string, path of native file, must have read permission
 * @return (Protein*), return pointer to rototraslated protein
 */
Protein* TmScore::TmImpose(string modelFile, string nativeFile) {

	Protein* prot = NULL;
	//if have read permission
	if (access(modelFile.c_str(), R_OK) == 0
			&& access(nativeFile.c_str(), R_OK) == 0) {

		//if binary have execute permission
		if (access(binary.c_str(), X_OK) == 0) {
			pid_t pid;
			//Pipeing TMS output
			int pipefd[2];
			pipe(pipefd);
			//Forking child for TMS
			pid = fork();
			if (pid < 0) {
				ERROR("Unable to fork child process", error);
			} else {
				if (pid == 0) {
					close(pipefd[0]);						//close input pipe
					if (dup2(pipefd[1], 1) < 0)
						ERROR("Unable to redirect output TMScore", error);//stdout
					if (dup2(pipefd[1], 2) < 0)
						ERROR("Unable to redirect error output TMScore", error);//stderr
					close(pipefd[1]);						//close output pipe
					if (execl(binary.c_str(), binary.c_str(), modelFile.c_str(),
							nativeFile.c_str(), "-o", outputTmScore.c_str(),
							NULL))
						ERROR("Unable to exec", error);
				} else {
					//Read output from child pipe

					//////////////////////////////////////////////////////////////////////////
					int status;

					wait(&status);

					remove(outputTmScore.c_str());

					if (access((outputTmScore + "_atm").c_str(), R_OK) != 0)
						ERROR("Cannot read pdb file to fix", exception);

					ifstream inFile((outputTmScore + "_atm").c_str());
					stringstream buffer;
					string line;

					while (inFile) {
						line = readLine(inFile);
						if (line.substr(0, 16) == "REMARK  TM-score") {
							buffer << line << endl;
							buffer << "MODEL        1" << endl;

						} else if (line.substr(0, 6) == "ATOM  ") {
							buffer << line << endl;
						} else if (line.substr(0, 3) == "TER") {
							buffer << line << endl;
							buffer << "ENDMDL" << endl;
							break;
						}
					}

					buffer.clear();

					//Load from memory buffer
					PdbLoader pl(buffer);

					pl.setNoVerbose();

					pl.setModel(1);
					pl.checkModel();

					prot = new Protein();

					pl.loadProtein(*prot);

				}
			}
		} else
			ERROR("No access to " + binary + "  binary!", exception);
	} else
		ERROR("No access to pdb files " + modelFile + " or " + nativeFile,
				exception);
	return prot;
}
