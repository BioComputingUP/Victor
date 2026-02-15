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
@Description */
#include <string>
#include <GetArg.h>
#include <CIFLoader.h>
#include <CIFSaver.h>
#include <vector3.h>
#include <matrix3.h>
#include <IntCoordConverter.h>

using namespace Victor;using namespace Victor::Biopool;

const double LAMBDA = 1.5; // minimum distance between neighbouring CAs

void sShowHelp() {
    cout << "CIF Shifter\n"
            << "Allows to move all residues in file by fixed offset.\n"
            << " Options: \n"
            << "\t-i <filename> \t\t Input CIF file\n"
            << "\t-o <filename> \t\t Output CIF file\n"
            << "\t[-r <repeat>] \t\t Repeat unit length (default = no rotation)\n"
            << "\t[-l <double>] \t\t Angle lambda factor (default = 0.1)\n"
            << "\t[-s <number>] \t\t Start residue of fragment (default = first)\n"
            << "\t[-e <number>] \t\t End residue of fragment (default = last)\n"
            << "\n";
}

void sAddLine() {
    cout << "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n";
}

int main(int nArgs, char* argv[]) {
    // --------------------------------------------------
    // 0. treat options
    // --------------------------------------------------

    if (getArg("h", nArgs, argv)) {
        sShowHelp();
        return 1;
    };

    string inputFile, outputFile;
    unsigned int startOffset, endOffset, repeatLength;
    double lambdaAngle;

    getArg("i", inputFile, nArgs, argv, "!");
    getArg("o", outputFile, nArgs, argv, "!");
    getArg("s", startOffset, nArgs, argv, 0);
    getArg("e", endOffset, nArgs, argv, 9999);
    getArg("r", repeatLength, nArgs, argv, 9999);
    getArg("l", lambdaAngle, nArgs, argv, 0.1);

    vgVector3<double> transOff;
    for (unsigned int i = 0; i < 3; i++)
        transOff[i] = 0.0;

    if ((inputFile == "!") || (outputFile == "!")) {
        cout << "Missing file specification. Aborting. (-h for help)" << endl;
        return -1;
    }

    // --------------------------------------------------
    // 1. read structure
    // --------------------------------------------------

    Spacer sp;

    ifstream inFile(inputFile.c_str());

    if (!inFile)
        ERROR("File does not exist.\n", exception);

    CIFLoader CIFL(inFile);

    sp.load(CIFL);
    inFile.close();


    endOffset = sp.getIndexFromPdbNumber(endOffset);
    if (startOffset > 0)
        startOffset = sp.getIndexFromPdbNumber(startOffset);

    // --------------------------------------------------
    // 2. rotate spacer
    // --------------------------------------------------

    if (repeatLength < 9999) {
        IntCoordConverter icc;

        // find offset
        vgVector3<double> firstA = sp.getAmino(endOffset
                - (3 * repeatLength / 4))[CA].getCoords()
                - sp.getAmino(endOffset)[CA].getCoords();

        vgVector3<double> firstB = sp.getAmino(endOffset
                - (1 * repeatLength / 4))[CA].getCoords()
                - sp.getAmino(endOffset)[CA].getCoords();

        vgVector3<double> firstNorm = (firstA.normalize()).cross(firstB.normalize());

        vgVector3<double> secondA = sp.getAmino(endOffset + repeatLength
                - (3 * repeatLength / 4))[CA].getCoords()
                - sp.getAmino(endOffset)[CA].getCoords();

        vgVector3<double> secondB = sp.getAmino(endOffset + repeatLength
                - (1 * repeatLength / 4))[CA].getCoords()
                - sp.getAmino(endOffset)[CA].getCoords();

        vgVector3<double> secondNorm = (secondA.normalize()).cross(secondB.normalize());

        firstNorm.normalize();
        secondNorm.normalize();

        double scalar = icc.getAngle(firstNorm, secondNorm) * lambdaAngle;

        cout << "Scalar = " << setw(5) << setprecision(3)
                << RAD2DEG * scalar / lambdaAngle << "\n";


        vgVector3<double> axis = (firstNorm.normalize()).cross(secondNorm.normalize());
        vgMatrix3<double> res = vgMatrix3<double>::createRotationMatrix(axis, scalar);

        sp.getAmino(startOffset)[N].addRot(res);
    }

    // --------------------------------------------------
    // 3. translate spacer
    // --------------------------------------------------

    if (endOffset < sp.sizeAmino() - 1) {
        sp.getAmino(endOffset).unbindOut(sp.getAmino(endOffset + 1));
        sp.getAmino(endOffset)[C].unbindOut(sp.getAmino(endOffset + 1)[N]);

        // find offset
        vgVector3<double> first = sp.getAmino(endOffset)[C].getCoords();
        vgVector3<double> second = sp.getAmino(endOffset + 1)[N].getCoords();

        double d = sp.getAmino(endOffset)[C].distance(
                sp.getAmino(endOffset + 1)[N]);

        double frac = (d - LAMBDA) / d;

        for (unsigned int i = 0; i < 3; i++)
            transOff[i] = (second[i] - first[i]) * frac;
    } else {
        if (startOffset != 0)
            ERROR("Both start and end offset are undefined.", exception);

        endOffset = sp.sizeAmino() - 1;

        // NB: unbind has to be reversed after moving the atoms if the model
        // is to be used further in the same program
        sp.getAmino(startOffset - 1).unbindOut(sp.getAmino(startOffset));
        sp.getAmino(startOffset - 1)[C].unbindOut(sp.getAmino(startOffset)[N]);

        // find offset
        vgVector3<double> first = sp.getAmino(startOffset + 1)[N].getCoords();
        vgVector3<double> second = sp.getAmino(startOffset)[C].getCoords();

        double d = sp.getAmino(startOffset + 1)[N].distance(
                sp.getAmino(startOffset)[C]);

        double frac = (d - LAMBDA) / d;

        for (unsigned int i = 0; i < 3; i++)
            transOff[i] = (second[i] - first[i]) * frac;
    }

    sp.getAmino(startOffset)[N].addTrans(transOff);

    // --------------------------------------------------
    // 4. write model to disk
    // --------------------------------------------------

    ofstream outFile(outputFile.c_str());
    if (!outFile)
        ERROR("File not found.", exception);
    CIFSaver ps(outFile);

    sp.save(ps);
    outFile.close();
}
