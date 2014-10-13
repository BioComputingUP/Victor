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

#ifndef _LOOPMODEL_H_
#define _LOOPMODEL_H_

// Includes:
#include <Debug.h>
#include <string.h>
#include <LoopTable.h>
#include <VectorTransformation.h>
#include <Spacer.h>
#include <SeqConstructor.h>
#include <set>
#include <ranking_helper.h>
#include <ranking_helper2.h>
#include <SolvationPotential.h>
#include <RapdfPotential.h>
#include <PhiPsi.h>
#include <ctype.h>
using namespace Victor;
using namespace Victor::Lobo;
using namespace Victor::Biopool;
using namespace Victor::Energy;
namespace Victor { namespace Lobo {

    // Global constants, typedefs, etc. (to avoid):

    class LoopModel {
    public:

        // CONSTRUCTORS/DESTRUCTOR:
        LoopModel();
        LoopModel(const LoopModel& orig);
        virtual ~LoopModel();

        // PREDICATES:
        vector<string>& getTableFileName();

        vector<int> vdwValues(Spacer& sp, unsigned int index1,
                unsigned int index2, vector<Spacer>& solVec);

        void rankRawScore(Spacer& sp, unsigned int index1, unsigned int index2,
                vector<Spacer>& solVec, unsigned int maxWrite = 9999);

        void doScatterPlot(Spacer& sp, unsigned int index1, unsigned int index2,
                vector<Spacer>& solVec, bool withOxygen = true);

        void refineModel(Spacer& sp, unsigned int index1, unsigned int index2,
                vector<Spacer>& solVec);

        void optimizeModel(Spacer& sp, unsigned int index1, unsigned int index2,
                vector<Spacer>& solVec, bool verbose = true);

        vector<double> rankRms2(Spacer& sp, unsigned int index1, unsigned int index2,
                vector<Spacer>& solVec);

        double calcOrigEnergy(const Spacer& sp, unsigned int index1,
                unsigned int index2);
        double calculateEndRms(const Spacer& sp, unsigned int index1,
                unsigned int index2, const Spacer& sp2);
        double calculateRms(const Spacer& sp, unsigned int index1,
                unsigned int index2, const Spacer& sp2,
                bool output = true, bool withOxygen = false);
        double calculateRms2(const Spacer& sp, unsigned int index1,
                unsigned int index2, const Spacer& sp2, bool output,
                bool withOxygen = false);
        void setStructure(Spacer& sp, Spacer& sp2, unsigned int index1,
                unsigned int index2);

        vector<int> consistencyValues(Spacer& sp, unsigned int index1,
                unsigned int index2, vector<Spacer>& solVec);

        double calculatePropensities(Spacer& sp);
        double calculatePropensities(const Spacer& sp, unsigned int index1,
                unsigned int index2, Spacer& sp2);

        double getENDRMS_WEIGHT(unsigned int len = 0);
        void saveENDRMS_WEIGHT(ostream& output);

        static string getSCWRLConservedSequence(const Spacer& sp,
                unsigned int index1, unsigned int index2);

        void defineLoopAnchors(Spacer& sp, unsigned int& index1,
                unsigned int& index2, bool isDeletion);

        // MODIFIERS:

        void setInterpolated() {
            pInter = true;
        }

        void setNotInterpolated() {
            pInter = false;
        }

        void setLimitedVerbose() {
            pVerbose = 1;
        }

        void setVerbose(unsigned int _verb = 10) {
            pVerbose = _verb;
        }

        void setQuiet() {
            pVerbose = 0;
            pPlot = false;
        }

        void setScatterPlot(ostream* _sc) {
            pPlot = true;
            pScatter = _sc;
        }

        void releaseTables(unsigned int index = 0); // releases memory occupied by loop tables
        void setTableFileName(string basename, string ending = ".lt");
        void setTableFileName(vector<string>& _tFN);
        vector<Spacer> createLoopModel(const AminoAcid& start,
                const vgVector3<double>& startN, const AminoAcid& end,
                const vgVector3<double>& endNAtom, unsigned int indexS,
                unsigned int indexE, unsigned int numLoops,
                unsigned int numLoops2, vector<string> typeVec);
        void clusterLoops(vector<Spacer>& vsp);

        void copy(const LoopModel& orig);

        int amino_amino_collision(AminoAcid& aa1, AminoAcid& aa2, int pos1,
                int pos2, double distance);
        double compactness(Spacer& loop, unsigned int index1, unsigned int index2,
                Spacer& proteine);

        double calculateEnergy(const Spacer& sp, unsigned int index1,
                unsigned int index2, Spacer& sp2);
        double calculateSecondaryPreference(const Spacer& sp, unsigned int index1,
                unsigned int index2, Spacer& sp2);
        double calculatePacking(const Spacer& sp, unsigned int index1,
                unsigned int index2, Spacer& sp2);
        double calculateSolvation(const Spacer& sp, unsigned int index1,
                unsigned int index2, Spacer& sp2);
        double calculateHydrogen(const Spacer& sp, unsigned int index1,
                unsigned int index2, Spacer& sp2);
        double calculateCompactness(const Spacer& sp, unsigned int index1,
                unsigned int index2, Spacer& sp2);
        double calculateFlankingConformation(const Spacer& sp, unsigned int index1,
                unsigned int index2, Spacer& sp2);
        double calculateConsistency(Spacer& sp, unsigned int index1,
                unsigned int index2, Spacer& sp2);

        void setENDRMS_WEIGHT(double val, unsigned int len = 0);
        void setAllENDRMS_WEIGHT(double val, unsigned int max = 13);
        void loadENDRMS_WEIGHT(istream& input);

        // OPERATORS:
        LoopModel& operator=(const LoopModel& orig);

        static unsigned int OPT_MAX1;
        static unsigned int OPT_MAX2;
        static unsigned int OPT_NUM;
        static double VDW_LIMIT;
        static double SIM_LIMIT;
        static double ENERGY_LIMIT;
        static double ENERGY_WEIGTH;
        static double SECPREF_WEIGTH;
        static double SECPREF_TOL;
        static double PACKING_WEIGTH;

        //ATTRIBUTES
        static unsigned int MAX_SPAN;
        static unsigned int MAX_ITER_SOL;
        static string TABLE_PARAM_FILE;

        static SolvationPotential solv;
        static RapdfPotential rapdf;
        static PhiPsi tor;

    protected:
        // HELPERS: 
        // converts the input coords to something processable:
        void convertCoords(AminoAcid start, vgVector3<float> startN, AminoAcid end,
                vgVector3<float> endN, LoopTableEntry& startEntry,
                LoopTableEntry& endEntry, VectorTransformation& vt,
                unsigned int nAmino);
        // tries to find a ring closure between to locations:
        bool ringClosureBase(const LoopTableEntry&, const LoopTableEntry&,
                unsigned int, unsigned int, double,
                VectorTransformation vt, unsigned int num,
                unsigned int depth,
                vector<vgVector3<float> >& partialSolution);
        bool ringClosure(const LoopTableEntry&, const LoopTableEntry&, unsigned int,
                unsigned int, double, VectorTransformation vt,
                vector<vgVector3<float> >& partialSolution,
                unsigned int currentSelection = 1);
        // claculates the actual coordinates for output:
        vector<Spacer> calculateLoop(const vgVector3<float>& startN,
                unsigned int length, vector<string> typeVec);

    private:
        // HELPERS: 
        void pAddRot(AminoAcid& start, vgVector3<float>& startN,
                AminoAcid& end, vgVector3<float>& endN,
                const vgMatrix3<float>& rotMat);
        void pAddTrans(AminoAcid& start, vgVector3<float>& startN,
                AminoAcid& end, vgVector3<float>& endN,
                const vgVector3<float>& trans);
        void loadTables(unsigned int nAmino);
        void pAminoAcidSetup(AminoAcid* aa, string type, double bfac = 90.0);
        void pSetSideChain(AminoAcid& aa);
        double pCalculateLoopRms(Spacer& sp1, Spacer& sp2);

        int loop_loop_vdw(Spacer& sp, unsigned int index1);
        int loop_spacer_vdw(Spacer& loop, unsigned int index1, unsigned int index2,
                Spacer& proteine);

        double sCalcEn(AminoAcid& aa, AminoAcid& aa2);
        double sICalcEn(AminoAcid& aa, AminoAcid& aa2);
        double sCalcEnInt(AminoAcid& aa, AminoAcid& aa2);
        double sICalcEnInt(AminoAcid& aa, AminoAcid& aa2);

        // ATTRIBUTES:
        bool pInter, pPlot;
        unsigned int pVerbose;
        ostream* pScatter;
        vector<LoopTable*> table;
        vector<string> tableFileName;
        vector<vgVector3<float> > solution;
        static unsigned int MAX_CHAIN_LENGTH;
        static double BOND_ANGLE_N_TO_CB;
        static double BOND_ANGLE_CB_TO_C;
        static double BOND_LENGTH_CA_TO_CB;

        vector<double> ENDRMS_WEIGTH;

    };

    // ---------------------------------------------------------------------------
    //                                  LoopModel
    // -----------------x-------------------x-------------------x-----------------


    // PREDICATES:

    inline vector<string>&
    LoopModel::getTableFileName() {
        return tableFileName;
    }

    // MODIFIERS:

    inline void
    LoopModel::setTableFileName(string basename, string ending) {
        //  PRECOND( _tFN.size() == MAX_CHAIN_LENGTH, exception);
        tableFileName.clear();
        for (unsigned int i = 0; i < MAX_CHAIN_LENGTH; i++)
            tableFileName.push_back(basename + (uitos(i)).c_str() + ending);
    }

    inline void
    LoopModel::setTableFileName(vector<string>& _tFN) {
        PRECOND(_tFN.size() == MAX_CHAIN_LENGTH, exception);
        tableFileName.clear();
        for (unsigned int i = 0; i < _tFN.size(); i++)
            tableFileName.push_back(_tFN[i]);
    }

    inline void
    LoopModel::pAddRot(AminoAcid& start, vgVector3<float>& startN,
            AminoAcid& end, vgVector3<float>& endN, const vgMatrix3<float>& rotMat) {
        for (unsigned int i = 0; i < start.sizeBackbone(); i++)
            start[i].setCoords(convert(rotMat) * start[i].getCoords());
        startN = rotMat * startN;
        for (unsigned int i = 0; i < end.sizeBackbone(); i++)
            end[i].setCoords(convert(rotMat) * end[i].getCoords());
        endN = rotMat * endN;
    }

    inline void
    LoopModel::pAddTrans(AminoAcid& start, vgVector3<float>& startN,
            AminoAcid& end, vgVector3<float>& endN, const vgVector3<float>& trans) {
        for (unsigned int i = 0; i < start.sizeBackbone(); i++)
            start[i].setCoords(start[i].getCoords() - convert(trans));
        for (unsigned int i = 0; i < end.sizeBackbone(); i++)
            end[i].setCoords(end[i].getCoords() - convert(trans));
        startN -= trans;
        endN -= trans;
    }

    inline double
    LoopModel::pCalculateLoopRms(Spacer& sp1, Spacer& sp2) {
        if (sp1.sizeAmino() != sp2.sizeAmino())
            ERROR("Arguments do not match.", exception);

        double res = 0.0;
        for (unsigned int i = 0; i < sp1.sizeAmino(); i++)
            res += sqr(sp1.getAmino(i)[N].distance(sp2.getAmino(i)[N]))
            + sqr(sp1.getAmino(i)[CA].distance(sp2.getAmino(i)[CA]))
            + sqr(sp1.getAmino(i)[C].distance(sp2.getAmino(i)[C]));

        return sqrt(res / (3 * sp1.sizeAmino()));
    }

    inline string
    LoopModel::getSCWRLConservedSequence(const Spacer& sp, unsigned int index1,
            unsigned int index2) {
        string tmp = "";

        for (unsigned int i = 0; i <= index1; i++)
            tmp += static_cast<char> (tolower(threeLetter2OneLetter(
                sp.getAmino(i).getType())));

        for (unsigned int i = index1 + 1; i <= index2; i++)
            tmp += threeLetter2OneLetter(sp.getAmino(i).getType());

        for (unsigned int i = index2 + 1; i < sp.sizeAmino(); i++)
            tmp += static_cast<char> (tolower(threeLetter2OneLetter(
                sp.getAmino(i).getType())));

        return tmp;
    }

}} // namespace

#endif //_LOOPMODEL_H_



