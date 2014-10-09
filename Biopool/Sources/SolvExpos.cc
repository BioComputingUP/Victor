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
 *@Description A module to determine the solvent exposure/accessibility of residues in a
 *                              protein fragment.
 */
#include <SolvExpos.h>
using namespace Victor; using namespace Victor::Biopool;
using namespace Victor::Biopool;

namespace Victor { namespace Biopool { 

    Atom& getReprAtom(AminoAcid &amino);

    unsigned int getNumNeighbours(Spacer &chain, const unsigned int tgt,
            const unsigned int start, const unsigned int end);

}} //namespace

/**
 *@Description Returns the atom corresponding to C, 
 * aminoacids have only a single possible, hard coded, open in-bond
 *@param AminoAcid reference 
 *@return the representative atom of that residue as to solvent exposure. Such atom
 *    is CA for Glycine and CB for all other amino acid types.
 */
Victor::Biopool::Atom& Victor::Biopool::getReprAtom(AminoAcid &amino) {
    if (amino.getCode() == XXX) throw "errore amino type XXX\n";

    // Glycine
    if (amino.getCode() == GLY) {
        if (amino.isMember(CA))
            return amino[CA];
        else
            ERROR("No CA defined for Glycine", exception);
    }
    // Other aminoacid
    if (amino.getSideChain().isMember(CB))
        return amino.getSideChain()[CB];
    else
        throw "No CB defined for residue. Missing atom\n";
}

/**
 *@Description in the case of a C-terminal fragment, 'end' is assumed to be the index of
 *    a fictitious residue after the entire protein chain. 
 *    a given residue is considered to be a neighbour of the target if and only
 *    if its representative atom is within a CUTOFF distance of the representative
 *    atom of the target.
 *    no assumptions are made on which is the representative atom of a given
 *    amino acid type.
 *    the target residue need not be internal to the fragment.
 *@param chain: a Spacer object representing an entire protein chain.
 *        tgt: index, in the chain, of a target residue.
 *       start: index, in the chain, of the first residue of a fragment internal
 *           to the chain.
 *       end: index, in the chain, of the first residue after the end of that
 *         fragment.
 *@return number of residues, in the fragment, that are neighbour to the target
    residue.
 */

unsigned int Victor::Biopool::getNumNeighbours(Spacer &chain, const unsigned int tgt,
        const unsigned int start, const unsigned int end) {
    const double CUTOFF = 10; // Angstrom

    AminoAcid& ta = chain.getAmino(tgt);
    Atom& tr = getReprAtom(ta);

    unsigned int tot = 0;
    for (unsigned int i = start; i < end; ++i) {
        AminoAcid &ia = chain.getAmino(i);
        try {
            Atom &ir = getReprAtom(ia);
            double dist = tr.distance(ir);
            if (dist <= CUTOFF)
                tot++;
        } catch (const char* err) {
            ;
        }

    }

    if ((tgt >= start) && (tgt < end))
        tot--;

    return tot;
}
SolvExpos::SolvExpos(){
}
SolvExpos::~SolvExpos(){
}
/**
 *@Description in the case of a C-terminal fragment, 'end' is assumed to be the index of
 *    a fictitious residue after the entire protein chain.
 *   the target residue need not be internal to the fragment.
 *   the solvent exposure of the target is CORE if its number of neighbours
 *   in the fragment is greater than NGB_MIN, and EXPOSED otherwise.
 *@param  chain: a Spacer object representing an entire protein chain.
 * tgt: index, in the chain, of a target residue.
 * start: index, in the chain, of the first residue of a fragment internal
 *           to the chain.
 * end: index, in the chain, of the first residue after the end of that
 *         fragment.
 *@return solvent exposure of the target residue with respect to the fragment.
 */
SolvExpos::SolvExposEnum SolvExpos::getSolvExpos(Spacer &chain, const unsigned int tgt,
        const unsigned int start, const unsigned int end) {
    const unsigned int NGB_MIN = 20;

    unsigned int ngb = getNumNeighbours(chain, tgt, start, end);

    return (ngb > NGB_MIN) ? CORE : EXPOSED;
}

/**
 *@Description when delimiting a C-terminal fragment, tgtE and envE are assumed to be
 *    the index of a fictitious residue after the entire protein chain.
 *    target residues need not be internal to the environment fragment. The
 *   target fragment can even be completely external to the environment
 *   fragment.
 * the solvent-exposure vector is allocated on the heap memory and is
 *    expected to be explicitly deleted by code outside this function.
 *@param 
 * chain: a Spacer object representing an entire protein chain.
 *tgtS: index, in the chain, of the first of a fragment of target residues,
 *          namely, residues whose solvent exposure is to be computed.
 *tgtE: index, in the chain, of the first residue after the end of the
 *          target fragment.
 * envS: index, in the chain, of the first of a fragment of environment
 *          residues, namely, residues on which the solvent exposure of
 *          target residues is to be computed.
 *   envE: index, in the chain, of the first residue after the end of the
 *          environment fragment.
 *@return value:
 * pointer to a vector containing the solvent exposure of each of the target
 *    residues with respect to the environment fragment. The i-th element of
 *    this vector is the solvent exposure of the i-th target residue
 *    (i = 0, ... , tgtNum-1).
 */
vector<SolvExpos::SolvExposEnum>* SolvExpos::getSolvExposVec(Spacer& chain,
        const unsigned int tgtS, const unsigned int tgtE,
        const unsigned int envS, const unsigned int envE) {
    const unsigned int tgtNum = tgtE - tgtS;
    vector<SolvExpos::SolvExposEnum>* seVec = new vector<SolvExpos::SolvExposEnum>(tgtNum);

    for (unsigned int t = 0; t < tgtNum; ++t)
        (*seVec)[t] = getSolvExpos(chain, tgtS + t, envS, envE);

    return seVec;
}

/**
 *@Description in the case of a C-terminal fragment, 'end' is assumed to be the index of
 *    a fictitious residue after the entire protein chain.
 *    the target residue need not be internal to the fragment.
 *    the solvent accessibility of the target residue is between 0 (if the
 *    number of its neighbours in the fragment is greater or equal than NGB_MAX)
 *    and 1 (otherwise).
 *  @param 
 *  chain: a Spacer object representing an entire protein chain.
 *    tgt: index, in the chain, of a target residue.
 *    start: index, in the chain, of the first residue of a fragment internal
 *           to the chain.
 *  end: index, in the chain, of the first residue after the end of that
 *         fragment.
 *@return 
 * 
 *    solvent accessibility of the target residue with respect to the fragment.
 */
double SolvExpos::getSolvAccess(Spacer &chain, unsigned int tgt,
        unsigned int start, unsigned int end) {
    const double NGB_MAX = 30;

    double ngb = (double) getNumNeighbours(chain, tgt, start, end);

    return (NGB_MAX - min(ngb, NGB_MAX)) / NGB_MAX;
}


// -*- C++ -*-----------------------------------------------------------------
//
// Function:        Biopool::getSolvAccessVec()
//
// Arguments:
//
//    chain: a Spacer object representing an entire protein chain.
//
//    tgtS: index, in the chain, of the first of a fragment of target residues,
//          namely, residues whose solvent exposure is to be computed.
//
//    tgtE: index, in the chain, of the first residue after the end of the
//          target fragment.
//
//    envS: index, in the chain, of the first of a fragment of environment
//          residues, namely, residues on which the solvent exposure of
//          target residues is to be computed.
//
//    envE: index, in the chain, of the first residue after the end of the
//          environment fragment.
//
// Return value:
//
//    a vector containing the solvent accessibility of each of the target
//    residues with respect to the environment fragment. The i-th element of
//    this vector is the solvent accessibility of the i-th target residue
//    (i = 0, ... , tgtNum-1).
//
// Notes:
//
//    when delimiting a C-terminal fragment, tgtE and envE are assumed to be
//    the index of a fictitious residue after the entire protein chain.
//
//    target residues need not be internal to the environment fragment. The
//    target fragment can even be completely external to the environment fragment.
//
// ---------------------------------------------------------------------------

vector<double> SolvExpos::getSolvAccessVec(Spacer &chain, unsigned int tgtS,
        unsigned int tgtE, unsigned int envS, unsigned int envE) {
    const unsigned int tgtNum = tgtE - tgtS;

    vector<double> seVec(tgtNum);

    for (unsigned int t = 0; t < tgtNum; ++t)
        seVec[t] = getSolvAccess(chain, tgtS + t, envS, envE);

    return seVec;
}
