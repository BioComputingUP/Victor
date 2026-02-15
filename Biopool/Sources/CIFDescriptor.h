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

/* 
 * File:   CIFDescriptor.h
 * Author: Francesco Menniti
 *
 * Created on 18 dicembre 2015, 16.01
 * This file contains some CIF names.
 */

#ifndef CIFDESCRIPTOR_H
#define CIFDESCRIPTOR_H

using namespace std;

const string AtomSite = "atom_site";
const string SheetSite = "struct_sheet_range";
const string HelixSite = "struct_conf";

const string LoopStart = "loop_";

const string AtomID = "id";
const string AtomPDBModelNumber = "pdbx_PDB_model_num";
const string AtomChainID = "auth_asym_id";
const string AtomGroupPDB = "group_PDB";
const string AtomSeqID = "auth_seq_id";
const string AtomInsCode = "pdbx_PDB_ins_code";
const string AtomX = "Cartn_x";
const string AtomY = "Cartn_y";
const string AtomZ = "Cartn_z";
const string AtomBFact = "B_iso_or_equiv";
const string AtomName = "auth_atom_id";
const string AtomResidueName = "auth_comp_id";
const string AtomOccupancy = "occupancy";

const string SheetStart = "beg_auth_seq_id";
const string SheetEnd = "end_auth_seq_id";
const string SheetCode = "beg_auth_asym_id";

const string HelixStart = "beg_auth_seq_id";
const string HelixEnd = "end_auth_seq_id";
const string HelixCode = "beg_auth_asym_id";

#endif /* CIFDESCRIPTOR_H */

