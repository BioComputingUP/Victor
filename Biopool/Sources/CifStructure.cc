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

#include <iostream>
#include <algorithm>

#include <IoTools.h>
#include <Debug.h>

#include "CifStructure.h"

using namespace Victor;
using namespace Victor::Biopool;
using namespace std;

CifStructure::CifStructure(istream& input, ostream& output) :
input(input), output(output) {
    initData();
}

CifStructure::CifStructure(ostream& output) :
output(output), input(cin) {
    initData();
}

CifStructure::~CifStructure() {
}

void CifStructure::initData() {
    header = "_entry.id";

    // atom group
    atom = "_atom_site.";
    atomId = "id ";
    chain = "auth_asym_id ";
    asymId = "label_asym_id ";
    entityId = "label_entity_id ";
    residueIns = "pdbx_PDB_ins_code ";
    x = "Cartn_x ";
    y = "Cartn_y ";
    z = "Cartn_z ";
    occupancy = "occupancy ";
    tempFactor = "B_iso_or_equiv ";
    residueNum = "auth_seq_id ";
    residueName = "auth_comp_id ";
    atomName = "auth_atom_id ";
    model = "pdbx_PDB_model_num ";
    
    // helix group
    helix = "_struct_conf.";
    helixStart = "beg_auth_seq_id ";
    helixEnd = "end_auth_seq_id ";
    helixChainId = "beg_auth_asym_id ";

    // sheet group
    sheet = "_struct_sheet.";
    sheetOrder = "_struct_sheet_order.";
    sheetRange = "_struct_sheet_range.";
    sheetHbond = "_pdbx_struct_sheet_hbond.";
    sheetStart = "beg_auth_seq_id ";
    sheetEnd = "end_auth_seq_id ";
    sheetChainId = "beg_auth_asym_id ";

    // flags
    atomGroupParsed = false;
    helixGroupParsed = false;
    sheetGroupParsed = false;
    sheetOrderGroupParsed = false;
    sheetRangeGroupParsed = false;
    sheetHboundgroupParsed = false;
}

vector<string>& CifStructure::getGroup(string name) {
    if (name == "atom") {
	return atomGroup;
    } else if (name == "helix") {
	return helixGroup;
    } else if (name == "sheet") {
	return sheetGroup;
    } else if (name == "sheet order") {
	return sheetOrderGroup;
    } else if (name == "sheet range") {
	return sheetRangeGroup;
    } else if (name == "sheet hbond") {
	return sheetHbondGroup;
    } else {
	ERROR("getGroup (CifStructure): invalid group name", exception);
    }
}

string CifStructure::getTag(string name) {
    if (name == "header") {
	return header;
    }// atom group
    else if (name == "atom") {
	return atom;
    } else if (name == "atom id") {
	return atomId;
    } else if (name == "chain") {
	return chain;
    } else if (name == "atom asym") {
	return asymId;
    } else if (name == "atom entity") {
	return entityId;
    } else if (name == "residue ins") {
	return residueIns;
    } else if (name == "x") {
	return x;
    } else if (name == "y") {
	return y;
    } else if (name == "z") {
	return z;
    } else if (name == "occupancy") {
	return occupancy;
    } else if (name == "bfac") {
	return tempFactor;
    } else if (name == "residue num") {
	return residueNum;
    } else if (name == "residue name") {
	return residueName;
    } else if (name == "atom name") {
	return atomName;
    } else if (name == "model") {
	return model;
    }// helix group
    else if (name == "helix") {
	return helix;
    } else if (name == "helix start") {
	return helixStart;
    } else if (name == "helix end") {
	return helixEnd;
    } else if (name == "helix chain") {
	return helixChainId;
    }// sheet group
    else if (name == "sheet") {
	return sheet;
    } else if (name == "sheet order") {
	return sheetOrder;
    } else if (name == "sheet range") {
	return sheetRange;
    } else if (name == "sheet hbond") {
	return sheetHbond;
    } else if (name == "sheet start") {
	return sheetStart;
    } else if (name == "sheet end") {
	return sheetEnd;
    } else if (name == "sheet chain") {
	return sheetChainId;
    } else {
	ERROR("getTag (CifStructure): invalid tag name", exception);
    }
}

int CifStructure::getGroupColumnNumber(string name, string field) {
    int col = -1;
    vector<string>& group = getGroup(name);
    string tag = getTag(field);
    for (vector<string>::iterator it = group.begin(); it != group.end(); it++) {
	if (*it == tag) {
	    col = it - group.begin();
	    break;
	}
    }
    return col;
}

string CifStructure::getGroupField(string name, string& line, int columnNum) {
    istringstream iss(line);
    vector<string>& group = getGroup(name);
    vector<string> fields;
    string field;
    for (unsigned int i = 0; i < group.size(); i++) {
	iss >> field;
	fields.push_back(field);
    }
    if (columnNum != -1) {
	return fields[columnNum];
    } else {
	return "?";
    }
}

string CifStructure::getInlineField(string& line) {
    istringstream iss(line);
    string tag, field;
    iss >> tag >> field;
    return field;
}

void CifStructure::parseGroup(string name, string& line) {
    bool found = false;
    vector<string>& group = getGroup(name);

    // exit the function if the group name is already parsed
    if (!isGroupParsed(name)) {
	while (input) {
	    string groupName(getTag(name));
	    size_t pos = line.find(groupName);
	    if (pos != string::npos) {
		group.push_back(line.substr(pos + groupName.size(),
			line.size() - groupName.size()));
		found = true;
	    } else {
		found = false;
	    }

	    // exit the loop when the research of the fields is completed
	    if (!found && group.size() > 1) {
		setParsedFlag(name);
		break;
	    }

	    line = readLine(input);
	}
    }
}

void CifStructure::setParsedFlag(string name) {
    if (name == "atom") {
	atomGroupParsed = true;
    } else if (name == "helix") {
	helixGroupParsed = true;
    } else if (name == "sheet") {
	sheetGroupParsed = true;
    } else if (name == "sheet hbond") {
	sheetHboundgroupParsed = true;
    } else if (name == "sheet order") {
	sheetOrderGroupParsed = true;
    } else if (name == "sheet range") {
	sheetRangeGroupParsed = true;
    } else {
	ERROR("setParsedFlag (CifStructure): invalid flag name", exception);
    }
}

bool CifStructure::isGroupParsed(string name) {
    if (name == "atom") {
	return atomGroupParsed;
    } else if (name == "helix") {
	return helixGroupParsed;
    } else if (name == "sheet") {
	return sheetGroupParsed;
    } else if (name == "sheet hbond") {
	return sheetHboundgroupParsed;
    } else if (name == "sheet order") {
	return sheetOrderGroupParsed;
    } else if (name == "sheet range") {
	return sheetRangeGroupParsed;
    } else {
	ERROR("isGroupParsed (CifStructure): invalid group name", exception);
    }
}

void CifStructure::printGroup(string name) {

    output << "loop_" << endl;

    if (name == "atom") {
	output << "_atom_site.group_PDB " << endl <<
		"_atom_site.id " << endl <<
		"_atom_site.type_symbol " << endl <<
		"_atom_site.label_atom_id " << endl <<
		"_atom_site.label_alt_id " << endl <<
		"_atom_site.label_comp_id " << endl <<
		"_atom_site.label_asym_id " << endl <<
		"_atom_site.label_entity_id " << endl <<
		"_atom_site.label_seq_id " << endl <<
		"_atom_site.pdbx_PDB_ins_code " << endl <<
		"_atom_site.Cartn_x " << endl <<
		"_atom_site.Cartn_y " << endl <<
		"_atom_site.Cartn_z " << endl <<
		"_atom_site.occupancy " << endl <<
		"_atom_site.B_iso_or_equiv " << endl <<
		"_atom_site.Cartn_x_esd " << endl <<
		"_atom_site.Cartn_y_esd " << endl <<
		"_atom_site.Cartn_z_esd " << endl <<
		"_atom_site.occupancy_esd " << endl <<
		"_atom_site.B_iso_or_equiv_esd " << endl <<
		"_atom_site.pdbx_formal_charge " << endl <<
		"_atom_site.auth_seq_id " << endl <<
		"_atom_site.auth_comp_id " << endl <<
		"_atom_site.auth_asym_id " << endl <<
		"_atom_site.auth_atom_id " << endl <<
		"_atom_site.pdbx_PDB_model_num " << endl;
    } else if (name == "helix") {
	output << "_struct_conf.conf_type_id " << endl <<
		"_struct_conf.id " << endl <<
		"_struct_conf.pdbx_PDB_helix_id " << endl <<
		"_struct_conf.beg_label_comp_id " << endl <<
		"_struct_conf.beg_label_asym_id " << endl <<
		"_struct_conf.beg_label_seq_id " << endl <<
		"_struct_conf.pdbx_beg_PDB_ins_code " << endl <<
		"_struct_conf.end_label_comp_id " << endl <<
		"_struct_conf.end_label_asym_id " << endl <<
		"_struct_conf.end_label_seq_id " << endl <<
		"_struct_conf.pdbx_end_PDB_ins_code " << endl <<
		"_struct_conf.beg_auth_comp_id " << endl <<
		"_struct_conf.beg_auth_asym_id " << endl <<
		"_struct_conf.beg_auth_seq_id " << endl <<
		"_struct_conf.end_auth_comp_id " << endl <<
		"_struct_conf.end_auth_asym_id " << endl <<
		"_struct_conf.end_auth_seq_id " << endl <<
		"_struct_conf.pdbx_PDB_helix_class " << endl <<
		"_struct_conf.details " << endl <<
		"_struct_conf.pdbx_PDB_helix_length " << endl;
    } else if (name == "sheet") {
	output << "_struct_sheet.id " << endl <<
		"_struct_sheet.type " << endl <<
		"_struct_sheet.number_strands " << endl <<
		"_struct_sheet.details " << endl;
    } else if (name == "sheet range") {
	output << "_struct_sheet_range.sheet_id " << endl <<
		"_struct_sheet_range.id " << endl <<
		"_struct_sheet_range.beg_label_comp_id " << endl <<
		"_struct_sheet_range.beg_label_asym_id " << endl <<
		"_struct_sheet_range.beg_label_seq_id " << endl <<
		"_struct_sheet_range.pdbx_beg_PDB_ins_code " << endl <<
		"_struct_sheet_range.end_label_comp_id " << endl <<
		"_struct_sheet_range.end_label_asym_id " << endl <<
		"_struct_sheet_range.end_label_seq_id " << endl <<
		"_struct_sheet_range.pdbx_end_PDB_ins_code " << endl <<
		"_struct_sheet_range.symmetry " << endl <<
		"_struct_sheet_range.beg_auth_comp_id " << endl <<
		"_struct_sheet_range.beg_auth_asym_id " << endl <<
		"_struct_sheet_range.beg_auth_seq_id " << endl <<
		"_struct_sheet_range.end_auth_comp_id " << endl <<
		"_struct_sheet_range.end_auth_asym_id " << endl <<
		"_struct_sheet_range.end_auth_seq_id " << endl;
    } else if (name == "sheet hbond") {
	output << "_pdbx_struct_sheet_hbond.sheet_id " << endl <<
		"_pdbx_struct_sheet_hbond.range_id_1 " << endl <<
		"_pdbx_struct_sheet_hbond.range_id_2 " << endl <<
		"_pdbx_struct_sheet_hbond.range_1_label_atom_id " << endl <<
		"_pdbx_struct_sheet_hbond.range_1_label_comp_id " << endl <<
		"_pdbx_struct_sheet_hbond.range_1_label_asym_id " << endl <<
		"_pdbx_struct_sheet_hbond.range_1_label_seq_id " << endl <<
		"_pdbx_struct_sheet_hbond.range_1_PDB_ins_code " << endl <<
		"_pdbx_struct_sheet_hbond.range_1_auth_atom_id " << endl <<
		"_pdbx_struct_sheet_hbond.range_1_auth_comp_id " << endl <<
		"_pdbx_struct_sheet_hbond.range_1_auth_asym_id " << endl <<
		"_pdbx_struct_sheet_hbond.range_1_auth_seq_id " << endl <<
		"_pdbx_struct_sheet_hbond.range_2_label_atom_id " << endl <<
		"_pdbx_struct_sheet_hbond.range_2_label_comp_id " << endl <<
		"_pdbx_struct_sheet_hbond.range_2_label_asym_id " << endl <<
		"_pdbx_struct_sheet_hbond.range_2_label_seq_id " << endl <<
		"_pdbx_struct_sheet_hbond.range_2_PDB_ins_code " << endl <<
		"_pdbx_struct_sheet_hbond.range_2_auth_atom_id " << endl <<
		"_pdbx_struct_sheet_hbond.range_2_auth_comp_id " << endl <<
		"_pdbx_struct_sheet_hbond.range_2_auth_asym_id " << endl <<
		"_pdbx_struct_sheet_hbond.range_2_auth_seq_id " << endl;
    } else if (name == "entity poly") {
	output << "_entity_poly_seq.entity_id " << endl <<
		"_entity_poly_seq.num " << endl <<
		"_entity_poly_seq.mon_id " << endl <<
		"_entity_poly_seq.hetero " << endl;
    } else if (name == "pdbx poly") {
	output << "_pdbx_poly_seq_scheme.asym_id " << endl <<
		"_pdbx_poly_seq_scheme.entity_id " << endl <<
		"_pdbx_poly_seq_scheme.seq_id " << endl <<
		"_pdbx_poly_seq_scheme.mon_id " << endl <<
		"_pdbx_poly_seq_scheme.ndb_seq_num " << endl <<
		"_pdbx_poly_seq_scheme.pdb_seq_num " << endl <<
		"_pdbx_poly_seq_scheme.auth_seq_num " << endl <<
		"_pdbx_poly_seq_scheme.pdb_mon_id " << endl <<
		"_pdbx_poly_seq_scheme.auth_mon_id " << endl <<
		"_pdbx_poly_seq_scheme.pdb_strand_id " << endl <<
		"_pdbx_poly_seq_scheme.pdb_ins_code " << endl <<
		"_pdbx_poly_seq_scheme.hetero " << endl;
    } else {
	ERROR("printGroup (CifStructure): invalid group name", exception);
    }
}
