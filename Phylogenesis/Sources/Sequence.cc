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

#include <Sequence.h>
#include <iostream>

namespace Victor { namespace Phylogenesis{

static unsigned int HYDROPHILIC_STRETCH_SIZE = 5;

Sequence::Sequence(){
	id = -1;
	divergent = false;
}

/**
*
* @param name of the RawSequence
* @param seq  of the RawSequence
* @param id   of the RawSequence
*/
Sequence::Sequence(std::string name, std::string seq, int id): name(name), sequence(seq), id(id), divergent(false){
}

bool Sequence::isGap(int pos){
	if(pos < 0 || pos >= (int)sequence.size())
		return false;
	return sequence[pos] == '-';
}

bool Sequence::isHydrophilic(int pos){
	if(	sequence[pos] == 'D' ||
		sequence[pos] == 'E' ||
		sequence[pos] == 'G' ||
		sequence[pos] == 'K' ||
		sequence[pos] == 'N' ||
		sequence[pos] == 'Q' ||
		sequence[pos] == 'P' ||
		sequence[pos] == 'R' ||
		sequence[pos] == 'S'){
		return true;
	}
	return false;
}

bool Sequence::isHydrophilicStretch(int pos){
	if(pos < 0 || pos >= (int)sequence.size())
		return false;

	bool goOn = true;
	for(unsigned int i = pos, j = 0; i<sequence.size() && j < HYDROPHILIC_STRETCH_SIZE && goOn; i++, j++){
		if(!isHydrophilic(i)){
			goOn = false;
		}
	}
	return goOn;
}

std::string Sequence::GetPureSequence(const std::string &s) {
	std::string result = "";

   for (unsigned int i = 0; i < s.size(); ++i)
	   if (s[i] != '-')
		   result.append(1, s[i]);

   return result;
}

std::vector<Sequence*> Sequence::LoadFasta(std::ifstream& inputFile){
	std::vector<Sequence*> seqs;
	int id=0;
	std::string line, name, content;
	while(std::getline( inputFile, line ).good()){
		if( line.empty() || line[0] == '>' ){
			if( !name.empty() ){
				seqs.push_back(new Sequence(name, GetPureSequence(content), id++));
				name.clear();
			}
			if( !line.empty() ){
				name = line.substr(1);
			}
			content.clear();
		} else if( !name.empty() ){
			if( line.find(' ') != std::string::npos ){
				name.clear();
				content.clear();
			} else {
				content += line;
			}
		}
	}
	if( !name.empty() ){
		seqs.push_back(new Sequence(name, content, id++));
	}
	return seqs;
}

double Sequence::GetAvgLength(std::vector<Sequence*> seq){
	double length = 0;
	for(unsigned int i=0;i<seq.size();i++){
		length += Sequence::GetPureSequence(seq[i]->sequence).size();
	}
	return length / seq.size();
}

}}
