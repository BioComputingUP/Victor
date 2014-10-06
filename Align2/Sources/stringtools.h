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
// stringtools.h

#ifndef __stringtools_H__
#define __stringtools_H__

#include <iostream>
#include <sstream>
#include <string>
#include <vector>


void strip(std::string &str);

std::string strip_and_return_string(std::string);

int strip_and_return_int(std::string);

unsigned int strip_and_return_unsigned_int(std::string);

float strip_and_return_float(std::string);

std::vector<std::string> split(std::string, char);

int seq_length(const std::string);

// template<class T>
// std::string to_string(const T&);

template<class T> std::string to_string(const T &x)
{
	std::ostringstream oss;
	oss << x;
	return oss.str();
}

std::string StringToLower(std::string);
std::string StringToUpper(std::string);

#endif
