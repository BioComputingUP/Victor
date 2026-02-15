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
 * File:   CIFFileManager.cc
 * Author: Francesco Menniti
 *
 * Created on 1 dicembre 2015, 9.52
 */

#include<CIFFileManager.h>

using namespace std;
using namespace Victor;
using namespace Victor::Biopool;

void CIFFileManager::Init() {

    //PrintOut("CIFFileManager::Init");

    input.clear(); // reset file to previous content 
    input.seekg(0);

    while (input) {
        string atomLine = readLine(input, false);

        if (StartsWith(atomLine, LoopStart))
        {
            CIFLoop temp = processLoop();
            if(temp.GetCategoryGroup() == AtomSite)
                AtomLoop = temp;
            else if(temp.GetCategoryGroup() == SheetSite)
                SheetLoop = temp;
            else if(temp.GetCategoryGroup() == HelixSite)
                HelixLoop = temp;
            else
                OtherLoops.push_back(temp);
        }
        else
        {
            if(StartsWith(atomLine, "_"))
                processLine(atomLine.substr(1));
        }
    }
}

void CIFFileManager::processLine(string line){
    string key = Trim(SplitOnSpace(line)[0]);
    string value = Trim(line.substr(key.size()));
    inLineValues[key] = value;
}

CIFLoop CIFFileManager::processLoop() {
    
    //PrintOut("CIFFileManager::processLoop");
    
    string atomLine = readLine(input, false);
    if (atomLine.at(0) != '_') {
        ERROR("Loop don't start with _[definition]", exception);
    }

    CIFLoop* result = new CIFLoop(CIFLoop::DecodeCategory(atomLine));

    while (input) {

        if (atomLine.at(0) == '_') {
            result->AddColumn(atomLine);
        } else {
            if (atomLine.at(0) == '#')
            {
                return *result;
            }
            else
            {
                result->AddRow(atomLine);
            }
        }
        
        atomLine = readLine(input, false);
    }
    
    return *result;
    
}

string CIFFileManager::GetInLineValue(string key){
    return inLineValues[key];
}

