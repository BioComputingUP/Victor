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

#include <CIFLoop.h>

using namespace std;
using namespace Victor;
using namespace Victor::Biopool;

CIFLoop::CIFLoop(string _category) {
    //PrintOut("CIFLoop::CIFLoop + " + _category );
    categoryGroup = _category;
}

void CIFLoop::AddColumn(string _line) {
    //PrintOut("CIFLoop::AddColumn + " + _line);
    columns.push_back(decodeColName(_line));
}

string CIFLoop::GetCellValue(string _colName, int _rowNumber) {
    return rows[_rowNumber][getColNumber(_colName)];
}

string CIFLoop::GetCellValueTrim(string _colName, int _rowNumber) {
    return Trim(rows[_rowNumber][getColNumber(_colName)]);
}

vector<string> CIFLoop::splitRow(string _line) {
    //PrintOut("CIFLoop::splitRow + " + _line);
    return SplitOnSpace(_line);
}

string CIFLoop::decodeColName(string _line) {
    return Trim(_line.substr(_line.find('.') + 1, _line.size() - _line.find('.')));
}

string CIFLoop::DecodeCategory(string _line) {
    return Trim(_line.substr(1, _line.find('.') - 1));
}

bool CIFLoop::HasColumn(string _colName) {
    for (int i = 0; i < ColumnsCount(); i++)
        if (columns[i] == _colName)
            return true;

    return false;
}

int CIFLoop::getColNumber(string _colName) {
    if (!HasColumn(_colName))
        ERROR("Loop hasn't column " + _colName, exception);

    for (int i = 0; i < ColumnsCount(); i++)
        if (columns[i] == _colName)
            return i;

    return -1;

}
