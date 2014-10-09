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
#include <vector3.h> 
#include <vector> 
#include <iostream>
#include <Atom.h>
#include <IoTools.h>

using namespace Victor;using namespace Victor::Biopool;

int main()
{ 
  cout << "Hello World" << endl;
  Atom at, at2, at3;
  at.setCoords(0,0,0);
  at2.setCoords(1,1,0);
  at.setType("C");  
  cout << "AtomType = " << at.getType() << "\t Distance = " << at.distance(at2)
       << endl;
  at3 = at;
  cout << "AT3: \t Type = " << at3.getType()  << "\t Dist = " << at.distance(at3) 
       << endl;
  at.setType("O");
  at.setCoords(1,0,0);
  cout << "AT3: \t Type = " << at3.getType()  << "\t Dist = " << at.distance(at3) 
       << endl;
  cout << "New: " << endl;
  at2.bindIn(at);
  cout << "at: \t" << at.sizeInBonds() << "\t" << at.sizeOutBonds() << endl;
  cout << "at2: \t" << at2.sizeInBonds() << "\t" << at2.sizeOutBonds() << endl;
  cout << "at3: \t" << at3.sizeInBonds() << "\t" << at3.sizeOutBonds() << endl;
  at2.setType("N");
  cout << "at: \t" << at.getOutBond(0).getType() << endl;

  vector<Atom> vec;

  vec.push_back(at);
  vec.push_back(at2);
  vec.push_back(at3);

  cout << "--------------------------------\n";
  cout << "vec= ";
  for (unsigned int i = 0; i < vec.size(); i++)
    cout << "\t" << vec[i].getType();
  cout << endl;

  return 0;
}
