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
#include <Group.h>
#include <vector3.h>
#include <XyzLoader.h>
#include <XyzSaver.h>
#include <vglMath.h>

using namespace Victor;using namespace Victor::Biopool;

vg_ieee64 degreesToRadian(const vg_ieee64& deg)
{
return (deg / 180.0) * M_PI;
}

int main()
{ 
  cout << "Start" << endl;

  Group g, g2;
  ifstream inFile("test-grp.xyz");
  if (!inFile)
    ERROR("File not found.", exception);

  XyzLoader il(inFile);
  g.load(il);

  ifstream inFile2("test-grp.xyz");
  if (!inFile2)
    ERROR("File not found.", exception);
  XyzLoader il2(inFile2);
  g2.load(il2);

  g2.setType("2-C3H8");

  cout << "-------------------------------------------------------------\n";

  XyzSaver is(cout, 0);
  g.save(is);

  cout << "-------------------------------------------------------------\n";
  
  g2.save(is);
  g2.bindIn(g2[0], g, g[6]);

  cout << "-------------------------------------------------------------\n";
  g2.save(is);
  
  vgVector3<double> tmp(0,0,1);
  vgMatrix3<double> rotationMatrix = 
          vgMatrix3<double>::createRotationMatrix(tmp, degreesToRadian(180));
  g2[0].addRot(rotationMatrix);

  cout << "-------------------------------------------------------------\n";
  g2.save(is);

  g2.addRot(rotationMatrix);

  cout << "-------------------------------------------------------------\n";
  g2.save(is);

  cout << "-------------------------------------------------------------\n";

  cout << "superior: \t g.a1= " << g[1].getSuperior().getType() 
       << "\t g2.a3= " << g2[3].getSuperior().getType() << endl;

  return 0;
}
