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
 * Author: Alberto Tilocca
 *
 */
#ifndef CEFUNCTIONS_H
#define	CEFUNCTIONS_H

typedef struct {
	int first;
	int second;
} afp;
namespace Victor {
    namespace Biopool {
        class CEFunctions {            
        public:
            //CONSTRUCTOR/DESTRUCTOR
            CEFunctions();
            ~CEFunctions();
            static void stampa();
            static void innerDistance(Spacer* ,vector<vector<double> >& );
            static void createSMatrix(Spacer* ,Spacer*, int ,int ,vector<vector<double> >&,vector<vector<double> >& ,vector<vector<double> >&);
            static void CESearch(Spacer*,Spacer*, vector<vector<afp> >&,vector<double>&,vector<vector<double> >&,int ,int,int,vector<vector<double> >&,vector<vector<double> >&,vector<int>&);
            static int finalSteps(Spacer*,Spacer*,vector<vector<afp> >&,int ,vector<int>&,int ,string ,string );
  

        protected:
            
        private:
            
        };
    }
}

#endif	/* CEFUNCTIONS_H */

