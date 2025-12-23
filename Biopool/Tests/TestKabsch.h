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
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>
#include <KabschAlgorithm.h>

#include <PdbLoader.h>

using namespace std;
using namespace Victor::Biopool;

class TestKabsch : public CppUnit::TestFixture 
{
    private:
        
    public:

        // SETUP METHOD
        void setUp() 
        {}

        /// TEARDOWN METHOD

        void tearDown() 
        {}

    

        void testRotoTranslation()
        {
            KabschAlgorithm kb= KabschAlgorithm();
            Eigen::Matrix<double,10,3> A, B;
            A(0,0)= 10.5;
            A(0,1)= 14.7;
            A(0,2)= 20.5;
            
            A(1,0)= -10.2;
            A(1,1)= 30.1;
            A(1,2)= -8.6;
            
            A(2,0)= -31;
            A(2,1)= -10;
            A(2,2)= 2;
            
            A(3,0)= 26.1;
            A(3,1)= 15.0;
            A(3,2)= -2.1;
            
            A(4,0)= 20.5 ;
            A(4,1)= 41.0;
            A(4,2)= 3.6;
            
            A(5,0)= 5.5 ;
            A(5,1)= 4.1;
            A(5,2)= -8.3;
            
            A(6,0)= 0.8 ;
            A(6,1)= 4.3;
            A(6,2)= 10.7;
            
            A(7,0)= 10.2 ;
            A(7,1)= 6.3;
            A(7,2)= -1.5;
            
            A(8,0)= 50.2 ;
            A(8,1)= 10.6 ;
            A(8,2)= -9.3 ;
            
            A(9,0)= 45.1;
            A(9,1)= -1.0;
            A(9,2)= 47.2;

            B(0,0)= 0.5;
            B(0,1)= 24.7;
            B(0,2)= 0.5;
            
            B(1,0)= 11.2;
            B(1,1)= 50.1;
            B(1,2)= -0.6;
            
            B(2,0)= 3.5;
            B(2,1)= 0.0;
            B(2,2)= -4;
            
            B(3,0)= 16.1;
            B(3,1)= 2.0;
            B(3,2)= -10.1;
            
            B(4,0)= 10.5 ;
            B(4,1)= 31.0;
            B(4,2)= -3.6;
            
            B(5,0)= 8.5 ;
            B(5,1)= 6.1;
            B(5,2)= -1.3;
            
            B(6,0)= -0.8 ;
            B(6,1)= 14.3;
            B(6,2)= 0.7;
            
            B(7,0)= 20.2 ;
            B(7,1)= 16.3;
            B(7,2)= -11.5;
            
            B(8,0)= 10.2 ;
            B(8,1)= -10.6 ;
            B(8,2)= -9.3 ;
            
            B(9,0)= 25.1;
            B(9,1)= -11.0;
            B(9,2)= 27.2;

            bool test= kb.kabschTest(A,B,0.001);
            CPPUNIT_ASSERT(test == true);
        }

        void testIdentity()
        {
            KabschAlgorithm kb= KabschAlgorithm();
            Eigen::Matrix<double,10,3> A, B;
            A(0,0)= 10.5;
            A(0,1)= 14.7;
            A(0,2)= 20.5;
            
            A(1,0)= -10.2;
            A(1,1)= 30.1;
            A(1,2)= -8.6;
            
            A(2,0)= -31;
            A(2,1)= -10;
            A(2,2)= 2;
            
            A(3,0)= 26.1;
            A(3,1)= 15.0;
            A(3,2)= -2.1;
            
            A(4,0)= 20.5 ;
            A(4,1)= 41.0;
            A(4,2)= 3.6;
            
            A(5,0)= 5.5 ;
            A(5,1)= 4.1;
            A(5,2)= -8.3;
            
            A(6,0)= 0.8 ;
            A(6,1)= 4.3;
            A(6,2)= 10.7;
            
            A(7,0)= 10.2 ;
            A(7,1)= 6.3;
            A(7,2)= -1.5;
            
            A(8,0)= 50.2 ;
            A(8,1)= 10.6 ;
            A(8,2)= -9.3 ;
            
            A(9,0)= 45.1;
            A(9,1)= -1.0;
            A(9,2)= 47.2;

            B=A;

            bool test= kb.kabschTest(A,B,0.001);
            CPPUNIT_ASSERT(test == true);

        }

         
        static CppUnit::Test *suite()
        {
            CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite( "TestKabsch" );
            suiteOfTests->addTest( new CppUnit::TestCaller<TestKabsch>("testRotoTranslation",&TestKabsch::testRotoTranslation ) );
            suiteOfTests->addTest( new CppUnit::TestCaller<TestKabsch>("testIdentity",&TestKabsch::testIdentity ) );
            return suiteOfTests;
        }
    };