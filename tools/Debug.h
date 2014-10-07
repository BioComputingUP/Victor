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
#ifndef __DEBUG_H__
#define __DEBUG_H__

#include <exception>
using namespace std;

extern "C" {
#include <stdlib.h>
}

#include <iostream>
#include <fstream>
using namespace std;

// The user may specify other output streams than cerr.
#ifndef OSTREAM
#define OSTREAM cerr
#endif

#if VERBOSE >= 3
#define V3OUT OSTREAM
#else
#define V3OUT ofstream("/dev/null")
#endif

#if VERBOSE >= 2
#define V2OUT OSTREAM
#else
#define V2OUT ofstream("/dev/null")
#endif

#if VERBOSE >= 1
#define V1OUT OSTREAM
#else
#define V1OUT ofstream("/dev/null")
#endif

// If internalization is used, please include <libintl.h> *before* <debug.h>
#ifndef gettext
#define gettext(aString) aString
#endif

// The user may redefine TERMINATION to his needs (e.g. "abort()").
#ifndef TERMINATION
#define TERMINATION exit(1)
#endif

#ifndef NEXCEPTIONS
#define TERMINATION throw(exception)
#endif

#define LOCATION SourceLocator(__FILE__, __LINE__, __FUNCTION__)

/** SourceLocator: class to provide information about
    C++-Source-Location.

    SourceLocator is used for error reports to "OSTREAM". */
class SourceLocator {
public:

    SourceLocator(const char* f, long l, const char* fu = 0)
    : file(f), func(fu), line(l) {
    }

    friend ostream& operator<<(ostream& os, const SourceLocator& loc) {
        if (loc.func != 0) {
            return os << loc.file << ":" << loc.line << ": "
                    << loc.func << ": ";
        }
        else {
            return os << loc.file << ":" << loc.line << ": ";
        }
    };
private:
    const char* file;
    const char* func;
    long line;
};

/** Class helping to trace number of instances of a class. */
class ObjectTrace {
public:

    ObjectTrace() : ct(++count) {
        V3OUT << "DEBUG_TRACE: Object " << ct << " constructed." << endl;
    }
    // commented to remove useless compiler warning - ST 23.2.99:
    // ObjectTrace(const char* n) :  ct(++count)
    //  { V3OUT << "DEBUG_TRACE: Object " << ct << " constructed." << endl; }
    //  ObjectTrace(const ObjectTrace& orig) : ct(++count)
    //  { V3OUT << "DEBUG_TRACE: Object " << ct << " constructed." << endl; }

    ~ObjectTrace() {
        V3OUT << "DEBUG_TRACE: Object " << ct << " destroyed." << endl;
    }
private:
    static unsigned count;
    unsigned ct;
};

#ifndef NDEBUG

#define PRINT_LOCATION V3OUT << "PRINT_LOCATION: " << LOCATION << endl  
#define PRINT_NAME V3OUT << "PRINT_NAME: " << __PRETTY_FUNCTION__ << endl

#define ASSERT(assertion, exception) {                              \
   V3OUT << "ASSERTION: " << #assertion << endl;                    \
   if (!(assertion)) {                                              \
      OSTREAM << (LOCATION)                                         \
           << gettext("Assert: assertion violated.") << endl        \
           << #assertion << endl                                    \
           << gettext("Aborting.") << endl << flush;                \
      TERMINATION;                                                  \
   } }

#if VERBOSE >= 1
#define DEBUG_BREAK(aCondition, aString)                            \
   if (aCondition) {                                                \
      V1OUT << "DEBUG_BRK:  " << aString << endl << flush;          \
      getchar(); } 
#else
#define DEBUG_BREAK(a, b)
#endif

#define DEBUG_TRACE ObjectTrace _debug_trace_

#define DEBUG_MSG(aString) V1OUT << "DEBUG_MSG:  " << (aString) << endl
#define DUMP(x) V2OUT << "DEBUG_DUMP: " << #x << " == " << (x) << endl
#define PRECOND(a, b)   ASSERT(a, b)
#define POSTCOND(a, b)  ASSERT(a, b)
#define INVARIANT(a, b) ASSERT(a, b)

#else /* NDEBUG */

#define ASSERT(a, b)
#define PRECOND(a, b) 
#define POSTCOND(a, b) 
#define INVARIANT(a, b) 
#define DEBUG_MSG(a)
#define DEBUG_BREAK(a, b)
#define DEBUG_TRACE 
#define DUMP(a)
#define PRINT_LOCATION
#define PRINT_NAME

#endif /* NDEBUG */
#define exception string("")
#define ERROR(message, exception) {                                 \
  OSTREAM << (LOCATION)                                             \
       << gettext("Error: an error occured.") << endl               \
       << gettext(message) << endl                                  \
       << gettext("Aborting.") << endl << flush;                    \
  TERMINATION; }

#define ERROR_IF(condition, message, exception) {                   \
  if (condition) { ERROR(message, exception); } }

#endif /* __DEBUG_H__ */


