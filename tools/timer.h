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
/**  
@Description Timer class for benchmarking and more.

    @version 0.1*/

#ifndef __TIMER
#define __TIMER

extern "C" {
#include <time.h>
}

#include <iostream>
using namespace std;

/** Timer provides an interface for timing.

 *@version 0.1 */
class Timer {
public:

    Timer() : startTime(0), endTime(0) {
    }

    // We use default destructor and copy-constructor.

    // Modifiers.

    inline void start() {
        time(&startTime);
    }

    inline void stop() {
        time(&endTime);
    }

    inline void reset() {
        startTime = endTime = 0;
    }

    // Selectors.
    int seconds() const;
    int minutes() const;
    int hours() const;

protected:
    time_t startTime, endTime;
};

// output operator

inline ostream& operator<<(ostream& os, const Timer& time) {
    os << time.hours() << " hours " << (time.minutes() % 60) << " minutes "
            << (time.seconds() % 60) << " seconds.";
    return os;
}

#endif
