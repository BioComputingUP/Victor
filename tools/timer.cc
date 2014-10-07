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
#define VERSIONSTRING "$Date: 2014-05-09  $" 
//  Timer class for benchmarking and more.
//
//  @version 0.1

#include "timer.h"

int Timer::seconds() const {
    if (endTime != 0) {
        return endTime - startTime;
    }
    else {
        time_t now = time(0);
        return now - startTime;
    }
}

int Timer::minutes() const {
    if (endTime != 0) {
        return (endTime - startTime) / 60;
    }
    else {
        time_t now = time(0);
        return (now - startTime) / 60;
    }
}

int Timer::hours() const {
    if (endTime != 0) {
        return (endTime - startTime) / 3600;
    }
    else {
        time_t now = time(0);
        return (now - startTime) / 3600;
    }
}




