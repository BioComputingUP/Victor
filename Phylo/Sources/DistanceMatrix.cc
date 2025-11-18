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

#include <DistanceMatrix.h>
#include <limits>
#include <iostream>

namespace Victor {
namespace Phylo {
    typedef DistanceMatrix::key key;

    const double DistanceMatrix::min_distance = 0.0;
    const double DistanceMatrix::max_distance = std::numeric_limits<double>::max();


    DistanceMatrix::~DistanceMatrix() {
    }



    double
    DistanceMatrix::operator()(const string &a, const string &b) const {
        return getElement(a, b);
    }


    DistanceMatrix &
    DistanceMatrix::operator()(const string &a, const string &b, const double distance) {
        return setDistance(a, b, distance);
    }



    bool
    DistanceMatrix::isEmpty() const {
        return getSize() <= 1;
    }



    bool
    DistanceMatrix::isSet(const string &a, const string &b) const {
        return isSetDirectional(a, b) || isSetDirectional(b, a);
    }



    bool
    DistanceMatrix::hasOTU(const string &OTU) const {
        for (auto &label : OTUs) {
            if (label == OTU) {
                return true;
            }
        }
        return false;
    }



    double
    DistanceMatrix::getElement(const string &a, const string &b) const {
        if (a == b) {
            return min_distance;
        } else if (isSetDirectional(a, b)) {
            return find(make_pair(a, b))->second;
        } else if (isSetDirectional(b, a)) {
            return find(make_pair(b, a))->second;
        } else {
            return max_distance;
        }
    }



    set<string>
    DistanceMatrix::getOTUs() const {
        return OTUs;
    }



    size_t
    DistanceMatrix::getSize() const {
        return OTUs.size();
    }



    double
    DistanceMatrix::getMinimum() const {
        double min = max_distance;

        for (auto &entry : *this) {
            if (entry.second <= min) {
                min = entry.second;
            }
        }

        return min;
    }



    key
    DistanceMatrix::getMinimumPosition() const {
        double min = max_distance;
        key pos;
        for (auto &entry : *this) {
            if (entry.second <= min) {
                min = entry.second;
                pos = entry.first;
            }
        }

        return pos;
    }



    double
    DistanceMatrix::getMaximum() const {
        double max = min_distance;

        for (auto &entry : *this) {
            if (entry.second >= max) {
                max = entry.second;
            }
        }

        return max;
    }



    key
    DistanceMatrix::getMaximumPosition() const {
        double max = min_distance;
        key pos;

        for (auto &entry : *this) {
            if (entry.second >= max) {
                max = entry.second;
                pos = entry.first;
            }
        }

        return pos;
    }



    DistanceMatrix &
    DistanceMatrix::setDistance(const string &a, const string &b, const double distance) {
        key pos = (isSetDirectional(a, b)) ? make_pair(a, b) : make_pair(b, a);
        (*this)[pos] = distance;
        return *this;
    }



    DistanceMatrix &
    DistanceMatrix::unsetDistance(const string &a, const string &b) {
        if (isSetDirectional(a, b))
            this->erase(find(make_pair(a, b)));
        if (isSetDirectional(b, a))
            this->erase(find(make_pair(b, a)));
        return *this;
    }



    DistanceMatrix &
    DistanceMatrix::unsetDistance(const string &OTU) {
        for (auto &other : getOTUs()) {
            unsetDistance(OTU, other);
        }
        return *this;
    }



    DistanceMatrix &
    DistanceMatrix::addOTU(const string &OTU) {
        if (!hasOTU(OTU)) {
            OTUs.insert(OTU);
        }
        return *this;
    }


    DistanceMatrix &
    DistanceMatrix::removeOTU(const string &OTU) {
        // Removes OTU from list
        OTUs.erase(OTU);

        // Removes OTU from distances
        for (auto &i : OTUs) {
            erase(make_pair(i, OTU));
            erase(make_pair(OTU, i));
        }

        return *this;
    }



    bool
    DistanceMatrix::isSetDirectional(const string &a, const string &b) const {
        return find(make_pair(a, b)) != end();
    }
}  // namespace Phylo
}  // namespace Victor