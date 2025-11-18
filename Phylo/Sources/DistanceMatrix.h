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

#ifndef _VICTOR_PHYLO_DISTANCEMATRIX_H_
#define _VICTOR_PHYLO_DISTANCEMATRIX_H_

#include <string>
#include <vector>
#include <set>
#include <map>

using std::string;
using std::pair;
using std::vector;
using std::set;
using std::map;

namespace Victor {
namespace Phylo {

/** @brief Implements a distance matrix.
 *
 * A distance matrix holds distances among OTUs.
 * This class uses Method Cascading (through Method Chaining).
 * 
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 */
class DistanceMatrix: public map<pair<string, string>, double> {
 public:
    /** Key to access an entry in the matrix. */
    typedef pair<string, string> key;


    /** @brief Destructor. */
    virtual ~DistanceMatrix();



    /** @brief Returns distance between given OTUs.
     *
     * If there is no recorded distance between given OTUs, returns
     * maximum representable distance.
     * Distance matrices are symmetric, hence order of indeces does not
     * matter.
     * @param[in]       a         First OTU
     * @param[in]       b         Second OTU
     * @return Distance between given OTUs
     */
    double operator()(const string &a, const string &b) const;

    /** @brief Sets distance between given OTUs.
     *
     * If a distance was previously recorder for given OTUs, it is
     * overwritten.
     * @param[in]       a         First OTU
     * @param[in]       b         Second OTU
     * @param[in]       distance  Distance between a and b
     * @return This matrix itself
     */
    DistanceMatrix &operator()(const string &a, const string &b, const double distance);



    /** @brief Tells whether this matrix is empty.
     *
     * A distance matrix is empty when it contains 1 or less OTUs.
     * @retval true     This distance matrix is empty
     * @retval false    This distance matrix is not empty
     */
    bool isEmpty() const;

    /** @brief Tells whether there is a distance recorder between given OTUs.
     *
     * Distance matrices are symmetric, hence order of indeces does not
     * matter.
     * @param[in]       a         First OTU
     * @param[in]       b         Second OTU
     * @retval true     There is a record for given OTUs
     * @retval false    There is no record for given OTUs
     */
    bool isSet(const string &a, const string &b) const;

    /** @brief Tells whether given OTU exists in this matrix.
     * @param[in]       OTU       Given OTU
     * @retval true     Given OTU exists in this matrix
     * @retval false    Given OTU does not exist in this matrix
     */
    bool hasOTU(const string &OTU) const;



    /** @brief Returns distance between given OTUs.
     *
     * If there is no recorded distance between given OTUs, returns
     * maximum representable distance.
     * Distance matrices are symmetric, hence order of indeces does not
     * matter.
     * @param[in]       a         First OTU
     * @param[in]       b         Second OTU
     * @return Distance between given OTUs
     */
    double getElement(const string &a, const string &b) const;

    /** @brief Returns OTUs in this matrix.
     * @return OTUs as vector<string>
     */
    set<string> getOTUs() const;

    /** @brief Returns number of OTUs in this matrix.
     * @return Number of OTUs in this matrix
     */
    size_t getSize() const;

    /** @brief Returns minimum distance in this matrix.
     * @return Minimum distance in this matrix
     */
    double getMinimum() const;

    /** @brief Returns position of minimum distance in this matrix.
     *
     * Returns position as a pair<string, string>, where strings are
     * the OTUs.
     * @return Position of minimum distance
     */
    key getMinimumPosition() const;

    /** @brief Returns maximum distance in this matrix.
     * @return Maximum distance in this matrix
     */
    double getMaximum() const;

    /** @brief Returns position of maximum distance in this matrix.
     *
     * Returns position as a pair<string, string>, where strings are
     * the OTUs.
     * @return Position of maximum distance
     */
    key getMaximumPosition() const;



    /** @brief Sets distance between given OTUs.
     *
     * If a distance was previously recorder for given OTUs, it is
     * overridden.
     * @param[in]       a         First OTU
     * @param[in]       b         Second OTU
     * @param[in]       distance  Distance between a and b
     * @return This matrix itself
     */
    DistanceMatrix &setDistance(const string &a, const string &b, const double distance);

    /** @brief Unsets the distance between given OTUs.
     *
     * If no distance was set, nothing happens.
     * @param[in]       a         First OTU
     * @param[in]       b         Second OTU
     * @return This matrix itself
     */
    DistanceMatrix &unsetDistance(const string &a, const string &b);

    /** @brief Unsets every distance involving given OTU.
     * @param[in]       OTU       OTU
     * @return This matrix itself
     */
    DistanceMatrix &unsetDistance(const string &OTU);

    /** @brief Adds a new OTU to this matrix.
     *
     * If given OTU already exists in this matrix, nothing happens.
     * @param[in]       OTU       Given OTU
     * @return This matrix itself
     */
    DistanceMatrix &addOTU(const string &OTU);

    /** @brief Removes an OTU from this matrix.
     *
     * Data about given OTU (i.e. distances) is not removed.
     * @param[in]       OTU       Given OTU
     * @return This matrix itself
     */
    DistanceMatrix &removeOTU(const string &OTU);



 private:
    static const double min_distance;  ///< Minimum representable distance
    static const double max_distance;  ///< Maximum representable distance

    set<string> OTUs;   ///< OTUs in this matrix


    /** @brief Tells whether there is a distance recorder between given OTUs.
     *
     * Order of OTUs matters.
     * @param[in]       a         First OTU
     * @param[in]       b         Second OTU
     * @retval true     There is a record for given OTUs
     * @retval false    There is no record for given OTUs
     */
    bool isSetDirectional(const string &a, const string &b) const;
};

}  // namespace Phylo
}  // namespace Victor

#endif  // _VICTOR_PHYLO_DISTANCEMATRIX_