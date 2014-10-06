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


#ifndef __Vec_H__
#define __Vec_H__

#include <Debug.h>
#include <iostream>
#include <vector>

using namespace std;

/** vector with range checking.

	Behaves like STL's vector. The [] Selector is range checked,
	unless "NDEBUG" is defined at compile time. This implementation
	differs from Stroustrup's book, since our STL doesn't provide an
	"at" Selector.

	This class offers such an "at" method which *always* performs range
	checking (even if NDEBUG is defined at compile time).

	@see    Stroustrup, "The C++ Programming Language", 3rd Edition, p.53 */

/** @brief    vector with range checking.
 * 
* @Description  
* @This 
 **/
template<class T>
class Vec : public vector<T>
{

public:

	typedef typename vector<T>::size_type size_type;


public:

	Vec(): vector<T>()
	{ }

	Vec(const vector<T> &orig) : vector<T>(orig)
	{ }

	Vec(size_type s) : vector<T>(s)
	{ }

	Vec(size_type s, const T &def) : vector<T>(s, def)
	{ }


// PREDICATES:

	inline T& operator [] (size_type i);
	inline const T& operator [] (size_type i) const;

	inline T& at(size_type i);
	inline const T& at(size_type i) const;

};

	template <class T> istream& operator >> (istream &is, Vec<T> &v);
	template <class T> ostream& operator << (ostream &os, const Vec<T> &v);

// -----------------------------------------------------------------------------
//                                     Vec
// -----------------------------------------------------------------------------

template <class T>
inline T&
Vec<T>::operator [] (size_type i)
{
#ifndef NDEBUG
	if (i >= capacity())
	{
		DUMP(i);
		DUMP(size());
		ERROR("Range error.", exception);
	}
	else
		return *(begin() + i);
#else
	return *(vector<T>::begin() + i);
#endif
}


template <class T>
inline const T&
Vec<T>::operator [] (size_type i) const
{
#ifndef NDEBUG
	if (i >= size())
	{
		DUMP(i);
		DUMP(size());
		ERROR("Range error.", exception);
	}
	else
		return *(begin() + i);
#else
	return *(vector<T>::begin() + i);
#endif
}


/** "at" always performs range checking.
	Unlike operator[] "at" always checks for range errors. */

template <class T>
inline T&
Vec<T>::at(size_type i)
{
	if (i >= vector<T>::size())
	{
		DUMP(i);
		DUMP(size());
		ERROR("Range error.", exception);
	}
	else
		return *(vector<T>::begin() + i);
}


/** "at" always performs range checking.
	Unlike operator[] "at" always checks for range errors. */

template <class T>
inline const T&
Vec<T>::at(size_type i) const
{
	if (i >= vector<T>::size())
	{
		DUMP(i);
		DUMP(size());
		ERROR("Range error.", exception);
	}
	else
		return *(vector<T>::begin() + i);
}


/** Sends a Vec to stream.
	The data format is very simple: the first element is the length of
	the vector followed by a space-seperated list of elements. */

template <class T>
//inline
istream&
operator >> (istream &is, Vec<T> &v)
{
	PRECOND(is, exception);
	typename Vec<T>::size_type size;
	v.clear(); // first remove all existing elements EB 02/2001
	is >> size;
	INVARIANT(is, exception);
	for (typename Vec<T>::size_type i = 0; i < size; i++)
	{
		T element;
		is >> element;
		INVARIANT(is, exception);
		v.push_back(element);
	}

//	v.reserve(size);
//	while (v.size() > size)
//	{
//		v.erase(v.begin());
//	}
//	for (Vec<T>::size_type i = 0; i < size; i++)
//	{
//		is >> v[i];
//	}

	POSTCOND(is, exception);
	return is;
}


/** Writes a Vec to a stream.
	The data format is very simple: the first element is the length of
	the vector followed by a space-seperated list of elements. */

template <class T>
//inline
ostream&
operator << (ostream &os, const Vec<T> &v)
{
	ERROR("Program is commented.", exception);
	/** commented - ST, 07/11/05
	PRECOND(os, exception);
	os << v.size() << "   ";
	copy(v.begin(), v.end(), ostream_iterator<T>(os, " "));
	os << endl;
	POSTCOND(os, exception);
	return os; */
	return os;
}


/** Writes a Vec to a stream.
	The data format is very simple: the first element is the length of
	the vector followed by a space-seperated list of elements. */

template <class T>
inline ostream&
outList(ostream &os, const Vec<T> &v)
{
	ERROR("Program is commented.", exception);
	/** commented - ST, 07/11/05
	PRECOND(os, exception);
	os << v.size() << endl;
	copy(v.begin(), v.end(), ostream_iterator<T>(os, "\n"));
	os << endl;
	POSTCOND(os, exception);
	return os; */
}

#endif
