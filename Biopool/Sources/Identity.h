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


#ifndef __IDENTITY_H__
#define __IDENTITY_H__

// Includes:

#include <iostream>
#include <string>
#include <Debug.h>


namespace Victor { namespace Biopool { 
    
    


    /**
    *    @brief Implements object id. With a "number" a "name" and a "counter".
    */
    class Identity {
    public:


        Identity(const string& name = "", long number = 0);
        Identity(const Identity& orig);

        /* OPERATORS */

        /** Assigment. */
        Identity& operator =(const Identity& orig);

        /** Comparison. */
        bool operator ==(const Identity& other) const;

        bool operator !=(const Identity& other) const {
            return !((*this) == other);
        };

        /** Cast to int. */
        operator long() const;

        /** Cast to string. */
        //  operator string() const;

        /** Cast to string&. */
        operator const string&() const;

        /** Set a new name. */
        void setName(string _name);

        /** Set a new number. */
        void setNumber(long _number);



        /* FRIENDS */

        friend ostream& operator <<(ostream& os, const Identity& id);

        friend istream& operator >>(istream& is, Identity& id);

    private:
        long number;
        string name;

        static long counter;
    };

    // ---------------------------------------------------------------------------
    //                                  Identity
    // -----------------x-------------------x-------------------x-----------------

    inline
    Identity::Identity(const string& s, long n) : name(s) {
        if (n != 0) {
            number = n;
        } else {
            number = ++counter;
        }
    }

    inline
    Identity::Identity(const Identity& orig)
    : name(orig.name) {
        number = ++counter;
    }

    inline
    Identity&
            Identity::operator =(const Identity& orig) {
        name = orig.name;

        return *this;
    }

    inline
    bool
    Identity::operator ==(const Identity& other) const {
        return number == other.number;
    }

    inline
    Identity::operator long() const {
        return number;
    }

    inline
    Identity::operator const string&() const {
        return name;
    }

    inline
    void
    Identity::setName(string _name) {
        name = _name;
    }

    inline
    void
    Identity::setNumber(long _number) {
        number = _number;
        if (_number >= counter)
            counter = _number;
    }

    // ---------------------------------------------------------------------------
    //                                   Friends
    // -----------------x-------------------x-------------------x-----------------

    inline
    ostream&
    operator <<(ostream& os, const Identity& id) {
        return os << "(Id " << id.name << " " << id.number << " ) ";
    }

    inline
    istream&
    operator >>(istream& is, Identity& id) {
        string tag;

        is >> tag;

        if (tag != "(Id") {
            DEBUG_MSG("Identity::operator >> : bad input.");

            is.clear(ios::badbit);
            return is;
        }

        return is >> id.name >> id.number >> tag;
    }

}} // namespace
#endif /* __IDENTITY_H__ */

