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
@Description  functions for sorting 

*/

#ifndef _TOOLS_SORT_H_
#define _TOOLS_SORT_H_

#include <stddef.h>
#include <stdlib.h>
#include <vglSwap.h>

/* elementary sorting functions (in-place) */

template <class T> inline void sort2(T& a, T& b)
{
  if (a > b) swap(b,a);
}

template <class T> inline  void sort2(T* f)
{
  if (f[0] > f[1]) swap(f[1], f[0]);
}

template <class T> inline  void sort3(T& a, T& b, T& c)
{
  sort2(a,b);
  if (c < b) {
    swap(c,b);
    sort2(a,b);
  }
}

template <class T> inline  void sort3(T* f)
{
  sort2(f);
  if (f[2] < f[1]) {
    swap(f[2],f[1]);
    sort2(f);
  }
}

template <class T> inline void sort4(T& a, T& b, T& c, T& d)
{
  sort2(a,b);
  sort2(c,d);
  sort2(a,c);
  sort2(b,d);
  sort2(b,c);
}

template <class T> inline void sort4(T* f)
{
  sort2(f[0],f[3]);
  sort2(f+1);
  sort2(f); 
  sort2(f+2);
  sort2(f+1);
}

template <class T> inline void sort5(T* f)
{
  sort3(f);	
  sort2(f+3);		
  sort2(f[0],f[3]);
  sort2(f[2],f[4]);
  sort3(f+1);
}

template <class T> inline void sort6(T* f)
{
  sort3(f);
  sort3(f+3);
  sort2(f[0],f[3]);
  sort2(f[3],f[5]);
  sort4(f+1);
}

template <class T> inline void sort7(T* f)
{
  sort3(f);
  sort4(f+3);
  sort2(f[0],f[3]);
  sort2(f[2],f[6]);
  sort5(f+1);
}

template <class T> inline void sort8(T* f)
{
  sort4(f);
  sort4(f+4);
  sort2(f[0],f[4]);
  sort2(f[3],f[7]);
  sort6(f+1);
}

template <class T> inline void sort9(T* f)
{
  sort3(f);
  sort3(f+3);
  sort3(f+6);
  sort3(f[0],f[3],f[6]);
  sort3(f[2],f[5],f[8]);
  sort7(f+1);
}

/* sorting helper functions */

template <class T> inline void vgSortIn(T* f1, ptrdiff_t len1, 
	 T* f2, ptrdiff_t len2, 
	 T* f)
{
  T* fend1 = f1 + len1;
  T* fend2 = f2 + len2;

  while ((f1 < fend1) && (f2 < fend2)) {
    if (*f2 < *f1) *f++ = *f2++;
    else *f++ = *f1++;
  }
  while (f2 < fend2) {
    *f++ = *f2++;
  }
  while (f1 < fend1) {
    *f++ = *f1++;
  }
}

template <class T> inline ptrdiff_t vgSortOut(T* fout, ptrdiff_t outlen, 
	  T* f, ptrdiff_t len, 
	  T* fres) 
{
  T* fend = f + len;
  T* foutend = fout + outlen;
  while ((fout < foutend) && (f < fend)) {
    if (*fout == *f) fout++;
    else *fres++ = *f;
    f++;
  }
  while (f < fend) {
    *fres++ = *f++;
  }
  return outlen - (foutend - fout);
}

/* sort function, variante #1 */

template <class T>  T* vgSort(T* f1, ptrdiff_t len, T* f2)
{
  if (len <= 9) {		// note: 9 should be the performance cutoff
    T* f = f1;
    switch (len) {
    case 9:
      sort3(f[0],f[1],f[2]); sort3(f[3],f[4],f[5]); sort3(f[6],f[7],f[8]); 
      sort3(f[0],f[3],f[6]); sort3(f[2],f[5],f[8]); f++; // continue with 7!
    case 7: sort3(f); sort4(f+3); sort2(f[0], f[3]); sort2(f[2], f[6]); f++; //continue with 5!
    case 5: sort3(f); sort2(f+3); sort2(f[0], f[3]); sort2(f[2], f[4]); f++; //continue with 3!
    case 3: sort3(f); break;
    case 8: sort4(f); sort4(f+4); sort2(f[0], f[4]); sort2(f[3], f[7]); f++; //continue with 6!
    case 6: sort3(f); sort3(f+3); sort2(f[0], f[3]); sort2(f[2], f[5]); f++; //continue with 4!
    case 4: sort2(f); sort2(f+2); sort2(f[0], f[2]); sort2(f[1], f[3]); f++; //continue with 2!
    case 2: sort2(f); break;
    default: break;
    }
    return f1;			
  }
  else {			// sort algorithm
    ptrdiff_t mid = (len>>1);
    T* ra = sort(f1, mid, f2); // sort 1st half
    T* rb = sort(f1+mid, len-mid, f2+mid); // sort 2nd half

    if ((ra == f1) && (rb == (f1+mid))) { // both sorts in place
      vgSortIn(ra, mid, rb, len-mid, f2); // sort in 2nd buffer
      return f2;		// result buffer is 2nd buffer
    }

    if ((ra == f2) && (rb == (f2+mid))) {	// both sorts in 2nd buffer
      vgSortIn(ra, mid, rb, len-mid, f1); // sort in 1st buffer
      return f1;		// result buffer is 1st buffer
    }

    // if the result buffer are mixed (one portion in "f1", the other
    // in "f2"), all parts are copied to the buffer "f1".
    T* org = (ra == f1) ? f2+mid : f2;
    T* copy = (ra == f1) ? f1+mid : f1;
    ptrdiff_t index = (ra == f1) ? len-mid : mid;
    while (index > 0) {		// copy 2nd buffer portion to 1st buffer
      index--;
      copy[index] = org[index];
    }
    vgSortIn(f1, mid, f1+mid, len-mid, f2); // sort in 2nd buffer
    return f2;
  }
}

/* sort function, variante #2 ("quicksort") */

template <class T> void vgQuicksort(T* f, ptrdiff_t len)
{
  if (len <= 9) {		
    switch (len) {
    case 9:
      sort3(f[0],f[1],f[2]); sort3(f[3],f[4],f[5]); sort3(f[6],f[7],f[8]); 
      sort3(f[0],f[3],f[6]); sort3(f[2],f[5],f[8]); f++; // continue with 7!
    case 7: sort3(f); sort4(f+3); sort2(f[0], f[3]); sort2(f[2], f[6]); f++; //continue with 5!
    case 5: sort3(f); sort2(f+3); sort2(f[0], f[3]); sort2(f[2], f[4]); f++; //continue with 3!
    case 3: sort3(f); break;
    case 8: sort4(f); sort4(f+4); sort2(f[0], f[4]); sort2(f[3], f[7]); f++; //continue with 6!
    case 6: sort3(f); sort3(f+3); sort2(f[0], f[3]); sort2(f[2], f[5]); f++; //continue with 4!
    case 4: sort2(f); sort2(f+2); sort2(f[0], f[2]); sort2(f[1], f[3]); f++; //continue with 2!
    case 2: sort2(f); break;
    default: break;
    }
  }
  else {			// quicksort algorithm (Sedgewick page 148)
    T* left = f;
    T* right = f + len;
    T sample = *f;

    while (1) {
      while (*--right > sample) ;
      while ((++left < right) && (*left < sample)) ;
      if (left >= right) break;
      swap(*left, *right);
    }
    swap(*f, *right);

    vgQuicksort(f, right-f);
    vgQuicksort(right+1, len-(right-f)-1);
  }
}

#endif /* _TOOLS_SORT_H_ */
