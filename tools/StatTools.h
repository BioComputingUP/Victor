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
 * @descriptipton: defines common statistical operations, e.g. AVG, SD
 */
#ifndef __STAT_TOOLS_H__
#define __STAT_TOOLS_H__

#include <math.h>
#include <vector>
#include <Debug.h>
#include <algorithm>

/**
 * @Description finds the min value in the data vector
 * @param vector<T> data
 * @return min value (template <class T> inline double )
 */
template <class T> inline double minimum(vector<T> data) {
    T min = data[0];
    for (unsigned int i = 1; i < data.size(); i++)
        if (min > data[i])
            min = data[i];
    return min;
}

/**
 * @Description finds the maximum value in the data vector
 * @param vector<T> data
 * @return max value (template <class T> inline double )
 */
template <class T> inline double maximum(vector<T> data) {
    T max = data[0];
    for (unsigned int i = 1; i < data.size(); i++)
        if (max < data[i])
            max = data[i];
    return max;
}

/**
 * @Description returns an integer vector with local maxima positions evalued from a data vector;
 * @param  vector<T> data
 * @return local max value (template <class T> inline vector<int>)
 */
template <class T> inline vector<int> localMaxima(vector<T> data) {
    vector<T> tmp;
    tmp.reserve(data.size());
    vector<int> maxima;

    for (unsigned int i = 0; i < data.size() - 1; i++) {
        tmp.push_back(data[i + 1] - data[i]);
    }

    for (unsigned int i = 0; i < tmp.size() - 1; i++) {
        if (tmp[i] > 0 && tmp[i + 1] < 0)
            maxima.push_back(i + 1);
    }

    return maxima;
}

/**
 * @Description returns an integer vector with local minima positions evalued from a data vector;
 * @param  vector<T> data
 * @return local min value (template <class T> inline vector<int>)
 */
template <class T> inline vector<int> localMinima(vector<T> data) {
    vector<T> tmp;
    tmp.reserve(data.size());
    vector<int> minima;

    for (unsigned int i = 0; i < data.size() - 1; i++) {
        tmp.push_back(data[i + 1] - data[i]);
    }

    for (unsigned int i = 0; i < tmp.size() - 1; i++) {
        if (tmp[i] < 0 && tmp[i + 1] > 0)
            minima.push_back(i + 1);
    }

    return minima;
}

/**
 * @Description returns the value for the average  evalued from a data vector;
 * @param  vector<T> data
 * @return average value (template <class T> inline double)
 */
template <class T> inline double average(vector<T> data) {
    T avg = 0;
    for (unsigned int i = 0; i < data.size(); i++)
        avg += data[i];
    avg = avg / data.size();
    return avg;
}

/**
 * @Description returns the value for the standardDeviation  evalued from a data vector;
 * @param  vector<T> data(in) , Avg value (out)
 * @return standard deviation value (template <class T> inline double )
 */
template <class T> inline double standardDeviation(vector<T> data, T avg = -9999) {
    if (avg == -9999)
        avg = average(data);

    T sd = 0;
    for (unsigned int i = 0; i < data.size(); i++)
        sd += (data[i] - avg) * (data[i] - avg);
    sd = sqrt(sd / data.size());
    return sd;
}

/**
 * @Description returns the value for the standardError  evalued from a data vector;
 * @param  vector<T> data(in), standard deviation(out)
 * @return standard error value (template <class T> inline double
 */
template <class T> inline double standardError(vector<T> data, T sd = -9999) {
    if (sd == -9999)
        sd = standardDeviation(data);

    return (sd / sqrt(data.size()));
}

/**
 * @Description returns the value for the Zscore  evalued from a data vector;
 * @param  vector<T> data (in), average (out) , standDesv (out)
 * @return template <class T> inline double
 */
template <class T> inline vector<double> Zscore(vector<T> data, T avg = -9999, T sd = -9999) {
    if (avg == -9999)
        avg = average(data);
    if (sd == -9999)
        sd = standardDeviation(data, avg);

    vector<double> tmp;
    tmp.reserve(data.size());

    for (unsigned int i = 0; i < data.size(); i++) {
        tmp.push_back((data[i] - avg) / sd);
    }

    return tmp;
}

/**
 * @Description returns the value for the windowAverage  evalued from a data vector;
 * @param  vector<T> data (in), window
 * @return template <class T> inline double
 */
template <class T> inline vector<double> windowAverage(vector<T> data, unsigned int w) {
    vector<double> tmp;
    tmp.reserve(data.size());

    for (int i = 0; i < static_cast<int> (data.size()); i++) {
        double val = 0.0;

        for (unsigned int j = i - w; j <= i + w; j++) {
            if (j < 0)
                val = data[0] + val;
            else if (j >= data.size())
                val = data[data.size() - 1] + val;
            else
                val = data[j] + val;
        }

        tmp.push_back(val / (2 * w + 1));
    }

    return tmp;
}

/**
 * @Description returns the value for the pearsonCorrelation  evalued from a data vector;
 * @param  vector<T> data (in), vector<T> data (in)
 * @return template <class T> inline double
 */
template <class T> inline double pearsonCorrelation(vector<T> data, vector<T> data2) {
    T pearson = 0;
    long double agv1 = average(data);
    long double agv2 = average(data2);
    long double numerator = 0;
    for (unsigned int i = 0; i < data.size(); i++) {
        numerator = numerator + ((data[i] - agv1) * (data2[i] - agv2));
    }
    long double tmpdenom1 = 0;
    long double denom1 = 0;
    for (unsigned int i = 0; i < data.size(); i++) {
        tmpdenom1 = tmpdenom1 + ((data[i] - agv1) * (data[i] - agv1));
    }
    denom1 = sqrt(tmpdenom1);

    long double tmpdenom2 = 0;
    long double denom2 = 0;
    for (unsigned int i = 0; i < data2.size(); i++) {
        tmpdenom2 = tmpdenom2 + ((data2[i] - agv2) * (data2[i] - agv2));
    }
    denom2 = sqrt(tmpdenom2);

    pearson = (numerator / (denom1 * denom2));
    return pearson;
}

/**
 * @Description returns the value for the sRankHelper  evalued from a data vector;
 * @param  n (in), vector<T> data (in)
 * @return template <class T> inline double
 */
template <class T> inline double sRankHelper(T n, vector<T> data) {
    double result;
    for (unsigned int i = 0; i < data.size(); i++) {
        if (n == data[i]) {
            result = static_cast<double> (i);
            return result;
        }
    }
    return -10000.00;
}

/**
 * @Description returns the value for the spearmanCorrelation  evalued from a data vector;
 * @param  vector<T> data (in), vector<T> data (in)
 * @return template <class T> inline double
 */
template <class T> inline double spearmanCorrelation(vector<T> data, vector<T> data2) {
    vector<T> ordData;
    vector<T> ordData2;

    for (unsigned int i = 0; i < data.size(); i++) {
        ordData.push_back(data[i]);
        ordData2.push_back(data2[i]);
    }

    sort(&ordData[0], &ordData[ordData.size() - 1]);
    sort(&ordData2[0], &ordData2[ordData2.size() - 1]);

    T spearman = 0;
    long double D = 0;
    for (unsigned int i = 0; i < data.size(); i++)
        D = ((sRankHelper(data[i], ordData) -
            (sRankHelper(data2[i], ordData2))) *
            (sRankHelper(data[i], ordData) -
            (sRankHelper(data2[i], ordData2))));

    spearman = 1 - ((6 * D) / ((data.size() * data.size() * data.size())
            - data.size()));
    return spearman;
}


// ---------------------------------------------------------------------------
// NAME: getRandomNumber()
// DESC: given an interval [a,b], we generate random number inside it.
// NOTE: before calling this function, please initialize a random seed by
//       using srand( (unsigned)time(0) ) ONCE in your program.
// ---------------------------------------------------------------------------
double getRandomNumber(double low, double high);

// ---------------------------------------------------------------------------
// NAME: getGaussianRandomNumber()
// DESC: given an interval [a,b], we generate a random number inside it based
//       on a gaussian distribution of probability.
// NOTE: before calling this function, please initialize a random seed by
//       using srand( (unsigned)time(0) ) ONCE in your program.
//       Be aware on using as seed the value 12: this is an experimental
//       result that gets almost the best distribution of probability.
// ---------------------------------------------------------------------------
double getGaussianRandomNumber(double low, double high, double seed = 12);


#endif // __STAT_TOOLS_H__
