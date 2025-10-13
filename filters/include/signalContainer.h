/*
    * filter.h
    *
    * This file is part of the filters library.
    * It contain the signal and all characteristics.
    *
    * Author: Vitaly Makarov
    * Date: 2025-10-04
*/

#pragma once
#include <vector>



/**
 * @brief A class for storing and processing a signal.
 * This class provides methods for filtering and analyzing the signal data.
 * @memberof timeStamp_ - is the timestamp of the signal - discrete time.
 * @memberof signal_ - is the raw signal data.
 * @memberof average_ - is the average value of the signal.
 * @memberof variance_ - is the variance of the signal.
 * @memberof stdDeviation_ - is the standard deviation of the signal.
 * @memberof max_ - is the maximum value of the signal.
 * @memberof min_ - is the minimum value of the signal.
 * @memberof sum_ - is the sum of the raw signal data.
 *
 * @tparam T The data type of the signal (e.g., float, double).
 * @param signal The raw signal data.
 */
template <typename T>
class signalContainer
{
private:
    std::vector<T> signal_;
    double timeStamp_;
    T average_;
    T variance_;
    T stdDeviation_;
    T max_;
    T min_;
    T sum_;
    void calculateCharacteristics();
public:
    signalContainer();
    signalContainer(const std::vector<T>& signal);
    ~signalContainer();
    void updateSignal(const std::vector<T>& newSignal);
    T getAverage() const;
    T getVariance() const;
    T getStdDeviation() const;
    T getMax() const;
    T getMin() const;
    T getSum() const;
    std::vector<T> getSignal() const;
};