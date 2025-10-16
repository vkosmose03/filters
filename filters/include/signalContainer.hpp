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
#include <numeric>
#include <cmath>
#include <algorithm>


namespace filters
{
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
    std::vector<T> getSignal();
    void setSignal(std::vector<T> signal);
};

/**
 * @brief Updates the signal with new data and recalculates its characteristics.
 * @param Signal The signal to be contained.
 */
template <typename T>
signalContainer<T>::signalContainer(const std::vector<T>& signal)
{
    size_t n = signal.size();
    if (n == 0)
        return;
    
    for (size_t i = 0; i < n; ++i)
    {
        signal_.push_back(signal[i]);
    }
    this->calculateCharacteristics();
}



/**
 * @brief Default constructor for the signalContainer class.
 */
template <typename T>
signalContainer<T>::signalContainer()
: timeStamp_(0.0)
, average_(0)
, variance_(0)
, stdDeviation_(0)
, max_(0)
, min_(0)
, sum_(0)
{
}



/**
 * @brief Default constructor for the signalContainer class.
 */
template <typename T>
signalContainer<T>::~signalContainer()
{
}




/**
 * @brief Calculates the characteristics of the signal such as average, variance, standard deviation, max, min, and sum.
 */
template <typename T>
void signalContainer<T>::calculateCharacteristics()
{
    if (this->signal_.empty())
        return;

    size_t n = this->signal_.size();
    this->sum_ = std::accumulate(this->signal_.begin(), this->signal_.end(), static_cast<T>(0));
    this->average_ = this->sum_ / n;
    this->variance_ = 0;
    for (const auto& value : this->signal_)
    {
        this->variance_ += std::pow((value - this->average_), 2);
    }
    this->variance_ /= n;
    this->stdDeviation_ = std::sqrt(this->variance_);
    this->max_ = *std::max_element(this->signal_.begin(), this->signal_.end());
    this->min_ = *std::min_element(this->signal_.begin(), this->signal_.end());
}



template <typename T>
void signalContainer<T>::setSignal(std::vector<T> signal)
{
    this->signal_ = signal;
}



template <typename T>
std::vector<T> signalContainer<T>::getSignal()
{
    return this->signal_;
}

}
