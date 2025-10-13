/*
    * signalContainer.h
    *
    * This file is part of the filters library.
    * It contain the signal and all characteristics.
    *
    * Author: Vitaly Makarov
    * Date: 2025-10-13
*/

#include "../include/signalContainer.h"
#include <vector>



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
