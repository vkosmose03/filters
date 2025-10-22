/*
    * filterMAF.hpp
    *
    * This file is part of the filters library.
    * It defines the filterMAF class for applying Moving Average Filters.
    *
    * Author: Vitaly Makarov
    * Date: 2025-10-18
*/



#pragma once
#include <vector>
#include <cmath>
#include "signalContainer.hpp"
#include "filter.hpp"



namespace filters
{
/**
 * @brief A class for applying Moving Average Filter (MAF) to signals.
 *        This filter smooths the signal by averaging over a specified window size.
 * 
 * @memberof windowSize_ - is the size of the moving window for averaging.
 * @memberof signalData_ - is the container for the signal data.
 *
 * @tparam T The data type of the signal (e.g., float, double).
 * @param signal - signal data to be filtered.
 * @param windowSize The size of the moving window for averaging.
 */
template <typename T>
class filterMAF : public filterBase<T>
{
private:
    int windowSize_;
    signalContainer<T> originSignal_;
    signalContainer<T> filteredSignal_;
public:
    filterMAF(signalContainer<T> signal, int windowSize);
    filterMAF(int windowSize);
    ~filterMAF();

    void applyFilter();

    signalContainer<T>& getOriginalSignalContainerReference();
    signalContainer<T>& getFilteredSignalContainerReference();
};


/**
 * @brief Constructor for the filterMAF class.
 * @param originSignal The original signal container.
 * @param windowSize The size of the moving window for averaging.
 */
template <typename T>
filterMAF<T>::filterMAF(signalContainer<T> signal, int windowSize)
: originSignal_(signal)
, filteredSignal_()
, windowSize_(windowSize)
{
}



/**
 * @brief Constructor for the filterMAF class.
 * @param windowSize The size of the moving window for averaging.
 */
template <typename T>
filterMAF<T>::filterMAF(int windowSize)
: windowSize_(windowSize)
{
}



/**
 * @brief Destructor for the filterMAF class.
 */
template <typename T>
filterMAF<T>::~filterMAF()
{
}



/**
 * @brief Applies the Moving Average Filter to the original signal.
 */
template <typename T>
void filterMAF<T>::applyFilter()
{
    if (this->windowSize_ < 1)
        return;

    int temp = this->windowSize_;
    if (temp > this->originSignal_.signal_.size())
        this->windowSize_ = this->originSignal_.signal_.size();
    size_t n = this->originSignal_.signal_.size();

    this->filteredSignal_.signal_.resize(n);
    this->filteredSignal_.signal_[0] = this->originSignal_.signal_[0];
    for (size_t i = 1; i < n; ++i)
    {
        while (i < this->windowSize_)
        {
            this->filteredSignal_.signal_[i] = this->filteredSignal_.signal_[i - 1] * (i - 1) / i 
                                                        + this->originSignal_.signal_[i] / i;
            ++i;
        }
        
        for (size_t j = 0; j < this->windowSize_; ++j)
        {
            if (i - j < 0)
                break;
            this->filteredSignal_.signal_[i] += this->originSignal_.signal_[i - j] / this->windowSize_;
        }
    }
    this->windowSize_ = temp;
}



/**
 * @brief Provides a reference to the original signal container.
 * @return A reference to the original signal container.
 */
template <typename T>
inline signalContainer<T>& filterMAF<T>::getOriginalSignalContainerReference()
{
    return this->originSignal_;
}



/**
 * @brief Provides a reference to the filtered signal container.
 * @return A reference to the filtered signal container.
 */
template <typename T>
inline signalContainer<T>& filterMAF<T>::getFilteredSignalContainerReference()
{
    return this->filteredSignal_;
}

} // namespace filters