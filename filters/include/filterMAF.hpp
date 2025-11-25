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
    signalContainer<T> originalSignal_;
    signalContainer<T> filteredSignal_;
public:
    filterMAF(int windowSize);
    ~filterMAF();

    void applyFilter();

    signalContainer<T>& getOriginalSignalContainerReference() { return this->originalSignal_; };
    signalContainer<T>& getFilteredSignalContainerReference() { return this->filteredSignal_; };

    void setSignal(const std::vector<T>signal) { this->originalSignal_.setSignal(signal); }
    std::vector<T> getSignal() const { return this->filteredSignal_.getSignal(); }
};



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
    if (temp > this->originalSignal_.signal_.size())
        this->windowSize_ = this->originalSignal_.signal_.size();
    size_t n = this->originalSignal_.signal_.size();

    this->filteredSignal_.signal_.resize(n);
    this->filteredSignal_.signal_[0] = this->originalSignal_.signal_[0];
    for (size_t i = 1; i < n; ++i)
    {
        while (i < this->windowSize_)
        {
            this->filteredSignal_.signal_[i] = this->filteredSignal_.signal_[i - 1] * (i - 1) / i 
                                                        + this->originalSignal_.signal_[i] / i;
            ++i;
        }
        
        for (size_t j = 0; j < this->windowSize_; ++j)
        {
            if (i - j < 0)
                break;
            this->filteredSignal_.signal_[i] += this->originalSignal_.signal_[i - j] / this->windowSize_;
        }
    }
    this->windowSize_ = temp;
}

} // namespace filters