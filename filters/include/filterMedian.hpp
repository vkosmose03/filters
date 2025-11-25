/*
    * filterMedian.hpp
    *
    * This file is part of the filters library.
    * It defines the filterMedian class for applying median filters.
    *
    * Author: Vitaly Makarov
    * Date: 2025-10-17
*/



#pragma once
#include <vector>
#include <cmath>
#include "signalContainer.hpp"
#include "helpers.hpp"
#include "filter.hpp"



namespace filters
{
/**
 * @brief A class for applying Median Filter (MF) to signals.
 *        This filter smooths the signal by replacing each value with the median of its neighbors.
 *
 * @memberof windowSize_ - is the size of the moving window for median filtering.
 * @memberof signalData_ - is the container for the signal data.
 *
 * @tparam T The data type of the signal (e.g., float, double).
 * @param windowSize The size of the moving window for median filtering.
 */
template <typename T>
class filterMedian : public filterBase<T>
{
private:
    int windowSize_;
    signalContainer<T> originalSignal_;
    signalContainer<T> filteredSignal_;
public:
    filterMedian(int windowSize);
    ~filterMedian();

    void applyFilter();

    signalContainer<T>& getOriginalSignalContainerReference() { return this->originalSignal_; };
    signalContainer<T>& getFilteredSignalContainerReference() { return this->filteredSignal_; };

    void setSignal(const std::vector<T>signal) { this->originalSignal_.setSignal(signal); }
    std::vector<T> getSignal() const { return this->filteredSignal_.getSignal(); }
};



/**
 * @brief Constructor for the filterMedian class.
 * @param windowSize The size of the moving window for median filtering.
 */
template <typename T>
filterMedian<T>::filterMedian(int windowSize)
: windowSize_(windowSize)
{
}



/**
 * @brief Destructor for the filterMedian class.
 */
template <typename T>
filterMedian<T>::~filterMedian()
{
}



/**
 * @brief Applies the median filter to the original signal.
 */
template <typename T>
void filterMedian<T>::applyFilter()
{
    std::vector<T> originalSignal(this->originalSignal_.getSignal());
    std::vector<T> filteredSignal;
    filteredSignal.resize(0);
    size_t n = originalSignal.size();

    for (size_t i = 0; i < n; ++i)
    {

        if ((i + this->windowSize_) > n)
        {
            std::vector<T> bufferVector(originalSignal.begin() + i, originalSignal.end());
            std::vector<T> sortedVector = helpers::timSort(bufferVector);
            size_t sN = sortedVector.size();
            size_t middle = std::ceil(sN / 2.0);
            if (sN == 1)
                middle = 0;
            if (n < this->windowSize_)
            {
                filteredSignal.push_back(sortedVector[middle]);
            }
            else
            {
                filteredSignal.push_back(filteredSignal.back());
            }
        }
        else if ((i + this->windowSize_) <= n)
        {
            std::vector<T> bufferVector(originalSignal.begin() + i, originalSignal.begin() + i + this->windowSize_);
            std::vector<T> sortedVector = helpers::timSort(bufferVector);
            size_t sN = sortedVector.size();
            size_t middle = std::ceil(sN / 2.0);
            if (sN == 1)
                middle = 0;
            filteredSignal.push_back(sortedVector[middle]);
        }
    }

    this->filteredSignal_.setSignal(filteredSignal);
}

} // namespace filters