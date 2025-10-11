/*
    * filter.h
    *
    * This file is part of the filters library.
    * It defines the filterHAAR class for applying Haar wavelet transforms.
    *
    * Author: Vitaly Makarov
    * Date: 2025-10-04
*/

#pragma once

#include "filterTypes.h"
#include "signalContainer.h"

namespace filters
{
/**
 * @brief A class for applying Haar wavelet transforms to signals.
 *        This filter can clear signal from noise: if u know noise level use hard thresholding,
 *        if u don't know noise level use soft thresholding.
 * @tparam T The data type of the signal (e.g., float, double).
 * @param thresholdType The type of thresholding to apply (soft or hard).
 * @param thresholdValue The threshold value for filtering.
 */
template <typename T>
class filterHAAR
{
private:
    HAARthreshold thresholdType_;
    double thresholdValue_;
    signalContainer<T> signalData_;
public:
    filterHAAR(const std::vector<T>& signal, HAARthreshold thresholdType, double thresholdValue = 0.0);
    ~filterHAAR();
    void applyFilter();
    std::vector<T> getFilteredSignal() const;
};


/**
 * @brief A class for applying Moving Average Filter (MAF) to signals.
 *        This filter smooths the signal by averaging over a specified window size.
 * @tparam T The data type of the signal (e.g., float, double).
 * @param windowSize The size of the moving window for averaging.
 */
template <typename T>
class filterMAF
{
private:
    int windowSize_;
    signalContainer<T> signalData_;
public:
    filterMAF(const std::vector<T>& signal, int windowSize); 
    ~filterMAF();
    void applyFilter();
    std::vector<T> getFilteredSignal() const;
};


}