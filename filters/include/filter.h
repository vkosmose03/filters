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
 * 
 * @memberof thresholdType_ - is the type of thresholding to apply (soft or hard).
 * @memberof thresholdValue_ - is the threshold value for filtering.
 * @memberof filteringwindow_ - is the size of the filtering window.
 * @memberof signalData_ - is the container for the signal data.
 * 
 * @tparam T The data type of the signal (e.g., float, double).
 * @param signal - signal data to be filtered.
 * @param thresholdType The type of thresholding to apply (soft or hard).
 * @param thresholdValue The threshold value for filtering.
 * @param filteringwindow The size of the filtering window.
 */
template <typename T>
class filterHAAR
{
private:
    HAARthreshold thresholdType_;
    double thresholdValue_;
    signalContainer<T> originSignal_;
    signalContainer<T> filteredSignal_;
    int filteringwindow_;
    void HAARdeconstruction(std::vector<double>& coeffs);
    void HAARreconstruction(const std::vector<double>& coeffs);
    T thresholding(T value, T threshold);
public:
    filterHAAR(const std::vector<T>& signal, HAARthreshold thresholdType, double thresholdValue = 0.0, int filteringwindow = 0);
    ~filterHAAR();
    void applyFilter();
    std::vector<T> getFilteredSignal() const;
};



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
class filterMAF
{
private:
    int windowSize_;
    signalContainer<T> originSignal_;
    signalContainer<T> filteredSignal_;
public:
    filterMAF(const std::vector<T>& signal, int windowSize);
    ~filterMAF();
    void applyFilter();
    std::vector<T> getFilteredSignal() const;
};



/**
 * @brief A class for applying Exponential Moving Filter (EMF) to signal.
 *        In contrast to MAF this filter will give to each new value certain mass.
 * 
 * @memberof windowSize_ - is the size of the moving window for averaging.
 * @memberof signalData_ - is the container for the signal data.
 * @memberof signalType_ - is the type of the signal (e.g., physical, mathematical).
 * @memberof k - is the smoothing factor for the filter.
 * @memberof std_k_ - is the standard smoothing factor for the filter.
 * @memberof max_k_ - is the maximum  smoothing factor for the filter.
 * @memberof threshold_ - is the threshold value for the filter
 *           is value from witch max_k will be used with signalType_ = mathematical or undefined.
 * 
 * @tparam T The data type of the signal (e.g., float, double).
 * @param signal - signal data to be filtered.
 * @param windowSize - the size of the moving window for averaging.
 * @param signalType - the type of the signal (e.g., physical, mathematical).
 * @param std_k - the standard smoothing factor for the filter.
 * @param max_k - the maximum smoothing factor for the filter.
 * @param threshold - the threshold value for the filter.
 */
template <typename T>
class filterEMF
{
private:
    int windowSize_;
    signalContainer<T> originSignal_;
    signalContainer<T> filteredSignal_;
    EMFenvironment signalType_;
    double k;
    double std_k_, max_k_;
    double threshold_;
public:
    filterEMF(const std::vector<T>& signal, int windowSize, EMFenvironment signalType
              , double std_k = 0.0, double max_k = 0.0, double threshold = 0.0);
    ~filterEMF();
    void applyFilter();
    std::vector<T> getFilteredSignal() const;
};



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
class filterMedian
{
private:
    int windowSize_;
    signalContainer<T> originSignal_;
    signalContainer<T> filteredSignal_;
public:
    filterMedian(const std::vector<T>& signal, int windowSize);
    ~filterMedian();
    void applyFilter();
    std::vector<T> getFilteredSignal() const;
};
} // namespace filters