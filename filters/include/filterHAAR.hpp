/*
    * filterHAAR.hpp
    *
    * This file is part of the filters library.
    * It defines the filterHAAR class for applying Haar wavelet transforms.
    *
    * Author: Vitaly Makarov
    * Date: 2025-10-18
*/



#pragma once
#include <cmath>
#include "filterTypes.hpp"
#include "signalContainer.hpp"
#include "helpers.hpp"
#include "filter.hpp"



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
class filterHAAR : public filterBase<T>
{
private:
    HAARthreshold thresholdType_;
    double thresholdValue_;
    signalContainer<T> originalSignal_;
    signalContainer<T> filteredSignal_;
    int filteringWindow_;
    int depth_;

    void HAARdeconstruction(const std::vector<T>& signal, std::vector<double>& coeffs, int depth);
    void HAARreconstruction(std::vector<T>& signal, const std::vector<double>& coeffs);
    T thresholding(T value);
public:
    filterHAAR(HAARfilterSettings settings);
    ~filterHAAR();

    void applyFilter();

    signalContainer<T>& getOriginalSignalContainerReference() { return this->originalSignal_; };
    signalContainer<T>& getFilteredSignalContainerReference() { return this->filteredSignal_; };

    void setSignal(const std::vector<T>signal) { this->originalSignal_.setSignal(signal); }
    std::vector<T> getSignal() const { return this->filteredSignal_.getSignal(); }
};



/**
 * @brief Constructor for the filterHAAR class.
 * @param thresholdType The type of thresholding to apply (soft or hard).
 * @param thresholdValue The threshold value for filtering.
 * @param filteringWindow The size of the filtering window.
 */
template <typename T>
filterHAAR<T>::filterHAAR(HAARfilterSettings settings)
: thresholdType_(settings.thresholdType)
, thresholdValue_(settings.thresholdValue)
, filteringWindow_(settings.filteringWindow)
, depth_(settings.depth)
{
}



/**
 * @brief Destructor for the filterHAAR class.
 */
template <typename T>
filterHAAR<T>::~filterHAAR()
{
}



/**
 * @brief Applies thresholding to a given value based on the selected threshold type.
 * 
 * @param value The value to be thresholded.
 * @return The thresholded value.
 */
template <typename T>
T filterHAAR<T>::thresholding(T value)
{
    if (this->thresholdType_ == HAARthreshold::SOFT)
    {
        if (std::abs(value) >= thresholdValue_)
            return (value > 0 ? value - thresholdValue_ : value + thresholdValue_);
        return 0;
    }
    else
    {
        return (std::abs(value) >= this->thresholdValue_ ? value : 0);
    }
}



/**
 * @brief Applies the Haar wavelet deconstruction to the signal.
 * 
 * @param signal The input signal to be deconstructed.
 * @param coeffs The output coefficients after deconstruction.
 */
template <typename T>
void filterHAAR<T>::HAARdeconstruction(const std::vector<T>& signal, std::vector<double>& coeffs, int depth)
{
    size_t n = signal.size();
    if (n < 2)
        return;
    
    size_t half = (n / std::pow(2.0, depth));

    for (size_t i = 0; i < half; ++i)
    {
        coeffs[2 * i] = (signal[2 * i] + signal[2 * i + 1]) / std::sqrt(2.0);
        coeffs[2 * i + 1] = (signal[2 * i] - signal[2 * i + 1]) / std::sqrt(2.0);
    }
}



/**
 * @brief Applies the Haar wavelet reconstruction to the signal.
 * 
 * @param signal The output signal after reconstruction.
 * @param coeffs The input coefficients to be reconstructed.
 */
template <typename T>
void filterHAAR<T>::HAARreconstruction(std::vector<T>& signal, const std::vector<double>& coeffs)
{
    size_t n = coeffs.size();
    if (n < 2)
        return;

    size_t half = n / 2;

    signal.resize(n);
    for (size_t i = 0; i < half; ++i)
    {
        signal[2 * i] = (coeffs[2 * i] + coeffs[2 * i + 1]) / std::sqrt(2.0);
        signal[2 * i + 1] = (coeffs[2 * i] - coeffs[2 * i + 1]) / std::sqrt(2.0);
    }
}



/**
 * @brief Applies the Haar wavelet transform filter to the original signal.
 */
template <typename T>
void filterHAAR<T>::applyFilter()
{
    std::vector<T> originalSignal(this->originalSignal_.getSignal());
    std::vector<T> filteredSignal;
    if (originalSignal.size() < 2)
        return;

    double p = std::pow(originalSignal.size() / std::pow(2.0, this->depth_));
    if (p != std::floor(p))
        return;

    filteredSignal = originalSignal;
    size_t n = filteredSignal.size();
    size_t pow2 = 1;
    while (pow2 < n) pow2 *= 2;
    if (pow2 != n)
        filteredSignal.resize(pow2, filteredSignal.back());
    
    if (this->thresholdType_ == HAARthreshold::SOFT)
    {
        double variance = helpers::calculateVariance(originalSignal);
        this->thresholdValue_ = std::sqrt(variance * 2 * std::log10(n));
    }


    std::vector<double> coeffs;
    coeffs.resize(n);
    int i = 0;
    for (i; i < this->depth_; ++i)
    {
        HAARdeconstruction(filteredSignal, coeffs, i);

        size_t half = coeffs.size() / std::pow(2.0, (i + 1));
        for (size_t i = half; i < half * 2; ++i)
        {
            coeffs[i] = thresholding(coeffs[i]);
        }

        filteredSignal = coeffs;
    }



    for (i; i >= 0; --i)
    {
        std::vector<T> coeffsBuffer;
        coeffsBuffer.resize(n / std::pow(2.0, (i + 1)));
        for (size_t j = 0; j < coeffsBuffer.size(); ++j)
        {
            coeffsBuffer[j] = coeffs[j];
        }

        HAARreconstruction(filteredSignal, coeffsBuffer);

        for (size_t j = 0; j < coeffsBuffer.size(); ++j)
        {
           filteredSignal[j] = coeffsBuffer[j];
        }
    }
}

} // namespace filters