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
#include <cmath>
#include "filterTypes.hpp"
#include "signalContainer.hpp"
#include "helpers.hpp"



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
    int filteringWindow_;
    int depth_;
    void HAARdeconstruction(std::vector<double>& coeffs, int level);
    void HAARreconstruction(const std::vector<double>& coeffs);
    T thresholding(T value);
public:
    filterHAAR(const std::vector<T>& signal, HAARthreshold thresholdType, double thresholdValue = 0.0, int filteringWindow = 0);
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
    std::vector<T> getFilteredSignal();
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
    signalContainer<T> originSignal_;
    signalContainer<T> filteredSignal_;
    EMFenvironment signalType_;
    double k_;
    double std_k_, max_k_;
    double threshold_;
public:
    filterEMF(EMFenvironment signalType, double k = 0.5
              , double std_k = 0.0, double max_k = 0.0, double threshold = 0.0);
    ~filterEMF();
    void applyFilter();
    void setSignal(std::vector<T> signal);
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
    filterMedian(int windowSize);
    ~filterMedian();
    void applyFilter();
    void setSignal(std::vector<T> signal);
    std::vector<T> getFilteredSignal();
};



template <typename T>
filterHAAR<T>::filterHAAR(const std::vector<T>& signal, HAARthreshold thresholdType, double thresholdValue, int filteringWindow)
: thresholdType_(thresholdType)
, thresholdValue_(thresholdValue)
, filteringWindow_(filteringWindow)
, originSignal_(signal)
, filteredSignal_()
{
}



template <typename T>
filterHAAR<T>::~filterHAAR()
{
}



template <typename T>
filterMAF<T>::filterMAF(const std::vector<T>& signal, int windowSize)
: originSignal_(signal)
, filteredSignal_()
, windowSize_(windowSize)
{
}



template <typename T>
filterMAF<T>::~filterMAF()
{
}



template <typename T>
filterEMF<T>::filterEMF(EMFenvironment signalType
        , double k, double std_k, double max_k, double threshold)
: filteredSignal_()
, signalType_(signalType)
, k_(k)
, std_k_(std_k)
, max_k_(max_k)
, threshold_(threshold)
{
}



template <typename T>
filterEMF<T>::~filterEMF()
{
}



template <typename T>
filterMedian<T>::filterMedian(int windowSize)
: filteredSignal_()
, windowSize_(windowSize)
{
}



template <typename T>
filterMedian<T>::~filterMedian()
{
}



/** 
 * @brief Applies thresholding to a value.
 * 
 * @param value The value to be thresholded.
 * @param threshold The threshold value.
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
void filterHAAR<T>::HAARdeconstruction(std::vector<double>& coeffs, int level)
{
    size_t n = this->originSignal_.signal_.size();
    if (n < 2)
        return;
    
    size_t half = (n / std::pow(2.0, level));

    for (size_t i = 0; i < half; ++i)
    {
        coeffs[2 * i] = (this->filteredSignal_.signal_[2 * i] + this->filteredSignal_.signal_[2 * i + 1]) / std::sqrt(2.0);
        coeffs[2 * i + 1] = (this->filteredSignal_.signal_[2 * i] - this->filteredSignal_.signal_[2 * i + 1]) / std::sqrt(2.0);
    }
}



/**
 * @brief Applies the Haar wavelet reconstruction to the signal.
 * 
 * @param signal The output signal after reconstruction.
 * @param coeffs The input coefficients to be reconstructed.
 */
template <typename T>
void filterHAAR<T>::HAARreconstruction(const std::vector<double>& coeffs)
{
    size_t n = coeffs.size();
    if (n < 2)
        return;

    size_t half = n / 2;

    this->filteredSignal_.signal_.resize(n);
    for (size_t i = 0; i < half; ++i)
    {
        this->filteredSignal_.signal_[2 * i] = (coeffs[2 * i] + coeffs[2 * i + 1]) / std::sqrt(2.0);
        this->filteredSignal_.signal_[2 * i + 1] = (coeffs[2 * i] - coeffs[2 * i + 1]) / std::sqrt(2.0);
    }
}


template <typename T>
void filterHAAR<T>::applyFilter()
{
    if (this->originSignal_.signal_.size() < 2)
        return;

    double p = std::pow(this->originSignal_.signal_.size() / std::pow(2.0, this->depth_));
    if (p != std::floor(p))
        return;

    this->filteredSignal_ = signalContainer<T>(this->originSignal_.signal_);
    size_t n = this->filteredSignal_.signal_.size();
    size_t pow2 = 1;
    while (pow2 < n) pow2 *= 2;
    if (pow2 != n)
        this->filteredSignal_.signal_.resize(pow2, this->filteredSignal_.signal_.back());
    
    if (this->thresholdType_ == HAARthreshold::SOFT)
    {
        double variance = helpers::calculateVariance(this->originSignal_.signal_);
        this->thresholdValue_ = std::sqrt(variance * 2 * std::log10(n));
    }


    std::vector<double> coeffs;
    coeffs.resize(n);
    int i = 0;
    for (i; i < this->depth_; ++i)
    {
        HAARdeconstruction(coeffs, i);

        size_t half = coeffs.size() / std::pow(2.0, (i + 1));
        for (size_t i = half; i < half * 2; ++i)
        {
            coeffs[i] = thresholding(coeffs[i]);
        }

        this->filteredSignal_ = coeffs;
    }



    for (i; i >= 0; --i)
    {
        std::vector<T> coeffsBuffer;
        coeffsBuffer.resize(n / std::pow(2.0, (i + 1)));
        for (size_t j = 0; j < coeffsBuffer.size(); ++j)
        {
            coeffsBuffer[j] = coeffs[j];
        }

        HAARreconstruction(coeffsBuffer);

        for (size_t j = 0; j < coeffsBuffer.size(); ++j)
        {
            this->filteredSignal_[j] = coeffsBuffer[j];
        }
    }
}

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



template <typename T>
void filterEMF<T>::setSignal(std::vector<T> signal)
{
    this->originSignal_.setSignal(signal);
}



template <typename T>
void filterEMF<T>::applyFilter()
{
    if (this->max_k_ + this->std_k_ != 1.0)
        return;

    size_t n = this->originSignal_.signal_.size();

    this->filteredSignal_.signal_.resize(n);
    this->filteredSignal_.signal_[0] = this->std_k_ * this->originSignal_.signal_[0];

    if (this->signalType_ == EMFenvironment::RADIOTECHNICAL)
    {
        for (size_t i = 1; i < n; ++i)
        {
            if (std::abs(this->originSignal_.signal_[i] - this->filteredSignal_.signal_[i - 1]) >= this->threshold_)
                this->filteredSignal_.signal_[i] =  (1 - this->max_k_) * this->filteredSignal_.signal_[i - 1]
                                                    + this->max_k_ * this->originSignal_.signal_[i];
            else
                this->filteredSignal_.signal_[i] =  (1 - this->std_k_) * this->filteredSignal_.signal_[i - 1]
                                                    + this->std_k_ * this->originSignal_.signal_[i];
        }
    }
    else if (this->signalType_ == EMFenvironment::PHYSICALS)
    {
        double variance = helpers::calculateVariance(this->originSignal_.signal_);
        for (size_t i = 1; i < n; ++i)
        {
            double delta = std::abs(this->originSignal_.signal_[i] - this->filteredSignal_.signal_[i - 1]);
            if (delta > variance)
                this->filteredSignal_.signal_[i] =  (1 - this->k_ * (variance / delta)) * this->filteredSignal_.signal_[i - 1]
                                                    + this->k_ * (variance / delta) * this->originSignal_.signal_[i];
            else
                this->filteredSignal_.signal_[i] =  (1 - this->k_) * this->filteredSignal_.signal_[i - 1]
                                                    + this->k_* this->originSignal_.signal_[i];
        }
    }
    else if (this->signalType_ == EMFenvironment::UNDEFINED)
    {
        double variance = helpers::calculateVariance(this->originSignal_.signal_);
        double threshold = 2 * variance;
        for (size_t i = 1; i < n; ++i)
        {
            double delta = std::abs(this->originSignal_.signal_[i] - this->filteredSignal_.signal_[i - 1]);
            if (delta > threshold)
                this->filteredSignal_.signal_[i] =  (1 - this->max_k_) * this->filteredSignal_.signal_[i - 1]
                                                    + this->max_k_ * this->originSignal_.signal_[i];
            else
                this->filteredSignal_.signal_[i] =  (1 - this->std_k_) * this->filteredSignal_.signal_[i - 1]
                                                    + this->std_k_* this->originSignal_.signal_[i];
        }
    }
}



template <typename T>
void filterMedian<T>::applyFilter()
{
    std::vector<T> originalSignal = this->originSignal_.getSignal();
    std::vector<T> filteredSignal;
    size_t n = originalSignal.size();

    for (size_t i = 0; i < n; ++i)
    {
        if (n == 1)
        {
            filteredSignal[0] = originalSignal[0];
        }
        else if ((i + this->windowSize_) > n)
        {
            std::vector<T> bufferVector(originalSignal.begin() + i, originalSignal.end());
            std::vector<T> sortedVector = helpers::timSort(bufferVector);
            size_t sN = sortedVector.size();
            size_t middle = std::ceil(sN / 2.0);
            filteredSignal.push_back(sortedVector[middle]);
        }
        else
        {
            std::vector<T> bufferVector(originalSignal.begin() + i, originalSignal.begin() + i + this->windowSize_);
            std::vector<T> sortedVector = helpers::timSort(bufferVector);
            size_t sN = sortedVector.size();
            size_t middle = std::ceil(sN / 2.0);
            filteredSignal.push_back(sortedVector[middle]);
        }
    }

    this->filteredSignal_.setSignal(filteredSignal);
}



template <typename T>
void filterMedian<T>::setSignal(std::vector<T> signal)
{
    this->originSignal_.setSignal(signal);
}



template <typename T>
std::vector<T> filterMedian<T>::getFilteredSignal()
{
    return this->filteredSignal_.getSignal();
}

} // namespace filters