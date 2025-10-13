/*
    * helpers.h
    *
    * This file is part of the filters library.
    * It defines helpers functions to write filters functions.
    *
    * Author: Vitaly Makarov
    * Date: 2025-10-12
*/

#pragma once
#include "../include/helpers.h"
#include "../include/filter.h"

namespace filters
{
    template <typename T>
    filterHAAR<T>::filterHAAR(const std::vector<T>& signal, HAARthreshold thresholdType, double thresholdValue = 0.0, int filteringwindow = 0)
    : thresholdType_(thresholdType)
    , thresholdValue_(thresholdValue)
    , filteringwindow_(filteringwindow)
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
    filterEMF<T>::filterEMF(const std::vector<T>& signal, int windowSize, EMFenvironment signalType
              , double std_k, double max_k, double threshold)
    : originSignal_(signal)
    , filteredSignal_()
    , windowSize_(windowSize)
    , signalType_(signalType)
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
    filterMedian<T>::filterMedian(const std::vector<T>& signal, int windowSize)
    : originSignal_(signal)
    , filteredSignal_()
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
    T filterHAAR<T>::thresholding(T value, T threshold)
    {
        if (std::abs(value) > threshold)
            return (value > 0 ? value - threshold : value + threshold);
        return value;
    }


    /**
     * @brief Applies the Haar wavelet deconstruction to the signal.
     * 
     * @param signal The input signal to be deconstructed.
     * @param coeffs The output coefficients after deconstruction.
     */
    template <typename T>
    void filterHAAR<T>::HAARdeconstruction(std::vector<double>& coeffs)
    {
        size_t n = this->originSignal_.signal_.size();
        if (n < 2)
            return;
        
        size_t half = n / 2;

        coeffs.resize(n);
        for (size_t i = 0; i < half; ++i)
        {
            coeffs[i] = (this->originSignal_.signal_[2 * i] + this->originSignal_.signal_[2 * i + 1]) / std::sqrt(2.0);
            coeffs[half + i] = (this->originSignal_.signal_[2 * i] - this->originSignal_.signal_[2 * i + 1]) / std::sqrt(2.0);
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
            this->filteredSignal_.signal_[i] = (coeffs[2 * i] + coeffs[2 * i + 1]) / std::sqrt(2.0);
            this->filteredSignal_.signal_[half + i] = (coeffs[2 * i] - coeffs[2 * i + 1]) / std::sqrt(2.0);
        }
    }


    template <typename T>
    void filterHAAR<T>::applyFilter()
    {
        if (this->originSignal_.signal_.size() < 2)
            return;

        this->filteredSignal_ = signalContainer<T>(this->originSignal_.signal_);
        size_t n = this->filteredSignal_.signal_.size();
        size_t pow2 = 1;
        while (pow2 < n) pow2 *= 2;
        if (pow2 != n)
            this->filteredSignal_.signal_.resize(pow2, this->filteredSignal_.signal_.back());
        
        if (this->thresholdType_ == HAARthreshold::SOFT)
        {
            double mean = std::accumulate(this->filteredSignal_.signal_.begin(), this->filteredSignal_.signal_.end(), 0.0) / this->filteredSignal_.signal_.size();
            double variance = 0.0;
            for (const auto& val : this->filteredSignal_.signal_)
                variance += std::pow((val - mean), 2);
            variance /=  this->filteredSignal_.signal_.size();
            this->thresholdValue_ = std::sqrt(variance) / 2.0;
        }

        std::vector<double> coeffs;
        HAARdeconstruction(coeffs);

        size_t half = coeffs.size() / 2;
        for (size_t i = half; i < coeffs.size(); ++i)
        {
            coeffs[i] = thresholding(coeffs[i], this->thresholdValue_);
        }

        HAARreconstruction(coeffs);
    }

    template <typename T>
    void filterMAF<T>::applyFilter()
    {
        if (this->originSignal_.signal_.size() < this->windowSize_ || this->windowSize_ < 1)
            return;

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
            
            this->filteredSignal_.signal_[i] = this->filteredSignal_.signal_[i - 1] * (this->windowSize_ - 1)/ this->windowSize_
                                                + this->originSignal_.signal_[i] / this->windowSize_;
        }
    }
}// namespace filters
