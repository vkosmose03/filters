/*
    * helpers.h
    *
    * This file is part of the filters library.
    * It defines helpers functions to write filters functions.
    *
    * Author: Vitaly Makarov
    * Date: 2025-10-12
*/

#include <cmath>
#include "../include/helpers.hpp"
#include "../include/filter.hpp"

namespace filters
{
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
        this->originSignal_.signal_ = signal;
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
        size_t n = this->originSignal_.signal_();

        for (size_t i = 0; i < n; ++i)
        {
            if ((i + this->windowSize_) > n)
            {
                std::vector<T> bufferVector(this->originSignal_.signal_.begin() + i, this->originSignal_.signal_.end());
                std::vector<T> sortedVector = helpers::timSort(bufferVector);
                size_t sN = sortedVector;
                size_t middle = std::ceil(sortedVector / 2.0);
                this->filteredSignal_.signal_.push_back(sortedVector[middle]);
            }
            else
            {
                std::vector<T> bufferVector(this->originSignal_.signal_.begin() + i, this->originSignal_.signal_.begin() + i + this->windowSize_);
                std::vector<T> sortedVector = helpers::timSort(bufferVector);
                                size_t sN = sortedVector;
                size_t middle = std::ceil(sortedVector / 2.0);
                this->filteredSignal_.signal_.push_back(sortedVector[middle]);
            }
        }
    }

    template <typename T>
    void filterMedian<T>::setSignal(std::vector<T> signal)
    {
        this->originSignal_.signal_ = signal;
    }

}// namespace filters
