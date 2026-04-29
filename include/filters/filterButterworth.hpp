/*
    * filterButterworth.hpp
    *
    * Realization of digital Butterworth filter (IIR).
    * Uses to smooth signals with save original wave form.
    *
    * Author: Vitaly Makarov
    * Date: 2026-04-29
*/

#pragma once
#include <vector>
#include <cmath>
#include "signalContainer.hpp"
#include "filter.hpp"
#include "filterTypes.hpp"

namespace filters
{
template <typename T>
class filterButterworth : public filterBase<T>
{
private:
    signalContainer<T> originalSignal_;
    signalContainer<T> filteredSignal_;
    ButterworthSettings settings_;
    
    std::vector<T> x_prev_;
    std::vector<T> y_prev_;
    std::vector<T> b_;
    std::vector<T> a_;

    void calculateCoefficients();
public:
    filterButterworth(ButterworthSettings settings);
    ~filterButterworth();

    void applyFilter();

    signalContainer<T>& getOriginalSignalContainerReference() { return this->originalSignal_; };
    signalContainer<T>& getFilteredSignalContainerReference() { return this->filteredSignal_; };

    void setSignal(const std::vector<T>signal) { this->originalSignal_.setSignal(signal); }
    std::vector<T> getSignal() const { return this->filteredSignal_.getSignal(); }
};

template <typename T>
filterButterworth<T>::filterButterworth(ButterworthSettings settings)
: originalSignal_()
, filteredSignal_()
, settings_(settings)
{
    calculateCoefficients();
}

template <typename T>
filterButterworth<T>::~filterButterworth()
{
}

template <typename T>
void filterButterworth<T>::calculateCoefficients()
{
    int n = settings_.order;
    double fc = settings_.cutoffFreq;
    double fs = settings_.sampleRate;
    
    double Wn = std::tan(M_PI * fc / fs);
    
    double k = std::tan(M_PI * fc / fs);
    double k2 = k * k;
    double V = std::pow(1.0 + std::sqrt(2.0) * k + k2, n);
    
    double sqrt2 = std::sqrt(2.0);
    double norm = 1.0 / (1.0 + sqrt2 * k + k2);
    
    b_.resize(3);
    a_.resize(3);
    
    b_[0] = k * k * norm;
    b_[1] = 2.0 * b_[0];
    b_[2] = b_[0];
    
    a_[0] = 1.0;
    a_[1] = 2.0 * (k * k - 1.0) * norm;
    a_[2] = (1.0 - sqrt2 * k + k2) * norm;
}

template <typename T>
void filterButterworth<T>::applyFilter()
{
    std::vector<T> signal = this->originalSignal_.getSignal();
    std::vector<T> filteredSignal(signal.size());
    
    if (signal.size() == 0)
        return;

    if (x_prev_.size() != signal.size()) {
        x_prev_.resize(signal.size(), 0.0);
        y_prev_.resize(signal.size(), 0.0);
    }

    for (size_t i = 0; i < signal.size(); ++i)
    {
        double xn = (i < 2) ? 0.0 : signal[i - 2];
        double yn = (i < 2) ? 0.0 : y_prev_[i - 2];

        filteredSignal[i] = (b_[0] * signal[i] + 
                             b_[1] * x_prev_[i] + 
                             b_[2] * xn - 
                             a_[1] * y_prev_[i] - 
                             a_[2] * yn) / a_[0];

        x_prev_[i] = signal[i];
        y_prev_[i] = filteredSignal[i];
    }

    this->filteredSignal_.setSignal(filteredSignal);
}
} // namespace filters
