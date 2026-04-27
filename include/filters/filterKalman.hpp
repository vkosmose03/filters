/*
    * filterKalman.hpp
    *
    * This file is part of the filters library.
    * It defines the filterKalman class for applying Kalman filtering to IMU signals.
 * 
 * Author: Vitaly Makarov
 * Date: 2025-10-20
*/

#pragma once
#include <vector>
#include <cmath>
#include "filterTypes.hpp"
#include "signalContainer.hpp"
#include "filter.hpp"

namespace filters
{
/**
 * @brief A class for applying Kalman Filter to signal.
 *        This filter is ideal for MEMS IMUs to reduce noise and drift.
 * 
 * @tparam T The data type of the signal (e.g., float, double).
 * @param signal - signal data to be filtered.
 * @param processNoise - Process noise covariance (Q).
 * @param measurementNoise - Measurement noise covariance (R).
 * @param estimationError - Initial estimation error covariance (P).
 */
template <typename T>
class filterKalman : public filterBase<T>
{
private:
    signalContainer<T> originalSignal_;
    signalContainer<T> filteredSignal_;
    KalmanSettings settings_;

    T x_;
    T p_;
    T q_;
    T r_;
    T k_;

public:
    filterKalman(KalmanSettings settings);
    ~filterKalman();

    void applyFilter();

    signalContainer<T>& getOriginalSignalContainerReference() { return this->originalSignal_; };
    signalContainer<T>& getFilteredSignalContainerReference() { return this->filteredSignal_; };

    void setSignal(const std::vector<T>signal) { this->originalSignal_.setSignal(signal); }
    std::vector<T> getSignal() const { return this->filteredSignal_.getSignal(); }
};

template <typename T>
filterKalman<T>::filterKalman(KalmanSettings settings)
: originalSignal_()
, filteredSignal_()
, settings_(settings)
, x_(0)
, p_(settings.estimationError)
, q_(settings.processNoise)
, r_(settings.measurementNoise)
{
}

template <typename T>
filterKalman<T>::~filterKalman()
{
}

template <typename T>
void filterKalman<T>::applyFilter()
{
    std::vector<T> originalSignal(this->originalSignal_.getSignal());
    std::vector<T> filteredSignal;
    filteredSignal.resize(originalSignal.size());

    if (originalSignal.size() == 0)
        return;

    for (size_t i = 0; i < originalSignal.size(); ++i)
    {
        T measurement = originalSignal[i];

        T predicted_x = x_; 
        T predicted_p = p_ + q_;

        T k = predicted_p / (predicted_p + r_);
        x_ = predicted_x + k * (measurement - predicted_x);
        p_ = (1 - k) * predicted_p;

        filteredSignal[i] = x_;
    }

    this->filteredSignal_.setSignal(filteredSignal);
}

} // namespace filters
