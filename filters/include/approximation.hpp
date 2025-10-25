/*
    * approximation.hpp
    *
    * This file is part of the filters library.
    * It defines the filterHAAR class for applying Haar wavelet transforms.
    *
    * Author: Vitaly Makarov
    * Date: 2025-10-20
*/


#pragma once
#include <vector>
#include <cmath>
#include "filter.hpp"
#include "signalContainer.hpp"
#include "filterTypes.hpp"

namespace filters
{
template <typename T>
class approximation : public filterBase<T>
{
private:
    signalContainer<T> originalSignal_;
    signalContainer<T> filteredSignal_;
    bool stabilize_;
    double stabilizeIncline_;
    double maxIncline_;
    int windowSize_;
    ErrorEstimate errorEstimate_;
    LinearizationType type_;
public:
    approximation(approximationSettings settings);
    ~approximation();

    void applyFilter();

    signalContainer<T>& getOriginalSignalContainerReference() { return this->originalSignal_; };
    signalContainer<T>& getFilteredSignalContainerReference() { return this->filteredSignal_; };

    void setSignal(const std::vector<T>signal) { this->originalSignal_.setSignal(signal); }
    std::vector<T> getSignal() const { return this->filteredSignal_.getSignal(); }
};



template <typename T>
approximation<T>::approximation(approximationSettings settings)
: originalSignal_()
, filteredSignal_()
, stabilize_(settings.useStabilization)
, stabilizeIncline_(settings.stabilizeIncline)
, maxIncline_(settings.maxIncline)
, windowSize_(settings.windowSize)
, errorEstimate_(settings.errorEstimate)
, type_(settings.type)
{
}



template <typename T>
approximation<T>::~approximation()
{
}



template <typename T>
void approximation<T>::applyFilter()
{
    std::vector<T> signal = this->originalSignal_.getSignal();
    std::vector<T> approxSignal(signal.size());
    double inclineSum = 0.0;
    int n = signal.size();
    if (n == 0)
        return;

    int windowSize = this->windowSize_;
    if (windowSize <= 0 || windowSize > n)
        windowSize = n;

    if (this->type_ == LinearizationType::LINEAR)
    {
        int startIndex = n - 1;
        while (startIndex >= 0)
        {
            int currentWindowSize = std::min(windowSize, startIndex + 1);
            int windowStart = startIndex - currentWindowSize + 1;
            double incline = 0.0;
            double intercept = 0.0;

            std::vector<T> windowSignal(currentWindowSize);
            for (int i = 0; i < currentWindowSize; ++i) {
                windowSignal[i] = signal[windowStart + i];
            }

            if (this->errorEstimate_ == ErrorEstimate::MSE || this->errorEstimate_ == ErrorEstimate::RMSE)
            {
                double sumX = 0.0, sumY = 0.0, sumXY = 0.0, sumX2 = 0.0;

                for (int i = 0; i < currentWindowSize; ++i)
                {
                    sumX += i;
                    sumY += windowSignal[i];
                    sumXY += i * windowSignal[i];
                    sumX2 += i * i;
                }

                incline = (currentWindowSize * sumXY - sumX * sumY) / (currentWindowSize * sumX2 - sumX * sumX);
                if (std::isnan(incline) || std::isinf(incline)) incline = 0.0;
                intercept = (sumY - incline * sumX) / currentWindowSize;

                for (int i = 0; i < currentWindowSize; ++i) 
                {
                    if (incline > this->maxIncline_)
                    {
                        incline = this->maxIncline_;
                    }
                    else if (incline < -this->maxIncline_)
                    {
                        incline = -this->maxIncline_;
                    }
                    approxSignal[windowStart + i] = incline * i + intercept;
                }
            }
            else if (this->errorEstimate_ == ErrorEstimate::MAE)
            {
                const int maxIterations = 10000;
                const double tolerance = 1e-6;
                double learningRate = 0.0001;

                for (int iter = 0; iter < maxIterations; ++iter) 
                {
                    double gradIncline = 0.0;
                    double gradIntercept = 0.0;

                    for (size_t i = 0; i < currentWindowSize; ++i)
                    {
                        double predicted = incline * i + intercept;
                        double error = windowSignal[i] - predicted;
                        double sign = (error >= 0) ? 1.0 : -1.0;
                        gradIncline += -sign * i;
                        gradIntercept += -sign;
                    }
                    gradIncline /= currentWindowSize;
                    gradIntercept /= currentWindowSize;

                    double newIncline = incline - learningRate * gradIncline;
                    double newIntercept = intercept - learningRate * gradIntercept;

                    if (std::abs(newIncline - incline) < tolerance && std::abs(newIntercept - intercept) < tolerance)
                    {
                        break;
                    }

                    incline = newIncline;
                    intercept = newIntercept;
                }

                for (size_t i = 0; i < currentWindowSize; ++i)
                {
                    if (incline > this->maxIncline_)
                    {
                        incline = this->maxIncline_;
                    }
                    else if (incline < -this->maxIncline_)
                    {
                        incline = -this->maxIncline_;
                    }
                    approxSignal[windowStart + i] = incline * i + intercept;
                }
            }

            inclineSum += incline;
            startIndex -= windowSize;
        }

        if (this->stabilize_)
        {
            double averageIncline = inclineSum / (n / windowSize);
            if (std::abs(averageIncline) < this->stabilizeIncline_)
            {
                double mean = std::accumulate(approxSignal.begin(), approxSignal.end(), 0.0) / n;
                for (size_t i = 0; i < n; ++i)
                {
                    approxSignal[i] = approxSignal[i] - mean;
                }
            }
        }
    }

    this->filteredSignal_.setSignal(approxSignal);
}
} // namespace filters