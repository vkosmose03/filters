/*
    * filterEMF.hpp
    *
    * This file is part of the filters library.
    * It defines the filterEMF class for applying exponential moving filters.
    *
    * Author: Vitaly Makarov
    * Date: 2025-10-17
*/



#pragma once
#include <vector>
#include <cmath>
#include "filterTypes.hpp"
#include "signalContainer.hpp"
#include "helpers.hpp"
#include "filter.hpp"



namespace filters
{
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
 * @tparam T The data type of the signal (e.g., float, double).
 * @param signal - signal data to be filtered.
 * @param windowSize - the size of the moving window for averaging.
 * @param signalType - the type of the signal (e.g., physical, mathematical).
 * @param std_k - the standard smoothing factor for the filter.
 * @param max_k - the maximum smoothing factor for the filter.
 * @param threshold - the threshold value for the filter.
 */
template <typename T>
class filterEMF : public filterBase<T>
{
private:
    signalContainer<T> originSignal_;
    signalContainer<T> filteredSignal_;
    EMFenvironment signalType_;
    double k_;
    double std_k_, max_k_;
    double threshold_;
public:
    filterEMF(signalContainer<T> originSignal, EMFenvironment signalType, double k = 0.5
              , double std_k = 0.5, double max_k = 0.5, double threshold = 0.0);
    filterEMF(EMFenvironment signalType, double k = 0.5
              , double std_k = 0.5, double max_k = 0.5, double threshold = 0.0);
    ~filterEMF();

    void applyFilter();

    signalContainer<T>& getOriginalSignalContainerReference();
    signalContainer<T>& getFilteredSignalContainerReference();
};



/**
 * @brief Constructor for the filterEMF class.
 * @param originSignal The original signal container.
 * @param signalType The type of the signal (e.g., physical, mathematical).
 * @param k The smoothing factor for the filter.
 * @param std_k The standard smoothing factor for the filter.
 * @param max_k The maximum smoothing factor for the filter.
 * @param threshold The threshold value for the filter.
 */
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



/**
 * @brief Destructor for the filterEMF class.
 */
template <typename T>
filterEMF<T>::~filterEMF()
{
}



/**
 * @brief Applies the Exponential Moving Filter (EMF) to the original signal.
 */
template <typename T>
void filterEMF<T>::applyFilter()
{
    if ((this->max_k_ + this->std_k_ != 1.0) && 
    (this->signalType_ == EMFenvironment::UNDEFINED || this->signalType_ == EMFenvironment::RADIOTECHNICAL))
        return;

    std::vector<T> originalSignal(this->originSignal_.getSignal());
    std::vector<T> filteredSignal;
    size_t n = originalSignal.size();

    filteredSignal.resize(n);
    filteredSignal[0] = originalSignal[0];

    if (this->signalType_ == EMFenvironment::RADIOTECHNICAL)
    {
        for (size_t i = 1; i < n; ++i)
        {
            if (std::abs(originalSignal[i] - filteredSignal[i - 1]) >= this->threshold_)
                filteredSignal[i] =  (1 - this->max_k_) * filteredSignal[i - 1]
                                                    + this->max_k_ * originalSignal[i];
            else
                filteredSignal[i] =  (1 - this->std_k_) * filteredSignal[i - 1]
                                                    + this->std_k_ * originalSignal[i];
        }
    }
    else if (this->signalType_ == EMFenvironment::PHYSICALS)
    {
        double variance = helpers::calculateVariance(originalSignal);
        for (size_t i = 1; i < n; ++i)
        {
            double delta = std::abs(originalSignal[i] - filteredSignal[i - 1]);
            if (delta > variance)
                filteredSignal[i] =  (1 - this->k_ * (variance / delta)) * filteredSignal[i - 1]
                                                    + this->k_ * (variance / delta) * originalSignal[i];
            else
                filteredSignal[i] =  (1 - this->k_) * filteredSignal[i - 1]
                                                    + this->k_ * originalSignal[i];
        }
    }
    else if (this->signalType_ == EMFenvironment::UNDEFINED)
    {
        double variance = helpers::calculateVariance(originalSignal);
        double threshold = 2 * variance;
        for (size_t i = 1; i < n; ++i)
        {
            double delta = std::abs(originalSignal[i] - filteredSignal[i - 1]);
            if (delta > threshold)
                filteredSignal[i] =  (1 - this->max_k_) * filteredSignal[i - 1]
                                                    + this->max_k_ * originalSignal[i];
            else
                filteredSignal[i] =  (1 - this->std_k_) * filteredSignal[i - 1]
                                                    + this->std_k_ * originalSignal[i];
        }
    }

    this->filteredSignal_.setSignal(filteredSignal);
}



/**
 * @brief Provides a reference to the original signal container.
 * @return A reference to the original signal container.
 */
template <typename T>
inline signalContainer<T>& filterEMF<T>::getOriginalSignalContainerReference()
{
    return this->originSignal_;
}



/**
 * @brief Provides a reference to the filtered signal container.
 * @return A reference to the filtered signal container.
 */
template <typename T>
inline signalContainer<T>& filterEMF<T>::getFilteredSignalContainerReference()
{
    return this->filteredSignal_;
}

} // namespace filters