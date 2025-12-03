/*
    * filter.hpp
    *
    * This file is part of the filters library.
    * It defines the filterHAAR class for applying Haar wavelet transforms.
    *
    * Author: Vitaly Makarov
    * Date: 2025-10-04
*/



#pragma once
#include <cmath>
#include <memory>
#include "signalContainer.hpp"



namespace filters
{
/**
 * @brief A base class for different types of filters.
 * 
 * @tparam T The data type of the signal (e.g., float, double).
 */
template <typename T>
class filterBase
{
public:
    virtual void applyFilter() = 0;
    virtual signalContainer<T>& getOriginalSignalContainerReference() = 0;
    virtual signalContainer<T>& getFilteredSignalContainerReference() = 0;
    virtual void setSignal(const std::vector<T> signal) = 0;
    virtual std::vector<T> getSignal() const = 0;

    virtual ~filterBase() = default;
};


/**
 * @brief A class for chaining multiple filters together.
 * 
 * @tparam T The data type of the signal (e.g., float, double).
 */
template <typename T>
class filterChain
{
private:
    std::vector<std::unique_ptr<filterBase<T>>> filters_;
    signalContainer<T> originalSignal_;
    signalContainer<T> filteredSignal_;
public:
    filterChain();
    ~filterChain();

    void appendFilter(std::unique_ptr<filterBase<T>>& filter);
    std::unique_ptr<filterBase<T>>& operator[](size_t index);
    void clearFilters();
    void removeFilter(size_t index);
    void applyFilters();

    void setSignal(const signalContainer<T> signal);
    signalContainer<T> getFilteredSignal() const;

    signalContainer<T>& getOriginalSignalReference();
    signalContainer<T>& getFilteredSignalReference();
};



template <typename T>
filterChain<T>::filterChain()
{
}



template <typename T>
filterChain<T>::~filterChain()
{
}



template <typename T>
inline void filterChain<T>::appendFilter(std::unique_ptr<filterBase<T>>& filter)
{
    this->filters_.push_back(std::move(filter));
}



template <typename T>
inline void filterChain<T>::clearFilters()
{
    this->filters_.clear();
}



template <typename T>
void filterChain<T>::removeFilter(size_t index)
{
    if (filters_.size() < index + 1)
        return;
    auto begin = this->filters_.begin();
    this->filters_.erase(begin + index);
}



template <typename T>
void filterChain<T>::applyFilters()
{
    if (this->filters_.size() == 0) return;

    this->filteredSignal_ = this->originalSignal_;

    for (std::unique_ptr<filterBase<T>>& filter : this->filters_)
    {
        filter->getOriginalSignalContainerReference() = this->filteredSignal_;
        filter->applyFilter();
        this->filteredSignal_ = filter->getFilteredSignalContainerReference();
    }
}



template <typename T>
inline void filterChain<T>::setSignal(const signalContainer<T> signal)
{
    this->originalSignal_ = signal;
}



template <typename T>
inline signalContainer<T>& filterChain<T>::getOriginalSignalReference()
{
    return this->originalSignal_;
}



template <typename T>
inline signalContainer<T>& filterChain<T>::getFilteredSignalReference()
{
    return this->filteredSignal_;
}

} // namespace filters

#include "filterTypes.hpp"
#include "filterHAAR.hpp"
#include "filterMAF.hpp"
#include "filterEMF.hpp"
#include "filterMedian.hpp"
#include "approximation.hpp"

