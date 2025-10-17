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
#include <memory>
#include "filterTypes.hpp"
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
    virtual signalContainer<T>& getOriginSignalContainerReference();
    virtual signalContainer<T>& getFilteredSignalContainerReference();
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
public:
    filterChain();
    ~filterChain();

    void addFilter(std::unique_ptr<filterBase<T>> filter);
    std::unique_ptr<filterBase<T>>& operator[](size_t index);
    void clearFilters();
    void removeFilter(size_t index);
    void applyFilters();
};




} // namespace filters