/*
    * helpers.cpp
    *
    * This file is part of the filters library.
    * It defines helpers functions to write filters functions.
    *
    * Author: Vitaly Makarov
    * Date: 2025-10-12
*/

#include <numeric>
#include <cmath>
#include "../include/helpers.hpp"

namespace filters
{
namespace helpers
{
    template <typename T>
    void insertionSort(std::vector<T>& arr, int left, int right);



    template <typename T>
    void merge(const std::vector<T>& arr, int left, int middle, int right);



    template <typename T>
    std::vector<T> timSort(std::vector<T>& arr, int runSize);



    template <typename T>
    double calculateVariance(const std::vector<T>& signal);
}

}