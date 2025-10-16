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
#include <vector>

namespace filters
{
namespace helpers
{
// /**
//  * @brief Sorts a portion of an array using insertion sort.
//  * 
//  * @tparam T The data type of the elements in the array.
//  * @param arr The array to sort.
//  * @param left The starting index of the portion to sort.
//  * @param right The ending index of the portion to sort.
//  */
// template <typename T>
// void insertionSort(std::vector<T>& arr, int left, int right);



// /**
//  * @brief Merges two sorted subarrays into a single sorted array.
//  * 
//  * @tparam T The data type of the elements in the array.
//  * @param arr The array containing the subarrays to merge.
//  * @param left The starting index of the left subarray.
//  * @param middle The ending index of the left subarray.
//  * @param right The ending index of the right subarray.
//  */
// template <typename T>
// void merge(std::vector<T>& arr, int left, int middle, int right);



// /**
//  * @brief Sorts an array using TimSort.
//  * 
//  * @tparam T The data type of the elements in the array.
//  * @param arr The array to sort.
//  * @param runSize The size of the runs to use for sorting.
//  */
// template <typename T>
// std::vector<T> timSort(std::vector<T>& arr, int runSize = 32);

// /**
//  * @brief Calculates the variance of a signal.
//  * 
//  * @tparam T The data type of the elements in the signal.
//  * @param signal The signal data.
//  * @return The variance of the signal.
//  */
// template <typename T>
// double calculateVariance(const std::vector<T>& signal);

template <typename T>
void insertionSort(std::vector<T>& arr, int left, int right)
{
    for (int i = left + 1; i <= right; ++i)
    {
        T key = arr[i];
        int j = i - 1;
        while (j >= left && arr[j] > key)
        {
            arr[j + 1] = arr[j];
            --j;
        }
        arr[j + 1] = key;
    }
}



template <typename T>
void merge(std::vector<T>& arr, int left, int middle, int right)
{

    int len1 = middle - left + 1, len2 = right - middle;
    std::vector<T> leftV, rightV;
    for (int i = 0; i < len1; i++)
        leftV.push_back(arr[left + i]);
    for (int i = 0; i < len2; i++)
        rightV.push_back(arr[middle + i + 1]);

    int i = 0;
    int j = 0;
    int k = left;

    while (i < len1 && j < len2)
    {
        if (leftV[i] <= rightV[j])
        {
            arr[k] = leftV[i];
            i++;
        }
        else
        {
            arr[k] = rightV[j];
            j++;
        }
        k++;
    }

    while (i < len1)
    {
        arr[k] = leftV[i];
        k++;
        i++;
    }

    while (j < len2)
    {
        arr[k] = rightV[j];
        k++;
        j++;
    }
}



template <typename T>
std::vector<T> timSort(std::vector<T>& arr, int runSize = 32)
{
    int n = arr.size();

    for (int i = 0; i < n; i += runSize)
    {
        insertionSort(arr, i, std::min(i + runSize - 1, n - 1));
    }

    for (int size = runSize; size < n; size *= 2)
    {
        for (int left = 0; left < n; left += 2 * size)
        {
            int middle = left + size - 1;
            int right = std::min(left + 2 * size - 1, n - 1);
            merge(arr, left, middle, right);
        }
    }
    return arr;
}



template <typename T>
double calculateVariance(const std::vector<T>& signal)
{
    double mean = std::accumulate(signal.begin(), signal.end(), 0.0) / signal.size();
    double variance = 0.0;
    for (const auto& val : signal)
        variance += std::pow((val - mean), 2);
    variance /=  signal.size();
    return variance;
}

}   // namespace helpers
}   // namespace filters