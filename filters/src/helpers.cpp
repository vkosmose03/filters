/*
    * helpers.cpp
    *
    * This file is part of the filters library.
    * It defines helpers functions to write filters functions.
    *
    * Author: Vitaly Makarov
    * Date: 2025-10-12
*/

#pragma once
#include "../include/helpers.h"

namespace helpers
{
    /**
     * @brief Sorts a portion of an array using insertion sort.
     * 
     * @tparam T The data type of the elements in the array.
     * @param arr The array to sort.
     * @param left The starting index of the portion to sort.
     * @param right The ending index of the portion to sort.
     */
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



    /**
     * @brief Merges two sorted subarrays into a single sorted array.
     * 
     * @tparam T The data type of the elements in the array.
     * @param arr The array containing the subarrays to merge.
     * @param left The starting index of the left subarray.
     * @param middle The ending index of the left subarray.
     * @param right The ending index of the right subarray.
     */
    template <typename T>
    void merge(const std::vector<T>& arr, int left, int middle, int right)
    {

        int len1 = middle - left + 1, len2 = right - middle;
        std::vector<T> left, right;
        for (int i = 0; i < len1; i++)
            left.push_back(arr[left + i]);
        for (int i = 0; i < len2; i++)
            right.push_back(arr[middle + i + 1]);

        int i = 0;
        int j = 0;
        int k = left;

        while (i < len1 && j < len2)
        {
            if (left[i] <= right[j])
            {
                arr[k] = left[i];
                i++;
            }
            else
            {
                arr[k] = right[j];
                j++;
            }
            k++;
        }

        while (i < len1)
        {
            arr[k] = left[i];
            k++;
            i++;
        }

        while (j < len2)
        {
            arr[k] = right[j];
            k++;
            j++;
        }
    }



    /**
     * @brief Sorts an array using TimSort.
     * 
     * @tparam T The data type of the elements in the array.
     * @param arr The array to sort.
     * @param runSize The size of the runs to use for sorting.
     */
    template <typename T>
    void timSort(std::vector<T>& arr, int runSize = 32)
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
    }
}