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
    void insertionSort(std::vector<T>& arr, int left, int right);



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
    void merge(std::vector<T>& arr, int left, int middle, int right);



    /**
     * @brief Sorts an array using TimSort.
     * 
     * @tparam T The data type of the elements in the array.
     * @param arr The array to sort.
     * @param runSize The size of the runs to use for sorting.
     */
    template <typename T>
    void timSort(std::vector<T>& arr, int runSize = 32);
}