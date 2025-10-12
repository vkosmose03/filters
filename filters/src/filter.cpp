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
#include "../include/helpers.h"
#include "../include/filter.h"

namespace filters
{
    template <typename T>
    filterHAAR<T>::filterHAAR(const std::vector<T>& signal, HAARthreshold thresholdType, double thresholdValue = 0.0, int filteringwindow = 0)
    : thresholdType_(thresholdType)
    , thresholdValue_(thresholdValue)
    , filteringwindow_(filteringwindow)
    , signalData_(signal)
    {
    }

    template <typename T>
    filterHAAR<T>::~filterHAAR()
    {
    }

    template <typename T>
    filterMAF<T>::filterMAF(const std::vector<T>& signal, int windowSize)
    : signalData_(signal)
    , windowSize_(windowSize)
    {
    }

    template <typename T>
    filterMAF<T>::~filterMAF()
    {
    }

    template <typename T>
    filterEMF<T>::filterEMF(const std::vector<T>& signal, int windowSize, EMFenvironment signalType
              , double std_k, double max_k, double threshold)
    : signalData_(signal)
    , windowSize_(windowSize)
    , signalType_(signalType)
    , std_k_(std_k)
    , max_k_(max_k)
    , threshold_(threshold)
    {
    }

    template <typename T>
    filterEMF<T>::~filterEMF()
    {
    }

    template <typename T>
    filterMedian<T>::filterMedian(const std::vector<T>& signal, int windowSize)
    : signalData_(signal)
    , windowSize_(windowSize)
    {
    }

    template <typename T>
    filterMedian<T>::~filterMedian()
    {
    }
}
