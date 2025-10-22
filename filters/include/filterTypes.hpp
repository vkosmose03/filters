/*
    * filterTypes.h
    *
    * This file is part of the filters library.
    * It defines minor filters settings (e.g. thresholds and types).
    *
    * Author: Vitaly Makarov
    * Date: 2025-10-09
*/

#pragma once



namespace filters
{
enum class HAARthreshold
{
    SOFT = 1,
    HARD = 2
};



enum class EMFenvironment
{
    PHYSICALS = 1,
    RADIOTECHNICAL = 2,
    UNDEFINED = 3
};



enum class ErrorEstimate
{
    MAE = 1,
    MSE = 2,
    RMSE = 3
};



enum class LinearizationType
{
    LINEAR = 1,
    PARABOLIC = 2
};

}// namespace filters