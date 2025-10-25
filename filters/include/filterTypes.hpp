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


/**
 * @brief This struct use to set parametrs to approximation filter.
 * 
 * @memberof useStabilization - set true if you want to stabilize
 * around zero.
 * @memberof stabilizeIncline - if incline less, then signal will 
 * stabilize if incline bigger, then signal wouldn stabilize.
 * @memberof maxIncline - maximum posible incline 
 * (upper limit to approximation incline).
 * @memberof windowSize - size of approximation step.
 * @memberof errorEstimate - which error type will be estimating.
 * @memberof type - type of line (e.g. linear, parabolic).
 */
struct approximationSettings
{
    bool useStabilization;
    double stabilizeIncline;
    double maxIncline;
    int windowSize;
    ErrorEstimate errorEstimate;
    LinearizationType type;
};
}// namespace filters