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
    double offset;
    ErrorEstimate errorEstimate;
    LinearizationType type;
};



/**
 * @brief This struct use to set parametrs to EMF filter.
 * 
 * @memberof signalType - type of signal (e.g. physical, radiotechnical).
 * @memberof physicalK - smoothing factor for physical signals.
 * @memberof standardK - standard smoothing factor for radiotechnical and undefined signals.
 * @memberof maximalK - maximum smoothing factor for radiotechnical and undefined signals.
 * @memberof threshold - threshold value for radiotechnical signals.
 */
struct EMFfilterSettings
{
    EMFenvironment signalType;
    double physicalK;
    double standardK;
    double maximalK;
    double threshold;
};



/**
 * @brief This struct use to set parametrs to HAAR filter.
 * 
 * @memberof thresholdType - type of thresholding (e.g. soft, hard).
 * @memberof thresholdValue - threshold value for filtering.
 * @memberof filteringWindow - size of filtering window.
 * @memberof depth - depth of the Haar transform.
 */
struct HAARfilterSettings
{
    HAARthreshold thresholdType;
    double thresholdValue;
    int filteringWindow;
    int depth;
};
}// namespace filters