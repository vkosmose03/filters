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

enum class HAARthreshold
{
    SOFT = 1,
    HARD = 2
};

enum class EMFenvironment
{
    PHYSICALS = 1,
    MATHEMATICAL = 2,
    UNDEFINED = 3
};