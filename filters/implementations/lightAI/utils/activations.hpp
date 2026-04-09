#pragma once
#include <cmath>
#include <functional>

namespace lightAI::utils {

using ActivFn = std::function<double(double)>;

inline ActivFn relu()  { return [](double v){ return v > 0.0 ? v : 0.0; }; }
inline ActivFn reluD() { return [](double v){ return v > 0.0 ? 1.0 : 0.0; }; }

inline ActivFn linear()  { return [](double v){ return v; }; }
inline ActivFn linearD() { return [](double)  { return 1.0; }; }

} // namespace lightAI::utils
