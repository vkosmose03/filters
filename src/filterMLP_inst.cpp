#include "filters/filterMLP.hpp"

namespace lightAI::core {
    template class filterMLP<double, 5>;
    template class filterMLP<double, 8>;
    template class filterMLP<double, 16>;
    template class filterMLP<double, 32>;
}
