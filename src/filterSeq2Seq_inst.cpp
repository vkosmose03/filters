#include "filters/filterSeq2Seq.hpp"

namespace lightAI::core {
    template class filterSeq2Seq<double, 5>;
    template class filterSeq2Seq<double, 8>;
    template class filterSeq2Seq<double, 16>;
    template class filterSeq2Seq<double, 32>;
}
