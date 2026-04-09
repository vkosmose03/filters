#include "filterMLP.hpp"
#include <fstream>

namespace lightAI::core {

template <typename T, int W>
filterMLP<T,W>::filterMLP(const std::vector<int>& topo,
                           int64_t warmupSteps, double lr0,
                           double lrMin, double decay,
                           const std::string& statePath)
    : net_(topo,
           // скрытые слои – ReLU, выходной – линейный
           [&]() {
               std::vector<utils::ActivFn> a, aD;
               for (size_t i = 0; i < topo.size()-2; ++i)
                   { a.push_back(utils::relu()); aD.push_back(utils::reluD()); }
               a.push_back(utils::linear()); aD.push_back(utils::linearD());
               return a;
           }(),
           [&]() {
               std::vector<utils::ActivFn> a, aD;
               for (size_t i = 0; i < topo.size()-2; ++i)
                   { a.push_back(utils::relu()); aD.push_back(utils::reluD()); }
               a.push_back(utils::linear()); aD.push_back(utils::linearD());
               return aD;
           }(),
           lr0, lrMin, decay)
    , window_(1)
    , warmupSteps_(warmupSteps)
    , statePath_(statePath)
{
    if (!statePath_.empty()) loadState();
}

template <typename T, int W>
void filterMLP<T,W>::setGnssLabel(T label) {
    gnssLabel_       = label;
    gnssUpdateReady_ = true;
}

template <typename T, int W>
void filterMLP<T,W>::resetBuffer() {
    window_.reset();
    gnssUpdateReady_ = false;
    // normIn_, normOut_ и веса сети не сбрасываются
}

template <typename T, int W>
void filterMLP<T,W>::applyFilter() {
    const std::vector<T>& sig = original_.getSignal();
    if (sig.empty()) { filtered_ = original_; return; }

    T rawVal = sig.back();

    // Обновляем скользящее окно
    Eigen::Matrix<T, Eigen::Dynamic, 1> meas(1);
    meas(0) = rawVal;
    window_.push(meas);

    // Обновляем нормализатор входных данных
    normIn_.update(static_cast<double>(rawVal));

    if (!window_.isFull()) { filtered_ = original_; return; }

    // Строим нормализованный входной вектор
    Eigen::VectorXd x = window_.buildEigenVector();
    Eigen::VectorXd xn(x.size());
    for (int i = 0; i < x.size(); ++i)
        xn(i) = normIn_.normalize(x(i));

    // Инференс
    Eigen::VectorXd corrNorm = net_.tick(xn);
    double corrVal = normOut_.denormalize(corrNorm(0));

    // В фазе прогрева сигнал не модифицируем
    std::vector<T> outSig = sig;
    if (trainCount_ >= warmupSteps_)
        outSig.back() = rawVal - static_cast<T>(corrVal);
    filtered_.setSignal(outSig);

    // Обучающий шаг – только при наличии новой RTK FIX метки
    if (gnssUpdateReady_) {
        double labelD = static_cast<double>(gnssLabel_);
        normOut_.update(labelD);
        double labelNorm = normOut_.normalize(labelD);
        Eigen::VectorXd target(1);
        target(0) = labelNorm;
        net_.learnStep(xn, target);
        ++trainCount_;
        gnssUpdateReady_ = false;
        if (!statePath_.empty() && (trainCount_ % 100 == 0))
            saveState();  // периодическое сохранение
    }
}

template <typename T, int W>
bool filterMLP<T,W>::saveState() const {
    if (statePath_.empty()) return false;
    if (!net_.saveWeights(statePath_ + ".weights")) return false;
    std::ofstream f(statePath_ + ".norm", std::ios::binary);
    if (!f) return false;
    // Сохраняем состояния обоих нормализаторов и счётчик
    f.write(reinterpret_cast<const char*>(&trainCount_), sizeof(trainCount_));
    for (const auto* nm : {&normIn_, &normOut_}) {
        f.write(reinterpret_cast<const char*>(&nm->n_),    sizeof(nm->n_));
        f.write(reinterpret_cast<const char*>(&nm->mean_), sizeof(nm->mean_));
        f.write(reinterpret_cast<const char*>(&nm->M_),    sizeof(nm->M_));
    }
    return f.good();
}

template <typename T, int W>
bool filterMLP<T,W>::loadState() {
    if (statePath_.empty()) return false;
    if (!net_.loadWeights(statePath_ + ".weights")) return false;
    std::ifstream f(statePath_ + ".norm", std::ios::binary);
    if (!f) return false;
    f.read(reinterpret_cast<char*>(&trainCount_), sizeof(trainCount_));
    for (auto* nm : {&normIn_, &normOut_}) {
        f.read(reinterpret_cast<char*>(&nm->n_),    sizeof(nm->n_));
        f.read(reinterpret_cast<char*>(&nm->mean_), sizeof(nm->mean_));
        f.read(reinterpret_cast<char*>(&nm->M_),    sizeof(nm->M_));
    }
    return f.good();
}

template class filterMLP<double, 5>;

} // namespace lightAI::core
