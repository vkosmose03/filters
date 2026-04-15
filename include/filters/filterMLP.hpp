/*
    * filterMAF.hpp
    *
    * This file is part of the filters library.
    * It defines the filterMAF class for applying Moving Average Filters.
    *
    * Author: Vitaly Makarov
    * Date: 2026-03-04
*/

#ifndef __FILTER_MLP_HPP
#define __FILTER_MLP_HPP

#include <core/lightAI.hpp>
#include <utils/onlineNormalizer.hpp>
#include <utils/windowBuffer.hpp>
#include "filter.hpp"
#include <string>
#include <cstdint>
#include <fstream>

namespace lightAI::core {

template <typename T, int W = 5>
class filterMLP : public filters::filterBase<T> {
 public:
    explicit filterMLP(const std::vector<int>& topology,
                       int64_t  warmupSteps = 500,
                       double   lr0         = 1e-3,
                       double   lrMin       = 1e-5,
                       double   decay       = 1e-4,
                       const std::string& statePath = "");

    void applyFilter() override;
    filters::signalContainer<T>&
    getOriginalSignalContainerReference() override { return original_; }
    filters::signalContainer<T>&
    getFilteredSignalContainerReference() override { return filtered_; }
    void setSignal(const std::vector<T> signal) override { original_.setSignal(signal); };
    std::vector<T> getSignal() const override { return filtered_.getSignal(); };

    void setGnssLabel(T label);

    void resetBuffer();
    bool saveState() const;
    bool loadState();

 private:
    lightAI                      net_;
    utils::WindowBuffer<T, W>    window_;
    utils::OnlineNormalizer      normIn_;
    utils::OnlineNormalizer      normOut_;
    filters::signalContainer<T>  original_;
    filters::signalContainer<T>  filtered_;
    bool    gnssUpdateReady_ = false;
    T       gnssLabel_       = T{};
    int64_t warmupSteps_;
    int64_t trainCount_ = 0;
    std::string statePath_;
    std::fstream corrValF;
};

#endif

template <typename T, int W>
inline filterMLP<T,W>::filterMLP(const std::vector<int>& topo,
                           int64_t warmupSteps, double lr0,
                           double lrMin, double decay,
                           const std::string& statePath)
    : net_(topo,
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
inline void filterMLP<T,W>::setGnssLabel(T label) {
    gnssLabel_       = label;
    gnssUpdateReady_ = true;
}

template <typename T, int W>
inline void filterMLP<T,W>::resetBuffer() {
    window_.reset();
    gnssUpdateReady_ = false;
}

template <typename T, int W>
inline void filterMLP<T,W>::applyFilter() {
    const std::vector<T>& sig = original_.getSignal();
    if (sig.empty()) { filtered_ = original_; return; }

    T rawVal = sig.back();

    Eigen::Matrix<T, Eigen::Dynamic, 1> meas(1);
    meas(0) = rawVal;
    window_.push(meas);

    normIn_.update(static_cast<double>(rawVal));

    if (!window_.isFull()) { filtered_ = original_; return; }

    Eigen::VectorXd x = window_.buildEigenVector();
    Eigen::VectorXd xn(x.size());

    double mean = 0.0,
           stdDev = 0.0;
    for (int i = 0; i < x.size(); ++i) {
        mean += x(i);
    }
    mean /= x.size();
    
    for (int i = 0; i < x.size(); ++i) {
        stdDev += (x(i) - mean) * (x(i) - mean);
    }

    stdDev = std::sqrt(stdDev / x.size());
    
    for (int i = 0; i < x.size(); ++i) {
        xn(i) = (x(i) - mean) / stdDev;
    }   

    Eigen::VectorXd corrNorm = net_.tick(xn);
    double corrVal = corrNorm(0) * stdDev + mean;

    // В фазе прогрева сигнал не модифицируем
    std::vector<T> outSig = sig;
    if (trainCount_ >= warmupSteps_ && !std::isnan(corrVal))
        outSig.back() = rawVal - static_cast<T>(corrVal);
    filtered_.setSignal(outSig);

    // Обучающий шаг
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
            saveState();
    }
}

template <typename T, int W>
inline bool filterMLP<T,W>::saveState() const {
    if (statePath_.empty()) return false;
    if (!net_.saveWeights(statePath_ + ".weights")) return false;
    std::ofstream f(statePath_ + ".norm", std::ios::binary);
    if (!f) return false;
    f.write(reinterpret_cast<const char*>(&trainCount_), sizeof(trainCount_));
    for (const auto* nm : {&normIn_, &normOut_}) {
        f.write(reinterpret_cast<const char*>(&nm->n_),    sizeof(nm->n_));
        f.write(reinterpret_cast<const char*>(&nm->mean_), sizeof(nm->mean_));
        f.write(reinterpret_cast<const char*>(&nm->M_),    sizeof(nm->M_));
    }
    return f.good();
}

template <typename T, int W>
inline bool filterMLP<T,W>::loadState() {
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
template class filterMLP<double, 8>;
} // namespace lightAI::core
