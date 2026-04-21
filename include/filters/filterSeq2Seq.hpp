/*
    * filterSeq2Seq.hpp
    *
    * This file is part of the filters library.
    * It defines the filterSeq2Seq class for applying sequence to sequence AI filtering.
    *
    * Author: Vitaly Makarov
    * Date: 2026-04-20
*/

#ifndef __FILTER_SEQ2SEQ_HPP
#define __FILTER_SEQ2SEQ_HPP

#include <core/lightAI.hpp>
#include <utils/onlineNormalizer.hpp>
#include <utils/windowBuffer.hpp>
#include "filter.hpp"
#include <deque>
#include <vector>
#include <string>
#include <cstdint>
#include <fstream>
#include <Eigen/Dense>

namespace lightAI::core {

template <typename T, int W = 10>
class filterSeq2Seq : public filters::filterBase<T> {
public:
    explicit filterSeq2Seq(const std::vector<int>& topology,
                           int64_t warmupSteps = 500,
                           double lr0 = 1e-3,
                           double lrMin = 1e-5,
                           double decay = 1e-4,
                           const std::string& statePath = "");

    void applyFilter() override;
    filters::signalContainer<T>& getOriginalSignalContainerReference() override { return original_; }
    filters::signalContainer<T>& getFilteredSignalContainerReference() override { return filtered_; }
    void setSignal(const std::vector<T> signal) override { original_.setSignal(signal); }
    std::vector<T> getSignal() const override { return filtered_.getSignal(); }

    void setGnssLabels(const std::vector<T>& labels);
    bool isWarmup() const { return trainCount_ >= (warmupSteps_ + 5); }
    void resetBuffer();

    bool saveState() const;
    bool loadState();

private:
    lightAI net_;
    std::deque<Eigen::Matrix<T, Eigen::Dynamic, 1>> imuWindow_;
    std::deque<T> gnssLabelWindow_;
    utils::OnlineNormalizer normIn_;
    utils::OnlineNormalizer normOut_;
    filters::signalContainer<T> original_;
    filters::signalContainer<T> filtered_;
    bool gnssUpdateReady_ = false;
    int64_t warmupSteps_;
    int64_t trainCount_ = 0;
    std::string statePath_;

    #ifdef _FILTER_MLP_DBG 
        std::fstream corrValF;
    #endif
};

extern template class filterSeq2Seq<double, 5>;
extern template class filterSeq2Seq<double, 8>;
extern template class filterSeq2Seq<double, 16>;
extern template class filterSeq2Seq<double, 32>;


template <typename T, int W>
filterSeq2Seq<T, W>::filterSeq2Seq(const std::vector<int>& topo,
                                   int64_t warmupSteps, double lr0,
                                   double lrMin, double decay,
                                   const std::string& statePath)
    : net_(topo,
           [&]() {
               std::vector<utils::ActivFn> a, aD;
               for (size_t i = 0; i < topo.size() - 2; ++i) {
                   a.push_back(utils::relu());
                   aD.push_back(utils::reluD());
               }
               a.push_back(utils::linear());
               aD.push_back(utils::linearD());
               return a;
           }(),
           [&]() {
               std::vector<utils::ActivFn> a, aD;
               for (size_t i = 0; i < topo.size() - 2; ++i) {
                   a.push_back(utils::relu());
                   aD.push_back(utils::reluD());
               }
               a.push_back(utils::linear());
               aD.push_back(utils::linearD());
               return aD;
           }(),
           lr0, lrMin, decay)
    , warmupSteps_(warmupSteps)
    , statePath_(statePath)
{
    if (!statePath_.empty()) loadState();

    #ifdef _FILTER_SEQ2SEQ_DBG
        corrValF.open("SEQ2SEQ_out.log", std::ios_base::out);
    #endif

}

template <typename T, int W>
void filterSeq2Seq<T, W>::setGnssLabels(const std::vector<T>& labels) {
    gnssLabelWindow_.clear();
    for (const auto& v : labels) gnssLabelWindow_.push_back(v);
    gnssUpdateReady_ = !gnssLabelWindow_.empty();
}

template <typename T, int W>
void filterSeq2Seq<T, W>::resetBuffer() {
    imuWindow_.clear();
    gnssLabelWindow_.clear();
    gnssUpdateReady_ = false;
}

template <typename T, int W>
void filterSeq2Seq<T, W>::applyFilter() {
    const std::vector<T>& sig = original_.getSignal();
    if (sig.empty()) { filtered_ = original_; return; }

    T rawVal = sig.back();
    Eigen::Matrix<T, Eigen::Dynamic, 1> meas(1);
    meas(0) = rawVal;
    imuWindow_.push_back(meas);
    if (static_cast<int>(imuWindow_.size()) > W) imuWindow_.pop_front();

    normIn_.update(static_cast<double>(rawVal));

    if (static_cast<int>(imuWindow_.size()) < W) { filtered_ = original_; return; }

    Eigen::VectorXd x(W);
    for (int i = 0; i < W; ++i) x(i) = imuWindow_[i](0);

    double mean = normIn_.mean();
    double stdDev = normIn_.stddev();
    Eigen::VectorXd xn(W);
    for (int i = 0; i < W; ++i) xn(i) = (x(i) - mean) / stdDev;

    Eigen::VectorXd corrNorm = net_.tick(xn);
    Eigen::VectorXd corrVal = (corrNorm.array() * stdDev + mean).matrix();

    std::vector<T> outSig = sig;
    if (trainCount_ >= warmupSteps_ && !std::isnan(corrVal(W - 1)))
        outSig.back() = rawVal - static_cast<T>(corrVal(W - 1));
    filtered_.setSignal(outSig);

    #ifdef _FILTER_SEQ2SEQ_DBG
        corrValF << "Raw: " << xn << std::endl << ",    label: ";
        for (T v : gnssLabelWindow_) {
            corrValF << v << std::endl;
        }    
        corrValF << ",    out: " << corrVal << std::endl;
    #endif

    if (gnssUpdateReady_ && static_cast<int>(gnssLabelWindow_.size()) >= W) {
        std::vector<double> labels(W);
        for (int i = 0; i < W; ++i) labels[i] = static_cast<double>(gnssLabelWindow_[i]);
        Eigen::VectorXd y(W);
        for (int i = 0; i < W; ++i) y(i) = labels[i];
        Eigen::VectorXd yNorm(W);
        for (int i = 0; i < W; ++i) {
            normOut_.update(y(i));
            yNorm(i) = normOut_.normalize(y(i));
        }
        net_.learnStep(xn, yNorm);
        ++trainCount_;
        gnssUpdateReady_ = false;
        if (!statePath_.empty() && (trainCount_ % 100 == 0)) saveState();
    }
}

template <typename T, int W>
bool filterSeq2Seq<T, W>::saveState() const {
    if (statePath_.empty()) return false;
    if (!net_.saveWeights(statePath_ + ".weights")) return false;
    std::ofstream f(statePath_ + ".norm", std::ios::binary);
    if (!f) return false;
    f.write(reinterpret_cast<const char*>(&trainCount_), sizeof(trainCount_));
    for (const auto* nm : {&normIn_, &normOut_}) {
        f.write(reinterpret_cast<const char*>(&nm->n_), sizeof(nm->n_));
        f.write(reinterpret_cast<const char*>(&nm->mean_), sizeof(nm->mean_));
        f.write(reinterpret_cast<const char*>(&nm->M_), sizeof(nm->M_));
    }
    return f.good();
}

template <typename T, int W>
bool filterSeq2Seq<T, W>::loadState() {
    if (statePath_.empty()) return false;
    if (!net_.loadWeights(statePath_ + ".weights")) return false;
    std::ifstream f(statePath_ + ".norm", std::ios::binary);
    if (!f) return false;
    f.read(reinterpret_cast<char*>(&trainCount_), sizeof(trainCount_));
    for (auto* nm : {&normIn_, &normOut_}) {
        f.read(reinterpret_cast<char*>(&nm->n_), sizeof(nm->n_));
        f.read(reinterpret_cast<char*>(&nm->mean_), sizeof(nm->mean_));
        f.read(reinterpret_cast<char*>(&nm->M_), sizeof(nm->M_));
    }
    return f.good();
}

} // namespace lightAI::core

#endif
