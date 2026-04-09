#include "lightAI.hpp"
#include <fstream>
#include <cmath>
#include <stdexcept>

namespace lightAI::core {

lightAI::lightAI(const std::vector<int>&            topo,
                 const std::vector<utils::ActivFn>& acts,
                 const std::vector<utils::ActivFn>& actsD,
                 double lr0, double lrMin, double decay)
    : lr_(lr0), lr0_(lr0), lrMin_(lrMin), decay_(decay)
{
    if (topo.size() < 2)
        throw std::invalid_argument("Topology must have >= 2 elements");
    layers_.reserve(topo.size() - 1);
    for (size_t i = 0; i < topo.size() - 1; ++i)
        layers_.emplace_back(topo[i+1], topo[i], acts[i], actsD[i]);
}

void lightAI::updateLr() {
    lr_ = std::max(lrMin_, lr0_ * std::exp(-decay_ * static_cast<double>(step_)));
}

Eigen::VectorXd lightAI::tick(const Eigen::VectorXd& input) {
    Eigen::VectorXd cur = input;
    for (auto& l : layers_) cur = l.forward(cur);
    return cur;
}

double lightAI::learnStep(const Eigen::VectorXd& input,
                           const Eigen::VectorXd& target) {
    Eigen::VectorXd out   = tick(input);
    Eigen::VectorXd error = out - target;
    double mse = error.squaredNorm() / static_cast<double>(error.size());
    backPropagation(error);
    ++step_;
    updateLr();
    return mse;
}

void lightAI::backPropagation(const Eigen::VectorXd& error) {
    Eigen::VectorXd grad = error;
    for (int i = static_cast<int>(layers_.size())-1; i >= 0; --i)
        grad = layers_[i].backward(grad, lr_);
}

bool lightAI::saveWeights(const std::string& path) const {
    std::ofstream f(path, std::ios::binary);
    if (!f) return false;
    // Сохраняем step_ и lr_ для корректного восстановления decay
    f.write(reinterpret_cast<const char*>(&step_), sizeof(step_));
    f.write(reinterpret_cast<const char*>(&lr_),   sizeof(float));
    int32_t n = static_cast<int32_t>(layers_.size());
    f.write(reinterpret_cast<const char*>(&n), sizeof(n));
    for (const auto& l : layers_) {
        const auto& W = l.getWeights();
        int32_t r = W.rows(), c = W.cols();
        f.write(reinterpret_cast<const char*>(&r), sizeof(r));
        f.write(reinterpret_cast<const char*>(&c), sizeof(c));
        Eigen::MatrixXf Wf = W.cast<float>();
        f.write(reinterpret_cast<const char*>(Wf.data()),
                Wf.size()*sizeof(float));
        const auto& b = l.getBiases();
        int32_t bs = b.size();
        f.write(reinterpret_cast<const char*>(&bs), sizeof(bs));
        Eigen::VectorXf bf = b.cast<float>();
        f.write(reinterpret_cast<const char*>(bf.data()),
                bf.size()*sizeof(float));
    }
    return f.good();
}

bool lightAI::loadWeights(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return false;
    f.read(reinterpret_cast<char*>(&step_), sizeof(step_));
    float lrf; f.read(reinterpret_cast<char*>(&lrf), sizeof(float));
    lr_ = static_cast<double>(lrf);
    int32_t n; f.read(reinterpret_cast<char*>(&n), sizeof(n));
    if (n != static_cast<int32_t>(layers_.size())) return false;
    for (auto& l : layers_) {
        int32_t r, c;
        f.read(reinterpret_cast<char*>(&r), sizeof(r));
        f.read(reinterpret_cast<char*>(&c), sizeof(c));
        Eigen::MatrixXf Wf(r, c);
        f.read(reinterpret_cast<char*>(Wf.data()), Wf.size()*sizeof(float));
        l.setWeights(Wf.cast<double>());
        int32_t bs; f.read(reinterpret_cast<char*>(&bs), sizeof(bs));
        Eigen::VectorXf bf(bs);
        f.read(reinterpret_cast<char*>(bf.data()), bf.size()*sizeof(float));
        l.setBiases(bf.cast<double>());
    }
    return f.good();
}

} // namespace lightAI::core
