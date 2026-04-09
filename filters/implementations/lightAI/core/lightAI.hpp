#pragma once
#include "layer.hpp"
#include <vector>
#include <string>
#include <cstdint>

namespace lightAI::core {

class lightAI {
 public:
    explicit lightAI(const std::vector<int>&            topology,
                     const std::vector<utils::ActivFn>& activations,
                     const std::vector<utils::ActivFn>& activationsD,
                     double lr0    = 1e-3,
                     double lrMin  = 1e-5,
                     double decay  = 1e-4);

    Eigen::VectorXd tick(const Eigen::VectorXd& input);

    // Один шаг online SGD с decay lr; возвращает MSE.
    double learnStep(const Eigen::VectorXd& input,
                     const Eigen::VectorXd& target);

    bool saveWeights(const std::string& path) const;
    bool loadWeights(const std::string& path);

    double currentLr()   const { return lr_; }
    int64_t stepCount()  const { return step_; }

 private:
    std::vector<layer> layers_;
    double lr_, lr0_, lrMin_, decay_;
    int64_t step_ = 0;
    void backPropagation(const Eigen::VectorXd& error);
    void updateLr();
};

} // namespace lightAI::core
