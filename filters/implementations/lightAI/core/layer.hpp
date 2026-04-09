#pragma once
#include <Eigen/Dense>
#include "../utils/activations.hpp"

namespace lightAI::core {

class layer {
 public:
    layer(int layerSize, int prevLayerSize,
          utils::ActivFn activation  = utils::relu(),
          utils::ActivFn activationD = utils::reluD());

    Eigen::VectorXd forward (const Eigen::VectorXd& input);
    Eigen::VectorXd backward(const Eigen::VectorXd& gradOut,
                             double learningRate);

    const Eigen::MatrixXd& getWeights() const { return W_; }
    const Eigen::VectorXd& getBiases()  const { return b_; }
    void setWeights(const Eigen::MatrixXd& W) { W_ = W; }
    void setBiases (const Eigen::VectorXd& b) { b_ = b; }

 private:
    Eigen::MatrixXd W_;
    Eigen::VectorXd b_, z_, output_, input_;
    utils::ActivFn  act_, actD_;
};

} // namespace lightAI::core
