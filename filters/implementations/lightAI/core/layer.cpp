#include "layer.hpp"

namespace lightAI::core {

layer::layer(int n, int m,
             utils::ActivFn activation, utils::ActivFn activationD)
    : act_(std::move(activation)), actD_(std::move(activationD))
{
    // Нулевая инициализация: прозрачный старт (нулевая поправка на выходе)
    W_ = Eigen::MatrixXd::Random(n, m);
    b_ = Eigen::VectorXd::Zero(n);
    z_.resize(n); output_.resize(n); input_.resize(m);
}

Eigen::VectorXd layer::forward(const Eigen::VectorXd& input) {
    input_  = input;
    z_      = W_ * input + b_;
    output_ = z_.unaryExpr(act_);
    return output_;
}

Eigen::VectorXd layer::backward(const Eigen::VectorXd& gradOut,
                                double lr) {
    Eigen::VectorXd delta = gradOut.cwiseProduct(z_.unaryExpr(actD_));
    Eigen::VectorXd gIn   = W_.transpose() * delta;
    W_ -= lr * (delta * input_.transpose());
    b_ -= lr * delta;
    return gIn;
}

} // namespace lightAI::core
