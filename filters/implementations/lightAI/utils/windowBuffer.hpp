#pragma once
#include <array>
#include <Eigen/Dense>

namespace lightAI::utils {

template <typename T, int W>
class WindowBuffer {
 public:
    explicit WindowBuffer(int dims) : dims_(dims), head_(0), count_(0) {
        for (auto& v : buf_)
            v = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(dims);
    }

    void push(const Eigen::Matrix<T, Eigen::Dynamic, 1>& v) {
        buf_[head_] = v;
        head_ = (head_ + 1) % W;
        if (count_ < W) ++count_;
    }

    bool    isFull() const { return count_ == W; }
    void    reset()        { head_ = 0; count_ = 0; }
    int     dims()   const { return dims_; }

    Eigen::VectorXd buildEigenVector() const {
        Eigen::VectorXd x(W * dims_);
        for (int i = 0; i < W; ++i) {
            int idx = (head_ + i) % W;
            x.segment(i * dims_, dims_) =
                buf_[idx].template cast<double>();
        }
        return x;
    }

 private:
    std::array<Eigen::Matrix<T, Eigen::Dynamic, 1>, W> buf_;
    int dims_, head_, count_;
};

} // namespace lightAI::utils
