#pragma once
#include <cmath>
#include <cstdint>

namespace lightAI::utils {

// Инкрементальное вычисление среднего и стандартного отклонения
// алгоритмом Уэлфорда. O(1) памяти и O(1) вычислений на отсчёт.
class OnlineNormalizer {
 public:
    OnlineNormalizer() : n_(0), mean_(0.0), M_(0.0) {}

    // Обновить оценки по новому значению x.
    void update(double x) {
        ++n_;
        double delta  = x - mean_;
        mean_ += delta / static_cast<double>(n_);
        double delta2 = x - mean_;
        M_    += delta * delta2;
    }

    double mean()   const { return mean_; }
    double stddev() const {
        if (n_ < 2) return 1.0;          // защита: единичное σ до накопления
        double s = std::sqrt(M_ / static_cast<double>(n_));
        return (s < 1e-9) ? 1.0 : s;
    }
    int64_t count() const { return n_; }

    // Нормализовать значение: (x - μ) / σ
    double normalize(double x)   const { return (x - mean()) / stddev(); }

    // Денормализовать значение: y * σ + μ
    double denormalize(double y) const { return y * stddev() + mean(); }

    // Сериализация состояния (для сохранения между сессиями)
    int64_t n_;
    double  mean_;
    double  M_;
};

} // namespace lightAI::utils
