#pragma once
#include <vector>

template <typename T>
class signalContainer
{
private:
    std::vector<T> signal_;
    std::vector<T> filteredSignal_;
    double frequency_;
    T average_;
    T variance_;
    T stdDeviation_;
    T max_;
    T min_ ;
    T rawSum_;
    T filteredSum_;
public:
    signalContainer();
    signalContainer(const std::vector<T>& signal);
    ~signalContainer();
    void updateSignal(const std::vector<T>& newSignal);
    T getAverage() const;
    T getVariance() const;
    T getStdDeviation() const;
    T getMax() const;
    T getMin() const;
    T getRawSum() const;
    T getFilteredSum() const;
    std::vector<T> getFilteredSignal() const;
};