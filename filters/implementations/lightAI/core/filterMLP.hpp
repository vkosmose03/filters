#pragma once
#include "lightAI.hpp"
#include "../utils/onlineNormalizer.hpp"
#include "../utils/windowBuffer.hpp"
#include "filter.hpp"
#include <string>
#include <cstdint>

namespace lightAI::core {

template <typename T, int W = 5>
class filterMLP : public filters::filterBase<T> {
 public:
    // topology      – {W, 32, 32, 1}
    // warmupSteps   – число RTK FIX эпох до включения поправки
    // statePath     – путь к файлу состояния (пуст = не сохранять/не загружать)
    explicit filterMLP(const std::vector<int>& topology,
                       int64_t  warmupSteps = 50,
                       double   lr0         = 1e-3,
                       double   lrMin       = 1e-5,
                       double   decay       = 1e-4,
                       const std::string& statePath = "");

    // ── Интерфейс filterBase<T> ────────────────────────────────────────
    void applyFilter() override;
    filters::signalContainer<T>&
    getOriginalSignalContainerReference() override { return original_; }
    filters::signalContainer<T>&
    getFilteredSignalContainerReference() override { return filtered_; }

    // ── Синхронизация с ГНСС ──────────────────────────────────────────
    // Вызывается DataAggregator при RTK FIX эпохе.
    // label – эталонная поправка из алгоритма обратной кинематики.
    void setGnssLabel(T label);

    // ── Управление состоянием ─────────────────────────────────────────
    void resetBuffer();   // сброс окна; веса и нормализатор сохраняются
    bool saveState() const;
    bool loadState();

 private:
    lightAI                      net_;
    utils::WindowBuffer<T, W>    window_;
    utils::OnlineNormalizer      normIn_;   // для входного вектора (окно)
    utils::OnlineNormalizer      normOut_;  // для целевых меток
    filters::signalContainer<T>  original_;
    filters::signalContainer<T>  filtered_;
    bool    gnssUpdateReady_ = false;
    T       gnssLabel_       = T{};
    int64_t warmupSteps_;
    int64_t trainCount_ = 0;  // число выполненных обучающих шагов
    std::string statePath_;
};

} // namespace lightAI::core
