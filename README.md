# Filters Library

A C++17 library for signal filtering, combining classical algorithms and AI-based methods.
It provides implementations of various filters and allows combining them into flexible processing chains.

## About

This repository contains a library designed for filtering time series and signals.
It offers a universal interface based on the template base class `filterBase<T>` and includes
a wide range of methods: from simple moving averages to neural network-based filters (MLP, Seq2Seq).

Key features:
- **Template architecture:** Works with different data types (`float`, `double`).
- **Filter chains:** The `filterChain<T>` class allows applying multiple filters sequentially.
- **AI filters:** Includes implementations of a single-layer perceptron (MLP) and Sequence-to-Sequence models.
- **Cross-platform:** Build via CMake.

## Available Filters

| Filter | Header file | Description |
| :--- | :--- | :--- |
| **Moving Average Filter (MAF)** | `filterMAF.hpp` | Classical moving average filter for signal smoothing. |
| **Exponential Moving Filter (EMF)** | `filterEMF.hpp` | Exponential moving filter that gives more weight to recent values. |
| **Median Filter** | `filterMedian.hpp` | Median filter, effective for suppressing impulse noise. |
| **Haar Wavelet Filter** | `filterHAAR.hpp` | Filter based on Haar wavelet transform with soft and hard thresholding. |
| **Approximation Filter** | `approximation.hpp` | Filter for piecewise linear approximation with optional slope stabilization. |
| **Kalman Filter (1D)** | `filterKalman.hpp` | One-dimensional Kalman filter, ideal for reducing noise and drift in IMU signals. |
| **AI MLP Filter** | `filterMLP.hpp` | Filter based on a multi-layer perceptron (MLP) to transform a sequence into a value. |
| **AI Seq2Seq Filter** | `filterSeq2Seq.hpp` | Filter based on a Sequence-to-Sequence model for processing sequences. |
| **Butterworth Filter** | `filterButterworth.hpp` | Base IIR low-pass filter without signal form changing.

## Requirements

- Compiler with C++17 support (GCC, Clang, MSVC)
- [CMake](https://cmake.org/) (version 3.14 or higher)
- [lightAI](https://github.com/vkosmose03/lightAI) library (fetched automatically during build)

## Build and Integration

### Building the project

```bash
# Clone the repository
git clone https://github.com/vkosmose03/filters.git
cd filters

# Create a build directory
mkdir build && cd build

# Configure and build
cmake ..
make
