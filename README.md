# Intelligent and Adaptive Control Systems Project: Non-Linear Systems Stabilization & Backstepping

[![MATLAB](https://img.shields.io/badge/MATLAB-Simulation-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![Control](https://img.shields.io/badge/Control-Adaptive_%26_Non--Linear-red.svg)]()
[![Stability](https://img.shields.io/badge/Stability-Lyapunov-green.svg)]()

## 📌 Project Overview
This repository focuses on the design and simulation of advanced control laws for non-linear dynamical systems. The primary challenge addressed is the stabilization of a second-order system characterized by non-linearities and **parametric uncertainty** (unknown parameter $\theta^*$).

The project compares traditional linearization techniques with state-of-the-art **Adaptive Control** methods, proving stability through Lyapunov analysis and validating performance via MATLAB simulations.

---

## 🛠️ Control Methodologies Implemented

The project is structured into four major control design phases:

### 1. Feedback Linearization (Known Parameters)
- Design of a control law $u$ that cancels non-linearities to achieve a linear closed-loop response.
- Implementation of pole-placement to ensure asymptotic stability $(x_1, x_2) \rightarrow (0,0)$.

### 2. Backstepping Control (Known Parameters)
- A recursive design methodology for stabilizing non-linear systems.
- Step-by-step construction of a **Control Lyapunov Function (CLF)** to guarantee system convergence without aggressive gain requirements.

### 3. Adaptive Feedback Linearization (Unknown $\theta^*$)
- Design of an **Update Law** (Adaptive Mechanism) to estimate the unknown parameter $\theta^*$ in real-time.
- Proof of boundedness for all closed-loop signals and convergence of the states to the origin.

### 4. Adaptive Backstepping (Unknown $\theta^*$)
- The most robust controller in the project.
- Combines recursive backstepping with an online parameter estimator.
- **Key Advantage:** Demonstrated superior performance in disturbance rejection and faster settling times compared to standard linearization.

---

## 📊 Performance & Robustness Analysis
A critical part of this work is the comparison between **Feedback Linearization** and **Backstepping** under external disturbances:
- **Settling Time:** Backstepping achieved faster recovery with lower control effort (Gains $k=2$ vs $k=70$ in linearization).
- **Disturbance Rejection:** The Adaptive Backstepping controller proved significantly more robust against abrupt changes in system dynamics.
- **Stability:** All designs were mathematically verified to ensure that signals remain bounded (BIBO stability) under all operating conditions.
