# Chaos Simulation: Gravitational Attractors

## Overview

This project simulates the chaotic behavior of particles in a gravitational field with multiple attractors. It visualizes the **basins of attraction** for particles starting from different positions in a 2D plane. The simulation was originally developed in **Python** and later extended in **C** for performance and flexibility.

The simulations output color-coded images where:

* Each **colored pixel** represents the attractor the particle eventually falls into
* **White pixels** represent particles that escape the system

---

## 🗃️ Code Files Summary (Located in `/Codes/`)

| File               | Language | Description                                                                                                                                                     |
| ------------------ | -------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `chaos.py`         | Python   | The original version. Simulates attraction zones using NumPy and Matplotlib. Slower but easy to read and modify.                                                |
| `chaos.c`          | C        | Optimized version of `chaos.py`. Reproduces the same results using low-level numerical integration and outputs `.ppm` images. Much faster.                      |
| `dynamic_random.c` | C        | An extended version that generates a **random number of attractors** with **random positions and masses** on each run. Useful for exploring new configurations. |

---

## 📂 Images

The output images for all three simulations are stored in the [/Images/](../Images/) folder. Examples include:

* `chaos_map.png` from `chaos.c`
* `random_rk4.png`, `random_symplectic.png`, `2random_rk4.png`, `2random_symplectic.png` from `dynamic_random.c`
* Matplotlib-rendered PNGs from `chaos.py`

---

## Project History

* **2024 - Python Version**: Developed as a university programming challenge. Functional but very slow on standard hardware (\~2 hours per image).
* **2025 - C Version**: Rewritten in C with help from AI. Dramatic performance gains (\~3 minutes total for both images).
* **2025 - Randomized Version**: Introduced randomized attractor generation to study statistical behavior and zone evolution over many simulations.

---

## 🔧 Features

* Two integrators: Runge-Kutta 4 (RK4) and Symplectic Euler
* Fast C-based `.ppm` image generation
* Randomized attractor mode
* Python version uses `matplotlib` for visualization
* Simple vector algebra framework for 2D dynamics

---

## ⏱️ Performance Comparison

| Version            | Runtime     | Language | Notes                      |
| ------------------ | ----------- | -------- | -------------------------- |
| `chaos.py`         | \~4 hours   | Python   | Slowest, uses Matplotlib   |
| `chaos.c`          | \~3 minutes | C        | Deterministic layout       |
| `dynamic_random.c` | \~3 minutes | C        | Random attractors each run |

---

## 🧪 How to Use

### C Versions (Recommended)

```bash
# For fixed attractors
gcc chaos.c -o chaos -lm
./chaos  # Output: rk4.ppm, symplectic.ppm

# For random attractor simulations
gcc dynamic_random.c -o dynamic_random -lm
./dynamic_random  # Output: random_rk4.ppm, random_symplectic.ppm
```

### Python Version

```bash
python chaos.py
```

> Requires: `numpy`, `matplotlib`, `tqdm`

---

## 🧭 Future Directions

* [ ] Add Verlet and other integrators
* [ ] Real-time interactive viewer (OpenGL or web-based)
* [ ] 3D gravitational simulations
* [ ] Emscripten-based browser version

---

## ⚙️ Technical Notes

The gravitational potential is modeled as:

> **V(r) = -k / r**

The force is calculated as:

> **F = -∇V = -k r̂ / r²**

Each pixel in the image represents a unique initial position, and the simulation traces the trajectory until it:

* Reaches a nearby attractor (`distance < ε`)
* Escapes the system (`distance > R_max`)

---

## 🚀 Why the C Version Is So Much Faster

* No Python object overhead
* Direct memory access (structs vs NumPy arrays)
* Inlined vector math

This project showcases how switching from high-level scripting to compiled languages can **massively improve performance**, especially in computational physics.
