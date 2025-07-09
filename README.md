# Chaos Simulation: Gravitational Attractors

## Overview
This project simulates the chaotic behavior of particles in a gravitational field with multiple attractors. It visualizes the basins of attraction for particles starting from different positions in a 2D plane. The simulation was originally developed in Python for a university programming course in 2024, and later optimized in C for performance.

## Project History
- **2024 Python Version**: Created as a challenge for a university programming course. While functional, it took approximately **4 hours** to generate two graphs (2 hours each) on university PC hardware.
- **2025 C Version**: Ported to C with AI assistance to learn C programming. The optimized version generates the same graphs in **under 3 minutes** - a 99% performance improvement.

## Key Features
- Two integration methods:
  - Runge-Kutta 4th order (high accuracy)
  - Symplectic Euler (energy-conserving)
- Customizable parameters:
  - Attractor positions, strengths and colors
  - Grid resolution
  - Time step and simulation duration
- PPM image output (C version)
- Matplotlib visualization (Python version)

## Performance Comparison
| Version | Execution Time | Language | Lines of Code |
|---------|----------------|----------|---------------|
| Python  | ~4 hours       | Python   | ~240          |
| C       | <3 minutes     | C        | ~175          |

## How to Use
### C Version (Recommended)
```bash
gcc -O3 -o chaos chaos.c -lm
./chaos
```
Output: `rk4.ppm` and `symplectic.ppm` images

### Python Version
```bash
python chaos.py
```
Requires: `numpy matplotlib tqdm`

## Future Development Plans
- [ ] Interactive real-time visualization
- [ ] Additional integrators (Verlet, Adams-Bashforth)
- [ ] 3D simulation support
- [ ] Web interface (Emscripten port)

## Technical Notes
The simulation models particles under gravitational potentials of the form V = -k/r. Each pixel in the output image represents:
- **Colored pixel**: Final attractor captured by the particle
- **White pixel**: Particle escaped the system

## Why the Performance Difference?
The C version achieves 80-100x speedup due to:
1. Native compilation vs Python interpretation
2. Lower-level memory management
3. Efficient inlining of vector operations
4. Removal of Python object overhead

This project demonstrates how algorithmic optimization and language choice can dramatically impact computational physics simulations.
