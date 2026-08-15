[README.md](https://github.com/user-attachments/files/31107448/README.md)
# Newton and Gauss-Newton Methods for GPS Positioning

Numerical implementation and comparison of two methods for solving the GPS positioning problem: the Multivariate Newton Method on an exactly-determined 4-satellite system, and the Gauss-Newton Method on an overdetermined 10-satellite system.

## Overview

GPS receivers determine their location by solving a system of nonlinear equations derived from pseudorange measurements (signal travel time to each satellite) and known satellite positions. In theory, 4 satellites are enough to solve for the receiver's 3D position and clock bias exactly. In practice, real-world pseudorange measurements are corrupted by atmospheric delays and other noise, so an exact solution using only 4 satellites locks in that noise rather than correcting for it.

This project implements both approaches on a real satellite dataset to see how each handles that tradeoff:

- **Multivariate Newton's Method** — solves the exactly-determined 4-satellite square system
- **Gauss-Newton Method** — solves the overdetermined 10-satellite system via nonlinear least squares

## Results

| Method | Satellites | Convergence | Iterations | Result |
|---|---|---|---|---|
| Newton | 4 | Quadratic | 7 | Exact solution to the system, but geographically inaccurate |
| Gauss-Newton | 10 | Linear - Superlinear | 11 | Least-squares solution, more geographically accurate |

The key finding: Newton's method achieves quadratic convergence to an exact solution of the 4-satellite square system, but since the system is exactly determined, the solution has no residual degrees of freedom to absorb measurement noise. Gauss-Newton exhibits slower, linear convergence, but by solving the overdetermined 10-satellite system in a least-squares sense, it distributes the residual error across redundant measurements, attenuating the effect of atmospheric noise on the final position estimate. The results indicate that measurement redundancy combined with least-squares estimation yields a more geographically accurate solution than an exact fit to an underdetermined, noise-corrupted system.

## Repository Structure

```
.
├── 1_Data/           # Satellite pseudorange and position data
├── 2_Code/           # Python implementations of Newton and Gauss-Newton solvers
├── 3_Report/         # Full technical writeup (PDF)
└── README.md
```

## Methodology

Both solvers linearize the nonlinear pseudorange equations at each iteration and update the position estimate using the Jacobian of the residuals:

- **Newton's Method** solves the square (4-equation, 4-unknown) system directly, using the full Jacobian inverse at each step.
- **Gauss-Newton** solves the overdetermined (10-equation, 4-unknown) system by minimizing the sum of squared residuals, using the pseudoinverse of the Jacobian at each step.

See the full writeup in `/3_Report` the linearization of the pseudorange equations and the least-squares formulation.

## Data

Real satellite pseudorange data, used to construct both the 4-satellite and 10-satellite systems.

## Author

Alexandru Tudor
