# Heat Equation Explorer

Interactive web-based visualizations for exploring numerical methods applied to the heat equation. Designed for educational purposes to help students understand how different numerical schemes behave under various parameter choices.

## Live Demo

**App**: [https://gmcavallazzi.github.io/heat_equation/](https://gmcavallazzi.github.io/heat_equation/)

Use the 1D/2D toggle button to switch between dimensions.

## Features

### 1D Heat Equation Explorer

Solves the 1D heat equation:

```
∂T/∂t = α ∂²T/∂x²
```

**Numerical Methods**:
- Forward Euler (Explicit)
- Backward Euler (Implicit)
- Crank-Nicolson (Implicit)

**Key Features**:
- Real-time visualization of numerical vs analytical solution
- Interactive parameter adjustment (time step, grid points, max time)
- Stability analysis with Fourier number monitoring
- Error metrics (L2 and max error)
- FLOPs counting for computational cost analysis
- Educational tabs explaining methods, stability conditions, and equations

**Initial Condition**: Sine wave `u(x,0) = sin(3πx/L)` with analytical solution for comparison

### 2D Heat Equation Explorer

Solves the 2D heat equation:

```
∂T/∂t = α (∂²T/∂x² + ∂²T/∂y²)
```

**Numerical Methods**:
- Forward Euler (Explicit)
- Backward Euler (Implicit, ADI)
- Crank-Nicolson (Implicit, ADI)

**Key Features**:
- 3D surface visualization with Plotly.js
- ADI (Alternating Direction Implicit) implementation for efficient implicit solving
- Separate Fourier numbers for x and y directions
- Max relative error tracking
- Dual visualization (solution surface + error chart)

**Initial Condition**: 2D sine wave `u(x,y,0) = sin(2πx) * sin(πy)`

## Educational Purpose

This application is designed to demonstrate:

1. **Stability Conditions**: See what happens when you violate the CFL condition for explicit methods
2. **Accuracy Trade-offs**: Compare first-order (Euler) vs second-order (Crank-Nicolson) accuracy
3. **Computational Cost**: Observe FLOPs differences between explicit and implicit methods
4. **Method Behavior**: Understand how different schemes handle diffusion

**Important**: The app does not hide failures. If you choose unstable parameters (e.g., Forward Euler with r > 0.5), the solution will diverge as expected. This is intentional for learning purposes.

## Technical Details

- **Pure JavaScript** - No build process required
- **Client-side only** - All computation in the browser
- **Dependencies**: Chart.js, Plotly.js (2D only), MathJax
- **Numerical Algorithms**: Thomas Algorithm for tridiagonal systems, ADI for 2D implicit methods

## Local Development

Simply open the HTML file in a web browser:

```bash
open index.html
```

Or serve with a local server:

```bash
python3 -m http.server 8000
# Then visit http://localhost:8000
```

Use the dimension toggle button at the top to switch between 1D and 2D modes.

## Mathematical Background

### Fourier Number

The Fourier number (or diffusion number) is defined as:

```
r = α Δt / Δx²
```

For explicit methods like Forward Euler, stability requires **r ≤ 0.5** (in 1D). Implicit methods are unconditionally stable.

### Discretization Schemes

**Forward Euler**:
```
T_i^{n+1} = T_i^n + r(T_{i+1}^n - 2T_i^n + T_{i-1}^n)
```

**Backward Euler**:
```
T_i^{n+1} - r(T_{i+1}^{n+1} - 2T_i^{n+1} + T_{i-1}^{n+1}) = T_i^n
```

**Crank-Nicolson** (θ = 0.5 scheme):
```
(T_i^{n+1} - T_i^n)/Δt = (α/2)[(∂²T/∂x²)^{n+1} + (∂²T/∂x²)^n]
```

## License

Educational use for MEM407 course students.

## Author

Giorgio Cavallazzi - City, University of London
