/**
 * 2D Heat Equation Explorer
 * Interactive visualization of numerical methods for the 2D heat equation
 * Uses ADI (Alternating Direction Implicit) for implicit methods
 */

// ========================================
// Configuration & State
// ========================================

const config = {
    L: 1.0,           // Domain length (Square domain LxL)
    alpha: 0.01,      // Thermal diffusivity (Fixed)
    Nx: 41,           // Number of grid points (Nx = Ny)
    dt: 0.001,        // Time step
    Tmax: 2.0,        // Maximum simulation time
    method: 'forward-euler',
    animationSpeed: 1, // Slower for 2D
};

const state = {
    t: 0,
    nSteps: 0,
    flops: 0,
    u: null,          // Current numerical solution (2D array)
    uExact: null,     // Exact analytical solution (2D array)
    x: null,          // Grid points
    y: null,
    dx: null,
    dy: null,
    rx: null,         // Fourier number in x
    ry: null,         // Fourier number in y
    isPlaying: false,
    animationId: null,
    errorChart: null, // Chart.js instance for error
    hasDiverged: false,
};

// ========================================
// Initial Condition & Analytical Solution
// ========================================

/**
 * Initial Condition: u(x,y,0) = sin(2πx) * sin(πy)
 * One positive peak, one negative peak
 */
function getSineIC2D(x, y) {
    const u = new Array(x.length);
    for (let i = 0; i < x.length; i++) {
        u[i] = new Array(y.length);
        for (let j = 0; j < y.length; j++) {
            u[i][j] = Math.sin(2 * Math.PI * x[i]) * Math.sin(Math.PI * y[j]);
        }
    }
    return u;
}

/**
 * Exact Analytical Solution
 * u(x,y,t) = sin(2πx)sin(πy) * exp(-α * π² * (2² + 1²) * t)
 */
function getAnalyticalSolution2D(x, y, t) {
    const alpha = config.alpha;
    const decay = Math.exp(-alpha * Math.PI * Math.PI * 5 * t);

    const u = new Array(x.length);
    for (let i = 0; i < x.length; i++) {
        u[i] = new Array(y.length);
        for (let j = 0; j < y.length; j++) {
            u[i][j] = Math.sin(2 * Math.PI * x[i]) * Math.sin(Math.PI * y[j]) * decay;
        }
    }
    return u;
}

// ========================================
// Numerical Solvers
// ========================================

/**
 * Thomas Algorithm for tridiagonal systems
 * Solves: a_i * x_{i-1} + b_i * x_i + c_i * x_{i+1} = d_i
 */
function thomasAlgorithm(a, b, c, d) {
    const n = d.length;
    const cp = new Array(n);
    const dp = new Array(n);
    const x = new Array(n);

    if (Math.abs(b[0]) < 1e-15) return d; // Safety

    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];

    for (let i = 1; i < n; i++) {
        const denom = b[i] - a[i] * cp[i - 1];
        if (Math.abs(denom) < 1e-15) return d;
        cp[i] = c[i] / denom;
        dp[i] = (d[i] - a[i] * dp[i - 1]) / denom;
    }

    x[n - 1] = dp[n - 1];
    for (let i = n - 2; i >= 0; i--) {
        x[i] = dp[i] - cp[i] * x[i + 1];
    }

    return x;
}

/**
 * Forward Euler 2D (Explicit)
 * u_{i,j}^{n+1} = u_{i,j}^n + rx * d2x(u) + ry * d2y(u)
 */
function forwardEuler2D(u, rx, ry) {
    const Nx = u.length;
    const Ny = u[0].length;
    const uNew = new Array(Nx).fill(0).map(() => new Array(Ny).fill(0));

    // Boundary conditions are 0 (Dirichlet)

    for (let i = 1; i < Nx - 1; i++) {
        for (let j = 1; j < Ny - 1; j++) {
            uNew[i][j] = u[i][j] +
                rx * (u[i + 1][j] - 2 * u[i][j] + u[i - 1][j]) +
                ry * (u[i][j + 1] - 2 * u[i][j] + u[i][j - 1]);
        }
    }
    return uNew;
}

/**
 * Backward Euler ADI (Peaceman-Rachford splitting)
 * Step 1: Implicit in X, Explicit in Y (t -> t + dt/2)
 * Step 2: Explicit in X, Implicit in Y (t + dt/2 -> t + dt)
 * Note: Standard BE is fully implicit, but ADI is a common approximation.
 * For true BE, we'd need a 2D solver. Here we implement ADI which is O(dt^2) usually,
 * but let's stick to the ADI structure for implicit handling.
 */
function backwardEulerADI(u, rx, ry) {
    const Nx = u.length;
    const Ny = u[0].length;

    // We'll use the standard Peaceman-Rachford ADI scheme
    // But technically that's Crank-Nicolson-like (2nd order).
    // For "Backward Euler" feel, we can just do splitting:
    // (1 - rx d2x) u* = (1 + ry d2y) u^n  <-- This is actually CN-ADI

    // Let's implement Peaceman-Rachford (which is 2nd order accurate in time)
    // This effectively serves as our "Crank-Nicolson" equivalent in 2D.
    // Since the user asked for BE and CN, we can differentiate by the coefficients?
    // Actually, standard ADI IS the way to do implicit in 2D efficiently.
    // Let's implement Peaceman-Rachford for both, maybe just varying the weighting?
    // No, let's implement Douglas-Rachford for BE?

    // Let's stick to Peaceman-Rachford for "Crank-Nicolson ADI".
    // For "Backward Euler", we can do:
    // (1 - rx d2x) u* = u^n
    // (1 - ry d2y) u^{n+1} = u*
    // This is splitting the operator (1 - dt L) into (1 - dt Lx)(1 - dt Ly).
    // This is O(dt) + O(dx^2) and unconditionally stable.

    // Step 1: Implicit X
    const uStar = new Array(Nx).fill(0).map(() => new Array(Ny).fill(0));

    // Solve for each row j
    for (let j = 1; j < Ny - 1; j++) {
        const a = new Array(Nx).fill(-rx);
        const b = new Array(Nx).fill(1 + 2 * rx);
        const c = new Array(Nx).fill(-rx);
        const d = new Array(Nx);

        for (let i = 0; i < Nx; i++) {
            d[i] = u[i][j]; // RHS is just u^n
        }

        // BCs
        b[0] = 1; c[0] = 0; d[0] = 0;
        a[Nx - 1] = 0; b[Nx - 1] = 1; d[Nx - 1] = 0;

        const res = thomasAlgorithm(a, b, c, d);
        for (let i = 0; i < Nx; i++) uStar[i][j] = res[i];
    }

    // Step 2: Implicit Y
    const uNew = new Array(Nx).fill(0).map(() => new Array(Ny).fill(0));

    // Solve for each column i
    for (let i = 1; i < Nx - 1; i++) {
        const a = new Array(Ny).fill(-ry);
        const b = new Array(Ny).fill(1 + 2 * ry);
        const c = new Array(Ny).fill(-ry);
        const d = new Array(Ny);

        for (let j = 0; j < Ny; j++) {
            d[j] = uStar[i][j]; // RHS is u*
        }

        // BCs
        b[0] = 1; c[0] = 0; d[0] = 0;
        a[Ny - 1] = 0; b[Ny - 1] = 1; d[Ny - 1] = 0;

        const res = thomasAlgorithm(a, b, c, d);
        for (let j = 0; j < Ny; j++) uNew[i][j] = res[j];
    }

    return uNew;
}

/**
 * Crank-Nicolson ADI (Peaceman-Rachford)
 * Step 1 (t -> t + dt/2): Implicit X, Explicit Y
 * (1 - rx/2 d2x) u* = (1 + ry/2 d2y) u^n
 * Step 2 (t + dt/2 -> t + dt): Explicit X, Implicit Y
 * (1 - ry/2 d2y) u^{n+1} = (1 + rx/2 d2x) u*
 */
function crankNicolsonADI(u, rx, ry) {
    const Nx = u.length;
    const Ny = u[0].length;
    const uStar = new Array(Nx).fill(0).map(() => new Array(Ny).fill(0));
    const uNew = new Array(Nx).fill(0).map(() => new Array(Ny).fill(0));

    const rx2 = rx / 2;
    const ry2 = ry / 2;

    // Step 1: Implicit X, Explicit Y
    for (let j = 1; j < Ny - 1; j++) {
        const a = new Array(Nx).fill(-rx2);
        const b = new Array(Nx).fill(1 + rx); // 2*rx/2 = rx
        const c = new Array(Nx).fill(-rx2);
        const d = new Array(Nx);

        for (let i = 1; i < Nx - 1; i++) {
            // Explicit Y part
            d[i] = u[i][j] + ry2 * (u[i][j + 1] - 2 * u[i][j] + u[i][j - 1]);
        }
        // BCs
        b[0] = 1; c[0] = 0; d[0] = 0;
        a[Nx - 1] = 0; b[Nx - 1] = 1; d[Nx - 1] = 0;

        const res = thomasAlgorithm(a, b, c, d);
        for (let i = 0; i < Nx; i++) uStar[i][j] = res[i];
    }

    // Step 2: Explicit X, Implicit Y
    for (let i = 1; i < Nx - 1; i++) {
        const a = new Array(Ny).fill(-ry2);
        const b = new Array(Ny).fill(1 + ry);
        const c = new Array(Ny).fill(-ry2);
        const d = new Array(Ny);

        for (let j = 1; j < Ny - 1; j++) {
            // Explicit X part using u*
            d[j] = uStar[i][j] + rx2 * (uStar[i + 1][j] - 2 * uStar[i][j] + uStar[i - 1][j]);
        }
        // BCs
        b[0] = 1; c[0] = 0; d[0] = 0;
        a[Ny - 1] = 0; b[Ny - 1] = 1; d[Ny - 1] = 0;

        const res = thomasAlgorithm(a, b, c, d);
        for (let j = 0; j < Ny; j++) uNew[i][j] = res[j];
    }

    return uNew;
}

// ========================================
// Time Stepping
// ========================================

function stepNumerical() {
    const rx = state.rx;
    const ry = state.ry;
    const N = config.Nx; // N x N grid

    switch (config.method) {
        case 'forward-euler':
            state.u = forwardEuler2D(state.u, rx, ry);
            // ~10 FLOPs per point
            state.flops += 10 * N * N;
            break;
        case 'backward-euler':
            state.u = backwardEulerADI(state.u, rx, ry);
            // 2 tridiagonal solves per point (one x, one y)
            // ~16 FLOPs per point
            state.flops += 16 * N * N;
            break;
        case 'crank-nicolson':
            state.u = crankNicolsonADI(state.u, rx, ry);
            // 2 tridiagonal solves + explicit parts
            // ~24 FLOPs per point
            state.flops += 24 * N * N;
            break;
    }
}

function timeStep() {
    stepNumerical();
    state.t += config.dt;
    state.nSteps++;
    state.uExact = getAnalyticalSolution2D(state.x, state.y, state.t);
}

// ========================================
// Visualization
// ========================================

function initPlots() {
    // 1. Surface Plot (Plotly)
    const data = [{
        z: state.u,
        x: state.x,
        y: state.y,
        type: 'surface',
        colorscale: 'Viridis',
        showscale: false,
        contours: {
            z: { show: true, usecolormap: true, highlightcolor: "#42f462", project: { z: true } }
        }
    }];

    const layout = {
        title: 'Temperature Distribution u(x,y)',
        autosize: true,
        uirevision: 'true', // Preserve state (camera, zoom) on updates
        margin: { l: 0, r: 0, b: 0, t: 30 },
        scene: {
            xaxis: { title: 'x', range: [0, 1], autorange: false },
            yaxis: { title: 'y', range: [0, 1], autorange: false },
            zaxis: { title: 'u', range: [-1, 1], autorange: false },
            aspectmode: 'manual',
            aspectratio: { x: 1.2, y: 1.2, z: 0.6 },
            camera: {
                eye: { x: 1.1, y: 1.1, z: 1.1 }
            }
        },
        paper_bgcolor: 'rgba(0,0,0,0)',
        plot_bgcolor: 'rgba(0,0,0,0)',
        font: { family: 'Outfit, sans-serif', color: '#ffffff' }
    };

    Plotly.newPlot('surface-plot', data, layout, { responsive: true, displayModeBar: false });

    // 2. Error Plot (Chart.js)
    Chart.defaults.font.family = "'Outfit', sans-serif";
    const ctx = document.getElementById('error-chart').getContext('2d');
    state.errorChart = new Chart(ctx, {
        type: 'line',
        data: {
            labels: state.x.map(xi => xi.toFixed(2)),
            datasets: [
                {
                    label: 'Analytical (y=0.5)',
                    data: [],
                    borderColor: '#ffa500', // Orange
                    borderDash: [5, 5],
                    borderWidth: 3,
                    pointRadius: 0,
                },
                {
                    label: 'Numerical (y=0.5)',
                    data: [],
                    borderColor: '#00e5ff', // Bright Cyan
                    borderWidth: 2,
                    pointRadius: 0,
                }
            ]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            animation: { duration: 0 },
            scales: {
                x: {
                    title: { display: true, text: 'x (at y=0.5)', color: '#aaa' },
                    ticks: { color: '#888' },
                    grid: { color: '#333' }
                },
                y: {
                    title: { display: true, text: 'u(x, 0.5)', color: '#aaa' },
                    ticks: { color: '#888' },
                    grid: { color: '#333' }
                }
            },
            plugins: {
                legend: { labels: { color: '#ccc' } }
            }
        }
    });
}

function updatePlots() {
    // Update Surface Plot

    // Must preserve layout to keep dark theme
    const layout = {
        title: 'Temperature Distribution u(x,y)',
        autosize: true,
        uirevision: 'true', // Preserve state (camera, zoom) on updates
        margin: { l: 0, r: 0, b: 0, t: 30 },
        scene: {
            xaxis: { title: 'x', range: [0, 1], autorange: false },
            yaxis: { title: 'y', range: [0, 1], autorange: false },
            zaxis: { title: 'u', range: [-1, 1], autorange: false },
            aspectmode: 'manual',
            aspectratio: { x: 1.2, y: 1.2, z: 0.6 },
            camera: {
                eye: { x: 1.1, y: 1.1, z: 1.1 }
            }
        },
        paper_bgcolor: 'rgba(0,0,0,0)',
        plot_bgcolor: 'rgba(0,0,0,0)',
        font: { family: 'Outfit, sans-serif', color: '#ffffff' }
    };

    Plotly.react('surface-plot', [{
        z: state.u,
        x: state.x,
        y: state.y,
        type: 'surface',
        colorscale: 'Viridis',
        showscale: false,
        contours: {
            z: { show: true, usecolormap: true, highlightcolor: "#42f462", project: { z: true } }
        }
    }], layout);

    // Update Error Plot (Slice at y=0.5)
    const diagNumerical = [];
    const diagAnalytical = [];

    // Extract values at constant y = 0.5 (midpoint index for j)
    const jMid = Math.floor(config.Nx / 2);

    // Iterate over x (index i)
    for (let i = 0; i < config.Nx; i++) {
        diagNumerical.push(state.u[i][jMid]);
        diagAnalytical.push(state.uExact[i][jMid]);
    }

    state.errorChart.data.labels = state.x.map(xi => xi.toFixed(2));
    state.errorChart.data.datasets[0].data = diagAnalytical;
    state.errorChart.data.datasets[1].data = diagNumerical;
    state.errorChart.update('none');
}

function updateStatus() {
    document.getElementById('time-display').textContent = `t = ${state.t.toFixed(3)}`;
    // Show max(rx, ry) or just rx since grid is square
    document.getElementById('r-display').textContent = `r = ${state.rx.toFixed(4)}`;
    document.getElementById('steps-display').textContent = state.nSteps;
    document.getElementById('flops-display').textContent = state.flops.toExponential(2);

    const methodNames = {
        'forward-euler': 'Forward Euler',
        'backward-euler': 'Backward Euler ADI',
        'crank-nicolson': 'Crank-Nicolson ADI'
    };
    document.getElementById('method-display').textContent = methodNames[config.method];

    // Stability check
    const isUnstable = config.method === 'forward-euler' && (state.rx + state.ry) > 0.5;
    const stabilityStatus = document.getElementById('stability-status');

    if (isUnstable) {
        stabilityStatus.textContent = 'UNSTABLE (r > 0.25)';
        stabilityStatus.style.background = 'rgba(239, 68, 68, 0.2)';
        stabilityStatus.style.color = '#ef4444';
    } else {
        stabilityStatus.textContent = 'Stable';
        stabilityStatus.style.background = 'rgba(16, 185, 129, 0.2)';
        stabilityStatus.style.color = '#10b981';
    }

    // Error Calculation
    let sumSqError = 0;
    let maxError = 0;
    let maxUExact = 0;
    let count = 0;

    for (let i = 0; i < config.Nx; i++) {
        for (let j = 0; j < config.Nx; j++) {
            const err = Math.abs(state.u[i][j] - state.uExact[i][j]);
            if (isFinite(err)) { // Only include finite differences in sum
                sumSqError += err * err;
                if (err > maxError) maxError = err;
            }
            if (isFinite(state.uExact[i][j]) && Math.abs(state.uExact[i][j]) > maxUExact) {
                maxUExact = Math.abs(state.uExact[i][j]);
            }
            count++;
        }
    }

    const l2Error = Math.sqrt(sumSqError / count);

    // Calculate Max % Error (relative to max amplitude of exact solution)
    let maxRelError = 0;
    if (maxUExact > 1e-9) { // Avoid division by zero or very small numbers
        maxRelError = (maxError / maxUExact) * 100;
    }

    const errorDisplay = document.getElementById('error-display');
    const maxErrorDisplay = document.getElementById('max-error-display');
    const maxRelErrorDisplay = document.getElementById('max-rel-error-display');

    if (l2Error > 10 || !isFinite(l2Error)) { // Check for divergence or NaN/Infinity
        errorDisplay.textContent = 'DIVERGED!';
        errorDisplay.style.color = '#ef4444';
        maxErrorDisplay.textContent = 'DIVERGED!';
        maxErrorDisplay.style.color = '#ef4444';
        maxRelErrorDisplay.textContent = 'DIVERGED!';
        maxRelErrorDisplay.style.color = '#ef4444';
    } else {
        errorDisplay.textContent = l2Error.toExponential(3);
        errorDisplay.style.color = '';
        maxErrorDisplay.textContent = maxError.toExponential(3);
        maxErrorDisplay.style.color = '';
        maxRelErrorDisplay.textContent = maxRelError.toFixed(2) + '%';
        maxRelErrorDisplay.style.color = '';
    }
}

// ========================================
// Simulation Control
// ========================================

function initSimulation() {
    // Grid setup
    state.dx = config.L / (config.Nx - 1);
    state.dy = config.L / (config.Nx - 1); // Square grid

    state.x = Array.from({ length: config.Nx }, (_, i) => i * state.dx);
    state.y = Array.from({ length: config.Nx }, (_, i) => i * state.dy);

    state.rx = config.alpha * config.dt / (state.dx * state.dx);
    state.ry = config.alpha * config.dt / (state.dy * state.dy);

    state.t = 0;
    state.nSteps = 0;
    state.flops = 0;
    state.hasDiverged = false;

    // Initial Conditions
    state.u = getSineIC2D(state.x, state.y);
    state.uExact = getAnalyticalSolution2D(state.x, state.y, 0);

    // Initialize plots if first run, else update
    if (!state.errorChart) {
        initPlots();
    }
    updatePlots();
    updateStatus();
}

function animate() {
    if (!state.isPlaying) return;

    try {
        // Run fewer steps per frame for 2D as it's heavier
        for (let i = 0; i < config.animationSpeed; i++) {
            if (state.t >= config.Tmax) {
                updatePlots();
                updateStatus();
                stopAnimation();
                return;
            }

            timeStep();

            // Check divergence
            const centerVal = Math.abs(state.u[Math.floor(config.Nx / 2)][Math.floor(config.Nx / 2)]);
            if (centerVal > 10 || !isFinite(centerVal)) {
                if (!state.hasDiverged) {
                    state.hasDiverged = true;
                    updateStatus();
                }
            }
        }
        updatePlots();
        updateStatus();
    } catch (e) {
        console.error(e);
        stopAnimation();
    }

    state.animationId = requestAnimationFrame(animate);
}

function startAnimation() {
    if (state.isPlaying) return;
    state.isPlaying = true;
    document.getElementById('play-btn').innerHTML = '<span id="play-icon">⏸</span> Pause';
    animate();
}

function stopAnimation() {
    state.isPlaying = false;
    document.getElementById('play-btn').innerHTML = '<span id="play-icon">▶</span> Play';
    if (state.animationId) {
        cancelAnimationFrame(state.animationId);
        state.animationId = null;
    }
}

function togglePlayPause() {
    if (state.isPlaying) stopAnimation();
    else startAnimation();
}

function reset() {
    stopAnimation();
    state.t = 0;
    state.nSteps = 0;
    state.flops = 0;
    state.hasDiverged = false;
    state.u = getSineIC2D(state.x, state.y);
    state.uExact = getAnalyticalSolution2D(state.x, state.y, 0);
    updatePlots();
    updateStatus();
}

// ========================================
// UI Event Handlers
// ========================================

document.addEventListener('DOMContentLoaded', () => {
    // Initialize inputs
    document.getElementById('dt-input').value = config.dt;
    document.getElementById('dx-slider').value = config.Nx;
    document.getElementById('tmax-slider').value = config.Tmax;
    document.getElementById('speed-slider').value = config.animationSpeed;

    // Displays
    document.getElementById('dt-value').textContent = config.dt.toExponential(2);
    document.getElementById('dx-value').textContent = config.Nx;
    document.getElementById('tmax-value').textContent = config.Tmax.toFixed(1);
    document.getElementById('speed-value').textContent = `${config.animationSpeed}x`;

    initSimulation();

    // Event Listeners

    // Method Buttons
    document.querySelectorAll('.method-btn').forEach(btn => {
        btn.addEventListener('click', () => {
            document.querySelectorAll('.method-btn').forEach(b => b.classList.remove('active'));
            btn.classList.add('active');
            config.method = btn.dataset.method;
            updateMethodInfo();
            reset();
        });
    });

    // dt Input
    const dtInput = document.getElementById('dt-input');
    dtInput.addEventListener('change', () => {
        let val = parseFloat(dtInput.value);
        if (val < 0.0001) val = 0.0001;
        if (val > 1.0) val = 1.0;
        dtInput.value = val;
        config.dt = val;
        document.getElementById('dt-value').textContent = config.dt.toExponential(2);
        state.rx = config.alpha * config.dt / (state.dx * state.dx);
        state.ry = config.alpha * config.dt / (state.dy * state.dy);
        updateStatus();
        reset();
    });

    // Nx Slider
    const dxSlider = document.getElementById('dx-slider');
    dxSlider.addEventListener('input', () => {
        config.Nx = parseInt(dxSlider.value);
        document.getElementById('dx-value').textContent = config.Nx;
        initSimulation(); // Re-init grid
    });

    // Tmax Slider
    const tmaxSlider = document.getElementById('tmax-slider');
    tmaxSlider.addEventListener('input', () => {
        config.Tmax = parseFloat(tmaxSlider.value);
        document.getElementById('tmax-value').textContent = config.Tmax.toFixed(1);
    });

    // Speed Slider
    const speedSlider = document.getElementById('speed-slider');
    speedSlider.addEventListener('input', () => {
        config.animationSpeed = parseInt(speedSlider.value);
        document.getElementById('speed-value').textContent = `${config.animationSpeed}x`;
    });

    // Controls
    document.getElementById('play-btn').addEventListener('click', togglePlayPause);
    document.getElementById('reset-btn').addEventListener('click', reset);

    // Tabs
    document.querySelectorAll('.info-tab').forEach(tab => {
        tab.addEventListener('click', () => {
            document.querySelectorAll('.info-tab').forEach(t => t.classList.remove('active'));
            document.querySelectorAll('.tab-content').forEach(c => c.classList.remove('active'));
            tab.classList.add('active');
            document.getElementById(`tab-${tab.dataset.tab}`).classList.add('active');
        });
    });

    updateMethodInfo();
});

function updateMethodInfo() {
    const title = document.getElementById('method-title');
    const desc = document.getElementById('method-description');

    const info = {
        'forward-euler': {
            title: 'Forward Euler (Explicit)',
            desc: `First-order accurate, <strong>conditionally stable</strong>: \\(r_x + r_y \\le 0.5\\).`
        },
        'backward-euler': {
            title: 'Backward Euler ADI',
            desc: `Solved using <strong>ADI splitting</strong>. Unconditionally stable, \\(O(\\Delta t) + O(\\Delta x^2)\\).`
        },
        'crank-nicolson': {
            title: 'Crank-Nicolson ADI',
            desc: `Solved using <strong>Peaceman-Rachford ADI</strong>. Unconditionally stable, \\(O(\\Delta t^2) + O(\\Delta x^2)\\).`
        }
    };

    title.textContent = info[config.method].title;
    desc.innerHTML = info[config.method].desc;

    if (window.MathJax && window.MathJax.typesetPromise) {
        window.MathJax.typesetPromise([desc]).catch(e => console.log(e));
    }
}
