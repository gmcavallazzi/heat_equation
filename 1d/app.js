/**
 * 1D Heat Equation Explorer
 * Interactive visualization of numerical methods
 * Sine Wave IC with Exact Analytical Solution
 */

// ========================================
// Configuration & State
// ========================================

const config = {
    L: 1.0,           // Domain length
    alpha: 0.01,      // Thermal diffusivity (Fixed)
    Nx: 41,           // Number of grid points (Default)
    dt: 0.001,        // Time step
    Tmax: 2.0,        // Maximum simulation time
    method: 'forward-euler',
    animationSpeed: 5,
};

const state = {
    t: 0,
    nSteps: 0,        // Number of time steps taken
    flops: 0,         // Total floating point operations
    u: null,          // Current numerical solution
    uExact: null,     // Exact analytical solution
    x: null,          // Grid points
    dx: null,
    r: null,          // Fourier number
    isPlaying: false,
    animationId: null,
    chart: null,
    hasDiverged: false,
};

// ========================================
// Initial Condition & Analytical Solution
// ========================================

/**
 * Initial Condition: u(x,0) = sin(3πx/L)
 * Higher frequency (k=3) to amplify numerical errors
 */
function getSineIC(x) {
    const L = config.L;
    const k = 3;
    return x.map(xi => Math.sin(k * Math.PI * xi / L));
}

/**
 * Exact Analytical Solution
 * u(x,t) = sin(3πx/L) * exp(-α * (3π)² * t / L²)
 */
function getAnalyticalSolution(x, t) {
    const L = config.L;
    const alpha = config.alpha;
    const k = 3;
    const decay = Math.exp(-alpha * Math.pow(k * Math.PI, 2) * t / (L * L));
    return x.map(xi => Math.sin(k * Math.PI * xi / L) * decay);
}

// ========================================
// Numerical Solvers
// ========================================

/**
 * Forward Euler (Explicit) - O(dt) + O(dx²)
 * Stability: r <= 0.5
 */
function forwardEuler(u, r) {
    const n = u.length;
    const uNew = new Array(n);

    uNew[0] = 0;
    uNew[n - 1] = 0;

    for (let i = 1; i < n - 1; i++) {
        uNew[i] = u[i] + r * (u[i + 1] - 2 * u[i] + u[i - 1]);
    }

    return uNew;
}

/**
 * Backward Euler (Implicit) - O(dt) + O(dx²)
 * Unconditionally stable - Uses Thomas Algorithm
 */
function backwardEuler(u, r) {
    const n = u.length;

    const a = new Array(n).fill(-r);
    const b = new Array(n).fill(1 + 2 * r);
    const c = new Array(n).fill(-r);
    const d = [...u];

    // Boundary conditions
    b[0] = 1; c[0] = 0; d[0] = 0;
    a[n - 1] = 0; b[n - 1] = 1; d[n - 1] = 0;

    return thomasAlgorithm(a, b, c, d);
}

/**
 * Crank-Nicolson - O(dt²) + O(dx²)
 * Unconditionally stable
 * LHS: (1 + r)*T_i - (r/2)*T_{i-1} - (r/2)*T_{i+1} (at n+1)
 * RHS: (1 - r)*T_i + (r/2)*T_{i-1} + (r/2)*T_{i+1} (at n)
 */
function crankNicolson(u, r) {
    const n = u.length;
    const r2 = r / 2;

    const a = new Array(n).fill(-r2);
    const b = new Array(n).fill(1 + r);
    const c = new Array(n).fill(-r2);
    const d = new Array(n);

    for (let i = 1; i < n - 1; i++) {
        d[i] = (1 - r) * u[i] + r2 * (u[i + 1] + u[i - 1]);
    }

    // Boundary conditions
    b[0] = 1; c[0] = 0; d[0] = 0;
    a[n - 1] = 0; b[n - 1] = 1; d[n - 1] = 0;

    return thomasAlgorithm(a, b, c, d);
}

/**
 * Thomas Algorithm for tridiagonal systems
 * Solves: a_i * x_{i-1} + b_i * x_i + c_i * x_{i+1} = d_i
 * Robust implementation
 */
function thomasAlgorithm(a, b, c, d) {
    const n = d.length;
    const cp = new Array(n); // c prime
    const dp = new Array(n); // d prime
    const x = new Array(n);

    // Forward sweep
    if (Math.abs(b[0]) < 1e-15) {
        console.error("Thomas Algorithm: Zero pivot at index 0");
        return d;
    }

    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];

    for (let i = 1; i < n; i++) {
        const denom = b[i] - a[i] * cp[i - 1];
        if (Math.abs(denom) < 1e-15) {
            console.error(`Thomas Algorithm: Zero pivot at index ${i}`);
            return d;
        }
        cp[i] = c[i] / denom;
        dp[i] = (d[i] - a[i] * dp[i - 1]) / denom;
    }

    // Back substitution
    x[n - 1] = dp[n - 1];
    for (let i = n - 2; i >= 0; i--) {
        x[i] = dp[i] - cp[i] * x[i + 1];
    }

    return x;
}

// ========================================
// Time Stepping
// ========================================

function stepNumerical() {
    const r = state.r;
    const N = config.Nx;

    if (state.nSteps === 0) {
        console.log(`Starting simulation with method: ${config.method}, r=${r.toFixed(4)}`);
    }

    switch (config.method) {
        case 'forward-euler':
            state.u = forwardEuler(state.u, r);
            // Estimate: 5 FLOPs per point (inner loop)
            state.flops += 5 * (N - 2);
            break;
        case 'backward-euler':
            state.u = backwardEuler(state.u, r);
            // Estimate: Thomas Algorithm ~8 FLOPs per point
            state.flops += 8 * N;
            break;
        case 'crank-nicolson':
            state.u = crankNicolson(state.u, r);
            // Estimate: RHS (~4 ops) + Thomas Algorithm (~8 ops) = ~12 FLOPs per point
            state.flops += 12 * N;
            break;
        default:
            console.error("Unknown method:", config.method);
    }
}

function timeStep() {
    stepNumerical();
    state.t += config.dt;
    state.nSteps++;

    // Update exact solution
    state.uExact = getAnalyticalSolution(state.x, state.t);
}

// ========================================
// Visualization
// ========================================

function initChart() {
    const ctx = document.getElementById('solution-chart').getContext('2d');

    state.chart = new Chart(ctx, {
        type: 'line',
        data: {
            labels: state.x.map(xi => xi.toFixed(2)),
            datasets: [
                {
                    label: 'Numerical Solution',
                    data: state.u,
                    borderColor: 'rgb(99, 102, 241)',
                    backgroundColor: 'rgba(99, 102, 241, 0.1)',
                    borderWidth: 3,
                    fill: true,
                    tension: 0.1,
                    pointRadius: 3,
                    pointBackgroundColor: 'rgb(99, 102, 241)',
                },
                {
                    label: 'Analytical Solution',
                    data: state.uExact,
                    borderColor: 'rgba(16, 185, 129, 0.9)',
                    borderWidth: 2,
                    borderDash: [5, 5],
                    fill: false,
                    tension: 0.4, // Smooth curve for analytical
                    pointRadius: 0,
                },
            ]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            animation: { duration: 0 },
            scales: {
                x: {
                    title: {
                        display: true,
                        text: 'Position x',
                        color: 'rgba(255, 255, 255, 0.7)',
                        font: { size: 14, weight: 500 }
                    },
                    ticks: { color: 'rgba(255, 255, 255, 0.5)' },
                    grid: { color: 'rgba(255, 255, 255, 0.1)' }
                },
                y: {
                    title: {
                        display: true,
                        text: 'Temperature T',
                        color: 'rgba(255, 255, 255, 0.7)',
                        font: { size: 14, weight: 500 }
                    },
                    ticks: { color: 'rgba(255, 255, 255, 0.5)' },
                    grid: { color: 'rgba(255, 255, 255, 0.1)' },
                    min: -1.0,
                    max: 1.0,
                }
            },
            plugins: {
                legend: {
                    labels: {
                        color: 'rgba(255, 255, 255, 0.7)',
                        font: { size: 12 }
                    }
                }
            }
        }
    });
}

function updateChart() {
    if (!state.chart) return;

    // Clamp values to prevent chart issues but allow visual divergence
    const clamp = (v) => {
        if (!isFinite(v)) return v > 0 ? 1000 : -1000;
        return Math.max(-1000, Math.min(1000, v));
    };

    state.chart.data.labels = state.x.map(xi => xi.toFixed(2));
    state.chart.data.datasets[0].data = state.u.map(clamp);
    state.chart.data.datasets[1].data = state.uExact;
    state.chart.update('none');
}

function updateStatus() {
    document.getElementById('time-display').textContent = `t = ${state.t.toFixed(3)}`;
    document.getElementById('r-display').textContent = `r = ${state.r.toFixed(4)}`;
    document.getElementById('steps-display').textContent = state.nSteps;
    document.getElementById('flops-display').textContent = state.flops.toExponential(2);

    const methodNames = {
        'forward-euler': 'Forward Euler',
        'backward-euler': 'Backward Euler',
        'crank-nicolson': 'Crank-Nicolson'
    };
    document.getElementById('method-display').textContent = methodNames[config.method];

    const isUnstable = config.method === 'forward-euler' && state.r > 0.5;
    const stabilityStatus = document.getElementById('stability-status');

    if (isUnstable) {
        stabilityStatus.textContent = 'UNSTABLE (r > 0.5)';
        stabilityStatus.style.background = 'rgba(239, 68, 68, 0.2)';
        stabilityStatus.style.color = '#ef4444';
    } else {
        stabilityStatus.textContent = 'Stable';
        stabilityStatus.style.background = 'rgba(16, 185, 129, 0.2)';
        stabilityStatus.style.color = '#10b981';
    }

    // Compute errors vs Analytical
    let sum = 0;
    let maxErr = 0;
    for (let i = 0; i < state.u.length; i++) {
        const diff = Math.abs(state.u[i] - state.uExact[i]);
        if (isFinite(diff)) {
            sum += diff * diff;
            if (diff > maxErr) maxErr = diff;
        }
    }
    const l2Error = Math.sqrt(sum / state.u.length);

    const errorDisplay = document.getElementById('error-display');
    const maxErrorDisplay = document.getElementById('max-error-display');

    if (l2Error > 1) {
        errorDisplay.textContent = 'DIVERGED!';
        errorDisplay.style.color = '#ef4444';
        maxErrorDisplay.textContent = 'DIVERGED!';
        maxErrorDisplay.style.color = '#ef4444';
    } else {
        errorDisplay.textContent = l2Error.toExponential(3);
        errorDisplay.style.color = '';
        maxErrorDisplay.textContent = maxErr.toExponential(3);
        maxErrorDisplay.style.color = '';
    }
}

// ========================================
// Simulation Control
// ========================================

function initSimulation() {
    // Numerical grid
    state.dx = config.L / (config.Nx - 1);
    state.x = Array.from({ length: config.Nx }, (_, i) => i * state.dx);
    state.r = config.alpha * config.dt / (state.dx * state.dx);
    state.t = 0;
    state.nSteps = 0;
    state.flops = 0;
    state.hasDiverged = false;

    document.getElementById('error-display').style.color = '';
    document.getElementById('max-error-display').style.color = '';

    // Initial condition
    state.u = getSineIC(state.x);
    state.uExact = getAnalyticalSolution(state.x, 0);

    if (state.chart) {
        updateChart();
    }
    updateStatus();
}

function animate() {
    if (!state.isPlaying) return;

    try {
        for (let i = 0; i < config.animationSpeed; i++) {
            if (state.t >= config.Tmax) {
                updateChart();
                updateStatus();
                stopAnimation();
                return;
            }

            timeStep();

            const maxVal = Math.max(...state.u.map(Math.abs));
            if (maxVal > 10 || state.u.some(v => !isFinite(v))) {
                if (!state.hasDiverged) {
                    state.hasDiverged = true;
                    updateStatus();
                }
            }
        }
        updateChart();
        updateStatus();
    } catch (e) {
        console.error('Animation error:', e);
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
    initSimulation();
    updateChart();
    updateStatus();
}

function resetToDefaults() {
    // Reset controls to minimum/default values
    const dtInput = document.getElementById('dt-input');
    const dxSlider = document.getElementById('dx-slider');

    dtInput.value = "0.001";
    dxSlider.value = "41"; // Default to 41

    // Update config
    config.dt = parseFloat(dtInput.value);
    config.Nx = parseInt(dxSlider.value);

    // Update displays
    document.getElementById('dt-value').textContent = config.dt.toExponential(2);
    document.getElementById('dx-value').textContent = config.Nx;

    // Re-initialize
    reset();
}

// ========================================
// UI Event Handlers
// ========================================

function setupEventListeners() {
    // Method buttons
    document.querySelectorAll('.method-btn').forEach(btn => {
        btn.addEventListener('click', () => {
            document.querySelectorAll('.method-btn').forEach(b => b.classList.remove('active'));
            btn.classList.add('active');
            config.method = btn.dataset.method;
            updateMethodInfo();
            updateStatus();
            reset();
        });
    });

    // dt input
    const dtInput = document.getElementById('dt-input');
    dtInput.addEventListener('change', () => {
        let val = parseFloat(dtInput.value);
        // Clamp value
        if (val < 0.0001) val = 0.0001;
        if (val > 1.0) val = 1.0;
        dtInput.value = val;

        config.dt = val;
        document.getElementById('dt-value').textContent = config.dt.toExponential(2);
        state.r = config.alpha * config.dt / (state.dx * state.dx);
        updateStatus();
        reset();
    });

    // Nx slider
    const dxSlider = document.getElementById('dx-slider');
    dxSlider.addEventListener('input', () => {
        config.Nx = parseInt(dxSlider.value);
        document.getElementById('dx-value').textContent = config.Nx;
        reset();
    });

    // Alpha slider removed - fixed at 1.0

    // Tmax slider
    const tmaxSlider = document.getElementById('tmax-slider');
    tmaxSlider.addEventListener('input', () => {
        config.Tmax = parseFloat(tmaxSlider.value);
        document.getElementById('tmax-value').textContent = config.Tmax.toFixed(1);
    });

    // Speed slider
    const speedSlider = document.getElementById('speed-slider');
    speedSlider.addEventListener('input', () => {
        config.animationSpeed = parseInt(speedSlider.value);
        document.getElementById('speed-value').textContent = `${config.animationSpeed}x`;
    });
    // Play/Pause & Reset
    document.getElementById('play-btn').addEventListener('click', togglePlayPause);
    document.getElementById('reset-btn').addEventListener('click', reset);

    document.querySelectorAll('.info-tab').forEach(tab => {
        tab.addEventListener('click', () => {
            document.querySelectorAll('.info-tab').forEach(t => t.classList.remove('active'));
            document.querySelectorAll('.tab-content').forEach(c => c.classList.remove('active'));
            tab.classList.add('active');
            document.getElementById(`tab-${tab.dataset.tab}`).classList.add('active');
        });
    });
}

function updateMethodInfo() {
    const title = document.getElementById('method-title');
    const desc = document.getElementById('method-description');

    const info = {
        'forward-euler': {
            title: 'Forward Euler (Explicit)',
            desc: `First-order accurate, <strong>conditionally stable</strong>: \\(r = \\alpha \\Delta t / \\Delta x^2\\) must satisfy \\(r \\le 0.5\\). If \\(r > 0.5\\), the solution will oscillate and diverge!`
        },
        'backward-euler': {
            title: 'Backward Euler (Implicit)',
            desc: `First-order accurate, <strong>unconditionally stable</strong>. Large time steps allowed but accuracy suffers. Strongly damps oscillations.`
        },
        'crank-nicolson': {
            title: 'Crank-Nicolson (Implicit)',
            desc: `<strong>Second-order accurate</strong>, unconditionally stable. Best balance of stability and accuracy for smooth solutions.`
        }
    };

    title.textContent = info[config.method].title;
    desc.innerHTML = info[config.method].desc;

    // Trigger MathJax re-render if available
    if (window.MathJax && window.MathJax.typesetPromise) {
        window.MathJax.typesetPromise([desc]).catch((err) => console.log('MathJax error:', err));
    }
}

// ========================================
// Initialization
// ========================================

document.addEventListener('DOMContentLoaded', () => {
    // Initialize controls with config values
    document.getElementById('dt-input').value = config.dt;
    document.getElementById('dx-slider').value = config.Nx;
    document.getElementById('tmax-slider').value = config.Tmax;
    document.getElementById('speed-slider').value = config.animationSpeed;

    // Update displays
    document.getElementById('dt-value').textContent = config.dt.toExponential(2);
    document.getElementById('dx-value').textContent = config.Nx;
    // document.getElementById('alpha-value').textContent = config.alpha.toFixed(4); // Removed
    document.getElementById('tmax-value').textContent = config.Tmax.toFixed(1);
    document.getElementById('speed-value').textContent = `${config.animationSpeed}x`;

    initSimulation();
    initChart();
    setupEventListeners();
    updateStatus();
    updateMethodInfo();

    console.log('1D Heat Equation Explorer initialized');
});
