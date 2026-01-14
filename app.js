/**
 * Unified Heat Equation Explorer
 * Switches between 1D and 2D modes
 */

// ========================================
// Global Dimension State
// ========================================

let currentDimension = '1d'; // Start with 1D

// ========================================
// Shared Configuration
// ========================================

const sharedConfig = {
    L: 1.0,
    alpha: 0.01,
    Nx: 41,
    dt: 0.001,
    Tmax: 2.0,
    method: 'forward-euler',
    animationSpeed: 5,
};

// ========================================
// 1D State and Functions
// ========================================

const state1d = {
    t: 0,
    nSteps: 0,
    flops: 0,
    u: null,
    uExact: null,
    x: null,
    dx: null,
    r: null,
    isPlaying: false,
    animationId: null,
    chart: null,
    hasDiverged: false,
};

function getSineIC(x) {
    const L = sharedConfig.L;
    const k = 3;
    return x.map(xi => Math.sin(k * Math.PI * xi / L));
}

function getAnalyticalSolution(x, t) {
    const L = sharedConfig.L;
    const alpha = sharedConfig.alpha;
    const k = 3;
    const decay = Math.exp(-alpha * Math.pow(k * Math.PI, 2) * t / (L * L));
    return x.map(xi => Math.sin(k * Math.PI * xi / L) * decay);
}

function forwardEuler1D(u, r) {
    const n = u.length;
    const uNew = new Array(n);
    uNew[0] = 0;
    uNew[n - 1] = 0;
    for (let i = 1; i < n - 1; i++) {
        uNew[i] = u[i] + r * (u[i + 1] - 2 * u[i] + u[i - 1]);
    }
    return uNew;
}

function thomasAlgorithm(a, b, c, d) {
    const n = d.length;
    const cp = new Array(n);
    const dp = new Array(n);
    const x = new Array(n);

    if (Math.abs(b[0]) < 1e-15) return d;
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

function backwardEuler1D(u, r) {
    const n = u.length;
    const a = new Array(n).fill(-r);
    const b = new Array(n).fill(1 + 2 * r);
    const c = new Array(n).fill(-r);
    const d = [...u];
    b[0] = 1; c[0] = 0; d[0] = 0;
    a[n - 1] = 0; b[n - 1] = 1; d[n - 1] = 0;
    return thomasAlgorithm(a, b, c, d);
}

function crankNicolson1D(u, r) {
    const n = u.length;
    const r2 = r / 2;
    const a = new Array(n).fill(-r2);
    const b = new Array(n).fill(1 + r);
    const c = new Array(n).fill(-r2);
    const d = new Array(n);

    for (let i = 1; i < n - 1; i++) {
        d[i] = (1 - r) * u[i] + r2 * (u[i + 1] + u[i - 1]);
    }
    b[0] = 1; c[0] = 0; d[0] = 0;
    a[n - 1] = 0; b[n - 1] = 1; d[n - 1] = 0;
    return thomasAlgorithm(a, b, c, d);
}

// ========================================
// 2D State and Functions
// ========================================

const state2d = {
    t: 0,
    nSteps: 0,
    flops: 0,
    u: null,
    uExact: null,
    x: null,
    y: null,
    dx: null,
    dy: null,
    rx: null,
    ry: null,
    isPlaying: false,
    animationId: null,
    errorChart: null,
    hasDiverged: false,
};

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

function getAnalyticalSolution2D(x, y, t) {
    const alpha = sharedConfig.alpha;
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

function forwardEuler2D(u, rx, ry) {
    const Nx = u.length;
    const Ny = u[0].length;
    const uNew = Array(Nx).fill(0).map(() => Array(Ny).fill(0));

    for (let i = 1; i < Nx - 1; i++) {
        for (let j = 1; j < Ny - 1; j++) {
            uNew[i][j] = u[i][j] +
                rx * (u[i + 1][j] - 2 * u[i][j] + u[i - 1][j]) +
                ry * (u[i][j + 1] - 2 * u[i][j] + u[i][j - 1]);
        }
    }
    return uNew;
}

function backwardEuler2D_ADI(u, rx, ry) {
    const Nx = u.length;
    const Ny = u[0].length;
    const uHalf = Array(Nx).fill(0).map(() => Array(Ny).fill(0));
    const uNew = Array(Nx).fill(0).map(() => Array(Ny).fill(0));

    // X-sweep
    for (let j = 0; j < Ny; j++) {
        const a = new Array(Nx).fill(-rx);
        const b = new Array(Nx).fill(1 + 2 * rx);
        const c = new Array(Nx).fill(-rx);
        const d = new Array(Nx);

        for (let i = 1; i < Nx - 1; i++) {
            d[i] = u[i][j];
        }
        b[0] = 1; c[0] = 0; d[0] = 0;
        a[Nx - 1] = 0; b[Nx - 1] = 1; d[Nx - 1] = 0;
        const sol = thomasAlgorithm(a, b, c, d);
        for (let i = 0; i < Nx; i++) {
            uHalf[i][j] = sol[i];
        }
    }

    // Y-sweep
    for (let i = 0; i < Nx; i++) {
        const a = new Array(Ny).fill(-ry);
        const b = new Array(Ny).fill(1 + 2 * ry);
        const c = new Array(Ny).fill(-ry);
        const d = new Array(Ny);

        for (let j = 1; j < Ny - 1; j++) {
            d[j] = uHalf[i][j];
        }
        b[0] = 1; c[0] = 0; d[0] = 0;
        a[Ny - 1] = 0; b[Ny - 1] = 1; d[Ny - 1] = 0;
        const sol = thomasAlgorithm(a, b, c, d);
        for (let j = 0; j < Ny; j++) {
            uNew[i][j] = sol[j];
        }
    }
    return uNew;
}

function crankNicolson2D_ADI(u, rx, ry) {
    const Nx = u.length;
    const Ny = u[0].length;
    const uHalf = Array(Nx).fill(0).map(() => Array(Ny).fill(0));
    const uNew = Array(Nx).fill(0).map(() => Array(Ny).fill(0));
    const rx2 = rx / 2;
    const ry2 = ry / 2;

    // X-sweep
    for (let j = 0; j < Ny; j++) {
        const a = new Array(Nx).fill(-rx2);
        const b = new Array(Nx).fill(1 + rx);
        const c = new Array(Nx).fill(-rx2);
        const d = new Array(Nx);

        for (let i = 1; i < Nx - 1; i++) {
            d[i] = (1 - rx) * u[i][j] + rx2 * (u[i + 1][j] + u[i - 1][j]);
        }
        b[0] = 1; c[0] = 0; d[0] = 0;
        a[Nx - 1] = 0; b[Nx - 1] = 1; d[Nx - 1] = 0;
        const sol = thomasAlgorithm(a, b, c, d);
        for (let i = 0; i < Nx; i++) {
            uHalf[i][j] = sol[i];
        }
    }

    // Y-sweep
    for (let i = 0; i < Nx; i++) {
        const a = new Array(Ny).fill(-ry2);
        const b = new Array(Ny).fill(1 + ry);
        const c = new Array(Ny).fill(-ry2);
        const d = new Array(Ny);

        for (let j = 1; j < Ny - 1; j++) {
            d[j] = (1 - ry) * uHalf[i][j] + ry2 * (uHalf[i][j + 1] + uHalf[i][j - 1]);
        }
        b[0] = 1; c[0] = 0; d[0] = 0;
        a[Ny - 1] = 0; b[Ny - 1] = 1; d[Ny - 1] = 0;
        const sol = thomasAlgorithm(a, b, c, d);
        for (let j = 0; j < Ny; j++) {
            uNew[i][j] = sol[j];
        }
    }
    return uNew;
}

// ========================================
// Dimension Switching
// ========================================

function switchDimension(dimension) {
    if (dimension === currentDimension) return;

    // Stop current animation
    if (currentDimension === '1d' && state1d.isPlaying) {
        stopAnimation1D();
    } else if (currentDimension === '2d' && state2d.isPlaying) {
        stopAnimation2D();
    }

    currentDimension = dimension;

    // Update UI
    document.getElementById('btn-1d').classList.toggle('active', dimension === '1d');
    document.getElementById('btn-2d').classList.toggle('active', dimension === '2d');
    document.getElementById('viz-1d').classList.toggle('active', dimension === '1d');
    document.getElementById('viz-2d').classList.toggle('active', dimension === '2d');
    document.getElementById('dimension-display').textContent = dimension.toUpperCase();

    // Show/hide 2D-specific status items
    document.getElementById('status-max-rel-error').style.display = dimension === '2d' ? '' : 'none';

    // Update equations display
    document.getElementById('eq-1d').style.display = dimension === '1d' ? '' : 'none';
    document.getElementById('eq-2d').style.display = dimension === '2d' ? '' : 'none';

    // Initialize the new dimension
    if (dimension === '1d') {
        initSimulation1D();
        updateChart1D();
        updateStatus1D();
    } else {
        initSimulation2D();
        if (!state2d.errorChart) {
            initErrorChart();
        }
        updateSurfacePlot();
        updateErrorChart();
        updateStatus2D();
    }

    // Update discretisation display for new dimension
    updateDiscretizationDisplay();

    // Re-render MathJax
    if (window.MathJax && window.MathJax.typesetPromise) {
        window.MathJax.typesetPromise().catch((err) => console.log('MathJax error:', err));
    }
}

// ========================================
// 1D Simulation Functions
// ========================================

function initSimulation1D() {
    state1d.dx = sharedConfig.L / (sharedConfig.Nx - 1);
    state1d.x = Array.from({ length: sharedConfig.Nx }, (_, i) => i * state1d.dx);
    state1d.r = sharedConfig.alpha * sharedConfig.dt / (state1d.dx * state1d.dx);
    state1d.t = 0;
    state1d.nSteps = 0;
    state1d.flops = 0;
    state1d.hasDiverged = false;
    state1d.u = getSineIC(state1d.x);
    state1d.uExact = getAnalyticalSolution(state1d.x, 0);
}

function stepNumerical1D() {
    const r = state1d.r;
    const N = sharedConfig.Nx;

    switch (sharedConfig.method) {
        case 'forward-euler':
            state1d.u = forwardEuler1D(state1d.u, r);
            state1d.flops += 5 * (N - 2);
            break;
        case 'backward-euler':
            state1d.u = backwardEuler1D(state1d.u, r);
            state1d.flops += 8 * N;
            break;
        case 'crank-nicolson':
            state1d.u = crankNicolson1D(state1d.u, r);
            state1d.flops += 12 * N;
            break;
    }
}

function timeStep1D() {
    stepNumerical1D();
    state1d.t += sharedConfig.dt;
    state1d.nSteps++;
    state1d.uExact = getAnalyticalSolution(state1d.x, state1d.t);
}

function initChart1D() {
    if (typeof Chart === 'undefined') {
        console.error('Chart.js not loaded yet!');
        return;
    }
    const ctx = document.getElementById('solution-chart').getContext('2d');
    state1d.chart = new Chart(ctx, {
        type: 'line',
        data: {
            labels: state1d.x.map(xi => xi.toFixed(2)),
            datasets: [
                {
                    label: 'Numerical Solution',
                    data: state1d.u,
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
                    data: state1d.uExact,
                    borderColor: 'rgba(16, 185, 129, 0.9)',
                    borderWidth: 2,
                    borderDash: [5, 5],
                    fill: false,
                    tension: 0.4,
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
                    title: { display: true, text: 'Position x', color: 'rgba(255, 255, 255, 0.7)', font: { size: 16, weight: 500, family: 'Outfit' } },
                    ticks: { color: 'rgba(255, 255, 255, 0.5)', font: { size: 14, family: 'Outfit' } },
                    grid: { color: 'rgba(255, 255, 255, 0.1)' }
                },
                y: {
                    title: { display: true, text: 'Temperature T', color: 'rgba(255, 255, 255, 0.7)', font: { size: 16, weight: 500, family: 'Outfit' } },
                    ticks: { color: 'rgba(255, 255, 255, 0.5)', font: { size: 14, family: 'Outfit' } },
                    grid: { color: 'rgba(255, 255, 255, 0.1)' },
                    min: -1.0,
                    max: 1.0,
                }
            },
            plugins: {
                legend: { labels: { color: 'rgba(255, 255, 255, 0.7)', font: { size: 14, family: 'Outfit' } } }
            }
        }
    });
}

function updateChart1D() {
    if (!state1d.chart) return;
    const clamp = (v) => {
        if (!isFinite(v)) return v > 0 ? 1000 : -1000;
        return Math.max(-1000, Math.min(1000, v));
    };
    state1d.chart.data.labels = state1d.x.map(xi => xi.toFixed(2));
    state1d.chart.data.datasets[0].data = state1d.u.map(clamp);
    state1d.chart.data.datasets[1].data = state1d.uExact;
    state1d.chart.update('none');
}

function updateStatus1D() {
    document.getElementById('time-display').textContent = state1d.t.toFixed(3);
    document.getElementById('r-display').textContent = state1d.r.toFixed(4);
    document.getElementById('steps-display').textContent = state1d.nSteps;
    document.getElementById('flops-display').textContent = state1d.flops.toExponential(2);

    const methodNames = { 'forward-euler': 'Forward Euler', 'backward-euler': 'Backward Euler', 'crank-nicolson': 'Crank-Nicolson' };
    document.getElementById('method-display').textContent = methodNames[sharedConfig.method];

    const isUnstable = sharedConfig.method === 'forward-euler' && state1d.r > 0.5;
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

    let sum = 0, maxErr = 0;
    for (let i = 0; i < state1d.u.length; i++) {
        const diff = Math.abs(state1d.u[i] - state1d.uExact[i]);
        if (isFinite(diff)) {
            sum += diff * diff;
            if (diff > maxErr) maxErr = diff;
        }
    }
    const l2Error = Math.sqrt(sum / state1d.u.length);
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

function animate1D() {
    if (!state1d.isPlaying) return;

    try {
        for (let i = 0; i < sharedConfig.animationSpeed; i++) {
            if (state1d.t >= sharedConfig.Tmax) {
                updateChart1D();
                updateStatus1D();
                stopAnimation1D();
                return;
            }
            timeStep1D();
            const maxVal = Math.max(...state1d.u.map(Math.abs));
            if (maxVal > 10 || state1d.u.some(v => !isFinite(v))) {
                if (!state1d.hasDiverged) {
                    state1d.hasDiverged = true;
                    updateStatus1D();
                }
            }
        }
        updateChart1D();
        updateStatus1D();
    } catch (e) {
        console.error('Animation error:', e);
        stopAnimation1D();
    }
    state1d.animationId = requestAnimationFrame(animate1D);
}

function startAnimation1D() {
    if (state1d.isPlaying) return;
    state1d.isPlaying = true;
    document.getElementById('play-btn').innerHTML = '<span id="play-icon">⏸</span> Pause';
    animate1D();
}

function stopAnimation1D() {
    state1d.isPlaying = false;
    document.getElementById('play-btn').innerHTML = '<span id="play-icon">▶</span> Play';
    if (state1d.animationId) {
        cancelAnimationFrame(state1d.animationId);
        state1d.animationId = null;
    }
}

// ========================================
// 2D Simulation Functions
// ========================================

function initSimulation2D() {
    state2d.dx = sharedConfig.L / (sharedConfig.Nx - 1);
    state2d.dy = state2d.dx;
    state2d.x = Array.from({ length: sharedConfig.Nx }, (_, i) => i * state2d.dx);
    state2d.y = state2d.x.slice();
    state2d.rx = sharedConfig.alpha * sharedConfig.dt / (state2d.dx * state2d.dx);
    state2d.ry = state2d.rx;
    state2d.t = 0;
    state2d.nSteps = 0;
    state2d.flops = 0;
    state2d.hasDiverged = false;
    state2d.u = getSineIC2D(state2d.x, state2d.y);
    state2d.uExact = getAnalyticalSolution2D(state2d.x, state2d.y, 0);
}

function stepNumerical2D() {
    const rx = state2d.rx;
    const ry = state2d.ry;
    const N = sharedConfig.Nx;

    switch (sharedConfig.method) {
        case 'forward-euler':
            state2d.u = forwardEuler2D(state2d.u, rx, ry);
            state2d.flops += 9 * (N - 2) * (N - 2);
            break;
        case 'backward-euler':
            state2d.u = backwardEuler2D_ADI(state2d.u, rx, ry);
            state2d.flops += 16 * N * N;
            break;
        case 'crank-nicolson':
            state2d.u = crankNicolson2D_ADI(state2d.u, rx, ry);
            state2d.flops += 24 * N * N;
            break;
    }
}

function timeStep2D() {
    stepNumerical2D();
    state2d.t += sharedConfig.dt;
    state2d.nSteps++;
    state2d.uExact = getAnalyticalSolution2D(state2d.x, state2d.y, state2d.t);
}

function updateSurfacePlot() {
    const layout = {
        title: {
            text: 'Temperature Distribution u(x,y)',
            font: { size: 18, family: 'Outfit, sans-serif' }
        },
        autosize: true,
        uirevision: 'true',
        margin: { l: 0, r: 0, b: 0, t: 40 },
        scene: {
            xaxis: {
                title: { text: 'x', font: { size: 14, family: 'Outfit' } },
                range: [0, 1],
                autorange: false,
                tickfont: { size: 12, family: 'Outfit' }
            },
            yaxis: {
                title: { text: 'y', font: { size: 14, family: 'Outfit' } },
                range: [0, 1],
                autorange: false,
                tickfont: { size: 12, family: 'Outfit' }
            },
            zaxis: {
                title: { text: 'u', font: { size: 14, family: 'Outfit' } },
                range: [-1, 1],
                autorange: false,
                tickfont: { size: 12, family: 'Outfit' }
            },
            aspectmode: 'manual',
            aspectratio: { x: 1.2, y: 1.2, z: 0.6 },
            camera: {
                eye: { x: 1.1, y: 1.1, z: 1.1 }
            }
        },
        paper_bgcolor: 'rgba(0,0,0,0)',
        plot_bgcolor: 'rgba(0,0,0,0)',
        font: { family: 'Outfit, sans-serif', color: '#ffffff', size: 14 }
    };

    Plotly.react('surface-plot', [{
        z: state2d.u,
        x: state2d.x,
        y: state2d.y,
        type: 'surface',
        colorscale: 'Viridis',
        showscale: false,
        contours: {
            z: { show: true, usecolormap: true, highlightcolor: "#42f462", project: { z: true } }
        }
    }], layout, { responsive: true, displayModeBar: false });
}

function initErrorChart() {
    if (state2d.errorChart) {
        state2d.errorChart.destroy();
    }

    Chart.defaults.font.family = "'Outfit', sans-serif";
    const ctx = document.getElementById('error-chart').getContext('2d');
    state2d.errorChart = new Chart(ctx, {
        type: 'line',
        data: {
            labels: [],
            datasets: [
                {
                    label: 'Analytical (y=0.5)',
                    data: [],
                    borderColor: '#ffa500',
                    borderDash: [5, 5],
                    borderWidth: 3,
                    pointRadius: 0,
                },
                {
                    label: 'Numerical (y=0.5)',
                    data: [],
                    borderColor: '#00e5ff',
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
                    title: { display: true, text: 'x (at y=0.5)', color: '#aaa', font: { size: 16, family: 'Outfit' } },
                    ticks: { color: '#888', font: { size: 14, family: 'Outfit' } },
                    grid: { color: '#333' }
                },
                y: {
                    title: { display: true, text: 'u(x, 0.5)', color: '#aaa', font: { size: 16, family: 'Outfit' } },
                    ticks: { color: '#888', font: { size: 14, family: 'Outfit' } },
                    grid: { color: '#333' },
                    min: -0.5,
                    max: 0.5
                }
            },
            plugins: {
                legend: { labels: { color: '#ccc', font: { size: 14, family: 'Outfit' } } }
            }
        }
    });
}

function updateErrorChart() {
    if (!state2d.errorChart) return;

    // Extract slice at y=0.5 (midpoint)
    const diagNumerical = [];
    const diagAnalytical = [];
    const jMid = Math.floor(sharedConfig.Nx / 2);

    for (let i = 0; i < sharedConfig.Nx; i++) {
        diagNumerical.push(state2d.u[i][jMid]);
        diagAnalytical.push(state2d.uExact[i][jMid]);
    }

    state2d.errorChart.data.labels = state2d.x.map(xi => xi.toFixed(2));
    state2d.errorChart.data.datasets[0].data = diagAnalytical;
    state2d.errorChart.data.datasets[1].data = diagNumerical;
    state2d.errorChart.update('none');
}

function updateStatus2D() {
    document.getElementById('time-display').textContent = state2d.t.toFixed(3);
    document.getElementById('r-display').textContent = `${state2d.rx.toFixed(4)}`;
    document.getElementById('steps-display').textContent = state2d.nSteps;
    document.getElementById('flops-display').textContent = state2d.flops.toExponential(2);

    const methodNames = { 'forward-euler': 'Forward Euler', 'backward-euler': 'Backward Euler (ADI)', 'crank-nicolson': 'Crank-Nicolson (ADI)' };
    document.getElementById('method-display').textContent = methodNames[sharedConfig.method];

    const isUnstable = sharedConfig.method === 'forward-euler' && (state2d.rx + state2d.ry) > 0.5;
    const stabilityStatus = document.getElementById('stability-status');
    if (isUnstable) {
        stabilityStatus.textContent = 'UNSTABLE';
        stabilityStatus.style.background = 'rgba(239, 68, 68, 0.2)';
        stabilityStatus.style.color = '#ef4444';
    } else {
        stabilityStatus.textContent = 'Stable';
        stabilityStatus.style.background = 'rgba(16, 185, 129, 0.2)';
        stabilityStatus.style.color = '#10b981';
    }

    let sum = 0, maxErr = 0, maxRelErr = 0;
    for (let i = 0; i < state2d.u.length; i++) {
        for (let j = 0; j < state2d.u[0].length; j++) {
            const diff = Math.abs(state2d.u[i][j] - state2d.uExact[i][j]);
            if (isFinite(diff)) {
                sum += diff * diff;
                if (diff > maxErr) maxErr = diff;
                if (Math.abs(state2d.uExact[i][j]) > 1e-10) {
                    const relErr = 100 * diff / Math.abs(state2d.uExact[i][j]);
                    if (relErr > maxRelErr) maxRelErr = relErr;
                }
            }
        }
    }
    const l2Error = Math.sqrt(sum / (state2d.u.length * state2d.u[0].length));

    document.getElementById('error-display').textContent = l2Error.toExponential(3);
    document.getElementById('max-error-display').textContent = maxErr.toExponential(3);
    document.getElementById('max-rel-error-display').textContent = maxRelErr.toFixed(2) + '%';
}

function animate2D() {
    if (!state2d.isPlaying) return;

    try {
        for (let i = 0; i < sharedConfig.animationSpeed; i++) {
            if (state2d.t >= sharedConfig.Tmax) {
                updateSurfacePlot();
                updateErrorChart();
                updateStatus2D();
                stopAnimation2D();
                return;
            }
            timeStep2D();
        }
        updateSurfacePlot();
        updateErrorChart();
        updateStatus2D();
    } catch (e) {
        console.error('Animation error:', e);
        stopAnimation2D();
    }
    state2d.animationId = requestAnimationFrame(animate2D);
}

function startAnimation2D() {
    if (state2d.isPlaying) return;
    state2d.isPlaying = true;
    document.getElementById('play-btn').innerHTML = '<span id="play-icon">⏸</span> Pause';
    animate2D();
}

function stopAnimation2D() {
    state2d.isPlaying = false;
    document.getElementById('play-btn').innerHTML = '<span id="play-icon">▶</span> Play';
    if (state2d.animationId) {
        cancelAnimationFrame(state2d.animationId);
        state2d.animationId = null;
    }
}

// ========================================
// Unified Controls
// ========================================

function togglePlayPause() {
    if (currentDimension === '1d') {
        if (state1d.isPlaying) stopAnimation1D();
        else startAnimation1D();
    } else {
        if (state2d.isPlaying) stopAnimation2D();
        else startAnimation2D();
    }
}

function reset() {
    if (currentDimension === '1d') {
        stopAnimation1D();
        initSimulation1D();
        updateChart1D();
        updateStatus1D();
    } else {
        stopAnimation2D();
        initSimulation2D();
        updateSurfacePlot();
        updateErrorChart();
        updateStatus2D();
    }
}

function updateMethodInfo() {
    const title = document.getElementById('method-title');
    const desc = document.getElementById('method-description');
    const info = {
        'forward-euler': {
            title: 'Forward Euler (Explicit)',
            desc: `First-order accurate, <strong>conditionally stable</strong>: \\(r = \\alpha \\Delta t / \\Delta x^2\\) must satisfy \\(r \\le 0.5\\) (1D) or \\(r_x + r_y \\le 0.5\\) (2D).`
        },
        'backward-euler': {
            title: 'Backward Euler (Implicit)',
            desc: `First-order accurate, <strong>unconditionally stable</strong>. Large time steps allowed but accuracy suffers.`
        },
        'crank-nicolson': {
            title: 'Crank-Nicolson (Implicit)',
            desc: `<strong>Second-order accurate</strong>, unconditionally stable. Best balance of stability and accuracy.`
        }
    };
    title.textContent = info[sharedConfig.method].title;
    desc.innerHTML = info[sharedConfig.method].desc;
    if (window.MathJax && window.MathJax.typesetPromise) {
        window.MathJax.typesetPromise([desc]).catch((err) => console.log('MathJax error:', err));
    }
}

// ========================================
// Method Selection
// ========================================

function updateDiscretizationDisplay() {
    const dim = currentDimension;
    const method = sharedConfig.method;

    // Hide all discretisation blocks
    document.querySelectorAll('.discretisation-block').forEach(block => {
        block.style.display = 'none';
    });

    // Show the selected method's discretisation for current dimension
    const blockId = `disc-${dim}-${method}`;
    const block = document.getElementById(blockId);
    if (block) {
        block.style.display = 'block';
    }
}

function setupEventListeners() {
    // Method buttons
    document.querySelectorAll('.method-btn').forEach(btn => {
        btn.addEventListener('click', () => {
            document.querySelectorAll('.method-btn').forEach(b => b.classList.remove('active'));
            btn.classList.add('active');
            sharedConfig.method = btn.dataset.method;
            updateMethodInfo();
            updateDiscretizationDisplay();
            reset();
        });
    });

    // dt slider
    const dtSlider = document.getElementById('dt-slider');
    dtSlider.addEventListener('input', () => {
        sharedConfig.dt = parseFloat(dtSlider.value);
        document.getElementById('dt-value').textContent = sharedConfig.dt.toExponential(2);
    });
    dtSlider.addEventListener('change', () => {
        reset();
    });

    // Nx slider
    const dxSlider = document.getElementById('dx-slider');
    dxSlider.addEventListener('input', () => {
        sharedConfig.Nx = parseInt(dxSlider.value);
        document.getElementById('dx-value').textContent = sharedConfig.Nx;
    });
    dxSlider.addEventListener('change', () => {
        reset();
    });

    // Tmax slider
    const tmaxSlider = document.getElementById('tmax-slider');
    tmaxSlider.addEventListener('input', () => {
        sharedConfig.Tmax = parseFloat(tmaxSlider.value);
        document.getElementById('tmax-value').textContent = sharedConfig.Tmax.toFixed(1);
    });

    // Speed slider
    const speedSlider = document.getElementById('speed-slider');
    speedSlider.addEventListener('input', () => {
        sharedConfig.animationSpeed = parseInt(speedSlider.value);
        document.getElementById('speed-value').textContent = `${sharedConfig.animationSpeed}x`;
    });

    // Play/Pause & Reset
    document.getElementById('play-btn').addEventListener('click', togglePlayPause);
    document.getElementById('reset-btn').addEventListener('click', reset);

    // Info tabs
    document.querySelectorAll('.info-tab').forEach(tab => {
        tab.addEventListener('click', () => {
            document.querySelectorAll('.info-tab').forEach(t => t.classList.remove('active'));
            document.querySelectorAll('.tab-content').forEach(c => c.classList.remove('active'));
            tab.classList.add('active');
            document.getElementById(`tab-${tab.dataset.tab}`).classList.add('active');
        });
    });
}

// ========================================
// Initialization
// ========================================

document.addEventListener('DOMContentLoaded', () => {
    try {
        // Initialize controls
        document.getElementById('dt-slider').value = sharedConfig.dt;
        document.getElementById('dx-slider').value = sharedConfig.Nx;
        document.getElementById('tmax-slider').value = sharedConfig.Tmax;
        document.getElementById('speed-slider').value = sharedConfig.animationSpeed;
        document.getElementById('dt-value').textContent = sharedConfig.dt.toExponential(2);
        document.getElementById('dx-value').textContent = sharedConfig.Nx;
        document.getElementById('tmax-value').textContent = sharedConfig.Tmax.toFixed(1);
        document.getElementById('speed-value').textContent = `${sharedConfig.animationSpeed}x`;

        // Initialize 1D (default)
        initSimulation1D();
        initChart1D();
        updateStatus1D();
        updateMethodInfo();
        setupEventListeners();

        // Render LaTeX in status panel
        if (window.MathJax && window.MathJax.typesetPromise) {
            window.MathJax.typesetPromise([document.querySelector('.status-panel')]).catch((err) => console.log('MathJax error:', err));
        }
    } catch (error) {
        console.error('Initialization error:', error);
        console.error('Stack:', error.stack);
    }
});
