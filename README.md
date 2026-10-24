<div align="center">

# DC Motor Control in MATLAB

### Speed and position control with P, PI, PD, PID and Model Predictive Control (MPC)

![MATLAB](https://img.shields.io/badge/MATLAB-R2019b%2B-orange)
![Control Systems](https://img.shields.io/badge/Control%20Systems-P%20%7C%20PI%20%7C%20PD%20%7C%20PID-blue)
![MPC](https://img.shields.io/badge/MPC-constrained%20QP-green)
![GUI](https://img.shields.io/badge/GUI-uifigure-lightgrey)

</div>

An interactive MATLAB application that simulates a linearized DC motor under classical and predictive control. You choose **Speed Control** or **Position Control**, pick one of five controllers, run the simulation, and read current, speed and position plots together with a panel of performance metrics.

> **Status.** The source has not been executed in MATLAB by the author of this README, and the GUI has not been tested there. See [Known Limitations](#known-limitations-and-warnings).

## Table of Contents

- [Key Features](#key-features)
- [Project Structure](#project-structure)
- [Requirements](#requirements)
- [Installation and Execution](#installation-and-execution)
- [Motor Model](#motor-model)
- [Transfer Functions](#transfer-functions)
- [Classical Controllers](#classical-controllers)
- [Model Predictive Control](#model-predictive-control)
- [Using the Application](#using-the-application)
- [Performance Metrics](#performance-metrics)
- [Reproducible Controller Comparison](#reproducible-controller-comparison)
- [Known Limitations and Warnings](#known-limitations-and-warnings)
- [License](#license)
- [Author](#author)

## Key Features

- Two operating modes: **Speed Control** (feedback on $\omega$) and **Position Control** (feedback on $\theta$).
- Five controllers: P, PI, PD, PID and MPC, selectable from a dropdown.
- Separate, independently editable gain tables for speed mode and position mode.
- Exact closed-loop transfer functions for the classical controllers, simulated with `lsim`.
- Receding-horizon MPC with zero-order-hold discretization; MPC uses only base MATLAB unless limits are enforced.
- Optional enforcement of voltage and current limits inside the MPC optimization (needs `quadprog`).
- Three live plots: armature current, angular position, angular speed.
- **Show All Controllers** overlays all five controllers with identical settings.
- Performance Metrics panel: steady-state error, ISE, peak current, control effort, RMS voltage, maximum voltage and current, and limit violations.
- Input validation, stability checks from closed-loop poles, and explicit error messages instead of silent fallbacks.

## Project Structure

```text
DC-Motor-Control-MATLAB-MPC/
├── DC_Motor_classical_and_MPC.m   % complete application (main function + local functions)
└── README.md
```

The whole application is one file. The main function is `DC_Motor_Speed_Control_GUI`. The other functions are local functions in the same file:

| Function | Purpose |
|---|---|
| `simulateController` | Closed-loop simulation of P, PI, PD, PID |
| `controllerPolynomials` | Builds $N_c(s)$ and $D_c(s)$ for each controller |
| `classicalVoltage` | Derives the applied armature voltage of the classical loop |
| `simulateMPC` | Receding-horizon MPC simulation |
| `buildPredictionMatrices` | Prediction matrices $\Phi$ and $\Gamma$ for any output |
| `computeMotorMetrics` | Single function that computes every performance metric |
| `metricsText` | Formats metrics for the panel |
| `validateInputs`, `validateGains`, `validateMPC`, `validateLimits` | Input checks |

## Requirements

| Item | Requirement |
|---|---|
| MATLAB | R2019b or newer (`uifigure`, `uigridlayout`, `uitable` with table data) |
| Control System Toolbox | **Required** (`tf`, `lsim`); the app stops with an error if it is missing |
| Optimization Toolbox | Only when MPC limit enforcement is ticked (`quadprog`) |

Without enforcement, MPC needs only base MATLAB (`expm`, `chol`), and the metrics need only `trapz`.

## Installation and Execution

1. Clone the repository (replace `<your-username>` with your GitHub account):

   ```bash
   git clone https://github.com/<your-username>/DC-Motor-Control-MATLAB-MPC.git
   cd DC-Motor-Control-MATLAB-MPC
   ```

2. Start MATLAB in that folder, or add the folder to the path:

   ```matlab
   cd('path/to/DC-Motor-Control-MATLAB-MPC')
   ```

3. Launch the application. MATLAB runs a function file by its **file name**, so use:

   ```matlab
   DC_Motor_classical_and_MPC
   ```

   You can also open the file in the editor and press **Run**.

4. The app opens and immediately runs the default P controller in speed mode.

> The function inside the file is called `DC_Motor_Speed_Control_GUI`. If you rename the file to `DC_Motor_Speed_Control_GUI.m`, call it by that name instead. File name and function name need not match, but the file name is what MATLAB uses.

## Motor Model

The motor is linear with constant field. The electrical, mechanical and kinematic equations are:

$$
V_a(t) = L_a \frac{di_a(t)}{dt} + R_a i_a(t) + K_e \omega(t)
$$

$$
J \frac{d\omega(t)}{dt} + B\omega(t) = K_t i_a(t)
$$

$$
\frac{d\theta(t)}{dt} = \omega(t)
$$

Initial conditions are $i_a(0) = 0$, $\omega(0) = 0$, $\theta(0) = 0$ in both modes.

### Default parameters

| Symbol | Meaning | Default | Unit |
|---|---|---|---|
| L_a | Armature inductance | 0.12 | H |
| R_a | Armature resistance | 7 | Ω |
| J | Rotor inertia | 2.12e-6 | kg·m² |
| B | Viscous friction | 6.04e-6 | N·m·s/rad |
| K_t | Torque constant | 0.0395 | N·m/A |
| K_e | Back-EMF constant | 0.0395 | V·s/rad |

> **Illustrative values.** The source marks $K_t$ and $K_e$ as placeholders. Replace all parameters with datasheet or identified values for a real motor. Nothing here has been validated against hardware.

## Transfer Functions

With $D(s)$ defined as

$$
D(s) = (L_a s + R_a)(J s + B) + K_t K_e
$$

the open-loop plant transfer functions from armature voltage are:

$$
P_i(s) = \frac{J s + B}{D(s)}, \qquad
P_\omega(s) = \frac{K_t}{D(s)}, \qquad
P_\theta(s) = \frac{K_t}{s\,D(s)}
$$

for armature current, angular speed and angular position respectively.

### Speed control

The error is $e = \omega_{ref} - \omega$ and the controller output is $V_a$. Writing the controller as $C(s) = N_c(s)/D_c(s)$, the closed-loop characteristic polynomial is

$$
D_{cl}(s) = D(s)\,D_c(s) + K_t N_c(s)
$$

and, since $D(s)$ cancels exactly, the closed-loop transfer functions from $\omega_{ref}$ are:

$$
\frac{I(s)}{\Omega_{ref}(s)} = \frac{(J s + B)\,N_c(s)}{D_{cl}(s)}, \qquad
\frac{\Omega(s)}{\Omega_{ref}(s)} = \frac{K_t N_c(s)}{D_{cl}(s)}, \qquad
\frac{\Theta(s)}{\Omega_{ref}(s)} = \frac{K_t N_c(s)}{s\,D_{cl}(s)}
$$

### Position control

The error is $e = \theta_{ref} - \theta$, with direct position feedback and no inner speed loop. The plant's extra integrator changes the characteristic polynomial:

$$
D_{cl}(s) = s\,D(s)\,D_c(s) + K_t N_c(s)
$$

$$
\frac{\Theta(s)}{\Theta_{ref}(s)} = \frac{K_t N_c(s)}{D_{cl}(s)}, \qquad
\frac{\Omega(s)}{\Theta_{ref}(s)} = \frac{s\,K_t N_c(s)}{D_{cl}(s)}, \qquad
\frac{I(s)}{\Theta_{ref}(s)} = \frac{s\,(J s + B)\,N_c(s)}{D_{cl}(s)}
$$

### Applied voltage

The applied armature voltage of the classical loop is derived from the same closed loop (function `classicalVoltage`):

$$
\frac{V_a(s)}{\Omega_{ref}(s)} = \frac{N_c(s)\,D(s)}{D_{cl}(s)}, \qquad
\frac{V_a(s)}{\Theta_{ref}(s)} = \frac{s\,D(s)\,N_c(s)}{D_{cl}(s)}
$$

Stability is decided from the roots of $D_{cl}(s)$: all must lie strictly in the left half-plane. The source notes that for position-mode P control with the default motor this requires $K_p$ below about 2.48.

## Classical Controllers

| Controller | $C(s)$ | Gains used |
|---|---|---|
| P | $K_p$ | Kp |
| PI | $K_p + K_i/s$ | Kp, Ki |
| PD | $K_p + K_d s$ | Kp, Kd |
| PID | $K_p + K_i/s + K_d s$ | Kp, Ki, Kd |

If $K_i = 0$ the integral term is dropped exactly (no cancelling pole and zero at the origin). The derivative is **ideal**: it acts on the full error and is not filtered.

### Default gains

| Controller | Speed mode (Kp / Ki / Kd) | Position mode (Kp / Ki / Kd) |
|---|---|---|
| P | 0.05 / 0 / 0 | 1 / 0 / 0 |
| PI | 0.05 / 0.20 / 0 | 1 / 2 / 0 |
| PD | 0.05 / 0 / 0.001 | 1 / 0 / 0.005 |
| PID | 0.05 / 0.20 / 0.001 | 1 / 2 / 0.005 |

Gains must be finite and non-negative. Defaults are illustrative and chosen to be stable with the default motor.

## Model Predictive Control

### State-space model

The MPC uses the same linearized motor with state $x = [i_a, \omega, \theta]^T$ and input $u = V_a$:

$$
\dot{x} = A x + B_c u, \qquad
A = \begin{bmatrix} -R_a/L_a & -K_e/L_a & 0 \\ K_t/J & -B/J & 0 \\ 0 & 1 & 0 \end{bmatrix}, \qquad
B_c = \begin{bmatrix} 1/L_a \\ 0 \\ 0 \end{bmatrix}
$$

It is discretized exactly with a zero-order hold and sampling interval $T_s$ using the matrix exponential:

$$
\exp\left(\begin{bmatrix} A & B_c \\ 0 & 0 \end{bmatrix} T_s\right) =
\begin{bmatrix} A_d & B_d \\ 0 & 1 \end{bmatrix}, \qquad
x_{k+1} = A_d x_k + B_d u_k
$$

The tracked output is $y = \omega$ in speed mode ($C = [0\ 1\ 0]$) and $y = \theta$ in position mode ($C = [0\ 0\ 1]$). The same code serves both modes.

### Horizons and sampling

| Setting | Value in the source |
|---|---|
| Prediction horizon $N_p$ | GUI field, default 20, integer 1 to 200 |
| Control horizon $N_c$ | $\min(N_p, 10)$, fixed internally, not a GUI field |
| Sampling interval $T_s$ | At most 5 ms; reduced for short runs so the grid has about 1001 samples; the grid ends exactly at $T$ |
| Run length limit | 20000 samples per run (100 s at 5 ms) |
| Inputs after $N_c$ | Held at the last free move $u_{k+N_c-1}$ |

### Cost function and QP

At each sample the controller minimizes, over $U = [u_k, \dots, u_{k+N_c-1}]^T$,

$$
J = Q \sum_{j=1}^{N_p} \left( y_{k+j|k} - r \right)^2 + R \sum_{j=0}^{N_c-1} u_{k+j|k}^2
$$

where $r$ is the reference ($\omega_{ref}$ or $\theta_{ref}$). With the prediction $y_{pred} = \Phi x_k + \Gamma U$ this is the quadratic program

$$
\min_U \; U^T H U + 2 f^T U, \qquad
H = Q\,\Gamma^T \Gamma + R\,I, \qquad
f = Q\,\Gamma^T \left( \Phi x_k - r\mathbf{1} \right)
$$

Only the first move $u_k$ is applied; the plant advances one sample and the problem is re-solved. Defaults are $Q = 1$ and $R = 1$, and both must be strictly positive.

### Constraints: monitored versus enforced

Voltage and current limits are **always monitored** for every controller. They are **enforced** only inside the MPC optimization, and only for the tick boxes that are ticked:

$$
V_{min} \le u_{k+j} \le V_{max}, \qquad
-I_{max} \le i_{a,k+j|k} \le I_{max}, \qquad j = 1,\dots,N_p
$$

- The voltage bound applies to the $N_c$ free inputs.
- The current bound is imposed on predictions at the sampling instants; the current between samples is not constrained.
- Without ticked boxes, MPC is unconstrained and solved by a Cholesky factorization.
- With a ticked box, `quadprog` is required. If the limits are infeasible, the run stops with an error and no controller is substituted.
- P, PI, PD and PID are linear transfer-function simulations without saturation: their voltage and current are reported, never limited.

## Using the Application

1. **Motor parameters (A).** Edit $L_a$, $R_a$, $J$, $B$, $K_t$, $K_e$ or keep the defaults.
2. **Control mode.** Choose Speed Control or Position Control. Only the matching reference field and gain table are shown.
3. **Reference (B).** Enter the reference speed [rad/s] or position [rad] and the simulation duration [s] (greater than 0, at most 1000).
4. **Controller (C).** Select P, PI, PD, PID or MPC.
5. **Gains (D).** Edit Kp, Ki, Kd in the table for the current mode. Invalid entries are reverted with a message.
6. **MPC (E).** Set $Q$, $R$ and $N_p$ (used only when MPC is selected).
7. **Limits (right column).** Set Vmin, Vmax and maximum $|i_a|$, and optionally tick **MPC enforces V limits** and **MPC enforces I limit**. Defaults are illustrative placeholders: −12 V, +12 V and 1.5 A.
8. **Run.** Click **Run Selected Controller**, or **Show All Controllers** to overlay all five.
9. **Reset.** **Reset Defaults** restores every default, unticks enforcement and reruns the default controller.

Behavior to know:

- Unstable closed loops are flagged with a warning; in the overlay they are listed but not drawn.
- Changing mode or controller clears the metrics panel; editing an input after a run marks the numbers "out of date" until you rerun.
- In the overlay, the metrics panel shows the controller selected in section C.

## Performance Metrics

All metrics are computed by one function, `computeMotorMetrics`, from the simulated signals. The tracking error is $e = \omega_{ref} - \omega$ in speed mode and $e = \theta_{ref} - \theta$ in position mode.

| Metric | Definition | Unit |
|---|---|---|
| Steady-state error | e_ss = ref − y_ss, signed and absolute | rad/s or rad |
| ISE | ∫ e² dt via `trapz` on the simulation time vector | (rad/s)²·s or rad²·s |
| Peak armature current | max abs(i_a) over all samples | A |
| Control effort J_u | ∫ V_a² dt | V²·s |
| RMS applied voltage | √(J_u / T) | V |
| Max observed abs(V_a), abs(i_a) | Maxima over the run | V, A |
| Voltage violations | Samples with V_a > Vmax or V_a < Vmin | samples |
| Current violations | Samples with abs(i_a) > Imax | samples |

Mathematically, the integral metrics and the steady-state error are:

$$
\mathrm{ISE} = \int_0^T e^2(t)\,dt, \qquad
J_u = \int_0^T V_a^2(t)\,dt, \qquad
V_{RMS} = \sqrt{\frac{J_u}{T}}, \qquad
e_{ss} = r - y_{ss}
$$

Implementation details:

- **Steady state.** $y_{ss}$ is the time-weighted mean over the final 10 % of the run. The value is labelled *settled* only if the output spread in that window is at most 0.1 % of the larger of the reference and the output magnitude. Otherwise it is labelled **NOT settled** and is an end-of-run error, not a true steady-state error.
- **Applied voltage.** For MPC it is the held input $u_k$, integrated exactly as a sum over each sample interval. For classical controllers it is derived from the closed loop as shown above.
- **Violations.** A sample exactly on a limit is not a violation. A relative tolerance of $10^{-6}$ absorbs solver round-off. The panel states per limit whether it is *ENFORCED by MPC* or *MONITORED only*.
- **Unavailable data.** Missing or invalid signals produce "n/a", never zero.

### Ideal derivative and infinite metrics

For PD and PID, the ideal derivative acts on a reference **step**, so the applied voltage contains a Dirac impulse of weight $K_d \cdot ref$ at $t = 0$. This is what makes the current jump by $K_d \cdot ref / L_a$. An impulse has infinite energy, so for PD and PID the panel reports **Inf** for $J_u$, $V_{RMS}$ and max abs($V_a$), flags a voltage violation, and gives the finite values for $t > 0$ separately in its note. No filter or clipping is added to hide this.

## Reproducible Controller Comparison

1. Start MATLAB and run `DC_Motor_classical_and_MPC`.
2. Click **Reset Defaults** so every run starts from the documented defaults.
3. Choose the mode and set one reference and one duration; leave them unchanged for all controllers.
4. Keep the limits at the same values for every run, with enforcement boxes unticked for a like-for-like comparison.
5. For each of P, PI, PD, PID and MPC: select it, click **Run Selected Controller**, and record the Performance Metrics panel.
6. Optionally click **Show All Controllers** to compare the plots on identical axes.
7. Repeat in the other mode. Record the gains and $Q$, $R$, $N_p$ used.
8. To study constraints, tick the MPC enforcement boxes and rerun MPC with the same limits; compare violation counts with the monitored classical runs.

## Known Limitations and Warnings

**Assumptions**

- Linear model with constant field; no saturation, friction nonlinearities, load torque, disturbances, sensor noise or measurement delay.
- All motor parameters and limits are illustrative until replaced with real values.
- A constant (step) reference is used for every run.
- MPC predicts with the same model that simulates the plant, so there is no model mismatch.

**Implementation limits**

- Classical controllers do not enforce limits; only MPC can, and only at sampling instants.
- Hard constraints can be infeasible; the run then fails with an error.
- Metrics are evaluated on the simulation time grid (1001 points for classical runs). Integrals can be off by a few percent and sampled peaks can sit slightly below the true peak when the time step is not small compared with the fastest closed-loop time constant; the panel warns in that case.
- MPC stability is judged from the unconstrained law; with voltage enforcement the response is bounded because the input is bounded.
- Duration is limited to 1000 s, and MPC to 20000 samples per run.

**Verification status**

- The app has not been run or GUI-tested in MATLAB by the author of this README.
- The formulas and metric logic were cross-checked by reading the source and with independent numerical checks outside MATLAB; this is not a substitute for running the app.
- No screenshots are included, and no results or performance claims are made beyond what the code computes.

## License

No license file is included in this repository yet. Until one is added, default copyright applies. Add a `LICENSE` file (for example MIT) if you want others to reuse the code.

## Author

AmarDeep Dwivedi, M.Tech (Communication, Signal Processing and Machine Learning), IIT Dharwad.

Built as a control-systems study project combining classical transfer-function design and model predictive control. Feedback and issues are welcome.
