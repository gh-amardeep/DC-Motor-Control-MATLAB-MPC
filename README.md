<div align="center">

# DC Motor Speed & Position Control

**Classical Controllers & Model Predictive Control (MPC) in MATLAB**

![MATLAB](https://img.shields.io/badge/MATLAB-Control%20Systems-orange?logo=mathworks)
![Controllers](https://img.shields.io/badge/Controllers-P%20%7C%20PI%20%7C%20PD%20%7C%20PID-blue)
![MPC](https://img.shields.io/badge/Advanced%20Control-MPC-purple)
![GUI](https://img.shields.io/badge/Interface-Interactive%20GUI-success)

An interactive MATLAB application for simulating and studying DC motor speed and angular-position control.

</div>

---

## Contents

- [Overview](#overview)
- [Features](#features)
- [Control Modes](#control-modes)
- [DC Motor Model](#dc-motor-model)
- [MPC Formulation](#mpc-formulation)
- [Performance Metrics](#performance-metrics)
- [Using the GUI](#using-the-gui)
- [Requirements](#requirements)
- [Run the Application](#run-the-application)
- [Limitations and Notes](#limitations-and-notes)

## Overview

This project simulates a linearized, constant-field DC motor using classical feedback controllers and Model Predictive Control. The GUI allows the user to change motor parameters, reference inputs, controller gains, and MPC tuning parameters, then inspect the resulting current, position, and speed plots.

> **Important:** The default \(K_t\) and \(K_e\) values are illustrative placeholders. Replace them with appropriate motor datasheet values before interpreting results as representative of a real motor.

## Features

- **Operating modes:** Speed Control and Position Control.
- **Controllers:** P, PI, PD, PID, and MPC.
- **Editable motor parameters:** armature inductance and resistance, inertia, viscous friction, torque constant, and back-EMF constant.
- **Separate classical-controller gain tables** for speed and position modes.
- **MPC tuning:** tracking weight \(Q\), control-effort weight \(R\), and prediction horizon \(N_p\).
- **Response plots:** armature current, angular position, and angular speed.
- **Performance metrics:** steady-state error, ISE, peak armature current, control effort, RMS voltage, observed signal maxima, and limit-violation counts when valid signals are available.
- **Voltage/current limits:** configurable monitoring limits and optional MPC constraint enforcement.
- **Simulation controls:** run the selected controller or show controller responses together.

## Control Modes

### Speed control

Tracks the requested angular speed \(\omega_{\mathrm{ref}}\).

\[
e_\omega(t)=\omega_{\mathrm{ref}}-\omega(t)
\]

### Position control

Tracks the requested angular position \(\theta_{\mathrm{ref}}\).

\[
e_\theta(t)=\theta_{\mathrm{ref}}-\theta(t)
\]

The tracking error used for performance metrics changes with the selected mode.

## DC Motor Model

The application uses the following linearized model:

\[
V_a(t)=L_a\frac{di_a(t)}{dt}+R_a i_a(t)+K_e\omega(t)
\]

\[
J\frac{d\omega(t)}{dt}+B\omega(t)=K_t i_a(t)
\]

\[
\frac{d\theta(t)}{dt}=\omega(t)
\]

### Symbols

| Symbol | Meaning | Unit |
|:--|:--|:--|
| \(V_a\) | Armature voltage | V |
| \(i_a\) | Armature current | A |
| \(\omega\) | Angular speed | rad/s |
| \(\theta\) | Angular position | rad |
| \(L_a\) | Armature inductance | H |
| \(R_a\) | Armature resistance | \(\Omega\) |
| \(J\) | Equivalent rotational inertia | kg·m² |
| \(B\) | Viscous friction coefficient | N·m·s/rad |
| \(K_t\) | Torque constant | N·m/A |
| \(K_e\) | Back-EMF constant | V·s/rad |

The transfer-function denominator is

\[
D(s)=(L_a s+R_a)(Js+B)+K_tK_e
\]

The current, speed, and position transfer functions are

\[
\frac{I_a(s)}{V_a(s)}=\frac{Js+B}{D(s)}
\]

\[
\frac{\Omega(s)}{V_a(s)}=\frac{K_t}{D(s)}
\]

\[
\frac{\Theta(s)}{V_a(s)}=\frac{K_t}{sD(s)}
\]

These equations describe the simulation model. A real motor may require additional effects such as load torque, saturation, nonlinear friction, and experimentally identified parameters.

## MPC Formulation

MPC predicts future motor behavior with a discrete-time state-space model:

\[
x_k=
\begin{bmatrix}
i_{a,k}\\
\omega_k\\
\theta_k
\end{bmatrix},
\qquad
x_{k+1}=A_d x_k+B_d u_k
\]

The controller optimizes tracking error and control effort over finite horizons. The objective can be expressed in general form as

\[
J=\sum_{j=1}^{N_p} Q e_{k+j}^{\,2}
+\sum_{j=0}^{N_c-1} R u_{k+j}^{\,2}
\]

where \(N_p\) is the prediction horizon, \(N_c\) is the control horizon, \(Q\) weights tracking error, and \(R\) weights control effort. The tracked output is speed or position depending on the selected mode.

The GUI exposes \(Q\), \(R\), and \(N_p\). In the current implementation, \(N_c=\min(N_p,10)\) is fixed internally, and the sampling time is no greater than 5 ms.

### Constraints

When MPC enforcement is enabled, the optimization can include voltage and predicted-current limits:

\[
V_{\min}\leq u_k\leq V_{\max},
\qquad |i_{a,k}|\leq I_{\max}
\]

**Monitoring is not the same as enforcement.** Classical controllers report violations but do not enforce limits unless saturation is implemented. MPC constraints apply at optimization sampling instants; inter-sample behavior should still be checked when limits are important.

## Performance Metrics

| Metric | Formula | What it indicates |
|:--|:--|:--|
| Steady-state error | \(e_{ss}=r-y_{ss}\) | Final tracking accuracy |
| Integral squared error (ISE) | \(\int_0^T e(t)^2\,dt\) | Accumulated tracking error |
| Peak current | \(\max_t |i_a(t)|\) | Maximum absolute armature current |
| Control effort | \(\int_0^T V_a(t)^2\,dt\) | Integrated squared voltage |
| RMS voltage | \(\sqrt{\frac{1}{T}\int_0^T V_a(t)^2\,dt}\) | Effective voltage magnitude |
| Signal maxima | \(\max |V_a|\), \(\max |i_a|\) | Maximum observed voltage and current |
| Limit violations | Count of samples outside configured bounds | Observed voltage/current limit exceedances |

For speed control, \(e(t)=\omega_{\mathrm{ref}}-\omega(t)\). For position control, \(e(t)=\theta_{\mathrm{ref}}-\theta(t)\). ISE units therefore depend on the selected mode.

### Why can a metric show `Inf`?

An ideal derivative controller responding to a step reference can produce a mathematical impulse at \(t=0\). The idealized model may consequently report `Inf` for control effort, RMS voltage, or maximum voltage. This does **not** mean a physical motor literally consumes infinite power. A filtered derivative or derivative-on-measurement can avoid derivative kick, but changes the controller implementation.

A missing or invalid metric should not be interpreted as zero. Check the status shown in the GUI.

## Using the GUI

1. Enter the motor parameters.
2. Choose **Speed Control** or **Position Control**.
3. Set the reference input and simulation duration.
4. Select P, PI, PD, PID, or MPC.
5. Adjust classical-controller gains or MPC parameters as appropriate.
6. Click **Run Selected Controller** to update plots and metrics.
7. Use **Show All Controllers** to inspect controller responses using the application's comparison workflow.
8. Set voltage/current limits. Enable MPC enforcement options only when you want MPC to apply those constraints.
9. Rerun after changing settings; metrics marked **out of date** should not be used for comparison.

## Requirements

- MATLAB.
- **Control System Toolbox** for the control-system functions used by the application.
- **Optimization Toolbox (`quadprog`)** when MPC voltage/current constraint enforcement is enabled, according to the implementation report.

The application has not been verified here by launching the GUI in MATLAB. Run it locally and confirm the required toolboxes are available.

## Run the Application

1. Download or clone the repository.
2. Place the `.m` file in a local working directory.
3. Open MATLAB and set the Current Folder to that directory.
4. Run:

   ```matlab
   DC_Motor_classical_and_MPC
   ```

The MATLAB file name must be compatible with its main function name. If you rename the file, check the function declaration too.

## Limitations and Notes

- Default motor constants and voltage/current limits are illustrative, not validated hardware specifications.
- MPC constraints are applied at the optimization sampling instants; check behavior between samples if strict limits matter.
- Classical-controller limits are monitored, not enforced, unless the implementation explicitly includes saturation.
- Ideal derivative action can cause derivative kick and non-finite voltage metrics.
- Metric accuracy depends on the simulation grid and numerical integration.
- Static checks are not a substitute for testing the actual GUI in MATLAB.
- This is a simulation and learning project, not a validated safety controller for physical hardware.

---

<div align="center">

**Built for learning, simulation, and controller-performance analysis.**

</div>
