function fig = DC_Motor_Speed_Control_GUI()
%DC_MOTOR_SPEED_CONTROL_GUI  Interactive P / PI / PD / PID speed control of a DC motor.
%
%   Run with:   DC_Motor_Speed_Control_GUI
%
%   Model (linearized, constant field):
%       Va = La*di/dt + Ra*i + Ke*w        J*dw/dt + B*w = Kt*i        dtheta/dt = w
%       D(s)     = (La*s + Ra)*(J*s + B) + Kt*Ke
%       P_i(s)   = (J*s + B)/D(s)          P_w(s) = Kt/D(s)          P_theta(s) = Kt/(s*D(s))
%
%   Unity negative feedback on SPEED, error e = w_ref - w, controller output = Va:
%       X(s)/Wref(s) = P_X(s)*C(s) / (1 + P_w(s)*C(s))
%
%   Writing C = Nc/Dc and clearing fractions, the common factor D(s) cancels
%   EXACTLY (done analytically, which is numerically safer than minreal):
%       Dcl(s)        = D(s)*Dc(s) + Kt*Nc(s)        (closed-loop characteristic polynomial)
%       I(s)/Wref(s)  = (J*s + B)*Nc(s) / Dcl(s)
%       W(s)/Wref(s)  =  Kt*Nc(s)       / Dcl(s)
%       Th(s)/Wref(s) =  Kt*Nc(s)       / (s*Dcl(s))
%   These equal P_X*C/(1 + P_w*C) identically.
%
%   MODEL PREDICTIVE CONTROL (MPC) - added controller option "MPC"
%   Same linearized motor, state x = [ia; w; theta], input u = Va:
%       dx/dt = A x + Bc u,  A = [-Ra/La -Ke/La 0; Kt/J -B/J 0; 0 1 0],  Bc = [1/La; 0; 0]
%   Exact zero-order-hold discretization with sampling interval Ts (computed from the
%   CURRENT motor parameters via expm):   x(k+1) = Ad x(k) + Bd u(k).
%   At every sampling instant k (receding horizon):
%       1. take the current simulated state x(k)
%       2. predict  w(k+j|k) = C Ad^j x(k) + sum_m C Ad^(j-1-m) Bd u(k+m),  j = 1..Np
%          (inputs after the control horizon Nc are held at u(k+Nc-1))
%       3. minimise  J = Q*sum_{j=1..Np}(w(k+j|k) - wref)^2 + R*sum_{j=0..Nc-1} u(k+j|k)^2
%          -> QP in U = [u(k)..u(k+Nc-1)]:  H = Q*G'*G + R*I,  f = Q*G'*(Phi*x(k) - wref)
%       4. apply ONLY u(k) = first element of the optimum, advance the plant one Ts
%       5. repeat at k+1
%   Default parameters: Q = 1, R = 1, Np = 20, Nc = min(Np,10) (fixed, internal),
%   Ts = 5 ms (reduced automatically only if needed to keep 1001 samples for short runs).
%   By default MPC applies NO constraints (the classical controllers are not limited
%   either, so the comparison is like-for-like). The two tick boxes in the "Voltage /
%   current limits" panel (right column) switch on ENFORCEMENT inside the MPC
%   optimisation: "MPC enforces V limits" adds  Vmin <= u(k+j) <= Vmax  and "MPC enforces
%   I limit" adds  -Imax <= ia(k+j|k) <= Imax, j = 1..Np (predicted at the sampling
%   instants). Either one makes the QP constrained and then requires quadprog
%   (Optimization Toolbox). Without enforcement the QP is solved with a Cholesky
%   factorisation of H (base MATLAB). MPC never silently falls back to another
%   controller: any failure (including infeasible limits) is reported in the Status panel.
%
%   CONTROL MODE (radio buttons "Speed Control" / "Position Control")
%   Speed mode (original): reference wref [rad/s], feedback = w, e = wref - w, output Va.
%   Position mode (added): reference thref [rad], feedback = THETA, e = thref - theta, output Va.
%   Direct position feedback (no inner speed loop). With C = Nc/Dc the exact closed loop is
%       Dcl(s)        = s*D(s)*Dc(s) + Kt*Nc(s)          (extra integrator of P_theta)
%       Th(s)/Thref(s) = Kt*Nc(s) / Dcl(s)
%       W(s)/Thref(s)  = s*Kt*Nc(s) / Dcl(s)
%       I(s)/Thref(s)  = s*(J*s + B)*Nc(s) / Dcl(s)
%   which equal P_X*C/(1 + P_theta*C). Stability = all roots of Dcl strictly in the left
%   half plane (for P control this needs Kp < about 2.48 with the default motor).
%   As in speed mode the derivative term acts on the full error, so a reference step gives
%   PD/PID an initial current jump of about Kd*thref/La (ideal derivative, no filtering).
%   MPC in position mode uses the SAME code: output matrix C = [0 0 1] and the cost
%       J = Q*sum_j (theta(k+j|k) - thref)^2 + R*sum_j u(k+j|k)^2.
%   The GUI fields Q, R, Np keep their values when the mode changes; Q weights the speed
%   error in speed mode and the position error in position mode (the label says which).
%   Initial conditions in both modes: ia(0) = 0, w(0) = 0, theta(0) = 0.
%   Position-mode defaults (illustrative, chosen to be stable with the default motor):
%   thref = 1 rad; Kp/Ki/Kd = P 1/0/0, PI 1/2/0, PD 1/0/0.005, PID 1/2/0.005.
%
%   PERFORMANCE METRICS (panel under the plots; recomputed after every simulation)
%   All metrics are computed by ONE function, computeMotorMetrics, from the signals the
%   simulation actually produced (time, w, theta, ia and the applied armature voltage Va).
%   Tracking error: e = wref - w (speed mode) or e = thref - theta (position mode).
%     1. Steady-state error  e_ss = ref - y_ss, y_ss = time-weighted mean of y over the
%        final DEF.ssWindow (10 %) of the run. The output counts as settled only if its
%        spread inside that window is <= DEF.settleTol (0.1 %) of max(|ref|,|y|); otherwise
%        the value is labelled "NOT settled" (end-of-run error, not a steady-state error).
%     2. ISE = trapz(t, e.^2) on the simulation time vector     [(rad/s)^2.s or rad^2.s]
%     3. Peak armature current  max|ia| over all samples                         [A]
%     4. Control effort  Ju = integral of Va^2 dt  and  V_RMS = sqrt(Ju/T)  [V^2.s, V]
%        Va is the voltage applied to the plant: MPC -> u(k) held over each Ts (exact
%        integral sum(u.^2.*dt)); classical -> derived from the SAME closed loop,
%        Va/Ref = Nc*D/Dcl (speed) or s*D*Nc/Dcl (position), simulated with lsim.
%        PD/PID: the ideal derivative acts on a reference STEP, so Va contains a Dirac
%        impulse Kd*ref*delta(t) (it is what makes ia jump by Kd*ref/La at t = 0).
%        Its energy is infinite, so Ju, V_RMS and max|Va| are reported as Inf and the
%        values for t > 0 are given separately. Nothing is invented or clipped.
%     5. Limit violations: number of samples with Va > Vmax, Va < Vmin or |ia| > Imax
%        (a sample exactly ON a limit is not a violation; a 1e-6 relative numerical
%        tolerance, DEF.limTolRel, absorbs solver round-off), plus max|Va| and max|ia|.
%   MONITORING vs ENFORCEMENT: the limits entered in the GUI are always MONITORED for every
%   controller. They are ENFORCED only by MPC, only for the tick boxes that are ticked, and
%   only at the sampling instants. P/PI/PD/PID are linear transfer-function simulations that
%   are NOT saturated: their voltage and current are reported, never limited.
%   Defaults of the limits (+-12 V, 1.5 A) are ILLUSTRATIVE placeholders, like Kt and Ke.
%
%   NOTE: Kt = Ke = 0.0395 are ILLUSTRATIVE PLACEHOLDER values. Replace them with
%   datasheet values for your motor.
%
%   Requires: MATLAB R2019b or newer (uifigure/uigridlayout/uitable with table data)
%             and the Control System Toolbox (tf, lsim) for P/PI/PD/PID.
%             MPC itself needs only base MATLAB (expm, chol); quadprog is needed
%             only if a voltage and/or current limit is ENFORCED in MPC (tick boxes).
%             The performance metrics need only base MATLAB (trapz).

% ------------------------------------------------------------------------
% Toolbox check
% ------------------------------------------------------------------------
if ~hasControlToolbox()
    msg = ['The Control System Toolbox is required (functions tf and lsim) ' ...
           'but is not available in this MATLAB installation or license.'];
    tmpFig = uifigure('Name','Control System Toolbox required');
    uialert(tmpFig, msg, 'Toolbox Missing', 'Icon', 'error');
    error('DCMotorGUI:NoToolbox', '%s', msg);
end

% ------------------------------------------------------------------------
% Defaults and constants
% ------------------------------------------------------------------------
DEF.La    = 0.12;       % H
DEF.Ra    = 7;          % Ohm
DEF.J     = 2.12e-6;    % kg.m^2
DEF.B     = 6.04e-6;    % N.m.s/rad
DEF.Kt    = 0.0395;     % N.m/A   (illustrative placeholder)
DEF.Ke    = 0.0395;     % V.s/rad (illustrative placeholder)
DEF.wref  = 10;         % rad/s
DEF.T     = 5;          % s
DEF.names = {'P','PI','PD','PID'};
DEF.gains = [0.05 0    0;       % P   : Kp Ki Kd
             0.05 0.20 0;       % PI
             0.05 0    0.001;   % PD
             0.05 0.20 0.001];  % PID

% Control-mode defaults
DEF.mode     = 'Speed Control';
DEF.thref    = 1;       % rad  (position-mode reference)
DEF.gainsPos = [1 0    0;       % P   : Kp Ki Kd   (position-mode gains, independent of speed mode)
                1 2    0;       % PI
                1 0    0.005;   % PD
                1 2    0.005];  % PID

% MPC defaults (documented in the header comment)
DEF.all      = [DEF.names {'MPC'}];  % dropdown entries
DEF.Q        = 1;       % tracking-error coefficient on (w - wref)^2   [1/(rad/s)^2]
DEF.R        = 1;       % control-effort coefficient on Va^2           [1/V^2]
DEF.Np       = 20;      % prediction horizon [samples] (20 x 5 ms = 100 ms)
DEF.NcMax    = 10;      % control horizon Nc = min(Np, NcMax)  (internal, not in GUI)
DEF.TsMax    = 5e-3;    % sampling interval [s]
DEF.MaxNp    = 200;     % practical upper limit for Np
DEF.MaxSteps = 20000;   % practical limit on MPC samples per run (=> T <= 100 s at 5 ms)

% Voltage / current limits (GUI fields). ILLUSTRATIVE placeholders - replace with your motor's
% ratings. They are always MONITORED; MPC enforces them only if the tick boxes are ticked.
DEF.Vmin     = -12;     % lower armature-voltage limit [V]
DEF.Vmax     =  12;     % upper armature-voltage limit [V]
DEF.Imax     = 1.5;     % maximum permitted |armature current| [A]

% Performance-metric settings (documented in the header comment)
DEF.ssWindow  = 0.10;   % steady-state window = final 10 % of the run
DEF.settleTol = 1e-3;   % window spread must be <= 0.1 % of max(|ref|,|y|) to call it settled
DEF.limTolRel = 1e-6;   % relative tolerance for "exceeds a limit" (solver / round-off)

NPTS = 1001;                    % number of time samples per simulation

COL.bg     = [0.94 0.94 0.94];  % window background
COL.panel  = [0.98 0.98 0.98];  % section background
COL.blue   = [0.00 0.32 0.60];  % blue section headers
COL.ok     = [0.00 0.40 0.10];
COL.warn   = [0.78 0.42 0.00];
COL.err    = [0.75 0.00 0.00];
COL.info   = [0.20 0.20 0.20];

% ------------------------------------------------------------------------
% Build the interface
% ------------------------------------------------------------------------
fig = uifigure('Name','DC Motor Speed / Position Control - P, PI, PD, PID, MPC', ...
               'Position',[40 10 1400 1040], 'Color',COL.bg);

mainGrid = uigridlayout(fig,[1 2]);
mainGrid.ColumnWidth   = {450,'1x'};
mainGrid.RowHeight     = {'1x'};
mainGrid.Padding       = [8 8 8 8];
mainGrid.ColumnSpacing = 8;
mainGrid.BackgroundColor = COL.bg;

% ---- Left column: simulation settings -----------------------------------
leftGrid = uigridlayout(mainGrid,[8 1]);
leftGrid.Layout.Column = 1;
leftGrid.RowHeight     = {180,74,104,76,246,104,120,'1x'};
leftGrid.Padding       = [0 0 0 0];
leftGrid.RowSpacing    = 4;
leftGrid.BackgroundColor = COL.bg;
leftGrid.Scrollable    = 'on';

% A. Motor parameters
motorGrid = createSection(leftGrid,1,'A. Motor parameters', ...
                          {26,26,26,40}, {110,'1x',110,'1x'}, COL);
addLabel(motorGrid,'La [H]',2,1);
edLa = addNumField(motorGrid,DEF.La,2,2,0,Inf,'off');
addLabel(motorGrid,'Ra [Ohm]',2,3);
edRa = addNumField(motorGrid,DEF.Ra,2,4,0,Inf,'off');
addLabel(motorGrid,'J [kg.m^2]',3,1);
edJ  = addNumField(motorGrid,DEF.J,3,2,0,Inf,'off');
addLabel(motorGrid,'B [N.m.s/rad]',3,3);
edB  = addNumField(motorGrid,DEF.B,3,4,0,Inf,'on');
addLabel(motorGrid,'Kt [N.m/A] *',4,1);
edKt = addNumField(motorGrid,DEF.Kt,4,2,0,Inf,'off');
addLabel(motorGrid,'Ke [V.s/rad] *',4,3);
edKe = addNumField(motorGrid,DEF.Ke,4,4,0,Inf,'off');
noteMotor = addLabel(motorGrid,['* Kt and Ke are illustrative placeholder values - ' ...
                                'replace them with your motor''s datasheet values.'],5,1);
noteMotor.Layout.Column = [1 4];
noteMotor.WordWrap = 'on';
noteMotor.FontAngle = 'italic';
noteMotor.FontSize = 11;
noteMotor.FontColor = COL.err;

% Control Mode (two mutually exclusive radio buttons)
modeGrid = createSection(leftGrid,2,'Control Mode',{28},{'1x'},COL);
bgMode = uibuttongroup(modeGrid,'BorderType','none','BackgroundColor',COL.panel, ...
                       'SelectionChangedFcn',@onModeChanged);
bgMode.Layout.Row = 2;  bgMode.Layout.Column = 1;
rbSpeed = uiradiobutton(bgMode,'Text','Speed Control','FontSize',12, ...
                        'Position',[10 3 190 22]);
rbPos   = uiradiobutton(bgMode,'Text','Position Control','FontSize',12, ...
                        'Position',[220 3 190 22]);
bgMode.SelectedObject = rbSpeed;

% B. Reference input (only the field for the selected mode is visible)
[refGrid,hdrRef] = createSection(leftGrid,3,'B. Reference input (speed)', ...
                        {26,26}, {170,'1x'}, COL);
lblRefSpd = addLabel(refGrid,'Reference speed [rad/s]',2,1);
edWref = addNumField(refGrid,DEF.wref,2,2,-Inf,Inf,'on');
lblRefPos = addLabel(refGrid,'Reference position [rad]',2,1);
edThref = addNumField(refGrid,DEF.thref,2,2,-Inf,Inf,'on');
lblRefPos.Visible = 'off';  edThref.Visible = 'off';
addLabel(refGrid,'Simulation duration [s]',3,1);
edTsim = addNumField(refGrid,DEF.T,3,2,0,1000,'off');

% C. Controller selection
ctrlGrid = createSection(leftGrid,4,'C. Controller selection', ...
                         {28}, {170,'1x'}, COL);
addLabel(ctrlGrid,'Selected controller',2,1);
ddCtrl = uidropdown(ctrlGrid,'Items',DEF.all,'Value','P','FontSize',12, ...
                    'ValueChangedFcn',@onControllerChanged);
ddCtrl.Layout.Row = 2;  ddCtrl.Layout.Column = 2;

% D. Independent controller gain table
[tblGrid,hdrGain] = createSection(leftGrid,5,'D. Controller gains - speed mode (editable)', ...
                        {150,48}, {'1x'}, COL);
tblGains = uitable(tblGrid, ...
                   'Data',makeGainTable(DEF), ...
                   'ColumnEditable',[false true true true], ...
                   'ColumnWidth',{105,100,100,100}, ...
                   'RowName',{}, ...
                   'ColumnSortable',false, ...
                   'FontSize',13, ...
                   'CellEditCallback',@onTableEdited);
tblGains.Layout.Row = 2;
% Separate gain table for position mode, in the same cell (only one is visible).
tblGainsPos = uitable(tblGrid, ...
                   'Data',makeGainTable(DEF,DEF.gainsPos), ...
                   'ColumnEditable',[false true true true], ...
                   'ColumnWidth',{105,100,100,100}, ...
                   'RowName',{}, ...
                   'ColumnSortable',false, ...
                   'FontSize',13, ...
                   'Visible','off', ...
                   'CellEditCallback',@onTableEdited);
tblGainsPos.Layout.Row = 2;
noteTbl = addLabel(tblGrid,['Double-click a cell to edit; press Enter or click another cell ' ...
                            'to commit. Each controller keeps its own gains, separately for speed and position mode. Unused gains are ' ...
                            'ignored: P uses Kp; PI uses Kp, Ki; PD uses Kp, Kd; PID uses all.'],3,1);
noteTbl.WordWrap = 'on';
noteTbl.FontSize = 11;
noteTbl.FontColor = [0.3 0.3 0.3];

% E. MPC parameters (read at every simulation run)
mpcGrid = createSection(leftGrid,6,'E. MPC parameters (used when MPC is selected)', ...
                        {26,26}, {110,'1x',110,'1x'}, COL);
lblQ = addLabel(mpcGrid,'Q (speed)',2,1);
edQ  = addNumField(mpcGrid,DEF.Q,2,2,0,Inf,'off');
addLabel(mpcGrid,'R (effort)',2,3);
edR  = addNumField(mpcGrid,DEF.R,2,4,0,Inf,'off');
addLabel(mpcGrid,'Np [samples]',3,1);
edNp = addNumField(mpcGrid,DEF.Np,3,2,1,DEF.MaxNp,'on');
edNp.RoundFractionalValues = 'on';
edNp.ValueDisplayFormat    = '%d';
noteMpc = addLabel(mpcGrid,sprintf('Nc = min(Np,%d), Ts <= %g ms (fixed)',DEF.NcMax,DEF.TsMax*1e3),3,3);
noteMpc.Layout.Column = [3 4];
noteMpc.FontSize = 11;
noteMpc.FontColor = [0.3 0.3 0.3];

% Actions
btnGrid = createSection(leftGrid,7,'Actions',{34,34},{'1x','1x'},COL);
btnRun = uibutton(btnGrid,'push','Text','Run Selected Controller', ...
                  'FontSize',13,'FontWeight','bold', ...
                  'BackgroundColor',COL.blue,'FontColor',[1 1 1], ...
                  'ButtonPushedFcn',@(~,~) runSelected(''));
btnRun.Layout.Row = 2;  btnRun.Layout.Column = [1 2];
btnAll = uibutton(btnGrid,'push','Text','Show All Controllers', ...
                  'FontSize',12,'ButtonPushedFcn',@(~,~) runAll());
btnAll.Layout.Row = 3;  btnAll.Layout.Column = 1;
btnReset = uibutton(btnGrid,'push','Text','Reset Defaults', ...
                    'FontSize',12,'ButtonPushedFcn',@onReset);
btnReset.Layout.Row = 3;  btnReset.Layout.Column = 2;

% Status
statGrid = createSection(leftGrid,8,'Status',{'1x'},{'1x'},COL);
statusLabel = uilabel(statGrid,'Text','Ready.','WordWrap','on', ...
                      'VerticalAlignment','top','FontSize',12,'FontColor',COL.info);
statusLabel.Layout.Row = 2;

% ---- Right column: three vertically stacked plots + metrics / limits panels ----
rightGrid = uigridlayout(mainGrid,[4 1]);
rightGrid.Layout.Column = 2;
rightGrid.RowHeight     = {'1x','1x','1x',276};
rightGrid.ColumnWidth   = {'1x'};
rightGrid.Padding       = [0 0 0 0];
rightGrid.RowSpacing    = 8;
rightGrid.BackgroundColor = COL.bg;

axI  = uiaxes(rightGrid);  axI.Layout.Row  = 1;
axTh = uiaxes(rightGrid);  axTh.Layout.Row = 2;
axW  = uiaxes(rightGrid);  axW.Layout.Row  = 3;
prepareAxes(axI, 'Armature Current i(t)',   'Current i(t) [A]',            DEF.T,[0 0 0]);
prepareAxes(axTh,'Angular Position \theta(t)','Position \theta(t) [rad]', DEF.T,[0 0 0]);
prepareAxes(axW, 'Angular Speed \omega(t)', 'Speed \omega(t) [rad/s]',    DEF.T,[0 0 0]);

% ---- Performance Metrics panel (left) and voltage / current limits (right) ----
metricsShown = false;   % true only while the panel shows the result of the CURRENT inputs
metCont = uigridlayout(rightGrid,[1 2]);
metCont.Layout.Row      = 4;
metCont.ColumnWidth     = {'1x',245};
metCont.RowHeight       = {'1x'};
metCont.Padding         = [0 0 0 0];
metCont.ColumnSpacing   = 8;
metCont.BackgroundColor = COL.bg;

mg = createSection(metCont,1,'Performance Metrics', ...
                   {22,22,22,22,22,22,22,56},{135,'1x',135,'1x'},COL);
mg.Parent.Layout.Column = 1;
mg.RowSpacing = 3;
addLabel(mg,'Result for',2,1);
MH.valCase = addLabel(mg,'',2,2);   MH.valCase.Layout.Column = [2 4];
MH.valCase.FontWeight = 'bold';
addLabel(mg,'Steady-state error',3,1);
MH.valSS   = addLabel(mg,'',3,2);   MH.valSS.Layout.Column = [2 4];
addLabel(mg,'ISE',4,1);                    MH.valISE  = addLabel(mg,'',4,2);
addLabel(mg,'Peak armature current',4,3);  MH.valIpk  = addLabel(mg,'',4,4);
addLabel(mg,'Control effort Ju',5,1);      MH.valJu   = addLabel(mg,'',5,2);
addLabel(mg,'RMS applied voltage',5,3);    MH.valVrms = addLabel(mg,'',5,4);
addLabel(mg,'Max |Va| observed',6,1);      MH.valVmax = addLabel(mg,'',6,2);
addLabel(mg,'Max |ia| observed',6,3);      MH.valImax = addLabel(mg,'',6,4);
addLabel(mg,'Voltage violations',7,1);
MH.valVviol = addLabel(mg,'',7,2);  MH.valVviol.Layout.Column = [2 4];
addLabel(mg,'Current violations',8,1);
MH.valIviol = addLabel(mg,'',8,2);  MH.valIviol.Layout.Column = [2 4];
MH.valNote = addLabel(mg,'',9,1);   MH.valNote.Layout.Column = [1 4];
MH.valNote.WordWrap = 'on';  MH.valNote.FontSize = 11;  MH.valNote.VerticalAlignment = 'top';

lg = createSection(metCont,1,'Voltage / current limits', ...
                   {22,22,22,22,22,46},{100,'1x'},COL);
lg.Parent.Layout.Column = 2;
lg.RowSpacing = 3;
addLabel(lg,'Vmin [V]',2,1);       edVmin = addNumField(lg,DEF.Vmin,2,2,-Inf,Inf,'on');
addLabel(lg,'Vmax [V]',3,1);       edVmax = addNumField(lg,DEF.Vmax,3,2,-Inf,Inf,'on');
addLabel(lg,'Max |ia| [A]',4,1);   edImax = addNumField(lg,DEF.Imax,4,2,0,Inf,'off');
cbEnfV = uicheckbox(lg,'Text','MPC enforces V limits','Value',false,'FontSize',11, ...
                    'Tooltip',['Unticked: voltage limits are only MONITORED. Ticked: MPC adds ' ...
                               'Vmin <= u <= Vmax to its optimisation (needs quadprog).']);
cbEnfV.Layout.Row = 5;  cbEnfV.Layout.Column = [1 2];
cbEnfI = uicheckbox(lg,'Text','MPC enforces I limit','Value',false,'FontSize',11, ...
                    'Tooltip',['Unticked: the current limit is only MONITORED. Ticked: MPC adds ' ...
                               '|ia| <= Imax at the predicted sampling instants (needs quadprog).']);
cbEnfI.Layout.Row = 6;  cbEnfI.Layout.Column = [1 2];
noteLim = addLabel(lg,['Illustrative placeholders: enter your motor''s ratings. Limits are always ' ...
                       'monitored; only ticked boxes are enforced (MPC only).'],7,1);
noteLim.Layout.Column = [1 2];  noteLim.WordWrap = 'on';  noteLim.FontSize = 10;
noteLim.FontAngle = 'italic';   noteLim.VerticalAlignment = 'top';
noteLim.FontColor = [0.3 0.3 0.3];

% Any input change makes the numbers on screen out of date -> flag them (never silently stale).
staleList = {edLa,edRa,edJ,edB,edKt,edKe,edWref,edThref,edTsim,edQ,edR,edNp, ...
             edVmin,edVmax,edImax,cbEnfV,cbEnfI};
for kStale = 1:numel(staleList)
    staleList{kStale}.ValueChangedFcn = @(~,~) markMetricsStale();
end

highlightSelectedRow();
applyModeUI();
resetMetricsPanel('Run a controller to see its performance metrics.');
drawnow;

% ---- Initial default simulation so all three axes are populated ---------
runSelected('');

% ========================================================================
%                          NESTED CALLBACKS / HELPERS
%   (nested so they share the UI handles defined above)
% ========================================================================

    function onModeChanged(~,~)
        % Switching mode only changes which reference field / gain table is shown and how Q
        % is labelled. Motor parameters, gains and MPC values are NOT modified.
        applyModeUI();
        resetMetricsPanel();      % definitions / units of the new mode, no stale numbers
        setStatus(sprintf(['%s selected. Motor parameters, controller gains and MPC settings are ' ...
                           'unchanged. Click "Run Selected Controller" or "Show All Controllers".'], ...
                          curMode()),'info');
    end

    function applyModeUI()
        isPos = strcmp(curMode(),'Position Control');
        lblRefSpd.Visible = onoff(~isPos);   edWref.Visible  = onoff(~isPos);
        lblRefPos.Visible = onoff(isPos);    edThref.Visible = onoff(isPos);
        tblGains.Visible    = onoff(~isPos);
        tblGainsPos.Visible = onoff(isPos);
        if isPos
            hdrRef.Text  = '  B. Reference input (position)';
            hdrGain.Text = '  D. Controller gains - position mode (editable)';
            lblQ.Text    = 'Q (position)';
        else
            hdrRef.Text  = '  B. Reference input (speed)';
            hdrGain.Text = '  D. Controller gains - speed mode (editable)';
            lblQ.Text    = 'Q (speed)';
        end
    end

    function md = curMode()
        md = bgMode.SelectedObject.Text;
    end

    function k = modeKey()
        if strcmp(curMode(),'Position Control'), k = 'position'; else, k = 'speed'; end
    end

    function t = modeTag()
        % '' in speed mode so the original speed-mode messages are unchanged.
        if strcmp(modeKey(),'position'), t = ' (position control)'; else, t = ''; end
    end

    function tb = activeTable()
        if strcmp(modeKey(),'position'), tb = tblGainsPos; else, tb = tblGains; end
    end

    function onControllerChanged(~,~)
        % Changing the selection NEVER touches the gain table.
        highlightSelectedRow();
        resetMetricsPanel();      % the metrics on screen belonged to the previous controller
        if strcmp(ddCtrl.Value,'MPC')
            setStatus(['MPC selected. Q, R and Np (section E) are read when you run; the PID gain ' ...
                       'table is unchanged. Click "Run Selected Controller" to simulate.'],'info');
        else
            setStatus(sprintf(['%s selected. Gains in the table are unchanged. ' ...
                               'Click "Run Selected Controller" to simulate.'], ddCtrl.Value),'info');
        end
    end

    function onTableEdited(src,evt)
        % Validate the single edited cell; revert it if invalid.
        editRow = evt.Indices(1);
        editCol = evt.Indices(2);
        newVal  = evt.NewData;
        isGood  = isnumeric(newVal) && isscalar(newVal) && isreal(newVal) && ...
                  isfinite(newVal) && newVal >= 0;
        if ~isGood
            curTbl = src.Data;
            curTbl{editRow,editCol} = evt.PreviousData;   % revert only this cell
            src.Data = curTbl;
            colNames = src.Data.Properties.VariableNames;
            msgTxt = sprintf(['Invalid entry for %s of the %s controller. Gains must be ' ...
                              'finite, non-negative numbers. The previous value was restored.'], ...
                             colNames{editCol}, curTbl.Controller{editRow});
            setStatus(msgTxt,'error');
            uialert(fig,msgTxt,'Invalid Gain','Icon','warning');
        else
            markMetricsStale();   % a valid gain change makes the shown metrics out of date
        end
        drawnow;
    end

    function runSelected(prefixTxt)
        % Simulate the controller chosen in the dropdown and update the plots.
        if nargin < 1, prefixTxt = ''; end
        try
            setStatus('Simulating...','info');
            drawnow;
            mot  = readMotorAndReference();
            lim  = readLimits();
            nm   = ddCtrl.Value;
            if strcmp(nm,'MPC')
                mpc = readMPC();
                res = simulateMPC(mot,mpc,NPTS);
                plotSingle(nm,mot,res);
                showMetrics(nm,mot,res,lim,'');
                if res.stable
                    setStatus(sprintf(['%sMPC simulated (%s). %s'], prefixTxt, mpcString(res), ...
                        mpcStabilityText(res)),'ok');
                else
                    wTxt = sprintf(['MPC (%s) gives an UNSTABLE closed loop (largest discrete-time pole ' ...
                        'magnitude %.4f >= 1). The plotted response is NOT a successful control result. ' ...
                        'Change Q, R or Np.'], mpcString(res), res.rho);
                    setStatus(['WARNING: ' wTxt],'warn');
                    uialert(fig,wTxt,'Unstable Closed Loop','Icon','warning');
                end
            else
                [kp,ki,kd] = readGains(nm);
                validateGains(nm,kp,ki,kd);
                res  = simulateController(nm,mot,kp,ki,kd,NPTS);
                plotSingle(nm,mot,res);
                showMetrics(nm,mot,res,lim,'');
                gTxt = gainString(nm,kp,ki,kd);
                if res.stable
                    setStatus(sprintf(['%s%s controller%s simulated (%s). Closed loop is stable ' ...
                        '(largest pole real part %.4g 1/s).'], prefixTxt, nm, modeTag(), gTxt, res.maxRe),'ok');
                else
                    wTxt = sprintf(['%s controller%s (%s) gives an UNSTABLE closed loop ' ...
                        '(largest pole real part %.4g 1/s). The plotted response is diverging and is ' ...
                        'NOT a successful control result. Reduce the gains.'], nm, modeTag(), gTxt, res.maxRe);
                    setStatus(['WARNING: ' wTxt],'warn');
                    uialert(fig,wTxt,'Unstable Closed Loop','Icon','warning');
                end
            end
        catch ME
            resetMetricsPanel('Simulation failed - there is no valid result, so no metrics are shown.');
            setStatus(['Error: ' ME.message],'error');
            uialert(fig,ME.message,'Simulation Error','Icon','error');
        end
        drawnow;
    end

    function runAll()
        % Overlay P, PI, PD, PID and MPC with identical motor/reference/duration.
        try
            setStatus('Simulating all controllers...','info');
            drawnow;
            mot = readMotorAndReference();
            lim = readLimits();
            nCtl = numel(DEF.all);
            results = cell(1,nCtl);
            failList = {};
            for idx = 1:nCtl
                nm = DEF.all{idx};
                if strcmp(nm,'MPC')
                    % An MPC failure is reported but does not hide the classical results.
                    try
                        results{idx} = simulateMPC(mot,readMPC(),NPTS);
                    catch ME2
                        failList{end+1} = sprintf('MPC: %s',ME2.message); %#ok<AGROW>
                    end
                    continue;
                end
                try
                    [kp,ki,kd] = readGains(nm);
                    validateGains(nm,kp,ki,kd);
                    results{idx} = simulateController(nm,mot,kp,ki,kd,NPTS);
                catch ME2
                    error('DCMotorGUI:SimFailed','%s controller: %s',nm,ME2.message);
                end
            end
            sfx = modeTag();
            prepareAxes(axI, ['Armature Current i(t) - All Controllers' sfx], 'Current i(t) [A]',            mot.T,[0 0 0]);
            prepareAxes(axTh,['Angular Position \theta(t) - All Controllers' sfx],'Position \theta(t) [rad]', mot.T,[0 0 0]);
            prepareAxes(axW, ['Angular Speed \omega(t) - All Controllers' sfx], 'Speed \omega(t) [rad/s]',   mot.T,[0 0 0]);
            unstableList = {};
            for idx = 1:nCtl
                nm = DEF.all{idx};
                r  = results{idx};
                if isempty(r), continue; end   % MPC failed: reported below
                if ~r.stable
                    unstableList{end+1} = nm; %#ok<AGROW>
                    continue;       % unstable responses are not drawn (they would hide the others)
                end
                cc = ctrlColor(nm);
                plot(axI, r.t,r.i, 'Color',cc,'LineWidth',1.8,'DisplayName',nm);
                plot(axTh,r.t,r.th,'Color',cc,'LineWidth',1.8,'DisplayName',nm);
                plot(axW, r.t,r.w, 'Color',cc,'LineWidth',1.8,'DisplayName',nm);
            end
            if strcmp(mot.mode,'position')
                plot(axTh,[0 mot.T],[mot.thref mot.thref],'--','Color',[0.15 0.15 0.15], ...
                     'LineWidth',1.4,'DisplayName','Reference position');
            else
                plot(axW,[0 mot.T],[mot.wref mot.wref],'--','Color',[0.15 0.15 0.15], ...
                     'LineWidth',1.4,'DisplayName','Reference speed');
            end
            legend(axI, 'Location','best');
            legend(axTh,'Location','best');
            legend(axW, 'Location','best');
            % Metrics are per controller: show those of the controller selected in section C.
            selIdx = find(strcmp(DEF.all,ddCtrl.Value),1);
            if ~isempty(selIdx) && ~isempty(results{selIdx})
                showMetrics(ddCtrl.Value,mot,results{selIdx},lim, ...
                            'from the all-controllers run (controller selected in section C)');
            else
                resetMetricsPanel(sprintf(['%s did not simulate in the all-controllers run - ' ...
                                           'no metrics available.'],ddCtrl.Value));
            end
            if isempty(unstableList) && isempty(failList)
                setStatus(sprintf('All five controllers simulated%s; all closed loops are stable.%s', ...
                                  modeTag(),mpcEnforceNote()),'ok');
            else
                wTxt = '';
                if ~isempty(unstableList)
                    wTxt = sprintf(['Unstable closed loop for: %s. These controllers are NOT plotted ' ...
                        'because their responses diverge; this is not a successful control result. ' ...
                        'Reduce their gains (or change Q, R, Np for MPC) and run again. '], ...
                        strjoin(unstableList,', '));
                end
                if ~isempty(failList)
                    wTxt = [wTxt 'Simulation failed and is NOT plotted - ' strjoin(failList,'; ')];
                end
                setStatus(['WARNING: ' wTxt],'warn');
                uialert(fig,wTxt,'Controller Warning','Icon','warning');
            end
        catch ME
            resetMetricsPanel('Simulation failed - there is no valid result, so no metrics are shown.');
            setStatus(['Error: ' ME.message],'error');
            uialert(fig,ME.message,'Simulation Error','Icon','error');
        end
        drawnow;
    end

    function onReset(~,~)
        % Restore every default (motor, reference, duration, ALL gains, MPC Q/R/Np, selection = P).
        edLa.Value   = DEF.La;
        edRa.Value   = DEF.Ra;
        edJ.Value    = DEF.J;
        edB.Value    = DEF.B;
        edKt.Value   = DEF.Kt;
        edKe.Value   = DEF.Ke;
        edWref.Value = DEF.wref;
        edThref.Value = DEF.thref;
        edTsim.Value = DEF.T;
        tblGains.Data    = makeGainTable(DEF);
        tblGainsPos.Data = makeGainTable(DEF,DEF.gainsPos);
        bgMode.SelectedObject = rbSpeed;            % default mode = Speed Control
        applyModeUI();
        edQ.Value    = DEF.Q;
        edR.Value    = DEF.R;
        edNp.Value   = DEF.Np;
        edVmin.Value = DEF.Vmin;
        edVmax.Value = DEF.Vmax;
        edImax.Value = DEF.Imax;
        cbEnfV.Value = false;                       % limits back to "monitor only"
        cbEnfI.Value = false;
        ddCtrl.Value  = 'P';
        highlightSelectedRow();
        resetMetricsPanel();
        drawnow;
        runSelected('Defaults restored. ');
    end

    function mot = readMotorAndReference()
        % Always read the CURRENT values of the edit fields.
        mot.La   = edLa.Value;
        mot.Ra   = edRa.Value;
        mot.J    = edJ.Value;
        mot.B    = edB.Value;
        mot.Kt   = edKt.Value;
        mot.Ke   = edKe.Value;
        mot.wref  = edWref.Value;
        mot.thref = edThref.Value;
        mot.mode  = modeKey();                  % 'speed' or 'position'
        if strcmp(mot.mode,'position'), mot.ref = mot.thref; else, mot.ref = mot.wref; end
        mot.T    = edTsim.Value;
        validateInputs(mot);
    end

    function [kp,ki,kd] = readGains(nm)
        % Read the LATEST committed gains of controller 'nm' from the current mode's table.
        tb  = activeTable();            % gain table of the CURRENT mode
        tbl = tb.Data;
        if ~istable(tbl)
            error('DCMotorGUI:BadTable','The gain table has an unexpected data type.');
        end
        row = find(strcmp(tbl.Controller,nm),1);
        if isempty(row)
            error('DCMotorGUI:BadTable','Controller "%s" was not found in the gain table.',nm);
        end
        kp = tbl.Kp(row);
        ki = tbl.Ki(row);
        kd = tbl.Kd(row);
    end

    function mpc = readMPC()
        % Always read the CURRENT MPC edit fields (Q, R, Np), then validate them.
        mpc.Q     = edQ.Value;
        mpc.R     = edR.Value;
        mpc.Np    = edNp.Value;
        mpc.Nc    = min(mpc.Np,DEF.NcMax);
        mpc.TsMax = DEF.TsMax;
        mpc.MaxSteps = DEF.MaxSteps;
        % Limits are ENFORCED inside the MPC optimisation only for ticked boxes; otherwise
        % the corresponding bound is infinite here (the limit is then only monitored).
        lim = readLimits();
        if cbEnfV.Value, mpc.Vmin = lim.Vmin; mpc.Vmax = lim.Vmax;
        else,            mpc.Vmin = -Inf;     mpc.Vmax = Inf;      end
        if cbEnfI.Value, mpc.Imax = lim.Imax;
        else,            mpc.Imax = Inf;      end
        validateMPC(mpc,DEF.MaxNp);
    end

    function lim = readLimits()
        % Current values of the voltage / current limit fields (always finite, validated).
        lim.Vmin = edVmin.Value;
        lim.Vmax = edVmax.Value;
        lim.Imax = edImax.Value;
        validateLimits(lim);
    end

    function txt = mpcEnforceNote()
        % Extra sentence for the all-controllers status when MPC runs with enforced limits.
        parts = {};
        if cbEnfV.Value, parts{end+1} = 'voltage'; end
        if cbEnfI.Value, parts{end+1} = 'current'; end
        if isempty(parts)
            txt = '';
        else
            txt = sprintf([' MPC enforced its %s limit(s); P, PI, PD and PID are NOT limited ' ...
                           '(their limits are only monitored).'],strjoin(parts,' and '));
        end
    end

    function M = collectMetrics(res,mot,lim)
        % Gather the REAL simulated signals and call the single metrics function.
        opts.ssWindow  = DEF.ssWindow;
        opts.settleTol = DEF.settleTol;
        opts.tolRel    = DEF.limTolRel;
        Va = [];
        if isfield(res,'u'), Va = res.u; end
        opts.vKind = 'sampled';
        if isfield(res,'uHold'), opts.vKind = res.uHold; end
        opts.impulseV = 0;
        if isfield(res,'uImpulse'), opts.impulseV = res.uImpulse; end
        opts.minTau = NaN;
        if isfield(res,'poles') && ~isempty(res.poles) && all(isfinite(res.poles))
            pm = max(abs(res.poles));
            if pm > 0, opts.minTau = 1/pm; end
        end
        M = computeMotorMetrics(res.t,res.w,res.th,res.i,Va,mot.ref,mot.mode,lim,opts);
        if ~M.vAvail && isfield(res,'uNote') && ~isempty(res.uNote)
            M.notes{end+1} = res.uNote;
        end
    end

    function showMetrics(nm,mot,res,lim,originTxt)
        % Compute and display the metrics of the simulation that has JUST been run.
        try
            M = collectMetrics(res,mot,lim);
            ctx.name     = nm;
            if strcmp(mot.mode,'position'), ctx.modeName = 'Position Control';
            else,                           ctx.modeName = 'Speed Control'; end
            ctx.enfV     = isfield(res,'enforceV') && res.enforceV;
            ctx.enfI     = isfield(res,'enforceI') && res.enforceI;
            ctx.unstable = isfield(res,'stable') && ~res.stable;
            ctx.origin   = originTxt;
            S = metricsText(M,ctx);
            MH.valCase.Text = S.caseTxt;   MH.valSS.Text   = S.ssTxt;
            MH.valISE.Text  = S.iseTxt;    MH.valIpk.Text  = S.ipkTxt;
            MH.valJu.Text   = S.juTxt;     MH.valVrms.Text = S.vrmsTxt;
            MH.valVmax.Text = S.vmaxTxt;   MH.valImax.Text = S.imaxTxt;
            MH.valVviol.Text = S.vviolTxt; MH.valIviol.Text = S.iviolTxt;
            MH.valVviol.FontColor = violColor(S.vviolState);
            MH.valIviol.FontColor = violColor(S.iviolState);
            MH.valNote.Text = S.noteTxt;
            if S.noteWarn, MH.valNote.FontColor = COL.warn; else, MH.valNote.FontColor = COL.info; end
            metricsShown = true;
        catch MEm
            resetMetricsPanel(['Metrics could not be computed: ' MEm.message]);
        end
    end

    function c = violColor(state)
        % state: 0 = none, 1 = violations, -1 = unavailable
        if state > 0, c = COL.err; elseif state == 0, c = COL.ok; else, c = COL.info; end
    end

    function resetMetricsPanel(msg)
        % Clear every number (no stale metrics) but show the definition and units of the
        % CURRENT control mode and the currently selected controller.
        if nargin < 1 || isempty(msg)
            msg = 'No metrics for this selection yet - press "Run Selected Controller".';
        end
        if strcmp(modeKey(),'position')
            eDef = 'e = thref - theta';   eu = 'rad';    iu = 'rad^2.s';
        else
            eDef = 'e = wref - omega';    eu = 'rad/s';  iu = '(rad/s)^2.s';
        end
        MH.valCase.Text = sprintf('%s | %s | %s | no valid result',ddCtrl.Value,curMode(),eDef);
        MH.valSS.Text   = sprintf('-- %s',eu);
        MH.valISE.Text  = sprintf('-- %s',iu);
        MH.valIpk.Text  = '-- A';
        MH.valJu.Text   = '-- V^2.s';
        MH.valVrms.Text = '-- V';
        MH.valVmax.Text = '-- V';
        MH.valImax.Text = '-- A';
        MH.valVviol.Text = '--';
        MH.valIviol.Text = '--';
        MH.valVviol.FontColor = COL.info;
        MH.valIviol.FontColor = COL.info;
        MH.valNote.Text = msg;
        MH.valNote.FontColor = COL.info;
        metricsShown = false;
    end

    function markMetricsStale()
        % An input changed after the last simulation: keep the numbers (they are still the
        % true result of the previous run) but say clearly that they are out of date.
        if ~metricsShown, return; end
        MH.valNote.Text = ['OUT OF DATE: an input (motor, reference, gains, MPC, limits) changed after ' ...
                           'this simulation. The values above belong to the previous run - run again to update.'];
        MH.valNote.FontColor = COL.warn;
        MH.valCase.Text = [MH.valCase.Text '  [out of date]'];
        metricsShown = false;
    end

    function plotSingle(nm,mot,res)
        cc = ctrlColor(nm);
        isPos = strcmp(mot.mode,'position');
        if res.stable
            suffix = sprintf(' - %s Controller%s',nm,modeTag());
            tc = [0 0 0];
        else
            suffix = sprintf(' - %s Controller%s (UNSTABLE)',nm,modeTag());
            tc = COL.err;
        end
        prepareAxes(axI, ['Armature Current i(t)' suffix],  'Current i(t) [A]',            mot.T,tc);
        prepareAxes(axTh,['Angular Position \theta(t)' suffix],'Position \theta(t) [rad]', mot.T,tc);
        prepareAxes(axW, ['Angular Speed \omega(t)' suffix],'Speed \omega(t) [rad/s]',      mot.T,tc);
        plot(axI, res.t,res.i, 'Color',cc,'LineWidth',1.8,'DisplayName',nm);
        if isPos
            plot(axTh,res.t,res.th,'Color',cc,'LineWidth',1.8,'DisplayName',[nm ' position \theta(t)']);
            plot(axTh,[0 mot.T],[mot.thref mot.thref],'--','Color',[0.15 0.15 0.15], ...
                 'LineWidth',1.4,'DisplayName','Reference position');
            plot(axW, res.t,res.w, 'Color',cc,'LineWidth',1.8,'DisplayName',[nm ' speed \omega(t)']);
            legend(axTh,'Location','best');
            legend(axW,'Location','best');
        else
            plot(axTh,res.t,res.th,'Color',cc,'LineWidth',1.8,'DisplayName',nm);
            plot(axW, res.t,res.w, 'Color',cc,'LineWidth',1.8,'DisplayName',[nm ' speed \omega(t)']);
            plot(axW,[0 mot.T],[mot.wref mot.wref],'--','Color',[0.15 0.15 0.15], ...
                 'LineWidth',1.4,'DisplayName','Reference speed');
            legend(axW,'Location','best');
        end
    end

    function highlightSelectedRow()
        % Visually mark the selected controller's row (cosmetic only).
        try
            tbs = {tblGains,tblGainsPos};
            for q = 1:2, removeStyle(tbs{q}); end
            selRow = find(strcmp(DEF.names,ddCtrl.Value),1);
            if isempty(selRow), return; end    % MPC has no row in the gain tables
            for q = 1:2
                addStyle(tbs{q},uistyle('BackgroundColor',[0.80 0.89 1.00], ...
                                        'FontWeight','bold'),'row',selRow);
            end
        catch
            % Styling is optional; ignore if unsupported.
        end
    end

    function setStatus(txt,kind)
        switch kind
            case 'ok',    statusLabel.FontColor = COL.ok;
            case 'warn',  statusLabel.FontColor = COL.warn;
            case 'error', statusLabel.FontColor = COL.err;
            otherwise,    statusLabel.FontColor = COL.info;
        end
        statusLabel.Text = txt;
        drawnow limitrate;
    end

end  % ======================= end of main function =======================


% ========================================================================
%                             LOCAL FUNCTIONS
% ========================================================================

function okFlag = hasControlToolbox()
% True if the Control System Toolbox functions we need are available.
okFlag = false;
try
    okFlag = license('test','Control_System_Toolbox') || license('test','Control_Toolbox');
    okFlag = okFlag && ~isempty(which('tf')) && ~isempty(which('lsim'));
catch
    okFlag = false;
end
end

function [g,hdr] = createSection(parent,row,titleText,rowHeights,colWidths,COL)
% Panel with a blue header bar; returns the inner grid (row 1 = header).
p = uipanel(parent,'BorderType','line','BackgroundColor',COL.panel);
p.Layout.Row = row;
nCols = numel(colWidths);
g = uigridlayout(p,[numel(rowHeights)+1, nCols]);
g.RowHeight       = [{26}, rowHeights];
g.ColumnWidth     = colWidths;
g.Padding         = [6 6 6 6];
g.RowSpacing      = 5;
g.ColumnSpacing   = 6;
g.BackgroundColor = COL.panel;
hdr = uilabel(g,'Text',['  ' titleText],'FontWeight','bold','FontSize',13, ...
              'FontColor',[1 1 1],'BackgroundColor',COL.blue, ...
              'HorizontalAlignment','left');
hdr.Layout.Row = 1;
% A one-column grid requires a scalar column index; [1 1] is invalid.
if nCols > 1
    hdr.Layout.Column = [1 nCols];
else
    hdr.Layout.Column = 1;
end
end

function h = addLabel(parent,txt,row,col)
h = uilabel(parent,'Text',txt,'FontSize',12,'HorizontalAlignment','left');
h.Layout.Row = row;
h.Layout.Column = col;
end

function h = addNumField(parent,val,row,col,lo,hi,loInclusive)
h = uieditfield(parent,'numeric','Value',val,'Limits',[lo hi], ...
                'LowerLimitInclusive',loInclusive, ...
                'ValueDisplayFormat','%.6g','HorizontalAlignment','right','FontSize',12);
h.Layout.Row = row;
h.Layout.Column = col;
end

function v = onoff(tf)
if tf, v = 'on'; else, v = 'off'; end
end

function T = makeGainTable(DEF,G)
% Table with columns Controller | Kp | Ki | Kd (one independent row per controller).
% G = gain matrix [Kp Ki Kd] per controller; defaults to the speed-mode gains.
if nargin < 2, G = DEF.gains; end
T = table(DEF.names(:), G(:,1), G(:,2), G(:,3), ...
          'VariableNames',{'Controller','Kp','Ki','Kd'});
end

function prepareAxes(ax,ttl,ylab,T,ttlColor)
% Clear an existing uiaxes and apply title / labels / grid (no new figures).
try, legend(ax,'off'); catch, end %#ok<NOCOM>
cla(ax);
hold(ax,'on');
grid(ax,'on');
box(ax,'on');
ax.FontSize = 12;
title(ax,ttl,'FontSize',13,'FontWeight','bold','Color',ttlColor);
xlabel(ax,'Time [s]','FontSize',12);
ylabel(ax,ylab,'FontSize',12);
xlim(ax,[0 T]);
ylim(ax,'auto');
end

function c = ctrlColor(nm)
switch nm
    case 'P',   c = [0.000 0.447 0.741];
    case 'PI',  c = [0.850 0.325 0.098];
    case 'PD',  c = [0.466 0.674 0.188];
    case 'PID', c = [0.494 0.184 0.556];
    case 'MPC', c = [0.635 0.078 0.184];
    otherwise,  c = [0 0 0];
end
end

function s = gainString(nm,kp,ki,kd)
switch nm
    case 'P',   s = sprintf('Kp=%.4g',kp);
    case 'PI',  s = sprintf('Kp=%.4g, Ki=%.4g',kp,ki);
    case 'PD',  s = sprintf('Kp=%.4g, Kd=%.4g',kp,kd);
    otherwise,  s = sprintf('Kp=%.4g, Ki=%.4g, Kd=%.4g',kp,ki,kd);
end
end

function s = mpcString(r)
% One-line description of the MPC configuration actually used.
s = sprintf('%s tracking: Q=%.4g, R=%.4g, Np=%d, Nc=%d, Ts=%.4g ms, %d samples', ...
            r.mode,r.Q,r.R,r.Np,r.Nc,r.Ts*1e3,r.N);
end

function s = mpcStabilityText(r)
if r.constrained
    parts = {};
    if r.enforceV
        parts{end+1} = sprintf('voltage limits [%.4g, %.4g] V',r.Vmin,r.Vmax);
    end
    if r.enforceI
        parts{end+1} = sprintf('armature-current limit |ia| <= %.4g A (at the sampling instants)',r.Imax);
    end
    s = sprintf(['Enforced inside the MPC optimisation with quadprog: %s (the classical controllers ' ...
                 'are NOT limited by these settings; their limits are only monitored). '], ...
                strjoin(parts,' and '));
    if r.enforceV
        s = [s 'The motor is open-loop stable and the input is bounded, so the response is bounded. '];
    else
        s = [s sprintf('Stability is judged from the unconstrained MPC law (largest discrete pole magnitude %.4f). ',r.rho)];
    end
    s = [s sprintf('Peak |Va| = %.4g V.',r.peakV)];
else
    s = sprintf(['Stable (largest discrete pole magnitude %.4f). ' ...
                 'No voltage or current limit enforced; peak |Va| = %.4g V.'],r.rho,r.peakV);
end
end

% ---------------------------- validation --------------------------------

function validateMPC(p,maxNp)
if ~(isnumeric(p.Q) && isscalar(p.Q) && isreal(p.Q) && isfinite(p.Q) && p.Q > 0)
    error('DCMotorGUI:InvalidInput','MPC coefficient Q must be a finite, real, strictly positive number.');
end
if ~(isnumeric(p.R) && isscalar(p.R) && isreal(p.R) && isfinite(p.R) && p.R > 0)
    error('DCMotorGUI:InvalidInput','MPC coefficient R must be a finite, real, strictly positive number.');
end
if ~(isnumeric(p.Np) && isscalar(p.Np) && isreal(p.Np) && isfinite(p.Np) && ...
     p.Np >= 1 && p.Np == round(p.Np) && p.Np <= maxNp)
    error('DCMotorGUI:InvalidInput', ...
          'MPC prediction horizon Np must be an integer between 1 and %d.',maxNp);
end
if ~(p.Nc >= 1 && p.Nc <= p.Np && p.Nc == round(p.Nc))
    error('DCMotorGUI:InvalidInput','Internal control horizon Nc must satisfy 1 <= Nc <= Np.');
end
if ~(p.Vmin < p.Vmax)   % also false for NaN
    error('DCMotorGUI:InvalidInput','Voltage limits must satisfy Vmin < Vmax.');
end
if ~(isnumeric(p.Imax) && isscalar(p.Imax) && ~isnan(p.Imax) && p.Imax > 0)   % Inf = not enforced
    error('DCMotorGUI:InvalidInput','The enforced current limit must be strictly positive.');
end
end

function validateInputs(m)
posFields = {'La','Ra','J','Kt','Ke'};
for k = 1:numel(posFields)
    v = m.(posFields{k});
    if ~(isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v) && v > 0)
        error('DCMotorGUI:InvalidInput', ...
              'Motor parameter %s must be a finite, real, strictly positive number.',posFields{k});
    end
end
if ~(isnumeric(m.B) && isscalar(m.B) && isreal(m.B) && isfinite(m.B) && m.B >= 0)
    error('DCMotorGUI:InvalidInput', ...
          'Motor parameter B must be a finite, real, non-negative number.');
end
if ~(isnumeric(m.wref) && isscalar(m.wref) && isreal(m.wref) && isfinite(m.wref))
    error('DCMotorGUI:InvalidInput','Reference speed must be a finite real number.');
end
if ~(isnumeric(m.thref) && isscalar(m.thref) && isreal(m.thref) && isfinite(m.thref))
    error('DCMotorGUI:InvalidInput','Reference position must be a finite real number.');
end
if ~(isnumeric(m.T) && isscalar(m.T) && isreal(m.T) && isfinite(m.T) && m.T > 0 && m.T <= 1000)
    error('DCMotorGUI:InvalidInput', ...
          'Simulation duration must be a finite number in the range (0, 1000] seconds.');
end
end

function validateGains(nm,kp,ki,kd)
gNames = {'Kp','Ki','Kd'};
vals   = {kp,ki,kd};
switch nm
    case 'P',   used = [true false false];
    case 'PI',  used = [true true  false];
    case 'PD',  used = [true false true];
    case 'PID', used = [true true  true];
    otherwise
        error('DCMotorGUI:InvalidInput','Unknown controller "%s".',nm);
end
for k = 1:3
    if used(k)
        v = vals{k};
        if ~(isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v) && v >= 0)
            error('DCMotorGUI:InvalidInput', ...
                  '%s gain %s must be a finite, non-negative number.',nm,gNames{k});
        end
    end
end
end

% ------------------------- performance metrics ---------------------------

function M = computeMotorMetrics(t,w,th,ia,Va,ref,mode,lim,opts)
%COMPUTEMOTORMETRICS  Performance metrics from SIMULATED signals (no GUI dependency).
%
%   M = computeMotorMetrics(t,w,th,ia,Va,ref,mode,lim,opts)
%     t    simulation time vector [s] (finite, strictly increasing, may be nonuniform)
%     w    angular speed omega(t) [rad/s]        (same length as t)
%     th   angular position theta(t) [rad]       (same length as t)
%     ia   armature current [A]                  (same length as t)
%     Va   applied armature voltage [V]; [] if unavailable. Length numel(t)
%          (opts.vKind = 'sampled', trapezoidal integration) or numel(t)-1
%          (opts.vKind = 'zoh': Va(k) is held on [t(k),t(k+1)), exact integration)
%     ref  reference (scalar or vector) in the unit of the tracked output
%     mode 'speed' (e = ref - w, rad/s) or 'position' (e = ref - theta, rad)
%     lim  struct with finite Vmin, Vmax [V] and Imax [A] (violation counts are
%          reported as unavailable if a limit is missing / not finite)
%     opts optional struct: ssWindow (0.10), settleTol (1e-3), tolRel (0), vKind
%          ('sampled'), impulseV (0; weight [V.s] of a Dirac impulse in Va at t(1)),
%          minTau (NaN; fastest closed-loop time constant, only used for a note)
%
%   Unavailable metrics are NaN with the matching *Avail flag false (never zero).
%   A sample ON a limit is not a violation: Va > Vmax+tol, Va < Vmin-tol and
%   |ia| > Imax+tol are counted, with tol = opts.tolRel*max(1,|limit|).

if nargin < 9 || isempty(opts), opts = struct(); end
isModeOk = (ischar(mode) && isrow(mode)) || (isstring(mode) && isscalar(mode));
if ~isModeOk || ~any(strcmpi(char(mode),{'speed','position'}))
    M = emptyMetrics('speed');
    M.notes{end+1} = 'Unknown operating mode (expected speed or position): no metrics computed.';
    return;
end
mode  = lower(char(mode));
isPos = strcmp(mode,'position');
M = emptyMetrics(mode);
[opts,optNotes] = fillMetricOpts(opts);
M.notes = [M.notes optNotes];
M.ssFrac = opts.ssWindow;

% ---- time vector --------------------------------------------------------
if ~(isnumeric(t) && isreal(t) && numel(t) >= 2 && all(isfinite(t(:))))
    M.notes{end+1} = 'Invalid time vector (need >= 2 finite real samples): no metrics computed.';
    return;
end
t  = t(:);
dt = diff(t);
if ~all(dt > 0)
    M.notes{end+1} = 'Time vector is not strictly increasing: no metrics computed.';
    return;
end
n = numel(t);
M.nSamples = n;
M.T = t(end) - t(1);

% ---- tracking signal and reference ---------------------------------------
if isPos, y = th; else, y = w; end
yOK   = okSignal(y,n);
refOK = isnumeric(ref) && isreal(ref) && ~isempty(ref) && all(isfinite(ref(:))) && ...
        (isscalar(ref) || numel(ref) == n);

if yOK && refOK
    y = y(:);
    if isscalar(ref), r = ref*ones(n,1); else, r = ref(:); end
    e = r - y;                                    % tracking error, ref - output

    % Metric 2: ISE = integral of e^2 dt (trapezoidal rule on the actual time vector)
    ise = trapz(t,e.^2);
    if isfinite(ise)
        M.ISE = ise;  M.iseAvail = true;
    else
        M.notes{end+1} = 'ISE is not finite (overflow of a diverging response).';
    end

    % Metric 1: steady-state error from the final window (time-weighted mean)
    tStart = t(end) - opts.ssWindow*M.T;
    idx = find(t >= tStart - 1e-12*M.T);
    if numel(idx) < 3, idx = (max(1,n-2):n).'; end      % always >= 2-3 samples
    tw = t(idx);  yw = y(idx);  rw = r(idx);
    span = tw(end) - tw(1);
    yss  = trapz(tw,yw)/span;
    rss  = trapz(tw,rw)/span;
    M.yss = yss;  M.refEnd = rss;
    M.ess = rss - yss;                                  % signed  e_ss = r - y_ss
    M.essAbs = abs(M.ess);
    M.ssWindowStart = tw(1);
    M.ssPeakToPeak  = max(yw) - min(yw);
    scale = max([abs(rss), max(abs(yw)), realmin]);
    M.settled = M.ssPeakToPeak <= opts.settleTol*scale;
    M.essAvail = isfinite(M.ess);
else
    M.notes{end+1} = 'Tracking signal or reference is missing / invalid: steady-state error and ISE are unavailable.';
end

% ---- Metric 3: peak armature current -------------------------------------
if okSignal(ia,n)
    ia = ia(:);
    M.Ipeak = max(abs(ia));
    M.ipeakAvail = true;
else
    ia = [];
    M.notes{end+1} = 'Armature-current signal is missing / invalid: peak current and current violations are unavailable.';
end

% ---- Metric 4: control effort from the ACTUAL applied voltage ----------------
isZoh = strcmp(opts.vKind,'zoh');
if isZoh, nv = n-1; else, nv = n; end
impOK = isnumeric(opts.impulseV) && isscalar(opts.impulseV) && isreal(opts.impulseV) && isfinite(opts.impulseV);
v = [];
if isempty(Va)
    M.notes{end+1} = 'Applied-voltage signal is not available: control effort and voltage metrics are unavailable (not assumed).';
elseif ~okSignal(Va,nv)
    M.notes{end+1} = sprintf('Applied-voltage signal is invalid or has the wrong length (expected %d samples).',nv);
elseif ~impOK
    M.notes{end+1} = 'Impulse weight of the applied voltage is invalid: voltage metrics are unavailable.';
else
    v = Va(:);
    if isZoh
        JuReg = sum(v.^2 .* dt);                       % exact for a held (piecewise-constant) input
    else
        JuReg = trapz(t,v.^2);
    end
    if isfinite(JuReg)
        M.JuReg = JuReg;
        M.VrmsReg = sqrt(JuReg/M.T);
        M.VmaxAbsReg = max(abs(v));
        M.vImpulse = opts.impulseV;
        M.vHasImpulse = (opts.impulseV ~= 0);
        M.vZOH = isZoh;
        if M.vHasImpulse        % energy of a Dirac impulse is infinite
            M.Ju = Inf;  M.Vrms = Inf;  M.VmaxAbs = Inf;
        else
            M.Ju = JuReg;  M.Vrms = M.VrmsReg;  M.VmaxAbs = M.VmaxAbsReg;
        end
        M.vAvail = true;
    else
        v = [];
        M.notes{end+1} = 'Control effort is not finite (overflow of a diverging response).';
    end
end

% ---- Metric 5: limit violations (monitoring) -----------------------------------
if ~isempty(v) && limitOK(lim,'Vmin') && limitOK(lim,'Vmax') && lim.Vmin < lim.Vmax
    M.Vmin = lim.Vmin;  M.Vmax = lim.Vmax;
    tolHi = opts.tolRel*max(1,abs(lim.Vmax));
    tolLo = opts.tolRel*max(1,abs(lim.Vmin));
    M.nVHigh = nnz(v > lim.Vmax + tolHi);
    M.nVLow  = nnz(v < lim.Vmin - tolLo);
    M.nVViol = M.nVHigh + M.nVLow;
    M.nVSamples = numel(v);
    M.vImpulseViol = M.vHasImpulse;     % a Dirac impulse exceeds ANY finite limit
    M.vViolAvail = true;
elseif ~isempty(v)
    M.notes{end+1} = 'Voltage limits are not configured (finite Vmin < Vmax): voltage violations are unavailable.';
end
if ~isempty(ia) && limitOK(lim,'Imax') && lim.Imax > 0
    M.Imax = lim.Imax;
    tolI = opts.tolRel*max(1,abs(lim.Imax));
    M.nIViol = nnz(abs(ia) > lim.Imax + tolI);
    M.nISamples = numel(ia);
    M.iViolAvail = true;
elseif ~isempty(ia)
    M.notes{end+1} = 'Current limit is not configured (finite Imax > 0): current violations are unavailable.';
end

% ---- resolution warning ------------------------------------------------------
if isfinite(opts.minTau) && opts.minTau > 0 && max(dt) > 0.5*opts.minTau
    M.notes{end+1} = sprintf(['The time step (%.3g ms) is not small compared with the fastest closed-loop time ' ...
        'constant (%.3g ms): ISE, effort and sampled peaks are only approximate (errors of a few %%).'], ...
        1e3*max(dt),1e3*opts.minTau);
end

M.valid = M.essAvail || M.iseAvail || M.ipeakAvail || M.vAvail;
end

function M = emptyMetrics(mode)
% All metrics start as UNAVAILABLE (NaN + false flag), never as zero.
if strcmp(mode,'position')
    M.mode = 'position';  M.trackName = 'theta';  M.errUnit = 'rad';    M.iseUnit = 'rad^2.s';
else
    M.mode = 'speed';     M.trackName = 'omega';  M.errUnit = 'rad/s';  M.iseUnit = '(rad/s)^2.s';
end
M.valid = false;   M.notes = {};   M.nSamples = 0;   M.T = NaN;
M.ssFrac = NaN;
M.essAvail = false;  M.ess = NaN;  M.essAbs = NaN;  M.yss = NaN;  M.refEnd = NaN;
M.settled = false;   M.ssWindowStart = NaN;  M.ssPeakToPeak = NaN;
M.iseAvail = false;  M.ISE = NaN;
M.ipeakAvail = false;  M.Ipeak = NaN;
M.vAvail = false;    M.vZOH = false;  M.vHasImpulse = false;  M.vImpulse = 0;
M.Ju = NaN;  M.Vrms = NaN;  M.VmaxAbs = NaN;
M.JuReg = NaN;  M.VrmsReg = NaN;  M.VmaxAbsReg = NaN;
M.vViolAvail = false;  M.Vmin = NaN;  M.Vmax = NaN;
M.nVHigh = NaN;  M.nVLow = NaN;  M.nVViol = NaN;  M.nVSamples = NaN;  M.vImpulseViol = false;
M.iViolAvail = false;  M.Imax = NaN;  M.nIViol = NaN;  M.nISamples = NaN;
end

function [o,notes] = fillMetricOpts(o)
notes = {};
if ~isstruct(o), o = struct(); end
def = struct('ssWindow',0.10,'settleTol',1e-3,'tolRel',0,'vKind','sampled','impulseV',0,'minTau',NaN);
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(o,fn{k}) || isempty(o.(fn{k})), o.(fn{k}) = def.(fn{k}); end
end
if ~(isnumeric(o.ssWindow) && isscalar(o.ssWindow) && isfinite(o.ssWindow) && o.ssWindow > 0 && o.ssWindow <= 1)
    o.ssWindow = def.ssWindow;
    notes{end+1} = 'Invalid steady-state window; the default 10 % of the run was used.';
end
if ~(isnumeric(o.settleTol) && isscalar(o.settleTol) && isfinite(o.settleTol) && o.settleTol > 0)
    o.settleTol = def.settleTol;
    notes{end+1} = 'Invalid settling tolerance; the default 0.1 % was used.';
end
if ~(isnumeric(o.tolRel) && isscalar(o.tolRel) && isfinite(o.tolRel) && o.tolRel >= 0)
    o.tolRel = def.tolRel;
    notes{end+1} = 'Invalid limit tolerance; strict comparison (0) was used.';
end
if ~(ischar(o.vKind) && any(strcmp(o.vKind,{'sampled','zoh'})))
    o.vKind = def.vKind;
    notes{end+1} = 'Invalid voltage-signal kind; "sampled" was assumed.';
end
end

function tf = okSignal(x,n)
% Numeric, real, finite vector of exactly n samples.
tf = isnumeric(x) && isreal(x) && isvector(x) && numel(x) == n && all(isfinite(x(:)));
end

function tf = limitOK(lim,f)
tf = isstruct(lim) && isfield(lim,f) && isnumeric(lim.(f)) && isscalar(lim.(f)) && ...
     isreal(lim.(f)) && isfinite(lim.(f));
end

function S = metricsText(M,ctx)
% Turn the metrics structure into the strings shown in the panel (unavailable -> "n/a").
if strcmp(M.mode,'position')
    eDef = 'e = thref - theta';
else
    eDef = 'e = wref - omega';
end
S.caseTxt = sprintf('%s | %s | %s | T = %.4g s',ctx.name,ctx.modeName,eDef,M.T);
if ~M.valid
    S.caseTxt = sprintf('%s | %s | %s | no valid result',ctx.name,ctx.modeName,eDef);
end
notes = {};
if ~isempty(ctx.origin), notes{end+1} = ['Metrics ' ctx.origin '.']; end
if ctx.unstable
    notes{end+1} = ['UNSTABLE closed loop: the metrics describe a diverging response, ' ...
                    'not a successful control result.'];
end

% 1. steady-state error
if M.essAvail
    if M.settled, tag = 'settled'; else, tag = 'NOT settled'; end
    S.ssTxt = sprintf('e_ss = %+.4e %s    |e_ss| = %.4e %s    [%s]',M.ess,M.errUnit,M.essAbs,M.errUnit,tag);
    if ~M.settled
        notes{end+1} = sprintf(['e_ss is the mean error over the last %.3g %% of the run, but the output ' ...
            'is still changing there (spread %.3g %s): it is not a true steady-state error.'], ...
            100*M.ssFrac,M.ssPeakToPeak,M.errUnit);
    end
else
    S.ssTxt = sprintf('n/a [%s]',M.errUnit);
end
% 2. ISE
if M.iseAvail, S.iseTxt = sprintf('%.4e %s',M.ISE,M.iseUnit); else, S.iseTxt = 'n/a'; end
% 3. peak current (also the maximum observed |ia|)
if M.ipeakAvail
    S.ipkTxt  = sprintf('%.4g A',M.Ipeak);
    S.imaxTxt = sprintf('%.4g A',M.Ipeak);
else
    S.ipkTxt = 'n/a';  S.imaxTxt = 'n/a';
end
% 4. control effort (and RMS voltage, max |Va|)
if M.vAvail
    if M.vHasImpulse
        S.juTxt   = 'Inf (impulse in Va)';
        S.vrmsTxt = 'Inf (impulse in Va)';
        S.vmaxTxt = 'Inf (impulse in Va)';
        notes{end+1} = sprintf(['Va contains a Dirac impulse of %.4g V.s at t = 0 (ideal derivative acting on ' ...
            'the reference step), so Ju, V_RMS and max|Va| are infinite. Excluding it (t > 0): ' ...
            'Ju = %.4g V^2.s, V_RMS = %.4g V, max|Va| = %.4g V.'], ...
            M.vImpulse,M.JuReg,M.VrmsReg,M.VmaxAbsReg);
    else
        S.juTxt   = sprintf('%.4e V^2.s',M.Ju);
        S.vrmsTxt = sprintf('%.4g V',M.Vrms);
        S.vmaxTxt = sprintf('%.4g V',M.VmaxAbs);
    end
else
    S.juTxt = 'n/a';  S.vrmsTxt = 'n/a';  S.vmaxTxt = 'n/a';
end
% 5. violations: counts + whether the limit is ENFORCED by the controller or only MONITORED
if M.vViolAvail
    if ctx.enfV
        how = 'ENFORCED by MPC (at sampling instants)';
    else
        how = 'MONITORED only (not enforced)';
    end
    S.vviolTxt = sprintf('%d of %d samples outside [%.4g, %.4g] V - %s',M.nVViol,M.nVSamples,M.Vmin,M.Vmax,how);
    if M.vImpulseViol
        S.vviolTxt = [S.vviolTxt '; + impulse at t = 0'];
    end
    if M.nVViol > 0 || M.vImpulseViol, S.vviolState = 1; else, S.vviolState = 0; end
    if ctx.enfV && M.nVViol > 0
        notes{end+1} = 'WARNING: the voltage limit is enforced by MPC but samples exceed it.';
    end
else
    S.vviolTxt = 'n/a'; S.vviolState = -1;
end
if M.iViolAvail
    if ctx.enfI
        how = 'ENFORCED by MPC (at sampling instants)';
    else
        how = 'MONITORED only (not enforced)';
    end
    S.iviolTxt = sprintf('%d of %d samples with |ia| > %.4g A - %s',M.nIViol,M.nISamples,M.Imax,how);
    if M.nIViol > 0, S.iviolState = 1; else, S.iviolState = 0; end
    if ctx.enfI && M.nIViol > 0
        notes{end+1} = 'WARNING: the current limit is enforced by MPC but samples exceed it.';
    end
else
    S.iviolTxt = 'n/a'; S.iviolState = -1;
end
notes = [notes M.notes];
if isempty(notes)
    S.noteTxt = ['Steady state = mean of the final ' sprintf('%.3g',100*M.ssFrac) ' % of the run. ' ...
                 'Peaks and integrals are evaluated on the simulation time grid. ' ...
                 'Limits are monitored for every controller; only MPC can enforce them.'];
    S.noteWarn = false;
else
    S.noteTxt = strjoin(notes,' ');
    S.noteWarn = true;
end
end

function validateLimits(lim)
% Finite, physically meaningful limits are required: the app never assumes "unlimited".
ok = @(v) isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v);
if ~(ok(lim.Vmin) && ok(lim.Vmax) && ok(lim.Imax))
    error('DCMotorGUI:InvalidInput', ...
          ['The voltage limits (Vmin, Vmax) and the current limit must be finite numbers; ' ...
           'the application does not assume unlimited voltage or current.']);
end
if ~(lim.Vmin < lim.Vmax)
    error('DCMotorGUI:InvalidInput','Voltage limits must satisfy Vmin < Vmax.');
end
if ~(lim.Imax > 0)
    error('DCMotorGUI:InvalidInput','The maximum armature current must be strictly positive.');
end
end

% ------------------------- polynomial helpers ----------------------------

function p = trimPoly(p)
% Remove leading zeros (keeps at least one coefficient).
p = p(:).';
k = find(p ~= 0,1,'first');
if isempty(k), p = 0; else, p = p(k:end); end
end

function c = polyAdd(a,b)
% Add polynomials of different degree (descending-power coefficient vectors).
a = a(:).';  b = b(:).';
n = max(numel(a),numel(b));
c = [zeros(1,n-numel(a)) a] + [zeros(1,n-numel(b)) b];
end

function [Nc,Dc] = controllerPolynomials(nm,kp,ki,kd)
% C(s) = Nc(s)/Dc(s) exactly as specified:
%   P  : Kp            PI : Kp + Ki/s = (Kp*s + Ki)/s
%   PD : Kp + Kd*s     PID: Kp + Ki/s + Kd*s = (Kd*s^2 + Kp*s + Ki)/s
% If Ki = 0 the integral term is exactly zero, so the 1/s factor is dropped
% (this is the same transfer function, without a cancelling pole/zero at s = 0).
useInt = ki > 0;
switch nm
    case 'P'
        Nc = kp;            Dc = 1;
    case 'PI'
        if useInt, Nc = [kp ki];     Dc = [1 0]; else, Nc = kp;        Dc = 1; end
    case 'PD'
        Nc = [kd kp];       Dc = 1;
    case 'PID'
        if useInt, Nc = [kd kp ki];  Dc = [1 0]; else, Nc = [kd kp];   Dc = 1; end
    otherwise
        error('DCMotorGUI:InvalidInput','Unknown controller "%s".',nm);
end
Nc = trimPoly(Nc);
end

% ----------------------------- simulation --------------------------------

function sim = simulateController(nm,m,kp,ki,kd,nPts)
% Closed-loop simulation of current, position and speed for one controller.

[Nc,Dc] = controllerPolynomials(nm,kp,ki,kd);

% D(s) = (La*s + Ra)(J*s + B) + Kt*Ke
Dp = conv([m.La m.Ra],[m.J m.B]);
Dp(end) = Dp(end) + m.Kt*m.Ke;

% Closed-loop characteristic polynomial.
%   speed    : Dcl = D*Dc + Kt*Nc      (from 1 + P_w*C = 0)
%   position : Dcl = s*D*Dc + Kt*Nc    (from 1 + P_theta*C = 0, extra integrator in the plant)
isPos = isfield(m,'mode') && strcmp(m.mode,'position');
if isPos
    Dcl = trimPoly(polyAdd(conv(conv([1 0],Dp),Dc), m.Kt*Nc));
else
    Dcl = trimPoly(polyAdd(conv(Dp,Dc), m.Kt*Nc));
end

% Stability from closed-loop poles (all must have strictly negative real part).
p = roots(Dcl);
if isempty(p) || any(~isfinite(p))
    error('DCMotorGUI:Numerical','Could not compute closed-loop poles (non-finite values).');
end
tol    = 1e-9*max(1,max(abs(p)));
maxRe  = max(real(p));
stable = all(real(p) < -tol);

% Closed-loop transfer functions  X/Wref = P_X*C/(1 + P_w*C)  (D(s) cancelled exactly)
if isPos
    % Closed loop from the position reference (P_X*C/(1 + P_theta*C), D(s) cancelled exactly)
    Ti  = tf(trimPoly(conv(conv([1 0],[m.J m.B]),Nc)), Dcl);  % current  I/Thref
    Tw  = tf(trimPoly(conv([1 0],m.Kt*Nc)),            Dcl);  % speed    W/Thref
    Tth = tf(trimPoly(m.Kt*Nc),                        Dcl);  % position Th/Thref
else
    Ti  = tf(trimPoly(conv([m.J m.B],Nc)), Dcl);      % current  I/Wref
    Tw  = tf(trimPoly(m.Kt*Nc),            Dcl);      % speed    W/Wref
    Tth = tf(trimPoly(m.Kt*Nc),            [Dcl 0]);  % position Theta/Wref  (extra 1/s)
end

% Time grid and constant reference input (speed in rad/s, or position in rad)
t = linspace(0,m.T,nPts).';
u = m.ref*ones(size(t));

try
    yi  = lsim(Ti, u,t,'zoh');
    yw  = lsim(Tw, u,t,'zoh');
    yth = lsim(Tth,u,t,'zoh');
catch ME
    error('DCMotorGUI:SimFailed','lsim failed: %s',ME.message);
end

if ~(all(isfinite(yi)) && all(isfinite(yw)) && all(isfinite(yth)))
    error('DCMotorGUI:Numerical', ...
          ['The simulation produced non-finite values (overflow or numerical failure). ' ...
           'Check the parameters and reduce the controller gains.']);
end

% Applied armature voltage Va(t) = controller output, from the SAME closed loop
[Va,vImp,vNote] = classicalVoltage(m,Nc,Dp,Dcl,isPos,t,u);

sim.t = t;   sim.i = yi;   sim.w = yw;   sim.th = yth;
sim.u = Va;                 % applied voltage [V] on the same time grid (proper part, t > 0)
sim.uImpulse = vImp;        % weight [V.s] of a Dirac impulse in Va at t = 0 (PD / PID), else 0
sim.uHold = 'sampled';      % Va is a continuous-time signal sampled on t
sim.uNote = vNote;          % non-empty if Va could not be derived (then it is UNAVAILABLE)
sim.enforceV = false;       % classical controllers are not saturated: limits are only monitored
sim.enforceI = false;
sim.poles = p;   sim.maxRe = maxRe;   sim.stable = stable;
end

function [Va,vImp,note] = classicalVoltage(m,Nc,Dp,Dcl,isPos,t,uRef)
% Applied armature voltage of the classical loop for the constant reference m.ref:
%   speed    : Va/Wref  = C/(1 + P_w*C)     = Nc*D / Dcl
%   position : Va/Thref = C/(1 + P_theta*C) = s*D*Nc / Dcl
% (C = Nc/Dc, D(s) = Dp, Dcl = closed-loop characteristic polynomial used above.)
% For PD / PID the numerator is one degree higher than Dcl (ideal derivative on a reference
% step), so   Tu(s) = q1*s + q0 + r(s)/Dcl(s).   q1*s*(ref/s) is a Dirac impulse of weight
% q1*ref (= Kd*ref) at t = 0; the remainder is proper and is simulated with lsim.
Va = [];  vImp = 0;  note = '';
try
    if isPos
        NumU = trimPoly(conv(conv([1 0],Dp),Nc));
    else
        NumU = trimPoly(conv(Dp,Nc));
    end
    d = Dcl(:).';
    degN = numel(NumU) - 1;
    degD = numel(d) - 1;
    if degN <= degD
        Tu = tf(NumU,d);
    elseif degN == degD + 1
        % explicit polynomial division NumU = (q1*s + q0)*d + r, deg(r) < deg(d)
        nU = NumU(:).';
        q1 = nU(1)/d(1);
        R1 = nU - q1*[d 0];
        q0 = R1(2)/d(1);
        R2 = R1(2:end) - q0*d;
        r  = R2(2:end);
        Tu = tf(trimPoly(q0*d + [0 r]),d);
        vImp = q1*m.ref;
    else
        error('DCMotorGUI:Numerical','Controller output has a derivative of an impulse (not supported).');
    end
    Va = lsim(Tu,uRef,t,'zoh');
    if ~all(isfinite(Va))
        error('DCMotorGUI:Numerical','The applied-voltage simulation produced non-finite values.');
    end
catch MEv
    Va = [];  vImp = 0;
    note = ['Applied voltage of the classical controller could not be derived: ' MEv.message];
end
end

% ------------------------------- MPC -------------------------------------

function sim = simulateMPC(m,mpc,nPts)
% Receding-horizon MPC of the SAME linearized motor used by the classical controllers.
% State x = [ia; w; theta], input u = Va, zero initial state, constant reference m.wref.

% ---- sampling interval: Ts <= TsMax; for short runs Ts = T/(nPts-1) so the
%      MPC grid has the same nPts samples as the classical simulations.
TsNom = min(mpc.TsMax, m.T/(nPts-1));
N     = max(1,ceil(m.T/TsNom - 1e-9));          % number of sampling intervals
if N > mpc.MaxSteps
    error('DCMotorGUI:MPCTooLong', ...
          ['MPC would need %d samples (Ts <= %.4g ms, T = %.4g s), which exceeds the practical ' ...
           'limit of %d samples. Reduce the simulation duration to at most %.4g s.'], ...
          N, mpc.TsMax*1e3, m.T, mpc.MaxSteps, mpc.MaxSteps*mpc.TsMax);
end
Ts = m.T/N;                                     % exact grid ending at T

% ---- continuous model from the CURRENT motor parameters
A  = [-m.Ra/m.La, -m.Ke/m.La, 0;
       m.Kt/m.J,  -m.B/m.J,   0;
       0,          1,         0];
Bc = [1/m.La; 0; 0];
isPos = strcmp(m.mode,'position');
if isPos
    Cw = [0 0 1];                               % tracked output = position theta [rad]
    sim.mode = 'position';
else
    Cw = [0 1 0];                               % tracked output = speed w [rad/s]
    sim.mode = 'speed';
end

% ---- exact zero-order-hold discretization: expm([A Bc; 0 0]*Ts) = [Ad Bd; 0 1]
E  = expm([A Bc; zeros(1,4)]*Ts);
Ad = E(1:3,1:3);
Bd = E(1:3,4);
if ~(all(isfinite(Ad(:))) && all(isfinite(Bd(:))))
    error('DCMotorGUI:Numerical','Discretization of the motor model produced non-finite values.');
end

% ---- prediction matrices:  w_pred = Phi*x(k) + Gam*U,  U = [u(k) ... u(k+Nc-1)]'
Np = mpc.Np;  Nc = mpc.Nc;
[Phi,Gam] = buildPredictionMatrices(Ad,Bd,Cw,Np,Nc);

% ---- QP:  min_U  Q*||Phi*x - wref + Gam*U||^2 + R*||U||^2
%      = U'*H*U + 2*f'*U + const,   H = Q*Gam'*Gam + R*I,   f = Q*Gam'*(Phi*x - wref)
H = mpc.Q*(Gam'*Gam) + mpc.R*eye(Nc);
H = (H+H')/2;
[Lc,pd] = chol(H,'lower');
if pd ~= 0 || any(~isfinite(Lc(:)))
    error('DCMotorGUI:MPCQP','The MPC Hessian is not positive definite; check Q, R and Np.');
end
% Constraint ENFORCEMENT (inside the optimisation):  voltage  Vmin <= u(k+j) <= Vmax  and
% current  -Imax <= ia(k+j|k) <= Imax,  j = 1..Np  (current predicted at the sampling instants).
enforceV = isfinite(mpc.Vmin) || isfinite(mpc.Vmax);
enforceI = isfinite(mpc.Imax);
constrained = enforceV || enforceI;
if enforceI
    [PhiI,GamI] = buildPredictionMatrices(Ad,Bd,[1 0 0],Np,Nc);   % output = armature current ia
end
if constrained
    if isempty(which('quadprog'))
        error('DCMotorGUI:NoToolbox', ...
              ['A voltage / current limit is set to be ENFORCED, so MPC needs quadprog (Optimization ' ...
               'Toolbox), which is not available. Untick "MPC enforces V limits" / "MPC enforces I ' ...
               'limit" to run the unconstrained MPC (the limits are then only monitored).']);
    end
    lb = mpc.Vmin*ones(Nc,1);  ub = mpc.Vmax*ones(Nc,1);   % +-Inf where not enforced
    qpOpts = optimoptions('quadprog','Display','off');
end

% ---- receding-horizon closed loop
X  = zeros(3,N+1);                              % X(:,k+1) = x(k)
Uv = zeros(1,N);                                % applied voltage u(k)
x  = zeros(3,1);                                % same zero initial conditions as classical runs
for k = 1:N
    wPred0 = Phi*x;                             % free response of the speed over the horizon
    f      = mpc.Q*(Gam'*(wPred0 - m.ref));      % wPred0 = free response of the tracked output
    if constrained
        if enforceI
            iPred0 = PhiI*x;                    % free response of ia over the horizon
            Aineq  = [GamI; -GamI];
            bineq  = [mpc.Imax - iPred0; mpc.Imax + iPred0];
        else
            Aineq = [];  bineq = [];
        end
        [U,~,flag,out] = quadprog(2*H,2*f,Aineq,bineq,[],[],lb,ub,[],qpOpts);
        if flag <= 0 || isempty(U) || any(~isfinite(U))
            hint = '';
            if flag == -2
                hint = ' The enforced voltage / current limits are infeasible from this state - relax them.';
            end
            error('DCMotorGUI:MPCQP', ...
                  'MPC optimization failed or is infeasible at sample %d (t = %.4g s): quadprog exitflag %d (%s).%s', ...
                  k, (k-1)*Ts, flag, strtrim(out.message), hint);
        end
    else
        U = -(Lc'\(Lc\f));                      % solves H*U = -f
        if any(~isfinite(U))
            error('DCMotorGUI:MPCQP', ...
                  'MPC optimization returned NaN/Inf at sample %d (t = %.4g s).',k,(k-1)*Ts);
        end
    end
    u = U(1);                                   % apply ONLY the first move
    Uv(k)  = u;
    x      = Ad*x + Bd*u;                       % advance the plant one sampling interval
    if any(~isfinite(x))
        error('DCMotorGUI:Numerical', ...
              'MPC simulation produced non-finite states at sample %d (t = %.4g s).',k,k*Ts);
    end
    X(:,k+1) = x;
end

% ---- closed-loop pole check (unconstrained law). theta does not influence
%      (ia,w), so the (ia,w) block decides stability; theta is a free integrator,
%      exactly as in the classical speed-loop check.
Kx  = -(H\(mpc.Q*(Gam'*Phi)));                  % u(k) = e1'*U*(x) = Kx(1,:)*x + const
Acl = Ad + Bd*Kx(1,:);
if isPos
    rho = max(abs(eig(Acl)));                   % theta is fed back: all three poles count
else
    rho = max(abs(eig(Acl(1:2,1:2))));
end
if ~isfinite(rho), error('DCMotorGUI:Numerical','Could not compute the MPC closed-loop poles.'); end

sim.t  = (0:N).'*Ts;
sim.i  = X(1,:).';
sim.w  = X(2,:).';
sim.th = X(3,:).';
sim.u  = Uv(:);
sim.Q = mpc.Q;  sim.R = mpc.R;  sim.Np = Np;  sim.Nc = Nc;  sim.Ts = Ts;  sim.N = N;
sim.rho = rho;  sim.constrained = constrained;  sim.Vmin = mpc.Vmin;  sim.Vmax = mpc.Vmax;
sim.Imax = mpc.Imax;  sim.enforceV = enforceV;  sim.enforceI = enforceI;
sim.uHold = 'zoh';      % u(k) is held (zero-order hold) on [t(k), t(k+1)): numel(u) = numel(t)-1
sim.uImpulse = 0;       % MPC applies a bounded sampled voltage: no impulse
sim.uNote = '';
sim.peakV  = max(abs(Uv));
% bounded voltage + open-loop stable plant => bounded response (needs the VOLTAGE limit enforced)
sim.stable = enforceV || (rho < 1 - 1e-9);
sim.maxRe  = NaN;
end

function [Phi,Gam] = buildPredictionMatrices(Ad,Bd,Cout,Np,Nc)
% y_pred(j) = Cout*x(k+j|k):   y_pred = Phi*x(k) + Gam*U,   U = [u(k) ... u(k+Nc-1)]'
% (inputs after the control horizon Nc are held at u(k+Nc-1)).
Phi = zeros(Np,3);
Gam = zeros(Np,Nc);
Ap  = eye(3);                                   % Ad^(j-1)
h   = zeros(Np,1);                              % h(j) = Cout*Ad^(j-1)*Bd  (impulse response)
for j = 1:Np
    h(j)     = Cout*Ap*Bd;
    Ap       = Ap*Ad;                           % now Ad^j
    Phi(j,:) = Cout*Ap;
end
for j = 1:Np
    for mm = 0:j-1                              % input u(k+mm) affects y(k+j|k)
        col = min(mm,Nc-1)+1;                   % inputs beyond Nc-1 are held = u(k+Nc-1)
        Gam(j,col) = Gam(j,col) + h(j-mm);
    end
end
end
