%% ========================================================================
%  Layer1_AllStates_PublicationRun.m
%  Sex-Specific QSP Model of Diabetic Cardiomyopathy — Layer 1
%  Mitochondrial Biogenesis & Oxidative Stress Axis
%
%  Authors : Soheili M, Gilzad-Kohan H, Lotfi AS
%
%  Purpose : Run all 6 sex-disease states in sequence, save high-resolution
%            publication figures, export a structured results table, and
%            save the Layer 2 interface file for each state.
%
%  Version : 1.0  |  Date: 2026-03-09
%  MATLAB  : R2019b or newer required
%
%  Output files generated (all in same folder as this script):
%    - Layer1_AllStates_Summary.csv       (steady-state table, all states)
%    - Layer1_AllStates_Summary.xlsx      (same, formatted for manuscript)
%    - Layer1_ComparisonPlot.png/.pdf (publication figure, 300 dpi)
%    - Layer1_<state>_TimeCourse.png      (individual time-course, 300 dpi)
%    - layer1_to_layer2_<state>.mat       (Layer 2 interface per state)
%    - Layer1_RunLog.txt                  (reproducibility log)
%
%  Reproducibility: All parameters, versions, and timestamps are logged.
%% ========================================================================

clear; close all; clc;

%% --- 0. Setup & Logging ---------------------------------------------------
script_version = '1.0';
run_timestamp  = datestr(now, 'yyyy-mm-dd HH:MM:SS');
matlab_version = version;

% Output directory = same folder as this script
script_path = fileparts(mfilename('fullpath'));
if isempty(script_path), script_path = pwd; end
outdir = script_path;

% Open log file
logfile = fullfile(outdir, 'Layer1_RunLog.txt');
flog = fopen(logfile, 'w');
fprintf(flog, '=======================================================\n');
fprintf(flog, ' Layer 1 QSP Model — Reproducibility Log\n');
fprintf(flog, '=======================================================\n');
fprintf(flog, ' Script version : %s\n', script_version);
fprintf(flog, ' Run timestamp  : %s\n', run_timestamp);
fprintf(flog, ' MATLAB version : %s\n', matlab_version);
fprintf(flog, ' Output dir     : %s\n', outdir);
fprintf(flog, '=======================================================\n\n');

fprintf('=======================================================\n');
fprintf(' Layer 1 QSP — All States Publication Run\n');
fprintf(' %s\n', run_timestamp);
fprintf('=======================================================\n\n');

%% --- 1. Define All Sex-Disease States ------------------------------------
all_states = {
    'female_pre_healthy';
    'female_pre_T2DM';
    'female_post_healthy';
    'female_post_T2DM';
    'male_healthy';
    'male_T2DM'
};

n_states = length(all_states);

% Display labels for figures/tables
state_labels = {
    'Female Pre Healthy';
    'Female Pre T2DM';
    'Female Post Healthy';
    'Female Post T2DM';
    'Male Healthy';
    'Male T2DM'
};

% Color scheme: blues for female, reds for male; solid=healthy, dashed=T2DM
state_colors = [
    0.20  0.45  0.80;   % female pre healthy   — blue
    0.20  0.45  0.80;   % female pre T2DM      — blue
    0.55  0.75  0.95;   % female post healthy  — light blue
    0.55  0.75  0.95;   % female post T2DM     — light blue
    0.85  0.20  0.20;   % male healthy         — red
    0.85  0.20  0.20;   % male T2DM            — red
];
state_linestyle = {'-', '--', '-', '--', '-', '--'};

%% --- 2. ODE Parameters (load_parameters function embedded below) ----------
% Results storage
SS = struct();   % steady states
T  = cell(n_states,1);
Y  = cell(n_states,1);

%% --- 3. ODE Solver Settings -----------------------------------------------
tspan   = [0 80];   % hours
options = odeset('RelTol',1e-8, 'AbsTol',1e-10, 'MaxStep',0.1);

%% --- 4. Run All States ----------------------------------------------------
fprintf('Running all %d states...\n\n', n_states);

for s = 1:n_states
    sex_state = all_states{s};
    fprintf('[%d/%d] %s ... ', s, n_states, sex_state);
    fprintf(flog, '[%d/%d] %s\n', s, n_states, sex_state);

    % Load parameters for this state
    p = load_parameters(sex_state);

    % Initial conditions [E2, ERb, PGC1a, Mito, ROS, MnSOD, DeltaPsi]
    y0 = [p.E2_0; p.ERb_0; p.PGC1a_0; p.Mito_0; p.ROS_0; p.MnSOD_0; p.DeltaPsi_0];

    % Solve ODEs
    [t_sol, y_sol] = ode15s(@(t,y) layer1_odes(t,y,p), tspan, y0, options);

    T{s} = t_sol;
    Y{s} = y_sol;

    % Extract steady states (last time point)
    SS(s).state      = sex_state;
    SS(s).label      = state_labels{s};
    SS(s).E2         = y_sol(end,1);
    SS(s).ERb        = y_sol(end,2);
    SS(s).PGC1a      = y_sol(end,3);
    SS(s).Mito       = y_sol(end,4);
    SS(s).ROS        = y_sol(end,5);
    SS(s).MnSOD      = y_sol(end,6);
    SS(s).DeltaPsi   = y_sol(end,7);

    % Save Layer 2 interface
    layer2_input.PGC1a_SS     = SS(s).PGC1a;
    layer2_input.Mito_dens_SS = SS(s).Mito;
    layer2_input.ROS_SS       = SS(s).ROS;
    layer2_input.delta_psi_SS = SS(s).DeltaPsi;
    layer2_input.MnSOD_SS     = SS(s).MnSOD;
    layer2_input.sex_state    = sex_state;
    save(fullfile(outdir, ['layer1_to_layer2_' sex_state '.mat']), 'layer2_input');

    fprintf('Done.\n');
    fprintf(flog, '  PGC1a=%.4f  Mito=%.4f  ROS=%.4f  DeltaPsi=%.4f  MnSOD=%.4f\n', ...
        SS(s).PGC1a, SS(s).Mito, SS(s).ROS, SS(s).DeltaPsi, SS(s).MnSOD);

    % --- Individual time-course figure (publication quality) ---
    plot_timecourse(t_sol, y_sol, sex_state, state_labels{s}, p, outdir);
end

fprintf('\nAll states complete.\n\n');

%% --- 5. Export Summary Table ----------------------------------------------
fprintf('Exporting summary table...\n');

% Build table
varNames = {'State','Label','E2_pM','ERb_fraction','PGC1a_fold', ...
            'Mito_pct','ROS_fold','MnSOD_UmL','DeltaPsi_JC1'};

T_out = table(...
    all_states, state_labels, ...
    [SS.E2]', [SS.ERb]', [SS.PGC1a]', ...
    [SS.Mito]', [SS.ROS]', [SS.MnSOD]', [SS.DeltaPsi]', ...
    'VariableNames', varNames);

% Save CSV
csv_file = fullfile(outdir, 'Layer1_AllStates_Summary.csv');
writetable(T_out, csv_file);
fprintf('  Saved: %s\n', csv_file);

% Save Excel with formatting
xlsx_file = fullfile(outdir, 'Layer1_AllStates_Summary.xlsx');
writetable(T_out, xlsx_file, 'Sheet', 'Steady States');

% Add a second sheet with metadata
meta = {'Field','Value';
        'Script Version', script_version;
        'Run Timestamp', run_timestamp;
        'MATLAB Version', matlab_version;
        'Model', 'Layer 1 — Mitochondrial Biogenesis & ROS';
        'ODE Solver', 'ode15s';
        'RelTol', '1e-8';
        'AbsTol', '1e-10';
        'Simulation Time (hours)', '80';
        'Authors', 'Soheili M, Gilzad-Kohan H, Lotfi AS'};
writecell(meta, xlsx_file, 'Sheet', 'Metadata');
fprintf('  Saved: %s\n', xlsx_file);

%% --- 6. Publication Comparison Figure ------------------------------------
fprintf('Generating publication comparison figure...\n');

fig_comp = figure('Name','Layer1 All States Comparison',...
    'Units','inches','Position',[0.5 0.5 14 10],...
    'Color','white');

variables   = {'PGC1a','Mito','ROS','MnSOD','DeltaPsi','ERb'};
var_labels  = {'PGC-1\alpha (fold)','Mito Density (%)','ROS (fold)',...
               'MnSOD (U/mL)','\Delta\Psi_m (JC-1)','ER\beta Active (fraction)'};
n_vars = length(variables);

col_map = containers.Map({'E2','ERb','PGC1a','Mito','ROS','MnSOD','DeltaPsi'}, {1,2,3,4,5,6,7});

for v = 1:n_vars
    ax = subplot(2,3,v);
    hold on; box on;
    cidx = col_map(variables{v});
    for s = 1:n_states
        plot(T{s}, Y{s}(:,cidx), ...
            'Color', state_colors(s,:), ...
            'LineStyle', state_linestyle{s}, ...
            'LineWidth', 1.8);
    end
    xlabel('Time (hours)', 'FontSize',10);
    ylabel(var_labels{v}, 'FontSize',10);
    title(var_labels{v}, 'FontSize',11, 'FontWeight','bold');
    xlim([0 80]);
    set(gca, 'FontSize',9, 'TickDir','out', 'LineWidth',0.8);
end

% Legend
lgd = legend(state_labels, 'Location','southoutside', ...
    'Orientation','horizontal', 'FontSize',8, 'NumColumns',3);
lgd.Position = [0.1 0.01 0.8 0.04];

% Title
sgtitle({'Layer 1: Mitochondrial Biogenesis — All Sex-Disease States', ...
    'QSP Model of Diabetic Cardiomyopathy'}, ...
    'FontSize',13, 'FontWeight','bold');

% Annotation
annotation('textbox',[0 0.96 1 0.04],'String', ...
    sprintf('QSP-DCM | %s | v%s', run_timestamp, script_version), ...
    'EdgeColor','none','HorizontalAlignment','center','FontSize',7,'Color',[0.5 0.5 0.5]);

% Save high-resolution
save_figure(fig_comp, fullfile(outdir,'Layer1_ComparisonPlot'));
fprintf('  Saved: Layer1_ComparisonPlot .png\n');

%% --- 7. Bar Chart — Steady State Comparison (Manuscript Figure) -----------
fig_bar = figure('Name','Steady State Bar Comparison',...
    'Units','inches','Position',[0.5 0.5 14 9],'Color','white');

bar_vars   = {'PGC1a','Mito','ROS','MnSOD','DeltaPsi'};
bar_labels = {'PGC-1\alpha (fold)','Mito Density (%)','ROS (fold)',...
              'MnSOD (U/mL)','\Delta\Psi_m (JC-1)'};
bar_vals   = [[SS.PGC1a]; [SS.Mito]; [SS.ROS]; [SS.MnSOD]; [SS.DeltaPsi]]';

for v = 1:5
    subplot(2,3,v);
    vals = bar_vals(:,v);
    b = bar(vals, 0.7);
    b.FaceColor = 'flat';
    for s = 1:n_states
        b.CData(s,:) = state_colors(s,:);
    end
    set(gca,'XTickLabel', {'FPH','FPT2','FPoH','FPoT2','MH','MT2'}, ...
        'FontSize',9,'TickDir','out','XTickLabelRotation',30,'LineWidth',0.8);
    ylabel(bar_labels{v},'FontSize',10);
    title(bar_labels{v},'FontSize',11,'FontWeight','bold');
    box on;
end

% Legend patch
subplot(2,3,6); axis off;
patch_handles = zeros(n_states,1);
for s = 1:n_states
    patch_handles(s) = patch(NaN,NaN, state_colors(s,:), ...
        'LineStyle', state_linestyle{s}, 'LineWidth',1.5);
end
legend(patch_handles, state_labels, 'Location','best', 'FontSize',9);
title('State Legend','FontSize',11,'FontWeight','bold');

sgtitle({'Layer 1: Steady-State Comparison — All Sex-Disease States', ...
    'QSP Model of Diabetic Cardiomyopathy'}, ...
    'FontSize',13,'FontWeight','bold');

annotation('textbox',[0 0.96 1 0.04],'String', ...
    sprintf('QSP-DCM | %s | v%s', run_timestamp, script_version), ...
    'EdgeColor','none','HorizontalAlignment','center','FontSize',7,'Color',[0.5 0.5 0.5]);

save_figure(fig_bar, fullfile(outdir,'Layer1_SteadyState_Bars'));
fprintf('  Saved: Layer1_SteadyState_Bars .png\n');

%% --- 8. Print Summary to Command Window & Log ----------------------------
fprintf('\n');
fprintf('=======================================================\n');
fprintf(' LAYER 1 STEADY STATE SUMMARY — ALL STATES\n');
fprintf('=======================================================\n');
fprintf('%-22s %8s %8s %8s %8s %8s\n', ...
    'State','PGC1a','Mito%','ROS','MnSOD','DeltaPsi');
fprintf('%s\n', repmat('-',1,66));
for s = 1:n_states
    fprintf('%-22s %8.3f %8.2f %8.4f %8.2f %8.3f\n', ...
        state_labels{s}, SS(s).PGC1a, SS(s).Mito, ...
        SS(s).ROS, SS(s).MnSOD, SS(s).DeltaPsi);
end
fprintf('=======================================================\n');
fprintf(' Files saved to: %s\n', outdir);
fprintf('=======================================================\n\n');

% Mirror to log
fprintf(flog,'\n=== STEADY STATE SUMMARY ===\n');
fprintf(flog,'%-22s %8s %8s %8s %8s %8s\n','State','PGC1a','Mito%','ROS','MnSOD','DeltaPsi');
for s = 1:n_states
    fprintf(flog,'%-22s %8.3f %8.2f %8.4f %8.2f %8.3f\n', ...
        state_labels{s}, SS(s).PGC1a, SS(s).Mito, ...
        SS(s).ROS, SS(s).MnSOD, SS(s).DeltaPsi);
end
fprintf(flog,'\nRun completed: %s\n', datestr(now));
fclose(flog);

fprintf('Log saved: Layer1_RunLog.txt\n');
fprintf('Run complete. All files saved.\n\n');

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

%% --- load_parameters ------------------------------------------------------
function p = load_parameters(sex_state)
% Returns parameter struct for the specified sex-disease state.
% All parameters sourced from Layer1_Parameters.xlsx (114 HIGH confidence).
% References logged in Layer1_RunLog.txt.

    % --- Shared baseline parameters ---
    p.k_ERb_on      = 0.05;    % ERb activation rate (pM^-1 h^-1)
    p.k_ERb_off     = 0.10;    % ERb deactivation rate (h^-1)
    p.k_PGC_basal   = 0.08;    % PGC-1a basal transcription (h^-1)
    p.k_PGC_ERb     = 0.25;    % ERb→PGC-1a coupling (h^-1)
    p.k_PGC_deg     = 0.15;    % PGC-1a degradation (h^-1)
    p.k_Mito_form   = 0.12;    % Mitochondrial biogenesis (h^-1)
    p.k_Mito_deg    = 0.05;    % Mitochondrial degradation (h^-1)
    p.k_ROS_prod    = 0.20;    % ROS basal production (fold h^-1)
    p.k_ROS_mito    = 0.15;    % Mito-driven ROS production (h^-1)
    p.k_ROS_scav    = 0.30;    % ROS scavenging by MnSOD (h^-1)
    p.k_MnSOD_basal = 2.00;    % MnSOD basal synthesis (U/mL/h)
    p.k_MnSOD_ROS   = 8.00;    % ROS→MnSOD induction (U/mL/h per fold)
    p.k_MnSOD_PGC   = 5.00;    % PGC1a→MnSOD induction
    p.k_MnSOD_deg   = 0.05;    % MnSOD degradation (h^-1)
    p.k_Psi_form    = 3.00;    % DeltaPsi formation rate
    p.k_Psi_ROS     = 1.50;    % ROS-driven Psi dissipation
    p.k_Psi_deg     = 0.08;    % DeltaPsi decay (h^-1)

    % --- State-specific parameter overrides ---
    switch sex_state

        case 'female_pre_healthy'
            p.E2_ss         = 150;     % pM — premenopausal estradiol (follicular)
            p.k_ERb_sens    = 1.00;    % full receptor sensitivity
            p.T2DM_factor   = 1.00;    % no diabetes
            p.mito_stress   = 1.00;    % no stress
            p.ROS_drive     = 1.00;
            p.E2_0          = 150;
            p.ERb_0         = 0.20;
            p.PGC1a_0       = 0.50;
            p.Mito_0        = 30;
            p.ROS_0         = 1.00;
            p.MnSOD_0       = 60;
            p.DeltaPsi_0    = 1.00;

        case 'female_pre_T2DM'
            p.E2_ss         = 150;     % pM
            p.k_ERb_sens    = 0.65;    % receptor downregulation by hyperglycemia
            p.T2DM_factor   = 1.60;    % diabetes amplifies ROS/mitostress
            p.mito_stress   = 1.40;
            p.ROS_drive     = 1.80;
            p.E2_0          = 150;
            p.ERb_0         = 0.15;
            p.PGC1a_0       = 0.40;
            p.Mito_0        = 25;
            p.ROS_0         = 1.50;
            p.MnSOD_0       = 70;
            p.DeltaPsi_0    = 0.80;

        case 'female_post_healthy'
            p.E2_ss         = 15;      % pM — postmenopausal
            p.k_ERb_sens    = 0.80;    % mild receptor reduction
            p.T2DM_factor   = 1.00;
            p.mito_stress   = 1.20;    % mild age-related stress
            p.ROS_drive     = 1.20;
            p.E2_0          = 15;
            p.ERb_0         = 0.15;
            p.PGC1a_0       = 0.40;
            p.Mito_0        = 25;
            p.ROS_0         = 1.20;
            p.MnSOD_0       = 65;
            p.DeltaPsi_0    = 0.90;

        case 'female_post_T2DM'
            p.E2_ss         = 5;       % pM — very low postmenopausal+T2DM
            p.k_ERb_sens    = 0.50;    % significant receptor downregulation
            p.T2DM_factor   = 1.80;
            p.mito_stress   = 1.70;
            p.ROS_drive     = 2.20;
            p.E2_0          = 5;
            p.ERb_0         = 0.10;
            p.PGC1a_0       = 0.35;
            p.Mito_0        = 20;
            p.ROS_0         = 1.80;
            p.MnSOD_0       = 75;
            p.DeltaPsi_0    = 0.70;

        case 'male_healthy'
            p.E2_ss         = 30;      % pM — male baseline estradiol
            p.k_ERb_sens    = 0.60;    % lower ERb expression in male heart
            p.T2DM_factor   = 1.00;
            p.mito_stress   = 1.00;
            p.ROS_drive     = 1.10;    % slightly higher male baseline ROS
            p.E2_0          = 30;
            p.ERb_0         = 0.12;
            p.PGC1a_0       = 0.45;
            p.Mito_0        = 28;
            p.ROS_0         = 1.10;
            p.MnSOD_0       = 55;
            p.DeltaPsi_0    = 0.95;

        case 'male_T2DM'
            p.E2_ss         = 30;      % pM
            p.k_ERb_sens    = 0.45;    % T2DM receptor downregulation
            p.T2DM_factor   = 1.70;
            p.mito_stress   = 1.60;
            p.ROS_drive     = 2.00;
            p.E2_0          = 30;
            p.ERb_0         = 0.10;
            p.PGC1a_0       = 0.35;
            p.Mito_0        = 22;
            p.ROS_0         = 1.70;
            p.MnSOD_0       = 72;
            p.DeltaPsi_0    = 0.75;

        otherwise
            error('Unknown sex_state: %s', sex_state);
    end
end

%% --- layer1_odes ----------------------------------------------------------
function dydt = layer1_odes(~, y, p)
% 7 coupled ODEs — Layer 1 Mitochondrial Biogenesis & ROS Axis
%
%  y(1) = E2       (pM)           — free estradiol
%  y(2) = ERb      (fraction)     — active ERbeta
%  y(3) = PGC1a    (fold)         — PGC-1alpha
%  y(4) = Mito     (%)            — mitochondrial density
%  y(5) = ROS      (fold)         — reactive oxygen species
%  y(6) = MnSOD    (U/mL)         — manganese superoxide dismutase
%  y(7) = DeltaPsi (JC-1 ratio)   — mitochondrial membrane potential

    E2      = max(y(1), 0);
    ERb     = max(y(2), 0);
    PGC1a   = max(y(3), 0);
    Mito    = max(y(4), 0);
    ROS     = max(y(5), 0);
    MnSOD   = max(y(6), 0);
    DeltaPsi= max(y(7), 0);

    % --- ODE 1: E2 — clamped to steady-state (exogenous/endocrine input)
    dE2 = p.k_ERb_on * (p.E2_ss - E2);

    % --- ODE 2: ERbeta activation
    dERb = p.k_ERb_on * p.k_ERb_sens * E2 * (1 - ERb) ...
         - p.k_ERb_off * ERb;

    % --- ODE 3: PGC-1alpha
    dPGC1a = p.k_PGC_basal ...
           + p.k_PGC_ERb * ERb ...
           - p.k_PGC_deg * p.T2DM_factor * PGC1a;

    % --- ODE 4: Mitochondrial density
    dMito = p.k_Mito_form * PGC1a * 30 ...            % PGC drives biogenesis
          - p.k_Mito_deg * p.mito_stress * Mito;       % stress drives fission/mitophagy

    % --- ODE 5: ROS (net — production minus scavenging)
    ROS_prod  = p.k_ROS_prod * p.ROS_drive ...
              + p.k_ROS_mito * p.T2DM_factor * (1 / (PGC1a + 0.1));
    ROS_scav  = p.k_ROS_scav * MnSOD / 100 * ROS;
    dROS = ROS_prod - ROS_scav;

    % --- ODE 6: MnSOD (antioxidant defense)
    dMnSOD = p.k_MnSOD_basal ...
           + p.k_MnSOD_ROS * ROS ...
           + p.k_MnSOD_PGC * PGC1a ...
           - p.k_MnSOD_deg * MnSOD;

    % --- ODE 7: Mitochondrial membrane potential
    dDeltaPsi = p.k_Psi_form * Mito / 30 * PGC1a ...
              - p.k_Psi_ROS  * ROS * DeltaPsi ...
              - p.k_Psi_deg  * DeltaPsi;

    dydt = [dE2; dERb; dPGC1a; dMito; dROS; dMnSOD; dDeltaPsi];
end

%% --- plot_timecourse ------------------------------------------------------
function plot_timecourse(t, y, sex_state, label, p, outdir)
% Generates individual publication-quality time-course figure

    fig = figure('Name', label, 'Units','inches', ...
        'Position',[0.5 0.5 13 9], 'Color','white', 'Visible','off');

    var_names = {'E_2 free (pM)','ER\beta active (fraction)', ...
                 'PGC-1\alpha (fold)','Mito Density (%)','ROS (fold)', ...
                 'MnSOD (U/mL)','\Delta\Psi_m (JC-1)'};

    % Reference lines (literature values)
    % Format: {subplot_idx, value, color, style, label}
    refs = {
        3, 0.62, [0.8 0.2 0.2], '--', 'db/db (38222788)';
        3, 0.50, [0.8 0.0 0.8], '--', 'OVX (24984145)';
        4, 45.0, [0.2 0.7 0.2], '--', 'Healthy F (24984145)';
        4, 30.0, [0.9 0.5 0.0], '--', 'OVX 30% (24984145)';
        6, 80.0, [0.8 0.2 0.2], '--', 'HG T2DM (38222788)';
        7, 2.0,  [0.2 0.7 0.2], '--', 'EMPA (38222788)';
        7, 2.0,  [0.8 0.2 0.2], ':',  'HG T2DM (38222788)';
    };

    for v = 1:7
        ax = subplot(3,3,v);
        plot(t, y(:,v), 'b-', 'LineWidth',2.0, 'Color',[0.15 0.40 0.75]);
        hold on;

        % Add reference lines for this subplot
        for r = 1:size(refs,1)
            if refs{r,1} == v
                yline(refs{r,2}, refs{r,4}, refs{r,5}, ...
                    'Color', refs{r,3}, 'LineWidth',1.2, 'FontSize',7, ...
                    'LabelHorizontalAlignment','right');
            end
        end

        xlabel('Time (hours)','FontSize',9);
        ylabel(var_names{v},'FontSize',9);
        title(var_names{v},'FontSize',10,'FontWeight','bold');
        xlim([0 80]);
        set(ax,'FontSize',8,'TickDir','out','LineWidth',0.8,'Box','on');
    end

    % Steady-state summary bar (subplot 8)
    subplot(3,3,8);
    ss_vals  = [y(end,2), y(end,3), y(end,4)/100, y(end,5), y(end,6)/200, y(end,7)/30];
    ss_names = {'ER\beta','PGC1\alpha','Mito/100','ROS','MnSOD/200','\Delta\Psi/30'};
    bar(ss_vals, 0.6, 'FaceColor',[0.15 0.40 0.75], 'EdgeColor','k');
    set(gca,'XTickLabel',ss_names,'FontSize',7,'XTickLabelRotation',25,'TickDir','out');
    ylabel('Normalized SS Value','FontSize',8);
    title('Normalized Steady State','FontSize',10,'FontWeight','bold');
    box on;

    % Parameter annotation (subplot 9)
    subplot(3,3,9); axis off;
    txt = sprintf(['State: %s\n\nE2_{ss} = %.0f pM\n'...
        'k_{ERb,sens} = %.2f\nT2DM factor = %.2f\n'...
        'Mito stress = %.2f\nROS drive = %.2f'], ...
        label, p.E2_ss, p.k_ERb_sens, p.T2DM_factor, ...
        p.mito_stress, p.ROS_drive);
    text(0.05, 0.95, txt, 'Units','normalized','VerticalAlignment','top', ...
        'FontSize',8, 'FontName','Courier');
    title('State Parameters','FontSize',10,'FontWeight','bold');

    sgtitle({['Layer 1: ' label], ...
        'Mitochondrial Biogenesis & Oxidative Stress'}, ...
        'FontSize',12,'FontWeight','bold');

    annotation('textbox',[0 0.01 1 0.03],'String', ...
        'QSP Diabetic Cardiomyopathy | Layer 1', ...
        'EdgeColor','none','HorizontalAlignment','center','FontSize',7, ...
        'Color',[0.5 0.5 0.5]);

    save_figure(fig, fullfile(outdir, ['Layer1_' sex_state '_TimeCourse']));
    close(fig);
end

%% --- save_figure ----------------------------------------------------------
function save_figure(fig, basepath)
% Saves figure as PNG (300 dpi) for publication

    % PNG — 300 dpi
    print(fig, basepath, '-dpng', '-r300');

end
