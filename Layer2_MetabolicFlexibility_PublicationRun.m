%% ========================================================================
%  Layer2_MetabolicFlexibility_PublicationRun.m
%  Sex-Specific QSP Model of Diabetic Cardiomyopathy — Layer 2
%  Metabolic Flexibility: AMPK / FAO / Ceramide / Diastolic Stiffness
%
%  Authors : Soheili M, Gilzad-Kohan H, Lotfi AS
%
%  Purpose : Reads Layer 1 steady-state outputs (.mat files) for all 6
%            sex-disease states. Solves 6 coupled ODEs modeling metabolic
%            flexibility and lipotoxicity. Outputs diastolic stiffness
%            (E/e' surrogate) as the primary clinical endpoint.
%
%  ODEs    :
%    y(1) = AMPK       (fold)      — AMP-activated protein kinase
%    y(2) = FAO        (fold)      — fatty acid oxidation rate
%    y(3) = GlucOx     (fold)      — glucose oxidation rate
%    y(4) = Ceramide   (µM)        — ceramide (lipotoxicity marker)
%    y(5) = CollagenX  (fold)      — myocardial collagen crosslinking
%    y(6) = EeRatio    (ratio)     — E/e' diastolic stiffness index
%
%  Layer 1 → Layer 2 Interface (from layer1_to_layer2_<state>.mat):
%    PGC1a_SS, Mito_dens_SS, ROS_SS, delta_psi_SS, MnSOD_SS, sex_state
%
%  Version : 1.0  |  Date: 2026-03-09
%  MATLAB  : R2019b or newer required
%
%  Output files:
%    - Layer2_AllStates_Summary.csv / .xlsx
%    - Layer2_ComparisonPlot.png       (300 dpi)
%    - Layer2_DiastolicStiffness.png   (primary clinical figure)
%    - Layer2_<state>_TimeCourse.png   (individual state figures)
%    - layer2_to_layer3_<state>.mat           (Layer 3 interface)
%    - Layer2_RunLog.txt                      (reproducibility log)
%% ========================================================================

clear; close all; clc;

%% --- 0. Setup & Logging --------------------------------------------------
script_version = '1.0';
run_timestamp  = datestr(now, 'yyyy-mm-dd HH:MM:SS');
matlab_version = version;

script_path = fileparts(mfilename('fullpath'));
if isempty(script_path), script_path = pwd; end
outdir = script_path;

logfile = fullfile(outdir, 'Layer2_RunLog.txt');
flog = fopen(logfile, 'w');
fprintf(flog, '=======================================================\n');
fprintf(flog, ' Layer 2 QSP Model — Reproducibility Log\n');
fprintf(flog, '=======================================================\n');
fprintf(flog, ' Script version : %s\n', script_version);
fprintf(flog, ' Run timestamp  : %s\n', run_timestamp);
fprintf(flog, ' MATLAB version : %s\n', matlab_version);
fprintf(flog, ' Output dir     : %s\n', outdir);
fprintf(flog, '=======================================================\n\n');

fprintf('=======================================================\n');
fprintf(' Layer 2 QSP — Metabolic Flexibility Publication Run\n');
fprintf(' %s\n', run_timestamp);
fprintf('=======================================================\n\n');

%% --- 1. Define States & Visual Style ------------------------------------
all_states = {
    'female_pre_healthy';
    'female_pre_T2DM';
    'female_post_healthy';
    'female_post_T2DM';
    'male_healthy';
    'male_T2DM'
};

state_labels = {
    'Female Pre Healthy';
    'Female Pre T2DM';
    'Female Post Healthy';
    'Female Post T2DM';
    'Male Healthy';
    'Male T2DM'
};

n_states = length(all_states);

state_colors = [
    0.20  0.45  0.80;
    0.20  0.45  0.80;
    0.55  0.75  0.95;
    0.55  0.75  0.95;
    0.85  0.20  0.20;
    0.85  0.20  0.20;
];
state_linestyle = {'-','--','-','--','-','--'};

%% --- 2. ODE Solver Settings ---------------------------------------------
tspan   = [0 120];   % hours — longer window for metabolic dynamics
options = odeset('RelTol',1e-8, 'AbsTol',1e-10, 'MaxStep',0.1);

%% --- 3. Load Layer 1 Outputs & Run All States ---------------------------
SS2 = struct();
T2  = cell(n_states,1);
Y2  = cell(n_states,1);
L1  = cell(n_states,1);   % store Layer 1 inputs for reference

fprintf('Loading Layer 1 outputs and running Layer 2...\n\n');

for s = 1:n_states
    sex_state = all_states{s};
    matfile   = fullfile(outdir, ['layer1_to_layer2_' sex_state '.mat']);

    % --- Load Layer 1 interface ---
    if ~exist(matfile, 'file')
        error(['Layer 1 .mat file not found: ' matfile ...
               '\nPlease run Layer1_AllStates_PublicationRun.m first ' ...
               'and ensure .mat files are in the same folder as this script.']);
    end

    tmp = load(matfile);
    L1{s} = tmp.layer2_input;
    fprintf('[%d/%d] %s\n', s, n_states, sex_state);
    fprintf(flog, '[%d/%d] %s\n', s, n_states, sex_state);
    fprintf('  Layer 1 inputs: PGC1a=%.3f  Mito=%.1f%%  ROS=%.4f  dPsi=%.2f\n', ...
        L1{s}.PGC1a_SS, L1{s}.Mito_dens_SS, L1{s}.ROS_SS, L1{s}.delta_psi_SS);

    % --- Load Layer 2 parameters ---
    p = load_L2_parameters(sex_state, L1{s});

    % --- Initial conditions ---
    % [AMPK, FAO, GlucOx, Ceramide, CollagenX, EeRatio]
    y0 = [p.AMPK_0; p.FAO_0; p.GlucOx_0; p.Cer_0; p.ColX_0; p.Ee_0];

    % --- Solve ---
    [t_sol, y_sol] = ode15s(@(t,y) layer2_odes(t,y,p), tspan, y0, options);

    T2{s} = t_sol;
    Y2{s} = y_sol;

    % --- Steady states ---
    SS2(s).state     = sex_state;
    SS2(s).label     = state_labels{s};
    SS2(s).AMPK      = y_sol(end,1);
    SS2(s).FAO       = y_sol(end,2);
    SS2(s).GlucOx    = y_sol(end,3);
    SS2(s).MetFlex   = y_sol(end,2) / (y_sol(end,2) + y_sol(end,3) + 0.001); % FAO fraction
    SS2(s).Ceramide  = y_sol(end,4);
    SS2(s).CollagenX = y_sol(end,5);
    SS2(s).EeRatio   = y_sol(end,6);

    fprintf('  SS: AMPK=%.3f  FAO=%.3f  GlucOx=%.3f  Ceramide=%.2f  E/e''=%.2f\n', ...
        SS2(s).AMPK, SS2(s).FAO, SS2(s).GlucOx, SS2(s).Ceramide, SS2(s).EeRatio);
    fprintf(flog, '  AMPK=%.3f FAO=%.3f GlucOx=%.3f Ceramide=%.2f ColX=%.3f Ee=%.2f\n', ...
        SS2(s).AMPK, SS2(s).FAO, SS2(s).GlucOx, SS2(s).Ceramide, ...
        SS2(s).CollagenX, SS2(s).EeRatio);

    % --- Save Layer 3 interface ---
    layer3_input.AMPK_SS      = SS2(s).AMPK;
    layer3_input.FAO_SS       = SS2(s).FAO;
    layer3_input.GlucOx_SS    = SS2(s).GlucOx;
    layer3_input.MetFlex_SS   = SS2(s).MetFlex;
    layer3_input.Ceramide_SS  = SS2(s).Ceramide;
    layer3_input.CollagenX_SS = SS2(s).CollagenX;
    layer3_input.EeRatio_SS   = SS2(s).EeRatio;
    % Pass through Layer 1 values
    layer3_input.PGC1a_SS     = L1{s}.PGC1a_SS;
    layer3_input.Mito_dens_SS = L1{s}.Mito_dens_SS;
    layer3_input.ROS_SS       = L1{s}.ROS_SS;
    layer3_input.delta_psi_SS = L1{s}.delta_psi_SS;
    layer3_input.sex_state    = sex_state;
    save(fullfile(outdir, ['layer2_to_layer3_' sex_state '.mat']), 'layer3_input');

    % --- Individual figure ---
    plot_L2_timecourse(t_sol, y_sol, sex_state, state_labels{s}, p, L1{s}, outdir);

    fprintf('  Done.\n\n');
end

%% --- 4. Export Summary Table --------------------------------------------
fprintf('Exporting summary table...\n');

varNames = {'State','Label','AMPK_fold','FAO_fold','GlucOx_fold', ...
            'MetFlex_FAOfraction','Ceramide_uM','CollagenX_fold','EeRatio'};

T_out = table(...
    all_states, state_labels, ...
    [SS2.AMPK]', [SS2.FAO]', [SS2.GlucOx]', [SS2.MetFlex]', ...
    [SS2.Ceramide]', [SS2.CollagenX]', [SS2.EeRatio]', ...
    'VariableNames', varNames);

csv_file  = fullfile(outdir, 'Layer2_AllStates_Summary.csv');
xlsx_file = fullfile(outdir, 'Layer2_AllStates_Summary.xlsx');
writetable(T_out, csv_file);
fprintf('  Saved: %s\n', csv_file);

writetable(T_out, xlsx_file, 'Sheet', 'Metabolic Flexibility');
meta = {'Field','Value';
        'Script Version', script_version;
        'Run Timestamp',  run_timestamp;
        'MATLAB Version', matlab_version;
        'Layer',          '2 — Metabolic Flexibility';
        'ODEs',           'AMPK, FAO, GlucOx, Ceramide, CollagenX, E/e''';
        'ODE Solver',     'ode15s';
        'RelTol',         '1e-8';
        'AbsTol',         '1e-10';
        'Simulation Time (hours)', '120';
        'Primary Endpoint', 'E/e'' diastolic stiffness index';
        'Authors', 'Soheili M, Gilzad-Kohan H, Lotfi AS'};
writecell(meta, xlsx_file, 'Sheet', 'Metadata');
fprintf('  Saved: %s\n', xlsx_file);

%% --- 5. Publication Comparison Figure — 6 Variables ---------------------
fprintf('Generating comparison figure...\n');

fig_comp = figure('Name','Layer2 All States Comparison', ...
    'Units','inches','Position',[0.5 0.5 14 10],'Color','white');

comp_vars  = {'AMPK','FAO','GlucOx','Ceramide','CollagenX','EeRatio'};
comp_cols  = [1, 2, 3, 4, 5, 6];
comp_ylabs = {'AMPK (fold)','FAO (fold)','Glucose Oxidation (fold)', ...
              'Ceramide (\muM)','Collagen X-linking (fold)','E/e'' Ratio'};

% Reference lines [subplot, value, color, style, label]
comp_refs = {
    1,  1.0,  [0.5 0.5 0.5], '--', 'Baseline';
    4, 12.0,  [0.8 0.2 0.2], '--', 'T2DM threshold (PMID 29661658)';
    4,  6.0,  [0.2 0.7 0.2], '--', 'Normal (<6 \muM)';
    6, 15.0,  [0.8 0.2 0.2], '--', 'Diastolic dysfunction (E/e''>15)';
    6,  8.0,  [0.2 0.7 0.2], '--', 'Normal E/e'' (<8)';
};

col_map2 = containers.Map(comp_vars, {1,2,3,4,5,6});

for v = 1:6
    ax = subplot(2,3,v);
    hold on; box on;
    cidx = col_map2(comp_vars{v});
    for s = 1:n_states
        plot(T2{s}, Y2{s}(:,cidx), ...
            'Color',     state_colors(s,:), ...
            'LineStyle', state_linestyle{s}, ...
            'LineWidth', 1.8);
    end
    % Reference lines
    for r = 1:size(comp_refs,1)
        if comp_refs{r,1} == v
            yline(comp_refs{r,2}, comp_refs{r,4}, comp_refs{r,5}, ...
                'Color',comp_refs{r,3},'LineWidth',1.1,'FontSize',7, ...
                'LabelHorizontalAlignment','right');
        end
    end
    xlabel('Time (hours)','FontSize',10);
    ylabel(comp_ylabs{v},'FontSize',10);
    title(comp_ylabs{v},'FontSize',11,'FontWeight','bold');
    xlim([0 120]);
    set(ax,'FontSize',9,'TickDir','out','LineWidth',0.8);
end

lgd = legend(state_labels,'Location','southoutside', ...
    'Orientation','horizontal','FontSize',8,'NumColumns',3);
lgd.Position = [0.10 0.01 0.80 0.04];

sgtitle({'Layer 2: Metabolic Flexibility — All Sex-Disease States', ...
    'QSP Model of Diabetic Cardiomyopathy'}, ...
    'FontSize',13,'FontWeight','bold');
annotation('textbox',[0 0.96 1 0.04],'String', ...
    sprintf('QSP-DCM | %s | v%s', run_timestamp, script_version), ...
    'EdgeColor','none','HorizontalAlignment','center', ...
    'FontSize',7,'Color',[0.5 0.5 0.5]);

save_figure(fig_comp, fullfile(outdir,'Layer2_ComparisonPlot'));
fprintf('  Saved: Layer2_ComparisonPlot .png\n');

%% --- 6. PRIMARY CLINICAL FIGURE: Diastolic Stiffness + Ceramide ---------
fprintf('Generating primary clinical figure...\n');

fig_clin = figure('Name','Diastolic Stiffness & Ceramide', ...
    'Units','inches','Position',[0.5 0.5 14 11],'Color','white');

%--- Panel A: E/e' time courses ---
subplot(2,3,[1 2]);
hold on; box on;
for s = 1:n_states
    plot(T2{s}, Y2{s}(:,6), ...
        'Color',state_colors(s,:),'LineStyle',state_linestyle{s},'LineWidth',2.0);
end
yline(15,'r--','Diastolic Dysfunction (E/e''>15)','LineWidth',1.2,'FontSize',8, ...
    'LabelHorizontalAlignment','right');
yline(8,'g--','Normal E/e'' (<8)','LineWidth',1.2,'FontSize',8, ...
    'LabelHorizontalAlignment','right');
xlabel('Time (hours)','FontSize',11);
ylabel('E/e'' Ratio','FontSize',11);
title('Diastolic Stiffness Index (E/e'') — Time Course','FontSize',12,'FontWeight','bold');
legend(state_labels,'Location','northwest','FontSize',8);
xlim([0 120]); set(gca,'FontSize',10,'TickDir','out','LineWidth',0.8);

%--- Panel B: Ceramide time courses ---
subplot(2,3,[4 5]);
hold on; box on;
for s = 1:n_states
    plot(T2{s}, Y2{s}(:,4), ...
        'Color',state_colors(s,:),'LineStyle',state_linestyle{s},'LineWidth',2.0);
end
yline(12,'r--','T2DM threshold (PMID 29661658)','LineWidth',1.2,'FontSize',8, ...
    'LabelHorizontalAlignment','right');
yline(6,'g--','Normal ceramide (<6 \muM)','LineWidth',1.2,'FontSize',8, ...
    'LabelHorizontalAlignment','right');
xlabel('Time (hours)','FontSize',11);
ylabel('Ceramide (\muM)','FontSize',11);
title('Ceramide Accumulation (Lipotoxicity)','FontSize',12,'FontWeight','bold');
legend(state_labels,'Location','northwest','FontSize',8);
xlim([0 120]); set(gca,'FontSize',10,'TickDir','out','LineWidth',0.8);

%--- Panel C: Steady-state E/e' bar chart ---
subplot(2,3,3);
ee_vals = [SS2.EeRatio]';
b = bar(ee_vals, 0.7);
b.FaceColor = 'flat';
for s = 1:n_states, b.CData(s,:) = state_colors(s,:); end
yline(15,'r--','LineWidth',1.5);
yline(8,'g--','LineWidth',1.5);
set(gca,'XTickLabel',{'FPH','FPT2','FPoH','FPoT2','MH','MT2'}, ...
    'FontSize',9,'XTickLabelRotation',30,'TickDir','out');
ylabel('E/e'' at Steady State','FontSize',10);
title('Steady-State E/e''','FontSize',11,'FontWeight','bold');
box on;

%--- Panel D: Steady-state Ceramide bar chart ---
subplot(2,3,6);
cer_vals = [SS2.Ceramide]';
b2 = bar(cer_vals, 0.7);
b2.FaceColor = 'flat';
for s = 1:n_states, b2.CData(s,:) = state_colors(s,:); end
yline(12,'r--','LineWidth',1.5);
yline(6,'g--','LineWidth',1.5);
set(gca,'XTickLabel',{'FPH','FPT2','FPoH','FPoT2','MH','MT2'}, ...
    'FontSize',9,'XTickLabelRotation',30,'TickDir','out');
ylabel('Ceramide \muM at Steady State','FontSize',10);
title('Steady-State Ceramide','FontSize',11,'FontWeight','bold');
box on;

sgtitle({'Layer 2: Primary Clinical Endpoints', ...
    'Diastolic Stiffness (E/e'') & Lipotoxicity (Ceramide)'}, ...
    'FontSize',13,'FontWeight','bold');
annotation('textbox',[0 0.96 1 0.04],'String', ...
    sprintf('QSP-DCM | %s | v%s', run_timestamp, script_version), ...
    'EdgeColor','none','HorizontalAlignment','center', ...
    'FontSize',7,'Color',[0.5 0.5 0.5]);

save_figure(fig_clin, fullfile(outdir,'Layer2_DiastolicStiffness_Ceramide'));
fprintf('  Saved: Layer2_DiastolicStiffness_Ceramide .png\n');

%% --- 7. Metabolic Flexibility Radar / Summary Bar -----------------------
fig_bar = figure('Name','Layer2 Steady State Bars', ...
    'Units','inches','Position',[0.5 0.5 14 9],'Color','white');

bar_vars  = {'AMPK','FAO','GlucOx','MetFlex','Ceramide','EeRatio'};
bar_ylabs = {'AMPK (fold)','FAO (fold)','Glucose Ox (fold)', ...
             'MetFlex (FAO fraction)','Ceramide (\muM)','E/e'' Ratio'};
bar_data  = [[SS2.AMPK];[SS2.FAO];[SS2.GlucOx];[SS2.MetFlex]; ...
             [SS2.Ceramide];[SS2.EeRatio]]';

for v = 1:6
    subplot(2,3,v);
    b = bar(bar_data(:,v), 0.7);
    b.FaceColor = 'flat';
    for s = 1:n_states, b.CData(s,:) = state_colors(s,:); end
    set(gca,'XTickLabel',{'FPH','FPT2','FPoH','FPoT2','MH','MT2'}, ...
        'FontSize',9,'XTickLabelRotation',30,'TickDir','out');
    ylabel(bar_ylabs{v},'FontSize',10);
    title(bar_ylabs{v},'FontSize',11,'FontWeight','bold');
    box on;
end

sgtitle({'Layer 2: Steady-State Metabolic Profile — All States', ...
    'QSP Model of Diabetic Cardiomyopathy'}, ...
    'FontSize',13,'FontWeight','bold');
annotation('textbox',[0 0.96 1 0.04],'String', ...
    sprintf('QSP-DCM | %s | v%s', run_timestamp, script_version), ...
    'EdgeColor','none','HorizontalAlignment','center', ...
    'FontSize',7,'Color',[0.5 0.5 0.5]);

save_figure(fig_bar, fullfile(outdir,'Layer2_SteadyState_Bars'));
fprintf('  Saved: Layer2_SteadyState_Bars .png\n');

%% --- 8. Print Summary to Command Window ---------------------------------
fprintf('\n');
fprintf('=======================================================\n');
fprintf(' LAYER 2 STEADY STATE SUMMARY — ALL STATES\n');
fprintf('=======================================================\n');
fprintf('%-22s %7s %7s %7s %7s %7s %7s\n', ...
    'State','AMPK','FAO','GlucOx','Cer(uM)','ColX','E/e''');
fprintf('%s\n', repmat('-',1,73));
for s = 1:n_states
    fprintf('%-22s %7.3f %7.3f %7.3f %7.2f %7.3f %7.2f\n', ...
        state_labels{s}, SS2(s).AMPK, SS2(s).FAO, SS2(s).GlucOx, ...
        SS2(s).Ceramide, SS2(s).CollagenX, SS2(s).EeRatio);
end
fprintf('=======================================================\n');
fprintf(' Files saved to: %s\n', outdir);
fprintf('=======================================================\n\n');

fprintf(flog,'\n=== LAYER 2 STEADY STATE SUMMARY ===\n');
fprintf(flog,'%-22s %7s %7s %7s %7s %7s %7s\n', ...
    'State','AMPK','FAO','GlucOx','Cer(uM)','ColX','E/e''');
for s = 1:n_states
    fprintf(flog,'%-22s %7.3f %7.3f %7.3f %7.2f %7.3f %7.2f\n', ...
        state_labels{s}, SS2(s).AMPK, SS2(s).FAO, SS2(s).GlucOx, ...
        SS2(s).Ceramide, SS2(s).CollagenX, SS2(s).EeRatio);
end
fprintf(flog,'\nRun completed: %s\n', datestr(now));
fclose(flog);

fprintf('Log saved: Layer2_RunLog.txt\n');
fprintf('Run complete. All Layer 2 files saved.\n\n');

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

%% --- load_L2_parameters --------------------------------------------------
function p = load_L2_parameters(sex_state, L1)
% Layer 2 parameters — metabolic flexibility axis
% L1 = layer1_to_layer2 struct (PGC1a_SS, Mito_dens_SS, ROS_SS, delta_psi_SS)
%
% Key parameter sources:
%  AMPK: Canto & Auwerx 2010 (PMID 20559966); Noppe et al 2019
%  FAO:  Lopaschuk et al 2010 (PMID 20959499); Wende & Abel 2010
%  Ceramide: Chaurasia & Summers 2015 (PMID 26228190); Turpin et al 2016
%  E/e': Redfield et al 2003; Obokata et al 2017 (PMID 28716855)
%  Sex differences: Regitz-Zagrosek & Kararigas 2017 (PMID 28228494)

    % Pass Layer 1 steady states as driving inputs
    p.PGC1a   = L1.PGC1a_SS;
    p.Mito    = L1.Mito_dens_SS;
    p.ROS     = L1.ROS_SS;
    p.dPsi    = L1.delta_psi_SS;
    p.MnSOD   = L1.MnSOD_SS;

    % --- Shared kinetic parameters ---
    % AMPK
    p.k_AMPK_basal   = 0.10;   % basal activation (h^-1)
    p.k_AMPK_dPsi    = 0.08;   % ΔΨm collapse drives AMPK (energy stress)
    p.k_AMPK_ROS     = 0.06;   % ROS-driven AMPK activation
    p.k_AMPK_deg     = 0.12;   % AMPK inactivation (h^-1)
    p.dPsi_ref       = 50.0;   % reference ΔΨm (healthy female pre)

    % FAO
    p.k_FAO_AMPK     = 0.20;   % AMPK drives FAO (CPT1 upregulation)
    p.k_FAO_PGC      = 0.15;   % PGC-1a drives FAO gene expression
    p.k_FAO_deg      = 0.10;   % FAO rate decay
    p.k_FAO_cer_inh  = 0.05;   % ceramide inhibits FAO (feedback)

    % Glucose oxidation
    p.k_Gluc_basal   = 0.12;   % basal glucose oxidation
    p.k_Gluc_ins_res = 0.08;   % insulin resistance term (T2DM)
    p.k_Gluc_deg     = 0.10;

    % Ceramide
    p.k_Cer_FA       = 0.08;   % FA overflow drives ceramide synthesis
    p.k_Cer_ROS      = 0.10;   % ROS drives ceramide (de novo)
    p.k_Cer_deg      = 0.04;   % ceramide degradation (ceramidase)
    p.Cer_ref        = 6.0;    % normal ceramide (µM)

    % Collagen crosslinking
    p.k_Col_ROS      = 0.06;   % ROS drives collagen oxidation/crosslinking
    p.k_Col_cer      = 0.04;   % ceramide drives fibrosis
    p.k_Col_deg      = 0.02;   % collagen turnover
    p.Col_ref        = 1.0;    % normalized baseline

    % E/e' (diastolic stiffness)
    p.k_Ee_Col       = 0.10;   % collagen → stiffness
    p.k_Ee_Cer       = 0.06;   % ceramide → direct stiffness
    p.k_Ee_dPsi      = 0.04;   % energy failure → diastolic dysfunction
    p.k_Ee_deg       = 0.03;   % recovery/remodeling

    % --- State-specific overrides ---
    switch sex_state

        case 'female_pre_healthy'
            p.T2DM_ins_res   = 1.00;
            p.FA_overflow    = 1.00;
            p.estrogen_prot  = 1.00;   % full estrogen protection on ceramide/fibrosis
            p.AMPK_0  = 1.0;
            p.FAO_0   = 1.0;
            p.GlucOx_0= 1.0;
            p.Cer_0   = 4.0;    % µM — normal
            p.ColX_0  = 1.0;
            p.Ee_0    = 7.0;    % normal E/e'

        case 'female_pre_T2DM'
            p.T2DM_ins_res   = 2.20;   % severe insulin resistance
            p.FA_overflow    = 1.80;   % lipid overflow despite estrogen
            p.estrogen_prot  = 0.65;   % partial — receptor downregulation
            p.AMPK_0  = 0.9;
            p.FAO_0   = 1.2;
            p.GlucOx_0= 0.7;
            p.Cer_0   = 7.0;
            p.ColX_0  = 1.2;
            p.Ee_0    = 10.0;

        case 'female_post_healthy'
            p.T2DM_ins_res   = 1.30;   % mild age-related insulin resistance
            p.FA_overflow    = 1.30;
            p.estrogen_prot  = 0.70;   % reduced — low E2
            p.AMPK_0  = 0.95;
            p.FAO_0   = 1.0;
            p.GlucOx_0= 0.9;
            p.Cer_0   = 5.5;
            p.ColX_0  = 1.1;
            p.Ee_0    = 9.0;

        case 'female_post_T2DM'
            p.T2DM_ins_res   = 2.60;   % worst insulin resistance
            p.FA_overflow    = 2.20;   % maximal lipid overflow
            p.estrogen_prot  = 0.35;   % minimal — low E2 + downregulated receptor
            p.AMPK_0  = 0.7;
            p.FAO_0   = 1.5;   % paradoxically high FAO but inefficient
            p.GlucOx_0= 0.5;
            p.Cer_0   = 10.0;
            p.ColX_0  = 1.5;
            p.Ee_0    = 13.0;

        case 'male_healthy'
            p.T2DM_ins_res   = 1.10;
            p.FA_overflow    = 1.10;
            p.estrogen_prot  = 0.60;   % lower ERb expression in males
            p.AMPK_0  = 1.0;
            p.FAO_0   = 1.0;
            p.GlucOx_0= 1.0;
            p.Cer_0   = 4.5;
            p.ColX_0  = 1.0;
            p.Ee_0    = 7.5;

        case 'male_T2DM'
            p.T2DM_ins_res   = 2.10;
            p.FA_overflow    = 1.70;
            p.estrogen_prot  = 0.45;
            p.AMPK_0  = 0.8;
            p.FAO_0   = 1.3;
            p.GlucOx_0= 0.65;
            p.Cer_0   = 8.5;
            p.ColX_0  = 1.3;
            p.Ee_0    = 11.5;

        otherwise
            error('Unknown sex_state: %s', sex_state);
    end
end

%% --- layer2_odes ---------------------------------------------------------
function dydt = layer2_odes(~, y, p)
% 6 coupled ODEs — Layer 2 Metabolic Flexibility
%
%  y(1) = AMPK       — energy sensor (fold)
%  y(2) = FAO        — fatty acid oxidation (fold)
%  y(3) = GlucOx     — glucose oxidation (fold)
%  y(4) = Ceramide   — lipotoxicity marker (µM)
%  y(5) = CollagenX  — myocardial collagen crosslinking (fold)
%  y(6) = EeRatio    — diastolic stiffness E/e' surrogate

    AMPK    = max(y(1), 0);
    FAO     = max(y(2), 0);
    GlucOx  = max(y(3), 0);
    Cer     = max(y(4), 0);
    ColX    = max(y(5), 0);
    Ee      = max(y(6), 0);

    % Normalized ΔΨm: collapse (low dPsi) = energy stress → activates AMPK
    dPsi_norm = max(1 - p.dPsi / p.dPsi_ref, 0);

    % --- ODE 1: AMPK ---
    % Activated by: energy stress (low ΔΨm), ROS, low PGC-1a
    dAMPK = p.k_AMPK_basal ...
           + p.k_AMPK_dPsi * dPsi_norm ...
           + p.k_AMPK_ROS  * p.ROS ...
           - p.k_AMPK_deg  * AMPK;

    % --- ODE 2: FAO ---
    % Driven by AMPK (CPT1) and PGC-1a; inhibited by ceramide-mediated
    % incomplete oxidation; upregulated in T2DM (FA overflow)
    dFAO = p.k_FAO_AMPK * AMPK * p.FA_overflow ...
          + p.k_FAO_PGC  * p.PGC1a ...
          - p.k_FAO_deg  * FAO ...
          - p.k_FAO_cer_inh * Cer * FAO;

    % --- ODE 3: Glucose oxidation ---
    % Suppressed by T2DM insulin resistance; competes with FAO (Randle cycle)
    dGlucOx = p.k_Gluc_basal ...
            - p.k_Gluc_ins_res * p.T2DM_ins_res * GlucOx ...
            - p.k_Gluc_deg * GlucOx * FAO / (FAO + 0.5);

    % --- ODE 4: Ceramide ---
    % De novo synthesis from FA overflow + ROS; degraded by ceramidase
    % Estrogen suppresses ceramide synthesis (Turpin 2016)
    dCer = p.k_Cer_FA  * p.FA_overflow  * (1 - p.estrogen_prot * 0.4) ...
          + p.k_Cer_ROS * p.ROS ...
          - p.k_Cer_deg * Cer;

    % --- ODE 5: Collagen crosslinking ---
    % ROS-driven oxidative crosslinking + ceramide-driven fibrotic signaling
    dColX = p.k_Col_ROS * p.ROS * (1 - p.estrogen_prot * 0.3) ...
           + p.k_Col_cer * Cer / 10 ...
           - p.k_Col_deg * ColX;

    % --- ODE 6: E/e' (diastolic stiffness) ---
    % Driven by collagen stiffness, ceramide-related myocyte dysfunction,
    % and bioenergetic failure (low ΔΨm impairs active relaxation/SERCA)
    dEe = p.k_Ee_Col  * ColX ...
        + p.k_Ee_Cer  * Cer / 10 ...
        + p.k_Ee_dPsi * dPsi_norm * 5 ...
        - p.k_Ee_deg  * Ee;

    dydt = [dAMPK; dFAO; dGlucOx; dCer; dColX; dEe];
end

%% --- plot_L2_timecourse --------------------------------------------------
function plot_L2_timecourse(t, y, sex_state, label, p, L1, outdir)

    fig = figure('Name', label, 'Units','inches', ...
        'Position',[0.5 0.5 14 10], 'Color','white', 'Visible','off');

    var_names = {'AMPK (fold)','FAO (fold)','Glucose Ox (fold)', ...
                 'Ceramide (\muM)','Collagen X-link (fold)','E/e'' Ratio'};

    ref_lines = {
        4, 12.0, [0.8 0.2 0.2], '--', 'T2DM threshold';
        4,  6.0, [0.2 0.7 0.2], '--', 'Normal (<6\muM)';
        6, 15.0, [0.8 0.2 0.2], '--', 'Diastolic dysfn (>15)';
        6,  8.0, [0.2 0.7 0.2], '--', 'Normal E/e'' (<8)';
        5,  1.5, [0.8 0.2 0.2], '--', 'Fibrosis threshold';
    };

    for v = 1:6
        ax = subplot(3,3,v);
        plot(t, y(:,v), 'b-', 'LineWidth',2.0, 'Color',[0.15 0.40 0.75]);
        hold on;
        for r = 1:size(ref_lines,1)
            if ref_lines{r,1} == v
                yline(ref_lines{r,2}, ref_lines{r,4}, ref_lines{r,5}, ...
                    'Color',ref_lines{r,3},'LineWidth',1.2,'FontSize',7, ...
                    'LabelHorizontalAlignment','right');
            end
        end
        xlabel('Time (hours)','FontSize',9);
        ylabel(var_names{v},'FontSize',9);
        title(var_names{v},'FontSize',10,'FontWeight','bold');
        xlim([0 120]);
        set(ax,'FontSize',8,'TickDir','out','LineWidth',0.8,'Box','on');
    end

    % Subplot 7: MetFlex ratio (FAO / total)
    subplot(3,3,7);
    metflex = y(:,2) ./ (y(:,2) + y(:,3) + 0.001);
    plot(t, metflex, 'Color',[0.15 0.40 0.75],'LineWidth',2.0);
    yline(0.7,'r--','FAO dominant (>70%)','FontSize',7,'LineWidth',1.1);
    yline(0.4,'g--','Balanced (40-70%)','FontSize',7,'LineWidth',1.1);
    xlabel('Time (hours)','FontSize',9);
    ylabel('FAO fraction','FontSize',9);
    title('Metabolic Flexibility Index','FontSize',10,'FontWeight','bold');
    xlim([0 120]); box on; set(gca,'TickDir','out');

    % Subplot 8: Layer 1 inputs summary bar
    subplot(3,3,8);
    l1_vals  = [L1.PGC1a_SS/2.2, L1.Mito_dens_SS/153, ...
                L1.ROS_SS, L1.delta_psi_SS/64];
    l1_names = {'PGC1\alpha/2.2','Mito/153','ROS','dPsi/64'};
    bar(l1_vals, 0.6, 'FaceColor',[0.4 0.6 0.8],'EdgeColor','k');
    set(gca,'XTickLabel',l1_names,'FontSize',8,'XTickLabelRotation',20,'TickDir','out');
    ylabel('Normalized L1 Input','FontSize',8);
    title('Layer 1 Inputs (normalized)','FontSize',10,'FontWeight','bold');
    yline(1,'k--','Reference','FontSize',7); box on;

    % Subplot 9: Parameter annotation
    subplot(3,3,9); axis off;
    txt = sprintf(['State: %s\n\nT2DM ins. res. = %.2f\n'...
        'FA overflow = %.2f\nEstrogen prot = %.2f\n\n'...
        'SS E/e'' = %.2f\nSS Ceramide = %.2f \muM\n'...
        'MetFlex = %.2f'], ...
        label, p.T2DM_ins_res, p.FA_overflow, p.estrogen_prot, ...
        y(end,6), y(end,4), y(end,2)/(y(end,2)+y(end,3)+0.001));
    text(0.05, 0.95, txt, 'Units','normalized','VerticalAlignment','top', ...
        'FontSize',8,'FontName','Courier');
    title('State Parameters','FontSize',10,'FontWeight','bold');

    sgtitle({['Layer 2: ' label], ...
        'Metabolic Flexibility & Diastolic Stiffness'}, ...
        'FontSize',12,'FontWeight','bold');
    annotation('textbox',[0 0.01 1 0.03],'String', ...
        'QSP Diabetic Cardiomyopathy | Layer 2', ...
        'EdgeColor','none','HorizontalAlignment','center', ...
        'FontSize',7,'Color',[0.5 0.5 0.5]);

    save_figure(fig, fullfile(outdir, ['Layer2_' sex_state '_TimeCourse']));
    close(fig);
end

%% --- save_figure ---------------------------------------------------------
function save_figure(fig, basepath)
    print(fig, basepath, '-dpng', '-r300');
end
