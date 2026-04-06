%% ========================================================================
%  Layer3_Microvascular_PublicationRun.m
%  Sex-Specific QSP Model of Diabetic Cardiomyopathy — Layer 3
%  Coronary Microvascular Dysfunction (CMD)
%
%  Authors : Soheili M, Gilzad-Kohan H, Lotfi AS
%
%  Purpose : Reads Layer 2 steady-state outputs (.mat files) for all 6
%            sex-disease states. Solves 6 coupled ODEs modeling coronary
%            microvascular dysfunction. Outputs CFR, IMR, and integrated
%            clinical phenotype as primary endpoints.
%
%  ODEs    :
%    y(1) = NO         (µM)     — nitric oxide (vasodilatory signal)
%    y(2) = ET1        (fold)   — endothelin-1 (vasoconstrictive)
%    y(3) = MicroRes   (fold)   — microvascular resistance
%    y(4) = CFR        (ratio)  — coronary flow reserve
%    y(5) = IMR        (U)      — index of microcirculatory resistance
%    y(6) = WallStress (fold)   — myocardial wall stress (integrates L1+L2)
%
%  Clinical Reference Thresholds (ACC/AHA/SCAI):
%    CFR  < 2.0  = microvascular dysfunction (normal ≥ 2.5)
%    IMR  > 25 U = microvascular dysfunction (normal < 25 U)
%    IMR  > 40 U = severe CMD
%
%  Layer 2 → Layer 3 Interface (from layer2_to_layer3_<state>.mat):
%    AMPK_SS, FAO_SS, GlucOx_SS, MetFlex_SS, Ceramide_SS,
%    CollagenX_SS, EeRatio_SS, PGC1a_SS, Mito_dens_SS,
%    ROS_SS, delta_psi_SS, sex_state
%
%  Version : 1.0  |  Date: 2026-03-09
%  MATLAB  : R2019b or newer required
%
%  Output files:
%    - Layer3_AllStates_Summary.csv / .xlsx
%    - Layer3_ComparisonPlot.png              (300 dpi)
%    - Layer3_CFR_IMR_Clinical.png            (primary clinical figure)
%    - Layer3_IntegratedPhenotype.png         (3-layer synthesis figure)
%    - Layer3_<state>_TimeCourse.png          (individual state figures)
%    - Layer3_RunLog.txt                          (reproducibility log)
%
%  Key references:
%    Pepine et al. JAMA 2011 (PMID 21791690) — CFR sex differences
%    Taqueti & Di Carli 2018 (JACC PMID 30213362) — CMD clinical thresholds
%    Obokata et al. 2017 (JACC PMID 28716855) — HFpEF microvascular
%    Bugiardini et al. 2005 (JAMA PMID 15741530) — sex-specific CMD
%    Regitz-Zagrosek & Kararigas 2017 (PMID 28228494) — sex differences
%    Lotfi et al. SCAI 2014/2018 consensus — CFR/IMR thresholds
%% ========================================================================

clear; close all; clc;

%% --- 0. Setup & Logging --------------------------------------------------
script_version = '1.0';
run_timestamp  = datestr(now, 'yyyy-mm-dd HH:MM:SS');
matlab_version = version;

script_path = fileparts(mfilename('fullpath'));
if isempty(script_path), script_path = pwd; end
outdir = script_path;

logfile = fullfile(outdir, 'Layer3_RunLog.txt');
flog = fopen(logfile, 'w');
fprintf(flog, '=======================================================\n');
fprintf(flog, ' Layer 3 QSP Model — Reproducibility Log\n');
fprintf(flog, '=======================================================\n');
fprintf(flog, ' Script version : %s\n', script_version);
fprintf(flog, ' Run timestamp  : %s\n', run_timestamp);
fprintf(flog, ' MATLAB version : %s\n', matlab_version);
fprintf(flog, ' Output dir     : %s\n', outdir);
fprintf(flog, '=======================================================\n\n');

fprintf('=======================================================\n');
fprintf(' Layer 3 QSP — Microvascular CMD Publication Run\n');
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
tspan   = [0 168];   % hours — 7 days for vascular remodeling dynamics
options = odeset('RelTol',1e-8, 'AbsTol',1e-10, 'MaxStep',0.1);

%% --- 3. Load Layer 2 Outputs & Run All States ---------------------------
SS3 = struct();
T3  = cell(n_states,1);
Y3  = cell(n_states,1);
L2  = cell(n_states,1);

fprintf('Loading Layer 2 outputs and running Layer 3...\n\n');

for s = 1:n_states
    sex_state = all_states{s};
    matfile   = fullfile(outdir, ['layer2_to_layer3_' sex_state '.mat']);

    if ~exist(matfile, 'file')
        error(['Layer 2 .mat file not found: ' matfile ...
               '\nPlease run Layer2_MetabolicFlexibility_PublicationRun.m first.']);
    end

    tmp = load(matfile);
    L2{s} = tmp.layer3_input;

    fprintf('[%d/%d] %s\n', s, n_states, sex_state);
    fprintf(flog, '[%d/%d] %s\n', s, n_states, sex_state);
    fprintf('  L2 inputs: Ceramide=%.2f  ColX=%.3f  E/e''=%.2f  ROS=%.4f  dPsi=%.2f\n', ...
        L2{s}.Ceramide_SS, L2{s}.CollagenX_SS, L2{s}.EeRatio_SS, ...
        L2{s}.ROS_SS, L2{s}.delta_psi_SS);

    % Load Layer 3 parameters
    p = load_L3_parameters(sex_state, L2{s});

    % Initial conditions [NO, ET1, MicroRes, CFR, IMR, WallStress]
    y0 = [p.NO_0; p.ET1_0; p.MRes_0; p.CFR_0; p.IMR_0; p.WS_0];

    % Solve ODEs
    [t_sol, y_sol] = ode15s(@(t,y) layer3_odes(t,y,p), tspan, y0, options);

    T3{s} = t_sol;
    Y3{s} = y_sol;

    % Steady states
    SS3(s).state      = sex_state;
    SS3(s).label      = state_labels{s};
    SS3(s).NO         = y_sol(end,1);
    SS3(s).ET1        = y_sol(end,2);
    SS3(s).MicroRes   = y_sol(end,3);
    SS3(s).CFR        = y_sol(end,4);
    SS3(s).IMR        = y_sol(end,5);
    SS3(s).WallStress = y_sol(end,6);
    % Derived: CMD severity score (composite)
    SS3(s).CMD_score  = (SS3(s).IMR/25) * (2.5/max(SS3(s).CFR,0.1));

    fprintf('  SS: NO=%.3f  ET1=%.3f  CFR=%.3f  IMR=%.1f  WallStress=%.3f\n', ...
        SS3(s).NO, SS3(s).ET1, SS3(s).CFR, SS3(s).IMR, SS3(s).WallStress);
    fprintf(flog, '  NO=%.3f ET1=%.3f MRes=%.3f CFR=%.3f IMR=%.1f WS=%.3f CMD=%.3f\n', ...
        SS3(s).NO, SS3(s).ET1, SS3(s).MicroRes, SS3(s).CFR, ...
        SS3(s).IMR, SS3(s).WallStress, SS3(s).CMD_score);

    % Individual figure
    plot_L3_timecourse(t_sol, y_sol, sex_state, state_labels{s}, p, L2{s}, outdir);

    fprintf('  Done.\n\n');
end

%% --- 4. Export Summary Table --------------------------------------------
fprintf('Exporting summary table...\n');

varNames = {'State','Label','NO_uM','ET1_fold','MicroRes_fold', ...
            'CFR','IMR_U','WallStress_fold','CMD_score'};

T_out = table(...
    all_states, state_labels, ...
    [SS3.NO]', [SS3.ET1]', [SS3.MicroRes]', ...
    [SS3.CFR]', [SS3.IMR]', [SS3.WallStress]', [SS3.CMD_score]', ...
    'VariableNames', varNames);

csv_file  = fullfile(outdir, 'Layer3_AllStates_Summary.csv');
xlsx_file = fullfile(outdir, 'Layer3_AllStates_Summary.xlsx');
writetable(T_out, csv_file);
fprintf('  Saved: %s\n', csv_file);

writetable(T_out, xlsx_file, 'Sheet', 'Microvascular CMD');

% Combined 3-layer table
varNames_all = {'State','Label',...
    'L1_PGC1a','L1_Mito_pct','L1_DeltaPsi',...
    'L2_Ceramide','L2_CollagenX','L2_Ee',...
    'L3_CFR','L3_IMR','L3_CMD_score'};
l1_pgc_v  = cellfun(@(x) x.PGC1a_SS,     L2);  l1_pgc_v  = l1_pgc_v(:);
l1_mito_v = cellfun(@(x) x.Mito_dens_SS, L2);  l1_mito_v = l1_mito_v(:);
l1_dpsi_v = cellfun(@(x) x.delta_psi_SS, L2);  l1_dpsi_v = l1_dpsi_v(:);
l2_cer_v  = cellfun(@(x) x.Ceramide_SS,  L2);  l2_cer_v  = l2_cer_v(:);
l2_colx_v = cellfun(@(x) x.CollagenX_SS, L2);  l2_colx_v = l2_colx_v(:);
l2_ee_v   = cellfun(@(x) x.EeRatio_SS,   L2);  l2_ee_v   = l2_ee_v(:);
T_combined = table(...
    all_states, state_labels, ...
    l1_pgc_v, l1_mito_v, l1_dpsi_v, ...
    l2_cer_v, l2_colx_v, l2_ee_v, ...
    [SS3.CFR]', [SS3.IMR]', [SS3.CMD_score]', ...
    'VariableNames', varNames_all);
writetable(T_combined, xlsx_file, 'Sheet', '3-Layer Combined');

meta = {'Field','Value';
        'Script Version',   script_version;
        'Run Timestamp',    run_timestamp;
        'MATLAB Version',   matlab_version;
        'Layer',            '3 — Coronary Microvascular Dysfunction';
        'ODEs',             'NO, ET1, MicroRes, CFR, IMR, WallStress';
        'ODE Solver',       'ode15s';
        'RelTol',           '1e-8';
        'AbsTol',           '1e-10';
        'Simulation Time (hours)', '168 (7 days)';
        'Primary Endpoints','CFR, IMR';
        'CFR normal',       '>= 2.5';
        'CFR dysfunction',  '< 2.0';
        'IMR normal',       '< 25 U';
        'IMR dysfunction',  '> 25 U';
        'IMR severe',       '> 40 U';
        'Authors', 'Soheili M, Gilzad-Kohan H, Lotfi AS'};
writecell(meta, xlsx_file, 'Sheet', 'Metadata');
fprintf('  Saved: %s\n', xlsx_file);

%% --- 5. Publication Comparison Figure -----------------------------------
fprintf('Generating comparison figure...\n');

fig_comp = figure('Name','Layer3 All States Comparison', ...
    'Units','inches','Position',[0.5 0.5 14 10],'Color','white');

comp_vars  = {'NO','ET1','MicroRes','CFR','IMR','WallStress'};
comp_cols  = {1, 2, 3, 4, 5, 6};
comp_ylabs = {'NO (\muM)','Endothelin-1 (fold)','Microvascular Resistance (fold)',...
              'CFR (ratio)','IMR (U)','Wall Stress (fold)'};

comp_refs = {
    4, 2.5, [0.2 0.7 0.2], '--', 'Normal CFR (\geq2.5)';
    4, 2.0, [0.8 0.2 0.2], '--', 'CMD threshold (<2.0)';
    5, 25,  [0.8 0.2 0.2], '--', 'IMR dysfunction (>25U)';
    5, 40,  [0.6 0.0 0.0], '--', 'Severe CMD (>40U)';
};

col_map3 = containers.Map(comp_vars, {1,2,3,4,5,6});

for v = 1:6
    ax = subplot(2,3,v);
    hold on; box on;
    cidx = col_map3(comp_vars{v});
    for s = 1:n_states
        plot(T3{s}, Y3{s}(:,cidx), ...
            'Color',     state_colors(s,:), ...
            'LineStyle', state_linestyle{s}, ...
            'LineWidth', 1.8);
    end
    for r = 1:size(comp_refs,1)
        if comp_refs{r,1} == v
            yline(comp_refs{r,2}, comp_refs{r,4}, comp_refs{r,5}, ...
                'Color',comp_refs{r,3},'LineWidth',1.2,'FontSize',7,...
                'LabelHorizontalAlignment','right');
        end
    end
    xlabel('Time (hours)','FontSize',10);
    ylabel(comp_ylabs{v},'FontSize',10);
    title(comp_ylabs{v},'FontSize',11,'FontWeight','bold');
    xlim([0 168]);
    set(ax,'FontSize',9,'TickDir','out','LineWidth',0.8);
end

lgd = legend(state_labels,'Location','southoutside', ...
    'Orientation','horizontal','FontSize',8,'NumColumns',3);
lgd.Position = [0.10 0.01 0.80 0.04];

sgtitle({'Layer 3: Coronary Microvascular Dysfunction — All Sex-Disease States',...
    'QSP Model of Diabetic Cardiomyopathy'},...
    'FontSize',13,'FontWeight','bold');
annotation('textbox',[0 0.96 1 0.04],'String',...
    sprintf('QSP-DCM | %s | v%s', run_timestamp, script_version),...
    'EdgeColor','none','HorizontalAlignment','center',...
    'FontSize',7,'Color',[0.5 0.5 0.5]);

save_figure(fig_comp, fullfile(outdir,'Layer3_ComparisonPlot'));
fprintf('  Saved: Layer3_ComparisonPlot .png\n');

%% --- 6. PRIMARY CLINICAL FIGURE: CFR + IMR --------------------------------
fprintf('Generating primary clinical figure (CFR + IMR)...\n');

fig_clin = figure('Name','CFR and IMR Clinical Figure',...
    'Units','inches','Position',[0.5 0.5 14 11],'Color','white');

% Panel A: CFR time courses
subplot(2,3,[1 2]);
hold on; box on;
for s = 1:n_states
    plot(T3{s}, Y3{s}(:,4),...
        'Color',state_colors(s,:),'LineStyle',state_linestyle{s},'LineWidth',2.0);
end
yline(2.5,'g--','Normal CFR (\geq2.5)','LineWidth',1.4,'FontSize',9,...
    'LabelHorizontalAlignment','right');
yline(2.0,'r--','CMD Threshold (<2.0)','LineWidth',1.4,'FontSize',9,...
    'LabelHorizontalAlignment','right');
xlabel('Time (hours)','FontSize',11);
ylabel('Coronary Flow Reserve (CFR)','FontSize',11);
title('Coronary Flow Reserve — Time Course','FontSize',12,'FontWeight','bold');
legend(state_labels,'Location','northeast','FontSize',8);
xlim([0 168]); set(gca,'FontSize',10,'TickDir','out','LineWidth',0.8);

% Panel B: IMR time courses
subplot(2,3,[4 5]);
hold on; box on;
for s = 1:n_states
    plot(T3{s}, Y3{s}(:,5),...
        'Color',state_colors(s,:),'LineStyle',state_linestyle{s},'LineWidth',2.0);
end
yline(25,'r--','IMR Dysfunction (>25U)','LineWidth',1.4,'FontSize',9,...
    'LabelHorizontalAlignment','right');
yline(40,'Color',[0.6 0 0],'LineStyle','--','LineWidth',1.4);
text(160, 41.5,'Severe CMD (>40U)','FontSize',9,'Color',[0.6 0 0],...
    'HorizontalAlignment','right');
xlabel('Time (hours)','FontSize',11);
ylabel('IMR (U)','FontSize',11);
title('Index of Microcirculatory Resistance (IMR)','FontSize',12,'FontWeight','bold');
legend(state_labels,'Location','northeast','FontSize',8);
xlim([0 168]); set(gca,'FontSize',10,'TickDir','out','LineWidth',0.8);

% Panel C: Steady-state CFR bar
subplot(2,3,3);
cfr_vals = [SS3.CFR]';
b1 = bar(cfr_vals, 0.7); b1.FaceColor = 'flat';
for s = 1:n_states, b1.CData(s,:) = state_colors(s,:); end
yline(2.5,'g--','LineWidth',1.5);
yline(2.0,'r--','LineWidth',1.5);
set(gca,'XTickLabel',{'FPH','FPT2','FPoH','FPoT2','MH','MT2'},...
    'FontSize',9,'XTickLabelRotation',30,'TickDir','out');
ylabel('CFR at Steady State','FontSize',10);
title('Steady-State CFR','FontSize',11,'FontWeight','bold');
box on;

% Panel D: Steady-state IMR bar
subplot(2,3,6);
imr_vals = [SS3.IMR]';
b2 = bar(imr_vals, 0.7); b2.FaceColor = 'flat';
for s = 1:n_states, b2.CData(s,:) = state_colors(s,:); end
yline(25,'r--','LineWidth',1.5);
yline(40,'Color',[0.6 0 0],'LineStyle','--','LineWidth',1.5);
set(gca,'XTickLabel',{'FPH','FPT2','FPoH','FPoT2','MH','MT2'},...
    'FontSize',9,'XTickLabelRotation',30,'TickDir','out');
ylabel('IMR (U) at Steady State','FontSize',10);
title('Steady-State IMR','FontSize',11,'FontWeight','bold');
box on;

sgtitle('Layer 3: Primary Clinical Endpoints — CFR & IMR',...
    'FontSize',13,'FontWeight','bold');
annotation('textbox',[0 0.96 1 0.04],'String',...
    sprintf('QSP-DCM | %s | v%s', run_timestamp, script_version),...
    'EdgeColor','none','HorizontalAlignment','center',...
    'FontSize',7,'Color',[0.5 0.5 0.5]);

save_figure(fig_clin, fullfile(outdir,'Layer3_CFR_IMR_Clinical'));
fprintf('  Saved: Layer3_CFR_IMR_Clinical .png\n');

%% --- 7. INTEGRATED 3-LAYER PHENOTYPE FIGURE -----------------------------
fprintf('Generating integrated 3-layer phenotype figure...\n');

fig_int = figure('Name','Integrated 3-Layer Phenotype',...
    'Units','inches','Position',[0.5 0.5 16 12],'Color','white');

% Data for all 3 layers
l1_pgc   = cellfun(@(x) x.PGC1a_SS,    L2);
l1_mito  = cellfun(@(x) x.Mito_dens_SS, L2);
l1_dpsi  = cellfun(@(x) x.delta_psi_SS, L2);
l2_cer   = cellfun(@(x) x.Ceramide_SS,  L2);
l2_colx  = cellfun(@(x) x.CollagenX_SS, L2);
l2_ee    = cellfun(@(x) x.EeRatio_SS,   L2);
l3_cfr   = [SS3.CFR]';
l3_imr   = [SS3.IMR]';
l3_cmd   = [SS3.CMD_score]';

x_pos = 1:n_states;

% Row 1: Layer 1 variables
subplot(3,4,1);
b = bar(l1_pgc,0.7); b.FaceColor='flat';
for s=1:n_states, b.CData(s,:)=state_colors(s,:); end
set(gca,'XTickLabel',{'FPH','FPT2','FPoH','FPoT2','MH','MT2'},...
    'FontSize',8,'XTickLabelRotation',30,'TickDir','out');
ylabel('PGC-1\alpha (fold)','FontSize',9); title('L1: PGC-1\alpha','FontSize',10,'FontWeight','bold'); box on;

subplot(3,4,2);
b = bar(l1_mito,0.7); b.FaceColor='flat';
for s=1:n_states, b.CData(s,:)=state_colors(s,:); end
set(gca,'XTickLabel',{'FPH','FPT2','FPoH','FPoT2','MH','MT2'},...
    'FontSize',8,'XTickLabelRotation',30,'TickDir','out');
ylabel('Mito Density (%)','FontSize',9); title('L1: Mito Density','FontSize',10,'FontWeight','bold'); box on;

subplot(3,4,3);
b = bar(l1_dpsi,0.7); b.FaceColor='flat';
for s=1:n_states, b.CData(s,:)=state_colors(s,:); end
set(gca,'XTickLabel',{'FPH','FPT2','FPoH','FPoT2','MH','MT2'},...
    'FontSize',8,'XTickLabelRotation',30,'TickDir','out');
ylabel('\Delta\Psi_m (JC-1)','FontSize',9); title('L1: \Delta\Psi_m','FontSize',10,'FontWeight','bold'); box on;

% Flow arrow annotation (L1→L2)
subplot(3,4,4); axis off;
text(0.5,0.6,'Layer 1 \rightarrow Layer 2','HorizontalAlignment','center',...
    'FontSize',12,'FontWeight','bold','Color',[0.3 0.3 0.7]);
text(0.5,0.4,'Mitochondrial failure','HorizontalAlignment','center',...
    'FontSize',9,'Color',[0.4 0.4 0.4]);
text(0.5,0.3,'drives metabolic','HorizontalAlignment','center',...
    'FontSize',9,'Color',[0.4 0.4 0.4]);
text(0.5,0.2,'inflexibility','HorizontalAlignment','center',...
    'FontSize',9,'Color',[0.4 0.4 0.4]);

% Row 2: Layer 2 variables
subplot(3,4,5);
b = bar(l2_cer,0.7); b.FaceColor='flat';
for s=1:n_states, b.CData(s,:)=state_colors(s,:); end
yline(6,'r--','LineWidth',1.2);
set(gca,'XTickLabel',{'FPH','FPT2','FPoH','FPoT2','MH','MT2'},...
    'FontSize',8,'XTickLabelRotation',30,'TickDir','out');
ylabel('Ceramide (\muM)','FontSize',9); title('L2: Ceramide','FontSize',10,'FontWeight','bold'); box on;

subplot(3,4,6);
b = bar(l2_colx,0.7); b.FaceColor='flat';
for s=1:n_states, b.CData(s,:)=state_colors(s,:); end
yline(1.5,'r--','LineWidth',1.2);
set(gca,'XTickLabel',{'FPH','FPT2','FPoH','FPoT2','MH','MT2'},...
    'FontSize',8,'XTickLabelRotation',30,'TickDir','out');
ylabel('Collagen X-link (fold)','FontSize',9); title('L2: Collagen X-link','FontSize',10,'FontWeight','bold'); box on;

subplot(3,4,7);
b = bar(l2_ee,0.7); b.FaceColor='flat';
for s=1:n_states, b.CData(s,:)=state_colors(s,:); end
yline(15,'r--','LineWidth',1.2); yline(8,'g--','LineWidth',1.2);
set(gca,'XTickLabel',{'FPH','FPT2','FPoH','FPoT2','MH','MT2'},...
    'FontSize',8,'XTickLabelRotation',30,'TickDir','out');
ylabel('E/e'' Ratio','FontSize',9); title('L2: E/e'' (Diastolic)','FontSize',10,'FontWeight','bold'); box on;

% Flow arrow annotation (L2→L3)
subplot(3,4,8); axis off;
text(0.5,0.6,'Layer 2 \rightarrow Layer 3','HorizontalAlignment','center',...
    'FontSize',12,'FontWeight','bold','Color',[0.3 0.5 0.3]);
text(0.5,0.4,'Lipotoxicity & fibrosis','HorizontalAlignment','center',...
    'FontSize',9,'Color',[0.4 0.4 0.4]);
text(0.5,0.3,'impair microvascular','HorizontalAlignment','center',...
    'FontSize',9,'Color',[0.4 0.4 0.4]);
text(0.5,0.2,'endothelial function','HorizontalAlignment','center',...
    'FontSize',9,'Color',[0.4 0.4 0.4]);

% Row 3: Layer 3 variables (clinical endpoints)
subplot(3,4,9);
b = bar(l3_cfr,0.7); b.FaceColor='flat';
for s=1:n_states, b.CData(s,:)=state_colors(s,:); end
yline(2.5,'g--','LineWidth',1.2); yline(2.0,'r--','LineWidth',1.2);
set(gca,'XTickLabel',{'FPH','FPT2','FPoH','FPoT2','MH','MT2'},...
    'FontSize',8,'XTickLabelRotation',30,'TickDir','out');
ylabel('CFR','FontSize',9); title('L3: Coronary Flow Reserve','FontSize',10,'FontWeight','bold'); box on;

subplot(3,4,10);
b = bar(l3_imr,0.7); b.FaceColor='flat';
for s=1:n_states, b.CData(s,:)=state_colors(s,:); end
yline(25,'r--','LineWidth',1.2); yline(40,'Color',[0.6 0 0],'LineStyle','--','LineWidth',1.2);
set(gca,'XTickLabel',{'FPH','FPT2','FPoH','FPoT2','MH','MT2'},...
    'FontSize',8,'XTickLabelRotation',30,'TickDir','out');
ylabel('IMR (U)','FontSize',9); title('L3: IMR','FontSize',10,'FontWeight','bold'); box on;

subplot(3,4,11);
b = bar(l3_cmd,0.7); b.FaceColor='flat';
for s=1:n_states, b.CData(s,:)=state_colors(s,:); end
yline(1.0,'r--','LineWidth',1.2,'Label','CMD threshold');
set(gca,'XTickLabel',{'FPH','FPT2','FPoH','FPoT2','MH','MT2'},...
    'FontSize',8,'XTickLabelRotation',30,'TickDir','out');
ylabel('CMD Score (composite)','FontSize',9); title('L3: CMD Composite Score','FontSize',10,'FontWeight','bold'); box on;

% Legend panel
subplot(3,4,12); axis off;
for s = 1:n_states
    patch_h(s) = patch(NaN,NaN,state_colors(s,:),...
        'LineStyle',state_linestyle{s},'LineWidth',1.5);
end
legend(patch_h, state_labels,'Location','best','FontSize',8);
title('State Legend','FontSize',10,'FontWeight','bold');

sgtitle({'Integrated 3-Layer QSP Phenotype — Sex-Specific Diabetic Cardiomyopathy',...
    'Layer 1: Mitochondria  |  Layer 2: Metabolic Flexibility  |  Layer 3: Microvascular CMD'},...
    'FontSize',12,'FontWeight','bold');
annotation('textbox',[0 0.96 1 0.04],'String',...
    sprintf('QSP-DCM | %s | v%s', run_timestamp, script_version),...
    'EdgeColor','none','HorizontalAlignment','center',...
    'FontSize',7,'Color',[0.5 0.5 0.5]);

save_figure(fig_int, fullfile(outdir,'Layer3_IntegratedPhenotype'));
fprintf('  Saved: Layer3_IntegratedPhenotype .png\n');

%% --- 8. Steady State Bars -----------------------------------------------
fig_bar = figure('Name','Layer3 Steady State Bars',...
    'Units','inches','Position',[0.5 0.5 14 9],'Color','white');

bar_data  = [[SS3.NO];[SS3.ET1];[SS3.MicroRes];[SS3.CFR];[SS3.IMR];[SS3.WallStress]]';
bar_ylabs = {'NO (\muM)','ET-1 (fold)','Micro Resistance (fold)',...
             'CFR','IMR (U)','Wall Stress (fold)'};

for v = 1:6
    subplot(2,3,v);
    b = bar(bar_data(:,v),0.7); b.FaceColor='flat';
    for s=1:n_states, b.CData(s,:)=state_colors(s,:); end
    set(gca,'XTickLabel',{'FPH','FPT2','FPoH','FPoT2','MH','MT2'},...
        'FontSize',9,'XTickLabelRotation',30,'TickDir','out');
    ylabel(bar_ylabs{v},'FontSize',10);
    title(bar_ylabs{v},'FontSize',11,'FontWeight','bold');
    box on;
end

sgtitle({'Layer 3: Steady-State Microvascular Profile — All States',...
    'QSP Model of Diabetic Cardiomyopathy'},...
    'FontSize',13,'FontWeight','bold');
annotation('textbox',[0 0.96 1 0.04],'String',...
    sprintf('QSP-DCM | %s | v%s', run_timestamp, script_version),...
    'EdgeColor','none','HorizontalAlignment','center',...
    'FontSize',7,'Color',[0.5 0.5 0.5]);

save_figure(fig_bar, fullfile(outdir,'Layer3_SteadyState_Bars'));
fprintf('  Saved: Layer3_SteadyState_Bars .png\n');

%% --- 9. Print Summary ---------------------------------------------------
fprintf('\n');
fprintf('=======================================================\n');
fprintf(' LAYER 3 STEADY STATE SUMMARY — ALL STATES\n');
fprintf('=======================================================\n');
fprintf('%-22s %7s %7s %7s %7s %8s\n','State','NO','ET1','CFR','IMR','CMD');
fprintf('%s\n', repmat('-',1,65));
for s = 1:n_states
    % Flag CMD severity
    if SS3(s).CFR < 2.0 || SS3(s).IMR > 25
        flag = ' ** CMD **';
    elseif SS3(s).CFR < 2.5 || SS3(s).IMR > 20
        flag = ' * borderline';
    else
        flag = ' normal';
    end
    fprintf('%-22s %7.3f %7.3f %7.3f %7.1f %8.3f  %s\n',...
        state_labels{s}, SS3(s).NO, SS3(s).ET1, ...
        SS3(s).CFR, SS3(s).IMR, SS3(s).CMD_score, flag);
end
fprintf('=======================================================\n');
fprintf(' CFR: normal>=2.5 | CMD<2.0\n');
fprintf(' IMR: normal<25U  | severe>40U\n');
fprintf('=======================================================\n');
fprintf(' Files saved to: %s\n', outdir);
fprintf('=======================================================\n\n');

fprintf(flog,'\n=== LAYER 3 STEADY STATE SUMMARY ===\n');
fprintf(flog,'%-22s %7s %7s %7s %7s %8s\n','State','NO','ET1','CFR','IMR','CMD');
for s = 1:n_states
    fprintf(flog,'%-22s %7.3f %7.3f %7.3f %7.1f %8.3f\n',...
        state_labels{s}, SS3(s).NO, SS3(s).ET1,...
        SS3(s).CFR, SS3(s).IMR, SS3(s).CMD_score);
end
fprintf(flog,'\nRun completed: %s\n', datestr(now));
fclose(flog);

fprintf('Log saved: Layer3_RunLog.txt\n');
fprintf('Run complete. All Layer 3 files saved.\n');
fprintf('The 3-layer QSP model is complete.\n\n');

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

%% --- load_L3_parameters --------------------------------------------------
function p = load_L3_parameters(sex_state, L2)
% Layer 3 parameters — coronary microvascular dysfunction
%
% Driving inputs from Layer 2:
%   Ceramide → endothelial dysfunction (eNOS uncoupling)
%   ROS      → NO scavenging, ET1 upregulation
%   CollagenX → perivascular fibrosis, increased baseline resistance
%   EeRatio  → elevated filling pressures compress microvasculature
%   delta_psi → energetic failure impairs vasodilatory reserve
%
% Parameter sources:
%   NO/eNOS: Förstermann & Münzel 2006 (PMID 16501158)
%   ET1:     Yanagisawa et al 1988; Dhaun & Webb 2019 (PMID 30538051)
%   CFR:     Taqueti & Di Carli 2018 (PMID 30213362)
%   IMR:     Fearon et al 2003 (PMID 14557361); Ng et al 2012
%   Sex/eNOS: Murphy & Steenbergen 2007 (PMID 17347480)

    % Pass through Layer 2 inputs
    p.Ceramide  = L2.Ceramide_SS;
    p.ROS       = L2.ROS_SS;
    p.ColX      = L2.CollagenX_SS;
    p.Ee        = L2.EeRatio_SS;
    p.dPsi      = L2.delta_psi_SS;
    p.PGC1a     = L2.PGC1a_SS;
    p.Mito      = L2.Mito_dens_SS;

    % --- Shared kinetic parameters ---
    % NO / eNOS axis
    p.k_NO_basal    = 0.50;    % basal eNOS NO production (µM/h)
    p.k_NO_cer_inh  = 0.08;   % ceramide inhibits eNOS (PMID 16501158)
    p.k_NO_ROS_scav = 0.20;   % ROS scavenges NO (peroxynitrite formation)
    p.k_NO_deg      = 0.15;   % NO degradation

    % ET-1 axis
    p.k_ET1_basal   = 0.05;   % basal ET-1 production
    p.k_ET1_ROS     = 0.12;   % ROS drives ET-1 transcription
    p.k_ET1_cer     = 0.08;   % ceramide upregulates ET-1
    p.k_ET1_NO_inh  = 0.10;   % NO suppresses ET-1 (negative feedback)
    p.k_ET1_deg     = 0.08;

    % Microvascular resistance
    p.k_MRes_ET1    = 0.15;   % ET-1 drives vasoconstriction
    p.k_MRes_colx   = 0.10;   % perivascular fibrosis → structural resistance
    p.k_MRes_NO_dil = 0.20;   % NO-driven vasodilation
    p.k_MRes_deg    = 0.05;

    % CFR — coronary flow reserve
    % CFR = hyperemic / resting flow; reduced by structural + functional resistance
    p.k_CFR_MRes    = 0.20;   % resistance reduces hyperemic capacity
    p.k_CFR_dPsi    = 0.10;   % energetic failure limits vasodilatory reserve
    p.k_CFR_deg     = 0.30;   % strong pull toward physiological ceiling (3.5)

    % IMR — index of microcirculatory resistance
    % IMR = distal pressure × hyperemic transit time; rises with resistance
    p.k_IMR_MRes    = 0.12;    % resistance drives IMR (scaled to U range 15-50)
    p.k_IMR_Ee      = 0.04;    % elevated filling pressure raises IMR
    p.k_IMR_deg     = 0.03;

    % Wall stress (LaPlace — integrates diastolic + microvascular load)
    p.k_WS_Ee       = 0.04;   % diastolic stiffness → wall stress
    p.k_WS_IMR      = 0.002;  % microvascular ischemia → wall stress
    p.k_WS_colx     = 0.06;   % fibrosis → wall stress
    p.k_WS_deg      = 0.03;

    % --- State-specific parameters ---
    switch sex_state

        case 'female_pre_healthy'
            p.eNOS_activity  = 1.00;   % full eNOS (estrogen-driven)
            p.vasc_estrogen  = 1.00;   % full vascular estrogen protection
            p.ET1_basal_mult = 1.00;
            p.NO_0    = 0.80;   p.ET1_0  = 1.00;
            p.MRes_0  = 1.00;   p.CFR_0  = 3.20;
            p.IMR_0   = 15.0;   p.WS_0   = 1.00;

        case 'female_pre_T2DM'
            p.eNOS_activity  = 0.65;   % eNOS uncoupling by hyperglycemia
            p.vasc_estrogen  = 0.65;
            p.ET1_basal_mult = 1.50;
            p.NO_0    = 0.55;   p.ET1_0  = 1.40;
            p.MRes_0  = 1.30;   p.CFR_0  = 2.40;
            p.IMR_0   = 28.0;   p.WS_0   = 1.30;

        case 'female_post_healthy'
            p.eNOS_activity  = 0.75;   % reduced without estrogen
            p.vasc_estrogen  = 0.70;
            p.ET1_basal_mult = 1.20;
            p.NO_0    = 0.65;   p.ET1_0  = 1.15;
            p.MRes_0  = 1.15;   p.CFR_0  = 2.80;
            p.IMR_0   = 20.0;   p.WS_0   = 1.10;

        case 'female_post_T2DM'
            p.eNOS_activity  = 0.35;   % severe eNOS uncoupling
            p.vasc_estrogen  = 0.30;   % minimal vascular protection
            p.ET1_basal_mult = 2.00;   % ET-1 dominance
            p.NO_0    = 0.30;   p.ET1_0  = 1.80;
            p.MRes_0  = 1.70;   p.CFR_0  = 1.80;
            p.IMR_0   = 42.0;   p.WS_0   = 1.70;

        case 'male_healthy'
            p.eNOS_activity  = 0.80;   % lower than female pre (less ERb-driven)
            p.vasc_estrogen  = 0.60;
            p.ET1_basal_mult = 1.05;
            p.NO_0    = 0.72;   p.ET1_0  = 1.05;
            p.MRes_0  = 1.05;   p.CFR_0  = 3.00;
            p.IMR_0   = 16.0;   p.WS_0   = 1.05;

        case 'male_T2DM'
            p.eNOS_activity  = 0.55;
            p.vasc_estrogen  = 0.45;
            p.ET1_basal_mult = 1.60;
            p.NO_0    = 0.48;   p.ET1_0  = 1.55;
            p.MRes_0  = 1.45;   p.CFR_0  = 2.10;
            p.IMR_0   = 32.0;   p.WS_0   = 1.45;

        otherwise
            error('Unknown sex_state: %s', sex_state);
    end
end

%% --- layer3_odes ---------------------------------------------------------
function dydt = layer3_odes(~, y, p)
% 6 coupled ODEs — Layer 3 Coronary Microvascular Dysfunction
%
%  y(1) = NO         (µM)   — nitric oxide
%  y(2) = ET1        (fold) — endothelin-1
%  y(3) = MicroRes   (fold) — microvascular resistance
%  y(4) = CFR        (ratio)— coronary flow reserve
%  y(5) = IMR        (U)    — index of microcirculatory resistance
%  y(6) = WallStress (fold) — myocardial wall stress

    NO      = max(y(1), 0.001);
    ET1     = max(y(2), 0);
    MRes    = max(y(3), 0);
    CFR     = max(y(4), 0);
    IMR     = max(y(5), 0);
    WS      = max(y(6), 0);

    % Normalized ΔΨm deficiency (0=normal, 1=complete failure)
    dPsi_def = max(1 - p.dPsi/50, 0);

    % --- ODE 1: NO ---
    % Produced by eNOS (estrogen/PGC-1a driven); scavenged by ROS + ceramide
    dNO = p.k_NO_basal * p.eNOS_activity * p.vasc_estrogen ...
        - p.k_NO_cer_inh  * p.Ceramide * NO ...
        - p.k_NO_ROS_scav * p.ROS * NO ...
        - p.k_NO_deg * NO;

    % --- ODE 2: ET-1 ---
    % Upregulated by ROS and ceramide; suppressed by NO
    dET1 = p.k_ET1_basal * p.ET1_basal_mult ...
         + p.k_ET1_ROS * p.ROS ...
         + p.k_ET1_cer * p.Ceramide / 10 ...
         - p.k_ET1_NO_inh * NO * ET1 ...
         - p.k_ET1_deg * ET1;

    % --- ODE 3: Microvascular resistance ---
    % ET-1 drives constriction; perivascular fibrosis adds structural resistance
    % NO provides vasodilation
    dMRes = p.k_MRes_ET1  * ET1 ...
          + p.k_MRes_colx * p.ColX ...
          - p.k_MRes_NO_dil * NO * MRes ...
          - p.k_MRes_deg * MRes;

    % --- ODE 4: CFR ---
    % Reduced by resistance and energetic failure (limits hyperemic response)
    % Reference CFR ~3.5 in healthy; floor at 0.5 (critical ischemia)
    % NOTE: max(MRes-1,0) ensures that vasodilation (MRes<1) does not
    %       spuriously push CFR above the physiological ceiling of 3.5.
    %       Only elevated resistance (MRes>1) penalises CFR.
    dCFR = - p.k_CFR_MRes * max(MRes - 1.0, 0) * CFR ...
           - p.k_CFR_dPsi * dPsi_def * CFR ...
           + p.k_CFR_deg  * (3.5 - CFR);  % pull toward physiological max

    % --- ODE 5: IMR ---
    % Rises with microvascular resistance and elevated filling pressures (E/e')
    dIMR = p.k_IMR_MRes * MRes ...
         + p.k_IMR_Ee   * p.Ee ...
         - p.k_IMR_deg  * IMR;

    % --- ODE 6: Wall stress ---
    % Integrates diastolic stiffness (L2), microvascular ischemia (L3),
    % and fibrotic collagen load (L2)
    dWS = p.k_WS_Ee   * p.Ee / 10 ...
        + p.k_WS_IMR  * IMR ...
        + p.k_WS_colx * p.ColX ...
        - p.k_WS_deg  * WS;

    dydt = [dNO; dET1; dMRes; dCFR; dIMR; dWS];
end

%% --- plot_L3_timecourse --------------------------------------------------
function plot_L3_timecourse(t, y, sex_state, label, p, L2, outdir)

    fig = figure('Name', label, 'Units','inches',...
        'Position',[0.5 0.5 14 10], 'Color','white', 'Visible','off');

    var_names = {'NO (\muM)','ET-1 (fold)','Micro Resistance (fold)',...
                 'CFR','IMR (U)','Wall Stress (fold)'};

    ref_lines = {
        4, 2.5, [0.2 0.7 0.2], '--', 'Normal CFR (\geq2.5)';
        4, 2.0, [0.8 0.2 0.2], '--', 'CMD threshold (<2.0)';
        5, 25,  [0.8 0.2 0.2], '--', 'Dysfunction (>25U)';
        5, 40,  [0.6 0.0 0.0], '--', 'Severe CMD (>40U)';
    };

    for v = 1:6
        ax = subplot(3,3,v);
        plot(t, y(:,v), 'Color',[0.15 0.40 0.75],'LineWidth',2.0);
        hold on;
        for r = 1:size(ref_lines,1)
            if ref_lines{r,1} == v
                yline(ref_lines{r,2}, ref_lines{r,4}, ref_lines{r,5},...
                    'Color',ref_lines{r,3},'LineWidth',1.2,'FontSize',7,...
                    'LabelHorizontalAlignment','right');
            end
        end
        xlabel('Time (hours)','FontSize',9);
        ylabel(var_names{v},'FontSize',9);
        title(var_names{v},'FontSize',10,'FontWeight','bold');
        xlim([0 168]);
        set(ax,'FontSize',8,'TickDir','out','LineWidth',0.8,'Box','on');
    end

    % Subplot 7: Layer 2 inputs bar
    subplot(3,3,7);
    l2_vals  = [L2.Ceramide_SS/6, L2.CollagenX_SS/1.5, ...
                L2.EeRatio_SS/15, L2.ROS_SS, L2.delta_psi_SS/50];
    l2_names = {'Cer/6',  'ColX/1.5', 'E/e''/15', 'ROS', 'd\Psi/50'};
    bar(l2_vals, 0.6,'FaceColor',[0.4 0.65 0.4],'EdgeColor','k');
    set(gca,'XTickLabel',l2_names,'FontSize',8,'XTickLabelRotation',20,'TickDir','out');
    yline(1,'k--','Reference','FontSize',7);
    ylabel('Normalized L2 Input','FontSize',8);
    title('Layer 2 Inputs (normalized)','FontSize',10,'FontWeight','bold'); box on;

    % Subplot 8: CMD severity gauge
    subplot(3,3,8);
    cmd_score = (y(end,5)/25) * (2.5/max(y(end,4),0.1));
    categories = {'Normal','Borderline','CMD','Severe CMD'};
    thresholds = [0.5, 1.0, 1.5, 3.0];
    colors_cat = [0.2 0.7 0.2; 1.0 0.8 0.0; 0.9 0.4 0.0; 0.8 0.1 0.1];
    bar_h = bar(cmd_score, 0.5);
    idx = find(cmd_score < thresholds, 1);
    if isempty(idx), idx = 4; end
    bar_h.FaceColor = colors_cat(idx,:);
    set(gca,'XTick',[],'YLim',[0 3],'TickDir','out');
    yline(0.5,'g--','Normal','FontSize',7);
    yline(1.0,'y--','Borderline','FontSize',7);
    yline(1.5,'r--','CMD','FontSize',7);
    ylabel('CMD Composite Score','FontSize',8);
    title(sprintf('CMD Score = %.2f', cmd_score),'FontSize',10,'FontWeight','bold');
    box on;

    % Subplot 9: State parameters
    subplot(3,3,9); axis off;
    txt = sprintf(['State: %s\n\neNOS activity = %.2f\n'...
        'Vasc estrogen = %.2f\nET-1 basal mult = %.2f\n\n'...
        'SS CFR  = %.2f\nSS IMR  = %.1f U\nSS NO   = %.3f uM'], ...
        label, p.eNOS_activity, p.vasc_estrogen, p.ET1_basal_mult,...
        y(end,4), y(end,5), y(end,1));
    text(0.05, 0.95, txt,'Units','normalized','VerticalAlignment','top',...
        'FontSize',8,'FontName','Courier');
    title('State Parameters','FontSize',10,'FontWeight','bold');

    sgtitle({['Layer 3: ' label],...
        'Coronary Microvascular Dysfunction'},...
        'FontSize',12,'FontWeight','bold');
    annotation('textbox',[0 0.01 1 0.03],'String',...
        'QSP Diabetic Cardiomyopathy | Layer 3',...
        'EdgeColor','none','HorizontalAlignment','center',...
        'FontSize',7,'Color',[0.5 0.5 0.5]);

    save_figure(fig, fullfile(outdir, ['Layer3_' sex_state '_TimeCourse']));
    close(fig);
end

%% --- save_figure ---------------------------------------------------------
function save_figure(fig, basepath)
    print(fig, basepath, '-dpng', '-r300');
end
