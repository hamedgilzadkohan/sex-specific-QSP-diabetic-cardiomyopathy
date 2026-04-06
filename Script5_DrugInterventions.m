%% ========================================================================
%  Script5_DrugInterventions.m
%  Sex-Specific QSP Model of Diabetic Cardiomyopathy
%  PHARMACOLOGICAL INTERVENTION SIMULATIONS
%
%  Authors : Soheili M, Gilzad-Kohan H, Lotfi AS
%
%  Drugs simulated (all at clinically relevant doses):
%    1. Empagliflozin  (SGLT2i)      — ROS↓, ΔΨm↑, AMPK↑, ceramide↓
%    2. Semaglutide    (GLP-1 RA)    — Insulin resistance↓, FA overflow↓
%    3. MHT            (Estrogen)    — eNOS↑, ERb↑, E2 restore, ET1↓
%    4. Sacubitril/Valsartan (ARNI)  — Wall stress↓, E/e'↓, collagen↓
%
%  For each drug: simulate on all 6 states, compare to vehicle (no drug)
%  Primary outputs: CFR, IMR, E/e', PGC-1α, Ceramide, CMD Score
%
%  Output files:
%    - Drug_ResponseProfiles.png       (primary intervention figure)
%    - Drug_SexComparison.png          (sex-stratified response)
%    - Drug_WaterfallPlot.png          (% change from baseline)
%    - Drug_AllStates_Results.xlsx          (full results table)
%    - Drug_RunLog.txt
%
%  Version : 1.0  |  Date: 2026-03-09
%  MATLAB  : R2019b or newer
%% ========================================================================

clear; close all; clc;

script_version = '1.0';
run_timestamp  = datestr(now, 'yyyy-mm-dd HH:MM:SS');
script_path = fileparts(mfilename('fullpath'));
if isempty(script_path), script_path = pwd; end
outdir = script_path;

flog = fopen(fullfile(outdir,'Drug_RunLog.txt'),'w');
fprintf(flog,'=======================================================\n');
fprintf(flog,' Drug Intervention Simulations — Reproducibility Log\n');
fprintf(flog,'=======================================================\n');
fprintf(flog,' Version   : %s\n', script_version);
fprintf(flog,' Timestamp : %s\n', run_timestamp);
fprintf(flog,' MATLAB    : %s\n', version);
fprintf(flog,'=======================================================\n\n');

fprintf('=======================================================\n');
fprintf(' Script 5: Pharmacological Intervention Simulations\n');
fprintf(' %s\n', run_timestamp);
fprintf('=======================================================\n\n');

%% --- Model settings ------------------------------------------------------
tspan   = [0 336];    % 14 days (longer horizon for drug effect)
options = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.1);

state_ids     = {'female_pre_healthy','female_pre_T2DM',...
                 'female_post_healthy','female_post_T2DM',...
                 'male_healthy','male_T2DM'};
state_labels  = {'Female Pre Healthy','Female Pre T2DM',...
                 'Female Post Healthy','Female Post T2DM',...
                 'Male Healthy','Male T2DM'};
n_states = length(state_ids);

%% --- Drug definitions ----------------------------------------------------
% Each drug modifies specific parameters in one or more layers
% Effect magnitudes calibrated to published clinical/preclinical data

drugs(1).name    = 'Empagliflozin';
drugs(1).abbrev  = 'SGLT2i';
drugs(1).color   = [0.20 0.60 0.85];
drugs(1).mech    = 'ROS reduction, mitochondrial protection, AMPK activation';
drugs(1).L1.k_ROS_prod   = 0.75;   % 25% ROS production reduction
drugs(1).L1.k_Psi_form   = 1.20;   % 20% improved ΔΨm formation
drugs(1).L1.T2DM_factor  = 0.85;   % 15% metabolic stress reduction
drugs(1).L2.T2DM_ins_res = 0.80;   % 20% insulin resistance improvement
drugs(1).L2.FA_overflow  = 0.85;   % 15% FA overflow reduction
drugs(1).L2.k_Cer_ROS    = 0.80;   % 20% ceramide synthesis reduction
drugs(1).L3.k_NO_ROS_scav= 0.80;   % 20% less NO scavenging by ROS
drugs(1).L3_eNOS_boost   = 1.10;   % 10% eNOS activity improvement

drugs(2).name    = 'Semaglutide';
drugs(2).abbrev  = 'GLP-1';
drugs(2).color   = [0.20 0.75 0.45];
drugs(2).mech    = 'Insulin sensitization, FA overflow reduction, ceramide lowering';
drugs(2).L1.T2DM_factor  = 0.80;   % 20% T2DM stress reduction
drugs(2).L1.ROS_drive    = 0.85;   % 15% ROS drive reduction
drugs(2).L2.T2DM_ins_res = 0.65;   % 35% insulin resistance improvement
drugs(2).L2.FA_overflow  = 0.70;   % 30% FA overflow reduction
drugs(2).L2.k_Cer_FA     = 0.75;   % 25% ceramide from FA reduction
drugs(2).L3.k_ET1_cer    = 0.85;   % 15% ET1 from ceramide reduction
drugs(2).L3_eNOS_boost   = 1.05;   % 5% eNOS improvement

drugs(3).name    = 'MHT (Estrogen)';
drugs(3).abbrev  = 'MHT';
drugs(3).color   = [0.90 0.35 0.65];
drugs(3).mech    = 'ERb restoration, eNOS activation, ET1 suppression, mitochondrial protection';
drugs(3).L1.k_ERb_sens   = 1.25;   % ERb sensitivity restoration
drugs(3).L1.E2_ss_mult   = 8.0;    % E2 level restoration (post→pre range)
drugs(3).L1.mito_stress  = 0.85;   % 15% mito stress reduction
drugs(3).L2.estrogen_prot= 1.40;   % 40% estrogen protection boost
drugs(3).L2.k_Cer_FA     = 0.80;   % 20% ceramide reduction
drugs(3).L3.eNOS_activity= 1.50;   % 50% eNOS activity restoration
drugs(3).L3.vasc_estrogen= 1.60;   % 60% vascular estrogen protection
drugs(3).L3.ET1_basal_mult=0.70;   % 30% ET1 reduction
drugs(3).L3_eNOS_boost   = 1.50;

drugs(4).name    = 'Sacubitril/Valsartan';
drugs(4).abbrev  = 'ARNI';
drugs(4).color   = [0.85 0.55 0.15];
drugs(4).mech    = 'Natriuretic peptide enhancement, wall stress reduction, fibrosis attenuation';
drugs(4).L1.T2DM_factor  = 0.92;   % modest mitochondrial benefit
drugs(4).L2.k_Col_ROS    = 0.75;   % 25% collagen crosslinking reduction
drugs(4).L2.k_Col_cer    = 0.75;   % 25% ceramide-driven fibrosis reduction
drugs(4).L2.k_Ee_Col     = 0.80;   % 20% E/e improvement via collagen
drugs(4).L3.k_MRes_colx  = 0.70;   % 30% structural resistance reduction
drugs(4).L3.k_WS_Ee      = 0.75;   % 25% wall stress reduction
drugs(4).L3_eNOS_boost   = 1.08;   % 8% modest eNOS improvement

n_drugs = length(drugs);

%% --- Run vehicle (no drug) baseline for all states ----------------------
fprintf('Running vehicle baselines...\n');
vehicle = struct();
for s = 1:n_states
    res = run_model(state_ids{s}, [], tspan, options);
    vehicle(s).CFR      = res.CFR;
    vehicle(s).IMR      = res.IMR;
    vehicle(s).Ee       = res.Ee;
    vehicle(s).PGC1a    = res.PGC1a;
    vehicle(s).Ceramide = res.Ceramide;
    vehicle(s).CMD      = res.CMD;
    vehicle(s).WS       = res.WS;
    fprintf('  Vehicle %s: CFR=%.3f  IMR=%.1f  E/e=%.2f  CMD=%.3f\n',...
        state_labels{s},res.CFR,res.IMR,res.Ee,res.CMD);
end

%% --- Run each drug on each state ----------------------------------------
fprintf('\nRunning drug simulations...\n');
drug_res = struct();
for d = 1:n_drugs
    fprintf('\n  Drug %d/%d: %s\n', d, n_drugs, drugs(d).name);
    for s = 1:n_states
        res = run_model(state_ids{s}, drugs(d), tspan, options);
        drug_res(d,s).CFR      = res.CFR;
        drug_res(d,s).IMR      = res.IMR;
        drug_res(d,s).Ee       = res.Ee;
        drug_res(d,s).PGC1a    = res.PGC1a;
        drug_res(d,s).Ceramide = res.Ceramide;
        drug_res(d,s).CMD      = res.CMD;
        drug_res(d,s).WS       = res.WS;

        % Percent change from vehicle
        drug_res(d,s).dCFR  = 100*(res.CFR      - vehicle(s).CFR)      / max(vehicle(s).CFR,0.01);
        drug_res(d,s).dIMR  = 100*(res.IMR      - vehicle(s).IMR)      / max(vehicle(s).IMR,0.01);
        drug_res(d,s).dEe   = 100*(res.Ee       - vehicle(s).Ee)       / max(vehicle(s).Ee,0.01);
        drug_res(d,s).dPGC  = 100*(res.PGC1a    - vehicle(s).PGC1a)    / max(vehicle(s).PGC1a,0.01);
        drug_res(d,s).dCer  = 100*(res.Ceramide - vehicle(s).Ceramide) / max(vehicle(s).Ceramide,0.01);
        drug_res(d,s).dCMD  = 100*(res.CMD      - vehicle(s).CMD)      / max(vehicle(s).CMD,0.01);

        fprintf('    %s: CFR=%+.1f%%  IMR=%+.1f%%  E/e=%+.1f%%  CMD=%+.1f%%\n',...
            state_labels{s},...
            drug_res(d,s).dCFR, drug_res(d,s).dIMR,...
            drug_res(d,s).dEe,  drug_res(d,s).dCMD);
    end
end

%% --- Figure 1: Drug Response Profiles -----------------------------------
fprintf('\nGenerating drug response profiles figure...\n');

state_colors = [0.20 0.45 0.80;   % FPH  - dark blue
                0.20 0.45 0.80;   % FPT2 - dark blue dashed
                0.55 0.75 0.95;   % FPoH - light blue
                0.55 0.75 0.95;   % FPoT2- light blue dashed
                0.85 0.20 0.20;   % MH   - dark red
                0.85 0.20 0.20];  % MT2  - dark red dashed

line_styles = {'-','--','-','--','-','--'};

fig1 = figure('Name','Drug Response Profiles',...
    'Units','inches','Position',[0.5 0.5 18 14],'Color','white');

metrics = {'CFR','IMR','Ee','PGC1a','Ceramide','CMD'};
metric_labels = {'CFR','IMR (U)','E/e''','PGC-1\alpha (fold)','Ceramide (uM)','CMD Score'};
thresholds_y  = {2.0, 25, 15, NaN, 6.0, 1.5};
thresh_colors = {'r','r','r','none','r','r'};
thresh_dir    = {'below','above','above','none','above','above'};

drug_abbrevs = {'SGLT2i','GLP-1','MHT','ARNI'};
drug_colors_list = {[0.20 0.60 0.85],[0.20 0.75 0.45],[0.90 0.35 0.65],[0.85 0.55 0.15]};

% Focus: female_post_T2DM (state 4) and male_T2DM (state 6)
focus_s = [4, 6];
focus_names = {'Female Post T2DM','Male T2DM'};

for mi = 1:6
    subplot(2,3,mi);
    hold on;

    metric = metrics{mi};
    x_pos = 1:n_drugs;
    bar_width = 0.35;

    % Get values for both focus states
    vals_fpo_t2 = zeros(1,n_drugs);
    vals_m_t2   = zeros(1,n_drugs);
    veh_fpo     = vehicle(4).(metric);
    veh_m       = vehicle(6).(metric);

    for d = 1:n_drugs
        vals_fpo_t2(d) = drug_res(d,4).(metric);
        vals_m_t2(d)   = drug_res(d,6).(metric);
    end

    % Add vehicle as position 0
    x_all = [0, x_pos];
    v_fpo = [veh_fpo, vals_fpo_t2];
    v_m   = [veh_m, vals_m_t2];

    % Grouped bar
    b1 = bar(x_all - bar_width/2, v_fpo, bar_width, ...
        'FaceColor',[0.55 0.75 0.95],'EdgeColor',[0.20 0.45 0.80],'LineWidth',1.0);
    b2 = bar(x_all + bar_width/2, v_m, bar_width, ...
        'FaceColor',[0.85 0.55 0.55],'EdgeColor',[0.70 0.15 0.15],'LineWidth',1.0);

    % Threshold line
    if ~isnan(thresholds_y{mi})
        yline(thresholds_y{mi},'r--','LineWidth',1.2);
    end

    set(gca,'XTick',x_all,'XTickLabel',['Vehicle',drug_abbrevs],...
        'FontSize',8,'TickDir','out');
    ylabel(metric_labels{mi},'FontSize',9);
    title(metric_labels{mi},'FontSize',10,'FontWeight','bold');
    if mi == 1
        legend([b1,b2],{'Female Post T2DM','Male T2DM'},...
            'Location','northwest','FontSize',8);
    end
    box on; grid on; grid minor;
    hold off;
end

sgtitle({'Drug Intervention Simulations — Primary Endpoints',...
    'Female Post T2DM vs Male T2DM | QSP Model of Diabetic Cardiomyopathy'},...
    'FontSize',12,'FontWeight','bold');
annotation('textbox',[0 0.96 1 0.04],'String',...
    sprintf('QSP-DCM | %s | v%s', run_timestamp, script_version),...
    'EdgeColor','none','HorizontalAlignment','center',...
    'FontSize',7,'Color',[0.5 0.5 0.5]);

save_figure(fig1, fullfile(outdir,'Drug_ResponseProfiles'));
fprintf('  Saved: Drug_ResponseProfiles\n');

%% --- Figure 2: Waterfall Plot (% change from vehicle) -------------------
fprintf('Generating waterfall plots...\n');

fig2 = figure('Name','Drug Waterfall - Percent Change',...
    'Units','inches','Position',[0.5 0.5 18 12],'Color','white');

delta_metrics = {'dCFR','dIMR','dEe','dPGC','dCer','dCMD'};
delta_labels  = {'\DeltaCFR (%)','  \DeltaIMR (%)','\DeltaE/e'' (%)',...
                 '\DeltaPGC-1\alpha (%)','  \DeltaCeramide (%)','  \DeltaCMD Score (%)'};
% For CFR: positive = good. For IMR,Ee,Cer,CMD: negative = good.
good_direction = [1,-1,-1,1,-1,-1];

for mi = 1:6
    subplot(2,3,mi);
    hold on;

    metric = delta_metrics{mi};
    good_dir = good_direction(mi);

    % Collect: 4 drugs x 6 states
    delta_mat = zeros(n_drugs, n_states);
    for d = 1:n_drugs
        for s = 1:n_states
            delta_mat(d,s) = drug_res(d,s).(metric);
        end
    end

    % Plot as grouped bar: x = states, groups = drugs
    x_pos = 1:n_states;
    offset = linspace(-0.35, 0.35, n_drugs);

    for d = 1:n_drugs
        b = bar(x_pos + offset(d), delta_mat(d,:), 0.18,...
            'FaceColor', drugs(d).color, 'EdgeColor','none');
    end

    yline(0,'k-','LineWidth',1.0);
    set(gca,'XTick',1:n_states,'XTickLabel',...
        {'FPH','FPT2','FPoH','FPoT2','MH','MT2'},...
        'FontSize',8,'TickDir','out');
    ylabel(delta_labels{mi},'FontSize',9);
    title(delta_labels{mi},'FontSize',10,'FontWeight','bold');

    if good_dir == 1
        title_str = ['\uparrow better'];
    else
        title_str = ['\downarrow better'];
    end
    text(0.98,0.95,title_str,'Units','normalized',...
        'HorizontalAlignment','right','FontSize',8,'Color',[0.4 0.4 0.4]);

    if mi == 1
        leg_h = zeros(1,n_drugs);
        for d = 1:n_drugs
            leg_h(d) = patch(NaN,NaN,drugs(d).color);
        end
        legend(leg_h,{drugs.abbrev},'Location','southwest','FontSize',8);
    end
    box on; grid on;
    hold off;
end

sgtitle({'Drug Response — Percent Change from Vehicle (All 6 States)',...
    'QSP Model of Sex-Specific Diabetic Cardiomyopathy'},...
    'FontSize',12,'FontWeight','bold');
annotation('textbox',[0 0.96 1 0.04],'String',...
    sprintf('QSP-DCM | %s | v%s', run_timestamp, script_version),...
    'EdgeColor','none','HorizontalAlignment','center',...
    'FontSize',7,'Color',[0.5 0.5 0.5]);

save_figure(fig2, fullfile(outdir,'Drug_WaterfallPlot'));
fprintf('  Saved: Drug_WaterfallPlot\n');

%% --- Figure 3: Sex-Stratified Drug Comparison ---------------------------
fprintf('Generating sex comparison figure...\n');

fig3 = figure('Name','Sex-Stratified Drug Comparison',...
    'Units','inches','Position',[0.5 0.5 16 10],'Color','white');

% Focus: Female Post T2DM vs Male T2DM — CFR and IMR response
subplot(1,2,1);  % CFR
hold on;
x = 1:n_drugs;
cfr_fpo = arrayfun(@(d) drug_res(d,4).CFR, 1:n_drugs);
cfr_mt2 = arrayfun(@(d) drug_res(d,6).CFR, 1:n_drugs);
cfr_veh_fpo = vehicle(4).CFR;
cfr_veh_mt2 = vehicle(6).CFR;

plot([0 x], [cfr_veh_fpo cfr_fpo], 'o-','Color',[0.20 0.45 0.80],...
    'LineWidth',2,'MarkerSize',8,'MarkerFaceColor',[0.55 0.75 0.95]);
plot([0 x], [cfr_veh_mt2 cfr_mt2], 's-','Color',[0.85 0.20 0.20],...
    'LineWidth',2,'MarkerSize',8,'MarkerFaceColor',[0.85 0.55 0.55]);
yline(2.5,'g--','LineWidth',1.5,'Label','Normal CFR','LabelHorizontalAlignment','right');
yline(2.0,'r--','LineWidth',1.5,'Label','CMD threshold','LabelHorizontalAlignment','right');
set(gca,'XTick',0:n_drugs,'XTickLabel',['Vehicle',{drugs.abbrev}],...
    'FontSize',9,'TickDir','out');
ylabel('Coronary Flow Reserve (CFR)','FontSize',10);
title('CFR Response to Treatment','FontSize',11,'FontWeight','bold');
legend({'Female Post T2DM','Male T2DM'},'Location','northwest','FontSize',9);
ylim([0 max([cfr_fpo cfr_mt2])*1.3]);
box on; grid on; hold off;

subplot(1,2,2);  % IMR
hold on;
imr_fpo = arrayfun(@(d) drug_res(d,4).IMR, 1:n_drugs);
imr_mt2 = arrayfun(@(d) drug_res(d,6).IMR, 1:n_drugs);
imr_veh_fpo = vehicle(4).IMR;
imr_veh_mt2 = vehicle(6).IMR;

plot([0 x], [imr_veh_fpo imr_fpo], 'o-','Color',[0.20 0.45 0.80],...
    'LineWidth',2,'MarkerSize',8,'MarkerFaceColor',[0.55 0.75 0.95]);
plot([0 x], [imr_veh_mt2 imr_mt2], 's-','Color',[0.85 0.20 0.20],...
    'LineWidth',2,'MarkerSize',8,'MarkerFaceColor',[0.85 0.55 0.55]);
yline(40,'r--','LineWidth',1.5,'Label','Severe CMD (>40U)','LabelHorizontalAlignment','right');
yline(25,'Color',[0.85 0.50 0.0],'LineStyle','--','LineWidth',1.5,...
    'Label','Dysfunction (>25U)','LabelHorizontalAlignment','right');
set(gca,'XTick',0:n_drugs,'XTickLabel',['Vehicle',{drugs.abbrev}],...
    'FontSize',9,'TickDir','out');
ylabel('Index of Microcirculatory Resistance (IMR, U)','FontSize',10);
title('IMR Response to Treatment','FontSize',11,'FontWeight','bold');
legend({'Female Post T2DM','Male T2DM'},'Location','northeast','FontSize',9);
box on; grid on; hold off;

sgtitle({'Sex-Stratified Treatment Response — CFR & IMR',...
    'The Same Drug Achieves Different Outcomes in Males vs Postmenopausal Women'},...
    'FontSize',12,'FontWeight','bold');
annotation('textbox',[0 0.96 1 0.04],'String',...
    sprintf('QSP-DCM | %s | v%s', run_timestamp, script_version),...
    'EdgeColor','none','HorizontalAlignment','center',...
    'FontSize',7,'Color',[0.5 0.5 0.5]);

save_figure(fig3, fullfile(outdir,'Drug_SexComparison'));
fprintf('  Saved: Drug_SexComparison\n');

%% --- Figure 4: MHT Spotlight (Estrogen — most sex-specific) -------------
fprintf('Generating MHT spotlight figure...\n');

fig4 = figure('Name','MHT Estrogen Therapy Spotlight',...
    'Units','inches','Position',[0.5 0.5 16 10],'Color','white');

% Run MHT time courses for Female Post T2DM vs vehicle
[tc_veh, t_veh] = run_timecourse('female_post_T2DM', [], tspan, options);
[tc_mht, t_mht] = run_timecourse('female_post_T2DM', drugs(3), tspan, options);
[tc_veh_m, t_vm]= run_timecourse('male_T2DM', [], tspan, options);
[tc_mht_m, t_mm]= run_timecourse('male_T2DM', drugs(3), tspan, options);

plot_vars = {'CFR','IMR','Ee','NO'};
plot_ylabels = {'CFR','IMR (U)','E/e''','NO (uM)'};
thresh_vals = {[2.0 2.5],[25 40],15,NaN};

for pi = 1:4
    subplot(2,2,pi); hold on;
    v_name = plot_vars{pi};

    % Female Post T2DM
    plot(t_veh, tc_veh.(v_name), '--','Color',[0.55 0.75 0.95],'LineWidth',2);
    plot(t_mht, tc_mht.(v_name), '-','Color',[0.20 0.45 0.80],'LineWidth',2.5);
    % Male T2DM
    plot(t_vm,  tc_veh_m.(v_name), '--','Color',[0.85 0.55 0.55],'LineWidth',2);
    plot(t_mm,  tc_mht_m.(v_name), '-','Color',[0.85 0.20 0.20],'LineWidth',2.5);

    if ~isnan(thresh_vals{pi})
        if isnumeric(thresh_vals{pi}) && length(thresh_vals{pi})==2
            yline(thresh_vals{pi}(1),'r--','LineWidth',1.0);
            yline(thresh_vals{pi}(2),'Color',[0.85 0.50 0],'LineStyle','--','LineWidth',1.0);
        else
            yline(thresh_vals{pi},'r--','LineWidth',1.0);
        end
    end

    xlabel('Time (hours)','FontSize',9);
    ylabel(plot_ylabels{pi},'FontSize',10);
    title(plot_ylabels{pi},'FontSize',10,'FontWeight','bold');

    if pi == 1
        legend({'FPoT2 Vehicle','FPoT2 + MHT','MT2 Vehicle','MT2 + MHT'},...
            'Location','best','FontSize',8);
    end
    box on; grid on; hold off;
end

sgtitle({'MHT (Estrogen Therapy) — Female-Specific Treatment Response',...
    'Female Post T2DM vs Male T2DM | Layer 3 Microvascular Endpoints'},...
    'FontSize',12,'FontWeight','bold');
annotation('textbox',[0 0.96 1 0.04],'String',...
    sprintf('QSP-DCM | %s | v%s', run_timestamp, script_version),...
    'EdgeColor','none','HorizontalAlignment','center',...
    'FontSize',7,'Color',[0.5 0.5 0.5]);

save_figure(fig4, fullfile(outdir,'Drug_MHT_Spotlight'));
fprintf('  Saved: Drug_MHT_Spotlight\n');

%% --- Export Results Table -----------------------------------------------
fprintf('Exporting drug results to Excel...\n');

xlsx_file = fullfile(outdir,'Drug_AllStates_Results.xlsx');

% Sheet per drug: rows = states, cols = vehicle vs drug vs delta
for d = 1:n_drugs
    rows = {};
    for s = 1:n_states
        rows(end+1,:) = {state_labels{s},...
            vehicle(s).CFR,    drug_res(d,s).CFR,    drug_res(d,s).dCFR,...
            vehicle(s).IMR,    drug_res(d,s).IMR,    drug_res(d,s).dIMR,...
            vehicle(s).Ee,     drug_res(d,s).Ee,     drug_res(d,s).dEe,...
            vehicle(s).PGC1a,  drug_res(d,s).PGC1a,  drug_res(d,s).dPGC,...
            vehicle(s).Ceramide,drug_res(d,s).Ceramide,drug_res(d,s).dCer,...
            vehicle(s).CMD,    drug_res(d,s).CMD,    drug_res(d,s).dCMD};
    end
    headers = {'State',...
        'CFR_veh','CFR_drug','CFR_delta%',...
        'IMR_veh','IMR_drug','IMR_delta%',...
        'Ee_veh','Ee_drug','Ee_delta%',...
        'PGC1a_veh','PGC1a_drug','PGC1a_delta%',...
        'Cer_veh','Cer_drug','Cer_delta%',...
        'CMD_veh','CMD_drug','CMD_delta%'};
    T = cell2table(rows,'VariableNames',headers);
    writetable(T, xlsx_file, 'Sheet', drugs(d).abbrev);
end

% Metadata sheet
meta = {'Field','Value';
    'Script','Script5_DrugInterventions.m';
    'Version',script_version;
    'Timestamp',run_timestamp;
    'Simulation horizon','336 hours (14 days)';
    'Drugs','Empagliflozin, Semaglutide, MHT, Sacubitril/Valsartan';
    'States','6 sex-disease states';
    'Reference female_post_T2DM vehicle CFR', num2str(vehicle(4).CFR,'%.3f');
    'Reference female_post_T2DM vehicle IMR', num2str(vehicle(4).IMR,'%.1f');
    'Reference male_T2DM vehicle CFR', num2str(vehicle(6).CFR,'%.3f');
    'Reference male_T2DM vehicle IMR', num2str(vehicle(6).IMR,'%.1f')};
writecell(meta, xlsx_file, 'Sheet','Metadata');
fprintf('  Saved: %s\n', xlsx_file);

%% --- Print Summary -------------------------------------------------------
fprintf('\n=======================================================\n');
fprintf(' DRUG INTERVENTION SUMMARY — KEY STATES\n');
fprintf('=======================================================\n');
fprintf('%-22s  %-8s  %6s  %6s  %6s  %6s\n',...
    'State','Drug','CFR','dCFR%','IMR','dIMR%');
fprintf('%s\n',repmat('-',1,72));
for d = 1:n_drugs
    for s = [4,6]
        fprintf('%-22s  %-8s  %6.3f  %+5.1f%%  %6.1f  %+5.1f%%\n',...
            state_labels{s}, drugs(d).abbrev,...
            drug_res(d,s).CFR, drug_res(d,s).dCFR,...
            drug_res(d,s).IMR, drug_res(d,s).dIMR);
    end
end
fprintf('\n');

fprintf(flog,'\nDrug simulations complete: %s\n', datestr(now));
fclose(flog);
fprintf('Script 5 complete. All drug intervention files saved.\n\n');

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function res = run_model(sex_state, drug, tspan, options)
    [p1,p2,p3] = get_all_params(sex_state, drug);
    opts = odeset('RelTol',1e-7,'AbsTol',1e-9,'MaxStep',0.5);

    y0_1 = [p1.E2_0;p1.ERb_0;p1.PGC1a_0;p1.Mito_0;p1.ROS_0;p1.MnSOD_0;p1.DeltaPsi_0];
    [~,y1] = ode15s(@(t,y) L1_odes(t,y,p1), tspan, y0_1, opts);
    p2.PGC1a=y1(end,3);p2.Mito=y1(end,4);p2.ROS=y1(end,5);
    p2.dPsi=y1(end,7);p2.MnSOD=y1(end,6);

    y0_2 = [p2.AMPK_0;p2.FAO_0;p2.GlucOx_0;p2.Cer_0;p2.ColX_0;p2.Ee_0];
    [~,y2] = ode15s(@(t,y) L2_odes(t,y,p2), tspan, y0_2, opts);
    p3.Ceramide=y2(end,4);p3.ColX=y2(end,5);p3.Ee=y2(end,6);
    p3.ROS=p2.ROS;p3.dPsi=p2.dPsi;p3.PGC1a=p2.PGC1a;p3.Mito=p2.Mito;

    y0_3 = [p3.NO_0;p3.ET1_0;p3.MRes_0;p3.CFR_0;p3.IMR_0;p3.WS_0];
    [~,y3] = ode15s(@(t,y) L3_odes(t,y,p3), tspan, y0_3, opts);

    res.CFR      = y3(end,4);
    res.IMR      = y3(end,5);
    res.WS       = y3(end,6);
    res.Ee       = y2(end,6);
    res.PGC1a    = y1(end,3);
    res.Ceramide = y2(end,4);
    cmd_cfr  = max(0, (2.5 - res.CFR)/2.5);
    cmd_imr  = max(0, (res.IMR - 10)/40);
    cmd_ee   = max(0, (res.Ee - 8)/20);
    res.CMD  = cmd_cfr + cmd_imr + cmd_ee;
end

function [tc, tvec] = run_timecourse(sex_state, drug, tspan, options)
    [p1,p2,p3] = get_all_params(sex_state, drug);
    opts = odeset('RelTol',1e-7,'AbsTol',1e-9,'MaxStep',0.5);

    y0_1 = [p1.E2_0;p1.ERb_0;p1.PGC1a_0;p1.Mito_0;p1.ROS_0;p1.MnSOD_0;p1.DeltaPsi_0];
    [~,y1] = ode15s(@(t,y) L1_odes(t,y,p1), tspan, y0_1, opts);
    p2.PGC1a=y1(end,3);p2.Mito=y1(end,4);p2.ROS=y1(end,5);
    p2.dPsi=y1(end,7);p2.MnSOD=y1(end,6);

    y0_2 = [p2.AMPK_0;p2.FAO_0;p2.GlucOx_0;p2.Cer_0;p2.ColX_0;p2.Ee_0];
    [~,y2] = ode15s(@(t,y) L2_odes(t,y,p2), tspan, y0_2, opts);
    p3.Ceramide=y2(end,4);p3.ColX=y2(end,5);p3.Ee=y2(end,6);
    p3.ROS=p2.ROS;p3.dPsi=p2.dPsi;p3.PGC1a=p2.PGC1a;p3.Mito=p2.Mito;

    y0_3 = [p3.NO_0;p3.ET1_0;p3.MRes_0;p3.CFR_0;p3.IMR_0;p3.WS_0];
    [t3,y3] = ode15s(@(t,y) L3_odes(t,y,p3), tspan, y0_3, opts);

    tvec = t3;
    tc.CFR = y3(:,4);
    tc.IMR = y3(:,5);
    tc.WS  = y3(:,6);
    tc.NO  = y3(:,1);
    tc.Ee  = repmat(y2(end,6), length(t3), 1);
end

function [p1,p2,p3] = get_all_params(sex_state, drug)
    p1 = get_L1_params(sex_state);
    if ~isempty(drug) && isfield(drug,'L1')
        fn = fieldnames(drug.L1);
        for k=1:length(fn)
            f = fn{k};
            if strcmp(f,'E2_ss_mult')
                p1.E2_ss = p1.E2_ss * drug.L1.(f);
                p1.E2_0  = p1.E2_ss;
            elseif isfield(p1,f)
                p1.(f) = p1.(f) * drug.L1.(f);
            end
        end
    end
    p2 = get_L2_params(sex_state, p1);
    if ~isempty(drug) && isfield(drug,'L2')
        fn = fieldnames(drug.L2);
        for k=1:length(fn)
            f = fn{k};
            if isfield(p2,f)
                p2.(f) = p2.(f) * drug.L2.(f);
            end
        end
    end
    p3 = get_L3_params(sex_state, p2);
    if ~isempty(drug) && isfield(drug,'L3')
        fn = fieldnames(drug.L3);
        for k=1:length(fn)
            f = fn{k};
            if isfield(p3,f)
                p3.(f) = p3.(f) * drug.L3.(f);
            end
        end
    end
    if ~isempty(drug) && isfield(drug,'L3_eNOS_boost')
        p3.eNOS_activity = min(1.0, p3.eNOS_activity * drug.L3_eNOS_boost);
    end
end

function p = get_L1_params(sex_state)
    p.k_ERb_on=0.05;p.k_ERb_off=0.10;p.k_PGC_basal=0.08;
    p.k_PGC_ERb=0.25;p.k_PGC_deg=0.15;p.k_Mito_form=0.12;
    p.k_Mito_deg=0.05;p.k_ROS_prod=0.20;p.k_ROS_mito=0.15;
    p.k_ROS_scav=0.30;p.k_MnSOD_basal=2.00;p.k_MnSOD_ROS=8.00;
    p.k_MnSOD_PGC=5.00;p.k_MnSOD_deg=0.05;p.k_Psi_form=3.00;
    p.k_Psi_ROS=1.50;p.k_Psi_deg=0.08;
    switch sex_state
        case 'female_pre_healthy'
            p.E2_ss=150;p.k_ERb_sens=1.00;p.T2DM_factor=1.00;
            p.mito_stress=1.00;p.ROS_drive=1.00;
            p.E2_0=150;p.ERb_0=0.20;p.PGC1a_0=0.50;p.Mito_0=30;
            p.ROS_0=1.00;p.MnSOD_0=60;p.DeltaPsi_0=1.00;
        case 'female_pre_T2DM'
            p.E2_ss=150;p.k_ERb_sens=0.65;p.T2DM_factor=1.60;
            p.mito_stress=1.40;p.ROS_drive=1.80;
            p.E2_0=150;p.ERb_0=0.15;p.PGC1a_0=0.40;p.Mito_0=25;
            p.ROS_0=1.50;p.MnSOD_0=70;p.DeltaPsi_0=0.80;
        case 'female_post_healthy'
            p.E2_ss=15;p.k_ERb_sens=0.80;p.T2DM_factor=1.00;
            p.mito_stress=1.20;p.ROS_drive=1.20;
            p.E2_0=15;p.ERb_0=0.15;p.PGC1a_0=0.40;p.Mito_0=25;
            p.ROS_0=1.20;p.MnSOD_0=65;p.DeltaPsi_0=0.90;
        case 'female_post_T2DM'
            p.E2_ss=5;p.k_ERb_sens=0.50;p.T2DM_factor=1.80;
            p.mito_stress=1.70;p.ROS_drive=2.20;
            p.E2_0=5;p.ERb_0=0.10;p.PGC1a_0=0.35;p.Mito_0=20;
            p.ROS_0=1.80;p.MnSOD_0=75;p.DeltaPsi_0=0.70;
        case 'male_healthy'
            p.E2_ss=30;p.k_ERb_sens=0.60;p.T2DM_factor=1.00;
            p.mito_stress=1.00;p.ROS_drive=1.10;
            p.E2_0=30;p.ERb_0=0.12;p.PGC1a_0=0.45;p.Mito_0=28;
            p.ROS_0=1.10;p.MnSOD_0=55;p.DeltaPsi_0=0.95;
        case 'male_T2DM'
            p.E2_ss=30;p.k_ERb_sens=0.45;p.T2DM_factor=1.70;
            p.mito_stress=1.60;p.ROS_drive=2.00;
            p.E2_0=30;p.ERb_0=0.10;p.PGC1a_0=0.35;p.Mito_0=22;
            p.ROS_0=1.70;p.MnSOD_0=72;p.DeltaPsi_0=0.75;
        otherwise, error('Unknown state: %s',sex_state);
    end
end

function dydt = L1_odes(~,y,p)
    E2=max(y(1),0);ERb=max(y(2),0);PGC1a=max(y(3),0);
    Mito=max(y(4),0);ROS=max(y(5),0);MnSOD=max(y(6),0);DeltaPsi=max(y(7),0);
    dE2=p.k_ERb_on*(p.E2_ss-E2);
    dERb=p.k_ERb_on*p.k_ERb_sens*E2*(1-ERb)-p.k_ERb_off*ERb;
    dPGC1a=p.k_PGC_basal+p.k_PGC_ERb*ERb-p.k_PGC_deg*p.T2DM_factor*PGC1a;
    dMito=p.k_Mito_form*PGC1a*30-p.k_Mito_deg*p.mito_stress*Mito;
    ROS_prod=p.k_ROS_prod*p.ROS_drive+p.k_ROS_mito*p.T2DM_factor*(1/(PGC1a+0.1));
    ROS_scav=p.k_ROS_scav*MnSOD/100*ROS;
    dROS=ROS_prod-ROS_scav;
    dMnSOD=p.k_MnSOD_basal+p.k_MnSOD_ROS*ROS+p.k_MnSOD_PGC*PGC1a-p.k_MnSOD_deg*MnSOD;
    dDeltaPsi=p.k_Psi_form*Mito/30*PGC1a-p.k_Psi_ROS*ROS*DeltaPsi-p.k_Psi_deg*DeltaPsi;
    dydt=[dE2;dERb;dPGC1a;dMito;dROS;dMnSOD;dDeltaPsi];
end

function p = get_L2_params(sex_state, L1)
    p.PGC1a=L1.PGC1a_0;p.Mito=L1.Mito_0;p.ROS=L1.ROS_0;
    p.dPsi=L1.DeltaPsi_0;p.MnSOD=L1.MnSOD_0;
    p.k_AMPK_basal=0.10;p.k_AMPK_dPsi=0.08;p.k_AMPK_ROS=0.06;p.k_AMPK_deg=0.12;
    p.dPsi_ref=50.0;p.k_FAO_AMPK=0.20;p.k_FAO_PGC=0.15;p.k_FAO_deg=0.10;
    p.k_FAO_cer_inh=0.05;p.k_Gluc_basal=0.12;p.k_Gluc_ins_res=0.08;
    p.k_Gluc_deg=0.10;p.k_Cer_FA=0.08;p.k_Cer_ROS=0.10;p.k_Cer_deg=0.04;
    p.Cer_ref=6.0;p.k_Col_ROS=0.06;p.k_Col_cer=0.04;p.k_Col_deg=0.02;
    p.Col_ref=1.0;p.k_Ee_Col=0.10;p.k_Ee_Cer=0.06;p.k_Ee_dPsi=0.04;p.k_Ee_deg=0.03;
    switch sex_state
        case 'female_pre_healthy'
            p.T2DM_ins_res=1.00;p.FA_overflow=1.00;p.estrogen_prot=1.00;
            p.AMPK_0=1.0;p.FAO_0=1.0;p.GlucOx_0=1.0;p.Cer_0=4.0;p.ColX_0=1.0;p.Ee_0=7.0;
        case 'female_pre_T2DM'
            p.T2DM_ins_res=2.20;p.FA_overflow=1.80;p.estrogen_prot=0.65;
            p.AMPK_0=0.9;p.FAO_0=1.2;p.GlucOx_0=0.7;p.Cer_0=7.0;p.ColX_0=1.2;p.Ee_0=10.0;
        case 'female_post_healthy'
            p.T2DM_ins_res=1.30;p.FA_overflow=1.30;p.estrogen_prot=0.70;
            p.AMPK_0=0.95;p.FAO_0=1.0;p.GlucOx_0=0.9;p.Cer_0=5.5;p.ColX_0=1.1;p.Ee_0=9.0;
        case 'female_post_T2DM'
            p.T2DM_ins_res=2.60;p.FA_overflow=2.20;p.estrogen_prot=0.35;
            p.AMPK_0=0.7;p.FAO_0=1.5;p.GlucOx_0=0.5;p.Cer_0=10.0;p.ColX_0=1.5;p.Ee_0=13.0;
        case 'male_healthy'
            p.T2DM_ins_res=1.10;p.FA_overflow=1.10;p.estrogen_prot=0.60;
            p.AMPK_0=1.0;p.FAO_0=1.0;p.GlucOx_0=1.0;p.Cer_0=4.5;p.ColX_0=1.0;p.Ee_0=7.5;
        case 'male_T2DM'
            p.T2DM_ins_res=2.10;p.FA_overflow=1.70;p.estrogen_prot=0.45;
            p.AMPK_0=0.8;p.FAO_0=1.3;p.GlucOx_0=0.65;p.Cer_0=8.5;p.ColX_0=1.3;p.Ee_0=11.5;
        otherwise, error('Unknown state: %s',sex_state);
    end
end

function dydt = L2_odes(~,y,p)
    AMPK=max(y(1),0);FAO=max(y(2),0);GlucOx=max(y(3),0);
    Cer=max(y(4),0);ColX=max(y(5),0);Ee=max(y(6),0);
    dPsi_norm=max(1-p.dPsi/p.dPsi_ref,0);
    dAMPK=p.k_AMPK_basal+p.k_AMPK_dPsi*dPsi_norm+p.k_AMPK_ROS*p.ROS-p.k_AMPK_deg*AMPK;
    dFAO=p.k_FAO_AMPK*AMPK*p.FA_overflow+p.k_FAO_PGC*p.PGC1a-p.k_FAO_deg*FAO-p.k_FAO_cer_inh*Cer*FAO;
    dGlucOx=p.k_Gluc_basal-p.k_Gluc_ins_res*p.T2DM_ins_res*GlucOx-p.k_Gluc_deg*GlucOx*FAO/(FAO+0.5);
    dCer=p.k_Cer_FA*p.FA_overflow*(1-p.estrogen_prot*0.4)+p.k_Cer_ROS*p.ROS-p.k_Cer_deg*Cer;
    dColX=p.k_Col_ROS*p.ROS*(1-p.estrogen_prot*0.3)+p.k_Col_cer*Cer/10-p.k_Col_deg*ColX;
    dEe=p.k_Ee_Col*ColX+p.k_Ee_Cer*Cer/10+p.k_Ee_dPsi*dPsi_norm*5-p.k_Ee_deg*Ee;
    dydt=[dAMPK;dFAO;dGlucOx;dCer;dColX;dEe];
end

function p = get_L3_params(sex_state, L2)
    p.Ceramide=L2.Cer_0;p.ROS=L2.ROS;p.ColX=L2.ColX_0;
    p.Ee=L2.Ee_0;p.dPsi=L2.dPsi;p.PGC1a=L2.PGC1a;p.Mito=L2.Mito;
    p.k_NO_basal=0.50;p.k_NO_cer_inh=0.08;p.k_NO_ROS_scav=0.20;p.k_NO_deg=0.15;
    p.k_ET1_basal=0.05;p.k_ET1_ROS=0.12;p.k_ET1_cer=0.08;p.k_ET1_NO_inh=0.10;p.k_ET1_deg=0.08;
    p.k_MRes_ET1=0.15;p.k_MRes_colx=0.10;p.k_MRes_NO_dil=0.20;p.k_MRes_deg=0.05;
    p.k_CFR_MRes=0.20;p.k_CFR_dPsi=0.10;p.k_CFR_deg=0.30;
    p.k_IMR_MRes=0.12;p.k_IMR_Ee=0.04;p.k_IMR_deg=0.03;
    p.k_WS_Ee=0.04;p.k_WS_IMR=0.002;p.k_WS_colx=0.06;p.k_WS_deg=0.03;
    switch sex_state
        case 'female_pre_healthy'
            p.eNOS_activity=1.00;p.vasc_estrogen=1.00;p.ET1_basal_mult=1.00;
            p.NO_0=0.80;p.ET1_0=1.00;p.MRes_0=1.00;p.CFR_0=3.20;p.IMR_0=15.0;p.WS_0=1.00;
        case 'female_pre_T2DM'
            p.eNOS_activity=0.65;p.vasc_estrogen=0.65;p.ET1_basal_mult=1.50;
            p.NO_0=0.55;p.ET1_0=1.40;p.MRes_0=1.30;p.CFR_0=2.40;p.IMR_0=28.0;p.WS_0=1.30;
        case 'female_post_healthy'
            p.eNOS_activity=0.75;p.vasc_estrogen=0.70;p.ET1_basal_mult=1.20;
            p.NO_0=0.65;p.ET1_0=1.15;p.MRes_0=1.15;p.CFR_0=2.80;p.IMR_0=20.0;p.WS_0=1.10;
        case 'female_post_T2DM'
            p.eNOS_activity=0.35;p.vasc_estrogen=0.30;p.ET1_basal_mult=2.00;
            p.NO_0=0.30;p.ET1_0=1.80;p.MRes_0=1.70;p.CFR_0=1.80;p.IMR_0=42.0;p.WS_0=1.70;
        case 'male_healthy'
            p.eNOS_activity=0.80;p.vasc_estrogen=0.60;p.ET1_basal_mult=1.05;
            p.NO_0=0.72;p.ET1_0=1.05;p.MRes_0=1.05;p.CFR_0=3.00;p.IMR_0=16.0;p.WS_0=1.05;
        case 'male_T2DM'
            p.eNOS_activity=0.55;p.vasc_estrogen=0.45;p.ET1_basal_mult=1.60;
            p.NO_0=0.48;p.ET1_0=1.55;p.MRes_0=1.45;p.CFR_0=2.10;p.IMR_0=32.0;p.WS_0=1.45;
        otherwise, error('Unknown state: %s',sex_state);
    end
end

function dydt = L3_odes(~,y,p)
    NO=max(y(1),0.001);ET1=max(y(2),0);MRes=max(y(3),0);
    CFR=max(y(4),0);IMR=max(y(5),0);WS=max(y(6),0);
    dPsi_def=max(1-p.dPsi/50,0);
    dNO=p.k_NO_basal*p.eNOS_activity*p.vasc_estrogen-p.k_NO_cer_inh*p.Ceramide*NO-p.k_NO_ROS_scav*p.ROS*NO-p.k_NO_deg*NO;
    dET1=p.k_ET1_basal*p.ET1_basal_mult+p.k_ET1_ROS*p.ROS+p.k_ET1_cer*p.Ceramide/10-p.k_ET1_NO_inh*NO*ET1-p.k_ET1_deg*ET1;
    dMRes=p.k_MRes_ET1*ET1+p.k_MRes_colx*p.ColX-p.k_MRes_NO_dil*NO*MRes-p.k_MRes_deg*MRes;
    dCFR=-p.k_CFR_MRes*max(MRes-1.0,0)*CFR-p.k_CFR_dPsi*dPsi_def*CFR+p.k_CFR_deg*(3.5-CFR);
    dIMR=p.k_IMR_MRes*MRes+p.k_IMR_Ee*p.Ee-p.k_IMR_deg*IMR;
    dWS=p.k_WS_Ee*p.Ee/10+p.k_WS_IMR*IMR+p.k_WS_colx*p.ColX-p.k_WS_deg*WS;
    dydt=[dNO;dET1;dMRes;dCFR;dIMR;dWS];
end

function save_figure(fig, basepath)
    print(fig, basepath, '-dpng', '-r300');
end
