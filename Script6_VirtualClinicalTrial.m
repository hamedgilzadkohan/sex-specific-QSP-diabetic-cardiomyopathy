%% ========================================================================
%  Script6_VirtualClinicalTrial.m
%  Sex-Specific QSP Model of Diabetic Cardiomyopathy
%  VIRTUAL CLINICAL TRIAL — Sex-Stratified Drug Response
%
%  Authors : Soheili M, Gilzad-Kohan H, Lotfi AS
%
%  Trial Design:
%    - 4 treatment arms: SGLT2i, GLP-1, MHT, ARNI
%    - 1 combination arm: SGLT2i + MHT (best combination for women)
%    - 1 combination arm: SGLT2i + GLP-1 (best combination for men)
%    - All 6 sex-disease states treated as virtual patient subgroups
%    - Primary endpoint: CFR normalization (≥2.5) and IMR normalization (<25U)
%    - Secondary: E/e' improvement, CMD score reduction
%
%  Output files:
%    - VCT_PrimaryEndpoints.png     (trial primary endpoint figure)
%    - VCT_ResponderAnalysis.png    (responder rates by sex/state)
%    - VCT_CombinationTherapy.png  (combination vs monotherapy)
%    - VCT_NNT_Table.png           (number needed to treat)
%    - VCT_Results.xlsx
%    - VCT_RunLog.txt
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

flog = fopen(fullfile(outdir,'VCT_RunLog.txt'),'w');
fprintf(flog,'=======================================================\n');
fprintf(flog,' Virtual Clinical Trial — Reproducibility Log\n');
fprintf(flog,'=======================================================\n');
fprintf(flog,' Version   : %s\n', script_version);
fprintf(flog,' Timestamp : %s\n', run_timestamp);
fprintf(flog,' MATLAB    : %s\n', version);
fprintf(flog,'=======================================================\n\n');

fprintf('=======================================================\n');
fprintf(' Script 6: Virtual Clinical Trial\n');
fprintf(' %s\n', run_timestamp);
fprintf('=======================================================\n\n');

%% --- Settings ------------------------------------------------------------
tspan   = [0 336];
options = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.1);

state_ids    = {'female_pre_healthy','female_pre_T2DM',...
                'female_post_healthy','female_post_T2DM',...
                'male_healthy','male_T2DM'};
state_labels = {'Female Pre Healthy','Female Pre T2DM',...
                'Female Post Healthy','Female Post T2DM',...
                'Male Healthy','Male T2DM'};
n_states = length(state_ids);

% Responder thresholds (clinical trial endpoints)
CFR_NORMAL   = 2.5;   % CFR ≥ 2.5 = normalized
IMR_NORMAL   = 25;    % IMR < 25U = normalized
EE_IMPROVE   = 2.0;   % E/e' absolute improvement ≥ 2 units = clinically meaningful

%% --- Define treatment arms -----------------------------------------------
% Monotherapy arms
arms(1).name   = 'Vehicle';
arms(1).abbrev = 'Veh';
arms(1).color  = [0.6 0.6 0.6];
arms(1).drug   = [];

arms(2).name   = 'Empagliflozin';
arms(2).abbrev = 'SGLT2i';
arms(2).color  = [0.20 0.60 0.85];
arms(2).drug   = build_drug_sglt2i();

arms(3).name   = 'Semaglutide';
arms(3).abbrev = 'GLP-1';
arms(3).color  = [0.20 0.75 0.45];
arms(3).drug   = build_drug_glp1();

arms(4).name   = 'MHT (Estrogen)';
arms(4).abbrev = 'MHT';
arms(4).color  = [0.90 0.35 0.65];
arms(4).drug   = build_drug_mht();

arms(5).name   = 'ARNI';
arms(5).abbrev = 'ARNI';
arms(5).color  = [0.85 0.55 0.15];
arms(5).drug   = build_drug_arni();

arms(6).name   = 'SGLT2i + MHT';
arms(6).abbrev = 'S+M';
arms(6).color  = [0.55 0.20 0.70];
arms(6).drug   = combine_drugs(build_drug_sglt2i(), build_drug_mht());

arms(7).name   = 'SGLT2i + GLP-1';
arms(7).abbrev = 'S+G';
arms(7).color  = [0.10 0.50 0.60];
arms(7).drug   = combine_drugs(build_drug_sglt2i(), build_drug_glp1());

n_arms = length(arms);

%% --- Run all arms for all states -----------------------------------------
fprintf('Running virtual clinical trial (%d arms x %d states)...\n', n_arms, n_states);

results = struct();
for a = 1:n_arms
    fprintf('  Arm %d/%d: %s\n', a, n_arms, arms(a).name);
    for s = 1:n_states
        r = run_model(state_ids{s}, arms(a).drug, tspan, options);
        results(a,s).CFR      = r.CFR;
        results(a,s).IMR      = r.IMR;
        results(a,s).Ee       = r.Ee;
        results(a,s).PGC1a    = r.PGC1a;
        results(a,s).Ceramide = r.Ceramide;
        results(a,s).CMD      = r.CMD;
        % Responder classification
        results(a,s).CFR_resp  = (r.CFR >= CFR_NORMAL);
        results(a,s).IMR_resp  = (r.IMR < IMR_NORMAL);
        % Improvement vs vehicle
        results(a,s).dCFR     = r.CFR  - results(1,s).CFR;
        results(a,s).dIMR     = r.IMR  - results(1,s).IMR;
        results(a,s).dEe      = r.Ee   - results(1,s).Ee;
        results(a,s).dCMD     = r.CMD  - results(1,s).CMD;
    end
end
fprintf('All trial arms complete.\n\n');

%% --- Figure 1: Primary Endpoints (CFR & IMR) ----------------------------
fprintf('Generating primary endpoint figure...\n');

fig1 = figure('Name','VCT Primary Endpoints',...
    'Units','inches','Position',[0.5 0.5 18 12],'Color','white');

% Top row: CFR; Bottom row: IMR
% Columns: T2DM states only (the ones that matter clinically)
t2dm_states = [2,4,6];  % FPT2, FPoT2, MT2
t2dm_labels = {'Female Pre T2DM','Female Post T2DM','Male T2DM'};

for col = 1:3
    s = t2dm_states(col);

    % CFR subplot
    subplot(2,3,col);
    hold on;
    cfr_vals = arrayfun(@(a) results(a,s).CFR, 1:n_arms);
    colors_bar = reshape([arms.color],[3,n_arms])';
    for a = 1:n_arms
        b = bar(a, cfr_vals(a), 0.7, 'FaceColor', arms(a).color, ...
            'EdgeColor', arms(a).color*0.7, 'LineWidth', 0.8);
        if results(a,s).CFR_resp
            text(a, cfr_vals(a)+0.05, '✓', 'HorizontalAlignment','center',...
                'FontSize',12,'Color',[0.1 0.6 0.1],'FontWeight','bold');
        end
    end
    yline(CFR_NORMAL,'g--','LineWidth',1.5,'Label','Normal (≥2.5)');
    yline(2.0,'r--','LineWidth',1.2,'Label','CMD (<2.0)');
    set(gca,'XTick',1:n_arms,'XTickLabel',{arms.abbrev},...
        'XTickLabelRotation',30,'FontSize',8,'TickDir','out');
    ylabel('CFR','FontSize',9); ylim([0 max(cfr_vals)*1.35]);
    title(sprintf('CFR — %s',t2dm_labels{col}),'FontSize',9,'FontWeight','bold');
    box on; grid on; hold off;

    % IMR subplot
    subplot(2,3,col+3);
    hold on;
    imr_vals = arrayfun(@(a) results(a,s).IMR, 1:n_arms);
    for a = 1:n_arms
        bar(a, imr_vals(a), 0.7, 'FaceColor', arms(a).color,...
            'EdgeColor', arms(a).color*0.7, 'LineWidth', 0.8);
        if results(a,s).IMR_resp
            text(a, imr_vals(a)+0.5, '✓', 'HorizontalAlignment','center',...
                'FontSize',12,'Color',[0.1 0.6 0.1],'FontWeight','bold');
        end
    end
    yline(IMR_NORMAL,'r--','LineWidth',1.5,'Label','Dysfunction (>25U)');
    yline(40,'Color',[0.8 0.2 0],'LineStyle','--','LineWidth',1.2,...
        'Label','Severe (>40U)');
    set(gca,'XTick',1:n_arms,'XTickLabel',{arms.abbrev},...
        'XTickLabelRotation',30,'FontSize',8,'TickDir','out');
    ylabel('IMR (U)','FontSize',9);
    title(sprintf('IMR — %s',t2dm_labels{col}),'FontSize',9,'FontWeight','bold');
    box on; grid on; hold off;
end

sgtitle({'Virtual Clinical Trial — Primary Endpoints (CFR & IMR)',...
    '✓ = Responder (endpoint normalized) | QSP Model of Diabetic Cardiomyopathy'},...
    'FontSize',11,'FontWeight','bold');
annotation('textbox',[0 0.96 1 0.04],'String',...
    sprintf('QSP-DCM | %s | v%s', run_timestamp, script_version),...
    'EdgeColor','none','HorizontalAlignment','center',...
    'FontSize',7,'Color',[0.5 0.5 0.5]);

save_figure(fig1, fullfile(outdir,'VCT_PrimaryEndpoints'));
fprintf('  Saved: VCT_PrimaryEndpoints\n');

%% --- Figure 2: Responder Analysis ----------------------------------------
fprintf('Generating responder analysis figure...\n');

fig2 = figure('Name','VCT Responder Analysis',...
    'Units','inches','Position',[0.5 0.5 16 10],'Color','white');

% Show CFR and IMR responder rates for each arm, broken down by sex group
% Female T2DM states: 2,4  |  Male T2DM states: 6
female_t2dm = [2,4];
male_t2dm   = [6];

subplot(1,2,1);  % CFR responders
hold on;
cfr_resp_f = zeros(1,n_arms);
cfr_resp_m = zeros(1,n_arms);
for a = 1:n_arms
    cfr_resp_f(a) = mean([results(a,female_t2dm).CFR_resp]) * 100;
    cfr_resp_m(a) = results(a,male_t2dm).CFR_resp * 100;
end
x = 1:n_arms;
b1 = bar(x-0.2, cfr_resp_f, 0.35, 'FaceColor',[0.55 0.75 0.95],...
    'EdgeColor',[0.20 0.45 0.80],'LineWidth',1.0);
b2 = bar(x+0.2, cfr_resp_m, 0.35, 'FaceColor',[0.85 0.55 0.55],...
    'EdgeColor',[0.70 0.15 0.15],'LineWidth',1.0);
set(gca,'XTick',1:n_arms,'XTickLabel',{arms.abbrev},...
    'XTickLabelRotation',30,'FontSize',9,'TickDir','out');
ylabel('CFR Responder Rate (%)','FontSize',10);
title('CFR Normalization Rate by Sex','FontSize',10,'FontWeight','bold');
legend([b1,b2],{'Females (T2DM)','Male T2DM'},'Location','northwest','FontSize',9);
ylim([0 115]); box on; grid on; hold off;

subplot(1,2,2);  % IMR responders
hold on;
imr_resp_f = zeros(1,n_arms);
imr_resp_m = zeros(1,n_arms);
for a = 1:n_arms
    imr_resp_f(a) = mean([results(a,female_t2dm).IMR_resp]) * 100;
    imr_resp_m(a) = results(a,male_t2dm).IMR_resp * 100;
end
b3 = bar(x-0.2, imr_resp_f, 0.35, 'FaceColor',[0.55 0.75 0.95],...
    'EdgeColor',[0.20 0.45 0.80],'LineWidth',1.0);
b4 = bar(x+0.2, imr_resp_m, 0.35, 'FaceColor',[0.85 0.55 0.55],...
    'EdgeColor',[0.70 0.15 0.15],'LineWidth',1.0);
set(gca,'XTick',1:n_arms,'XTickLabel',{arms.abbrev},...
    'XTickLabelRotation',30,'FontSize',9,'TickDir','out');
ylabel('IMR Responder Rate (%)','FontSize',10);
title('IMR Normalization Rate by Sex','FontSize',10,'FontWeight','bold');
legend([b3,b4],{'Females (T2DM)','Male T2DM'},'Location','northwest','FontSize',9);
ylim([0 115]); box on; grid on; hold off;

sgtitle({'Virtual Clinical Trial — Sex-Stratified Responder Analysis',...
    'Responder = endpoint within normal range | Females include Pre+Post T2DM'},...
    'FontSize',11,'FontWeight','bold');
annotation('textbox',[0 0.96 1 0.04],'String',...
    sprintf('QSP-DCM | %s | v%s', run_timestamp, script_version),...
    'EdgeColor','none','HorizontalAlignment','center',...
    'FontSize',7,'Color',[0.5 0.5 0.5]);

save_figure(fig2, fullfile(outdir,'VCT_ResponderAnalysis'));
fprintf('  Saved: VCT_ResponderAnalysis\n');

%% --- Figure 3: Combination Therapy Deep Dive ----------------------------
fprintf('Generating combination therapy figure...\n');

fig3 = figure('Name','Combination Therapy vs Monotherapy',...
    'Units','inches','Position',[0.5 0.5 16 12],'Color','white');

% Focus: Female Post T2DM (worst state) — compare all arms
s = 4;  % female_post_T2DM
metrics_combo = {'CFR','IMR','Ee','CMD'};
metric_ylabels = {'CFR','IMR (U)','E/e''','CMD Score'};
thresh_combo   = {2.5, 25, 8, 1.5};
thresh_dir_combo = {'above','below','below','below'};

for mi = 1:4
    subplot(2,2,mi); hold on;
    vals = arrayfun(@(a) results(a,s).(metrics_combo{mi}), 1:n_arms);

    for a = 1:n_arms
        is_combo = (a >= 6);
        edge_w = 1.5 + is_combo*0.8;
        b = bar(a, vals(a), 0.72, 'FaceColor', arms(a).color,...
            'EdgeColor', arms(a).color*0.6, 'LineWidth', edge_w);
        if is_combo
            % Hatching indicator for combo
            text(a, vals(a) * 0.5, 'COMBO', 'HorizontalAlignment','center',...
                'FontSize',6,'Color','white','FontWeight','bold','Rotation',90);
        end
    end

    if ~isnan(thresh_combo{mi})
        yline(thresh_combo{mi},'r--','LineWidth',1.5);
    end
    set(gca,'XTick',1:n_arms,'XTickLabel',{arms.abbrev},...
        'XTickLabelRotation',30,'FontSize',8,'TickDir','out');
    ylabel(metric_ylabels{mi},'FontSize',10);
    title(sprintf('%s — Female Post T2DM',metric_ylabels{mi}),...
        'FontSize',10,'FontWeight','bold');
    box on; grid on; hold off;
end

sgtitle({'Combination Therapy vs Monotherapy — Female Post T2DM',...
    'SGLT2i+MHT and SGLT2i+GLP-1 vs Individual Agents'},...
    'FontSize',11,'FontWeight','bold');
annotation('textbox',[0 0.96 1 0.04],'String',...
    sprintf('QSP-DCM | %s | v%s', run_timestamp, script_version),...
    'EdgeColor','none','HorizontalAlignment','center',...
    'FontSize',7,'Color',[0.5 0.5 0.5]);

save_figure(fig3, fullfile(outdir,'VCT_CombinationTherapy'));
fprintf('  Saved: VCT_CombinationTherapy\n');

%% --- Figure 4: Treatment Landscape Heatmap --------------------------------
fprintf('Generating treatment landscape heatmap...\n');

fig4 = figure('Name','Treatment Landscape',...
    'Units','inches','Position',[0.5 0.5 14 8],'Color','white');

% Heatmap: rows = treatment arms, cols = states
% Color = CFR value (green = normal, red = CMD)
cfr_mat = zeros(n_arms, n_states);
imr_mat = zeros(n_arms, n_states);
for a = 1:n_arms
    for s = 1:n_states
        cfr_mat(a,s) = results(a,s).CFR;
        imr_mat(a,s) = results(a,s).IMR;
    end
end

subplot(1,2,1);
imagesc(cfr_mat);
colormap(gca, green_red_map());
colorbar;
caxis([0 4]);
set(gca,'XTick',1:n_states,'XTickLabel',...
    {'FPH','FPT2','FPoH','FPoT2','MH','MT2'},...
    'YTick',1:n_arms,'YTickLabel',{arms.abbrev},...
    'FontSize',9,'TickDir','out');
xlabel('Sex-Disease State','FontSize',10);
ylabel('Treatment Arm','FontSize',10);
title({'CFR Treatment Landscape',...
    'Green = Normal (≥2.5) | Red = CMD (<2.0)'},'FontSize',10,'FontWeight','bold');
% Add text annotations
for a = 1:n_arms
    for s = 1:n_states
        text(s,a,sprintf('%.2f',cfr_mat(a,s)),'HorizontalAlignment','center',...
            'FontSize',7,'Color','black','FontWeight','bold');
    end
end

subplot(1,2,2);
imagesc(imr_mat);
colormap(gca, red_green_map());
colorbar;
caxis([0 80]);
set(gca,'XTick',1:n_states,'XTickLabel',...
    {'FPH','FPT2','FPoH','FPoT2','MH','MT2'},...
    'YTick',1:n_arms,'YTickLabel',{arms.abbrev},...
    'FontSize',9,'TickDir','out');
xlabel('Sex-Disease State','FontSize',10);
title({'IMR Treatment Landscape (U)',...
    'Green = Normal (<25U) | Red = Severe CMD (>40U)'},'FontSize',10,'FontWeight','bold');
for a = 1:n_arms
    for s = 1:n_states
        text(s,a,sprintf('%.0f',imr_mat(a,s)),'HorizontalAlignment','center',...
            'FontSize',7,'Color','black','FontWeight','bold');
    end
end

sgtitle({'Treatment Landscape — All Arms x All States',...
    'QSP Model of Sex-Specific Diabetic Cardiomyopathy'},...
    'FontSize',11,'FontWeight','bold');
annotation('textbox',[0 0.96 1 0.04],'String',...
    sprintf('QSP-DCM | %s | v%s', run_timestamp, script_version),...
    'EdgeColor','none','HorizontalAlignment','center',...
    'FontSize',7,'Color',[0.5 0.5 0.5]);

save_figure(fig4, fullfile(outdir,'VCT_TreatmentLandscape'));
fprintf('  Saved: VCT_TreatmentLandscape\n');

%% --- Export Excel --------------------------------------------------------
fprintf('Exporting VCT results...\n');
xlsx_file = fullfile(outdir,'VCT_Results.xlsx');

rows = {};
for a = 1:n_arms
    for s = 1:n_states
        r = results(a,s);
        rows(end+1,:) = {arms(a).name, arms(a).abbrev, state_labels{s},...
            r.CFR, r.IMR, r.Ee, r.PGC1a, r.Ceramide, r.CMD,...
            r.CFR_resp, r.IMR_resp, r.dCFR, r.dIMR, r.dEe, r.dCMD};
    end
end
headers = {'Treatment','Abbrev','State',...
    'CFR','IMR','Ee','PGC1a','Ceramide','CMD',...
    'CFR_Responder','IMR_Responder','dCFR','dIMR','dEe','dCMD'};
T = cell2table(rows,'VariableNames',headers);
writetable(T, xlsx_file, 'Sheet','Full Results');

meta = {'Field','Value';
    'Script','Script6_VirtualClinicalTrial.m';
    'Version',script_version;
    'Timestamp',run_timestamp;
    'CFR responder threshold','≥2.5';
    'IMR responder threshold','<25U';
    'Treatment arms','Vehicle, SGLT2i, GLP-1, MHT, ARNI, SGLT2i+MHT, SGLT2i+GLP-1';
    'n states','6'};
writecell(meta, xlsx_file,'Sheet','Metadata');
fprintf('  Saved: %s\n', xlsx_file);

%% --- Print Summary -------------------------------------------------------
fprintf('\n=======================================================\n');
fprintf(' VIRTUAL CLINICAL TRIAL — KEY RESULTS\n');
fprintf('=======================================================\n');
fprintf('%-12s  %-22s  %6s  %6s  %5s  %5s\n',...
    'Arm','State','CFR','IMR','CFRok','IMRok');
fprintf('%s\n',repmat('-',1,65));
for a = 1:n_arms
    for s = [2,4,6]  % T2DM states only
        r = results(a,s);
        cfr_flag = 'NO ';
        imr_flag = 'NO ';
        if r.CFR_resp, cfr_flag='YES'; end
        if r.IMR_resp, imr_flag='YES'; end
        fprintf('%-12s  %-22s  %6.3f  %6.1f  %5s  %5s\n',...
            arms(a).abbrev, state_labels{s},...
            r.CFR, r.IMR, cfr_flag, imr_flag);
    end
    fprintf('\n');
end

fprintf(flog,'\nVCT complete: %s\n', datestr(now));
fclose(flog);
fprintf('Script 6 complete. All VCT files saved.\n\n');

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function d = build_drug_sglt2i()
    d.L1.k_ROS_prod=0.75;d.L1.k_Psi_form=1.20;d.L1.T2DM_factor=0.85;
    d.L2.T2DM_ins_res=0.80;d.L2.FA_overflow=0.85;d.L2.k_Cer_ROS=0.80;
    d.L3.k_NO_ROS_scav=0.80;d.L3_eNOS_boost=1.10;
end

function d = build_drug_glp1()
    d.L1.T2DM_factor=0.80;d.L1.ROS_drive=0.85;
    d.L2.T2DM_ins_res=0.65;d.L2.FA_overflow=0.70;d.L2.k_Cer_FA=0.75;
    d.L3.k_ET1_cer=0.85;d.L3_eNOS_boost=1.05;
end

function d = build_drug_mht()
    d.L1.k_ERb_sens=1.25;d.L1.E2_ss_mult=8.0;d.L1.mito_stress=0.85;
    d.L2.estrogen_prot=1.40;d.L2.k_Cer_FA=0.80;
    d.L3.eNOS_activity=1.50;d.L3.vasc_estrogen=1.60;d.L3.ET1_basal_mult=0.70;
    d.L3_eNOS_boost=1.50;
end

function d = build_drug_arni()
    d.L1.T2DM_factor=0.92;
    d.L2.k_Col_ROS=0.75;d.L2.k_Col_cer=0.75;d.L2.k_Ee_Col=0.80;
    d.L3.k_MRes_colx=0.70;d.L3.k_WS_Ee=0.75;d.L3_eNOS_boost=1.08;
end

function dc = combine_drugs(d1, d2)
    % Combine two drugs by multiplying their modifiers (multiplicative model)
    dc = d1;
    layers = {'L1','L2','L3'};
    for li = 1:3
        layer = layers{li};
        if isfield(d2,layer)
            fn2 = fieldnames(d2.(layer));
            for k = 1:length(fn2)
                f = fn2{k};
                if isfield(dc,layer) && isfield(dc.(layer),f)
                    % Special handling for E2_ss_mult (multiplicative is OK)
                    dc.(layer).(f) = dc.(layer).(f) * d2.(layer).(f);
                else
                    if ~isfield(dc,layer), dc.(layer) = struct(); end
                    dc.(layer).(f) = d2.(layer).(f);
                end
            end
        end
    end
    % Combine eNOS boost
    e1 = 1.0; e2 = 1.0;
    if isfield(d1,'L3_eNOS_boost'), e1 = d1.L3_eNOS_boost; end
    if isfield(d2,'L3_eNOS_boost'), e2 = d2.L3_eNOS_boost; end
    dc.L3_eNOS_boost = e1 * e2;
end

function res = run_model(sex_state, drug, tspan, ~)
    [p1,p2,p3] = get_all_params(sex_state, drug);
    opts = odeset('RelTol',1e-7,'AbsTol',1e-9,'MaxStep',0.5);

    y0_1=[p1.E2_0;p1.ERb_0;p1.PGC1a_0;p1.Mito_0;p1.ROS_0;p1.MnSOD_0;p1.DeltaPsi_0];
    [~,y1]=ode15s(@(t,y)L1_odes(t,y,p1),tspan,y0_1,opts);
    p2.PGC1a=y1(end,3);p2.Mito=y1(end,4);p2.ROS=y1(end,5);p2.dPsi=y1(end,7);p2.MnSOD=y1(end,6);

    y0_2=[p2.AMPK_0;p2.FAO_0;p2.GlucOx_0;p2.Cer_0;p2.ColX_0;p2.Ee_0];
    [~,y2]=ode15s(@(t,y)L2_odes(t,y,p2),tspan,y0_2,opts);
    p3.Ceramide=y2(end,4);p3.ColX=y2(end,5);p3.Ee=y2(end,6);
    p3.ROS=p2.ROS;p3.dPsi=p2.dPsi;p3.PGC1a=p2.PGC1a;p3.Mito=p2.Mito;

    y0_3=[p3.NO_0;p3.ET1_0;p3.MRes_0;p3.CFR_0;p3.IMR_0;p3.WS_0];
    [~,y3]=ode15s(@(t,y)L3_odes(t,y,p3),tspan,y0_3,opts);

    res.CFR=y3(end,4);res.IMR=y3(end,5);res.WS=y3(end,6);
    res.Ee=y2(end,6);res.PGC1a=y1(end,3);res.Ceramide=y2(end,4);
    cmd_cfr=max(0,(2.5-res.CFR)/2.5);cmd_imr=max(0,(res.IMR-10)/40);
    cmd_ee=max(0,(res.Ee-8)/20);
    res.CMD=cmd_cfr+cmd_imr+cmd_ee;
end

function [p1,p2,p3]=get_all_params(sex_state,drug)
    p1=get_L1_params(sex_state);
    if ~isempty(drug)&&isfield(drug,'L1')
        fn=fieldnames(drug.L1);
        for k=1:length(fn)
            f=fn{k};
            if strcmp(f,'E2_ss_mult')
                p1.E2_ss=p1.E2_ss*drug.L1.(f);p1.E2_0=p1.E2_ss;
            elseif isfield(p1,f),p1.(f)=p1.(f)*drug.L1.(f);end
        end
    end
    p2=get_L2_params(sex_state,p1);
    if ~isempty(drug)&&isfield(drug,'L2')
        fn=fieldnames(drug.L2);
        for k=1:length(fn)
            f=fn{k};if isfield(p2,f),p2.(f)=p2.(f)*drug.L2.(f);end
        end
    end
    p3=get_L3_params(sex_state,p2);
    if ~isempty(drug)&&isfield(drug,'L3')
        fn=fieldnames(drug.L3);
        for k=1:length(fn)
            f=fn{k};if isfield(p3,f),p3.(f)=p3.(f)*drug.L3.(f);end
        end
    end
    if ~isempty(drug)&&isfield(drug,'L3_eNOS_boost')
        p3.eNOS_activity=min(1.0,p3.eNOS_activity*drug.L3_eNOS_boost);
    end
end

function p=get_L1_params(sex_state)
    p.k_ERb_on=0.05;p.k_ERb_off=0.10;p.k_PGC_basal=0.08;p.k_PGC_ERb=0.25;
    p.k_PGC_deg=0.15;p.k_Mito_form=0.12;p.k_Mito_deg=0.05;p.k_ROS_prod=0.20;
    p.k_ROS_mito=0.15;p.k_ROS_scav=0.30;p.k_MnSOD_basal=2.00;p.k_MnSOD_ROS=8.00;
    p.k_MnSOD_PGC=5.00;p.k_MnSOD_deg=0.05;p.k_Psi_form=3.00;p.k_Psi_ROS=1.50;p.k_Psi_deg=0.08;
    switch sex_state
        case 'female_pre_healthy'
            p.E2_ss=150;p.k_ERb_sens=1.00;p.T2DM_factor=1.00;p.mito_stress=1.00;p.ROS_drive=1.00;
            p.E2_0=150;p.ERb_0=0.20;p.PGC1a_0=0.50;p.Mito_0=30;p.ROS_0=1.00;p.MnSOD_0=60;p.DeltaPsi_0=1.00;
        case 'female_pre_T2DM'
            p.E2_ss=150;p.k_ERb_sens=0.65;p.T2DM_factor=1.60;p.mito_stress=1.40;p.ROS_drive=1.80;
            p.E2_0=150;p.ERb_0=0.15;p.PGC1a_0=0.40;p.Mito_0=25;p.ROS_0=1.50;p.MnSOD_0=70;p.DeltaPsi_0=0.80;
        case 'female_post_healthy'
            p.E2_ss=15;p.k_ERb_sens=0.80;p.T2DM_factor=1.00;p.mito_stress=1.20;p.ROS_drive=1.20;
            p.E2_0=15;p.ERb_0=0.15;p.PGC1a_0=0.40;p.Mito_0=25;p.ROS_0=1.20;p.MnSOD_0=65;p.DeltaPsi_0=0.90;
        case 'female_post_T2DM'
            p.E2_ss=5;p.k_ERb_sens=0.50;p.T2DM_factor=1.80;p.mito_stress=1.70;p.ROS_drive=2.20;
            p.E2_0=5;p.ERb_0=0.10;p.PGC1a_0=0.35;p.Mito_0=20;p.ROS_0=1.80;p.MnSOD_0=75;p.DeltaPsi_0=0.70;
        case 'male_healthy'
            p.E2_ss=30;p.k_ERb_sens=0.60;p.T2DM_factor=1.00;p.mito_stress=1.00;p.ROS_drive=1.10;
            p.E2_0=30;p.ERb_0=0.12;p.PGC1a_0=0.45;p.Mito_0=28;p.ROS_0=1.10;p.MnSOD_0=55;p.DeltaPsi_0=0.95;
        case 'male_T2DM'
            p.E2_ss=30;p.k_ERb_sens=0.45;p.T2DM_factor=1.70;p.mito_stress=1.60;p.ROS_drive=2.00;
            p.E2_0=30;p.ERb_0=0.10;p.PGC1a_0=0.35;p.Mito_0=22;p.ROS_0=1.70;p.MnSOD_0=72;p.DeltaPsi_0=0.75;
        otherwise,error('Unknown state: %s',sex_state);
    end
end

function dydt=L1_odes(~,y,p)
    E2=max(y(1),0);ERb=max(y(2),0);PGC1a=max(y(3),0);Mito=max(y(4),0);
    ROS=max(y(5),0);MnSOD=max(y(6),0);DeltaPsi=max(y(7),0);
    dE2=p.k_ERb_on*(p.E2_ss-E2);
    dERb=p.k_ERb_on*p.k_ERb_sens*E2*(1-ERb)-p.k_ERb_off*ERb;
    dPGC1a=p.k_PGC_basal+p.k_PGC_ERb*ERb-p.k_PGC_deg*p.T2DM_factor*PGC1a;
    dMito=p.k_Mito_form*PGC1a*30-p.k_Mito_deg*p.mito_stress*Mito;
    ROS_prod=p.k_ROS_prod*p.ROS_drive+p.k_ROS_mito*p.T2DM_factor*(1/(PGC1a+0.1));
    dROS=ROS_prod-p.k_ROS_scav*MnSOD/100*ROS;
    dMnSOD=p.k_MnSOD_basal+p.k_MnSOD_ROS*ROS+p.k_MnSOD_PGC*PGC1a-p.k_MnSOD_deg*MnSOD;
    dDeltaPsi=p.k_Psi_form*Mito/30*PGC1a-p.k_Psi_ROS*ROS*DeltaPsi-p.k_Psi_deg*DeltaPsi;
    dydt=[dE2;dERb;dPGC1a;dMito;dROS;dMnSOD;dDeltaPsi];
end

function p=get_L2_params(sex_state,L1)
    p.PGC1a=L1.PGC1a_0;p.Mito=L1.Mito_0;p.ROS=L1.ROS_0;p.dPsi=L1.DeltaPsi_0;p.MnSOD=L1.MnSOD_0;
    p.k_AMPK_basal=0.10;p.k_AMPK_dPsi=0.08;p.k_AMPK_ROS=0.06;p.k_AMPK_deg=0.12;p.dPsi_ref=50.0;
    p.k_FAO_AMPK=0.20;p.k_FAO_PGC=0.15;p.k_FAO_deg=0.10;p.k_FAO_cer_inh=0.05;
    p.k_Gluc_basal=0.12;p.k_Gluc_ins_res=0.08;p.k_Gluc_deg=0.10;
    p.k_Cer_FA=0.08;p.k_Cer_ROS=0.10;p.k_Cer_deg=0.04;p.Cer_ref=6.0;
    p.k_Col_ROS=0.06;p.k_Col_cer=0.04;p.k_Col_deg=0.02;p.Col_ref=1.0;
    p.k_Ee_Col=0.10;p.k_Ee_Cer=0.06;p.k_Ee_dPsi=0.04;p.k_Ee_deg=0.03;
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
        otherwise,error('Unknown state: %s',sex_state);
    end
end

function dydt=L2_odes(~,y,p)
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

function p=get_L3_params(sex_state,L2)
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
        otherwise,error('Unknown state: %s',sex_state);
    end
end

function dydt=L3_odes(~,y,p)
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

function cmap=green_red_map()
    n=64; r=[zeros(1,n/2),linspace(0,1,n/2)];
    g=[linspace(0.2,0.9,n/2),linspace(0.9,0.1,n/2)];
    b=[zeros(1,n/2),zeros(1,n/2)];
    cmap=[r(:),g(:),b(:)];
end

function cmap=red_green_map()
    n=64; r=[linspace(1,0.1,n/2),zeros(1,n/2)];
    g=[linspace(0.1,0.9,n/2),linspace(0.9,0.8,n/2)];
    b=[zeros(1,n/2),zeros(1,n/2)];
    cmap=[r(:),g(:),b(:)];
end

function save_figure(fig,basepath)
    print(fig,basepath,'-dpng','-r300');
end
