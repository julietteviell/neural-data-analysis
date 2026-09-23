close all
clear all

%% Percentage of activation per condition (comparaison deux à deux)
% Paramètres communs

files = {'processed_data_classification_script_50msbin'};
load(files{1})

   for j=1:2
    act_sig = 1.95;
    inh_sig = -1.95;
    if j==1 varstring={'heatstart','pressurestart'}; limi = [-2.5 7]; edges = linspace(-2.5, 7, 96); period= [125:220]; end
    if j==2 varstring={'heatW','PressureW'}; limi = [-2.5 2.5]; edges = linspace(-2.5, 2.5, 51); period= [125:175];end

    for v = 1:length(varstring)  

         %% useful for the segregation of the subpopulations
    conditions =  {
    struct('All',sensorycells,'type', 'sensory', 'axisZscore', [-1 15], 'axisNorm', [0 0.5]),...
    struct('All',unresponsive,'type', 'UR', 'axisZscore', [-1 15], 'axisNorm', [0 0.5]),...
    struct('All',PNLAspe,'type', 'PNLA', 'axisZscore', [-1 15], 'axisNorm', [0 0.5]),... 
    struct('All',MNLAspe,'type', 'MNLA', 'axisZscore', [-1 15], 'axisNorm', [0 0.5]),...
    struct('All',TLA,'type', 'TLA', 'axisZscore', [-1 15], 'axisNorm', [0 0.5]),...
    struct('All',MLA,'type', 'MLA', 'axisZscore', [-1 15], 'axisNorm', [0 0.5]),...
    struct('All',PLA,'type', 'PLA', 'axisZscore', [-1 15], 'axisNorm', [0 0.5]),...
    };

        % Initialisation des variables
        perc_actc_all_cond = cell(length(conditions), 1);
        perc_actr_all_cond = cell(length(conditions), 1);

        for c = 1:length(conditions)
            cond = conditions{c};

            %%% for pawW ipsi
            zscored_allc = [zscored_all_SC.LDH.L.(varstring{v}) zscored_all_SC.RDH.R.(varstring{v})];
            zscored_allc =zscored_allc(1:2:end, :); % 0.1ms bin
            zscored_allc =zscored_allc(period,:); %  cond.All to separate the different subpopulations

            perc_actc = zeros(size(zscored_allc, 1), 1);
            for k = 1:size(zscored_allc, 1)
                act = find(zscored_allc(k, :) >= act_sig);
                perc_actc(k) = numel(act) / size(zscored_allc, 2) * 100;
            end

               %%% for paW contra
            zscored_allr = [zscored_all_SC.RDH.L.(varstring{v}) zscored_all_SC.LDH.R.(varstring{v})]; %(:, cond.RDH)  (:, cond.LDH)
            zscored_allr =zscored_allr(1:2:end, :); % 0.1ms bin
            zscored_allr =zscored_allr(period, :); %  cond.All to separate the different subpopulations


            perc_actr = zeros(size(zscored_allr, 1), 1);
            for k = 1:size(zscored_allr, 1)
                act = find(zscored_allr(k, :) >= act_sig);
                perc_actr(k) = numel(act) / size(zscored_allr, 2) * 100;
            end

            % Stockage des résultats
            perc_actc_all{v,c} = perc_actc;
            perc_actr_all{v,c} = -1 .* perc_actr;
            zscore_actc_all{v,c} = zscored_allc;
            zscore_actr_all{v,c} = zscored_allr;
            [meanc, semc] = computeMeanSEM(zscored_allc');
            [meanr, semr] = computeMeanSEM(zscored_allc');
            meansemc_all{v,c}={meanc,semc};
            meansemr_all{v,c}={meanr,semr};

        end
    end

    % Création de la figure
    figure('Name', 'Comparison of activation percentages per condition', 'NumberTitle', 'off');%, 'Position', [100 100 1400 1000]);

    % Couleurs pour les conditions
    colors = [
        0.8 0.2 0.0;  % Rouge
        0.0 0.5 0.8;  % Bleu
        0.2 0.7 0.3;  % Vert
        0.8 0.3 0.7;  % Violet
        0.9 0.6 0.2;  % Orange
        0.1 0.8 0.5;  % Turquoise
        0.5 0.2 0.7;  % Mauve
        0.3 0.7 0.9;  % Bleu clair
        0.6 0.3 0.1   % Marron
        0.6 0.3 0.2   % Marron
        0.6 0.3 0.5   % Marron
        0.6 0.3 0.1   % Marron
        0.6 0.1 0.1   % Marron
        0.1 0.3 0.1   % Marron
        ];

    % Comparaison deux à deux des conditions
    num_conditions = length(conditions);
    num_varstring = length(varstring);
    subplot_idx = 1;

    for v = 1:num_varstring
        for i = 1 % :num_conditions  %% if subpopulation segregation is needed 
            cond1 = conditions{i};

            subplot_idx = (v-1) * num_conditions + i;
            subplot(num_varstring, num_conditions, subplot_idx);
            hold on;

            % Tracé pour la condition 1
            bar(edges(1:end), perc_actc_all{v,i}, 'LineWidth', 1.5, 'FaceAlpha', 0.3, 'EdgeAlpha', 0.2, 'FaceColor', colors(i, :), 'EdgeColor', colors(i, :));
            bar(edges(1:end), perc_actr_all{v,i}, 'LineWidth', 1.5, 'FaceAlpha', 0.3, 'EdgeAlpha', 0.2, 'FaceColor', colors(i, :) * 0.7, 'EdgeColor', colors(i, :) * 0.7);

            % Légende et labels
            title(sprintf('%s ', cond1.type), 'FontSize', 12); %,
            xlim(limi);
            ylim([-100 100]);
            if subplot_idx <= length(varstring) * (length(conditions) - 1)
                ylabel('% of modulated cells', 'FontSize', 12);
            end
            xlabel('time (sec)');
            box off;

            %Légende unique pour la première ligne
            if v == 1 && i == 1
                legend('Ipsi lateral Paw', 'Contralateral Paw' , 'Location', 'best', 'FontSize', 8);
            end

            subplot_idx = subplot_idx + 1;
        end
    end
    end

% Fonction pour calculer la moyenne et SEM
function [meanVal, semVal] = computeMeanSEM(y)
    if size(y, 1) == 1
        meanVal = y;
        semVal = y;
    else
        meanVal = mean(y, 1);
        N = size(y, 1);
        semVal = std(y) / sqrt(N);
    end
end