clear all;
close all;

%% ========================================================================
%  PARAMETRE PRINCIPAL : type de stimulus utilisé pour la synchronisation
%  'pressure'    -> applied force slope (dynamic Von Frey)
%  'temperature' -> Heat gradient (theraml hargreave ramping)
%  ========================================================================
modality = ['pressure'];   % <-- change here depending on mechanical or thermal trials

% mice data list 
files = {'m31'};%,'m40','m41','m43', 'm93', 'm96'};

switch modality
    case 'pressure'
        fileramp = {'m31ramptrials'};%,'m40ramptrials','m41ramptrials','m43ramptrials', 'm93ramptrials', 'm96ramptrials'};
    case 'temperature'
        fileramp = {'m31ramptrials'}%,'m40ramptrials','m41ramptrials','m43ramptrials', 'm93ramptrials', 'm96ramptrials'};
    otherwise
        error('modality is either ''pressure'' or ''temperature''');
end

%% Paramètres spécifiques à la modalité (seuils, axes, libellés, titres)
switch modality
    case 'pressure'
        params.corr_threshold  = 0.5;    % seuil |r| pour neurones "sélectionnés"
        params.interval_offset = 0;      % pas de décalage sur le début des intervalles
        params.plot_all_cumsum = true;   % graphe 6 : cumule TOUTES les colonnes
        params.axis2_max = 0.2;
        params.axis6     = [0 0.95 0 5];
        params.axis7_max = 0.7;
        params.varname   = 'force';
        params.titles = struct( ...
            'g1','normalized force cumulative  occurence', ...
            'g2','mean neural activity (± SEM) all neurons', ...
            'g3','Scatter plots with regression line', ...
            'g4','correlation histogram', ...
            'g5','mean activity (± SEM) of selected neurons', ...
            'g6','Cumulative activity (all neurons)', ...
            'g7','neural activity Heatmap');
    case 'temperature'
        params.corr_threshold  = 0.2;
        params.interval_offset = 2;      % décalage de 2s (évite l'artefact de début de rampe)
        params.plot_all_cumsum = false;  % graphe 6 : uniquement les neurones corrélés
        params.axis2_max = 0.8;
        params.axis6     = [0 0.95 0 10];
        params.axis7_max = 0.2;
        params.varname   = 'temperature';
        params.titles = struct( ...
            'g1','normalized thermal cumulative occurence', ...
            'g2','mean neural activity (± SEM) all neurons', ...
            'g3','Scatter plots with regression line', ...
            'g4','correlation histogram', ...
            'g5','mean activity (± SEM) of selected neurons', ...
            'g6','Cumulative activity (all neurons)', ...
            'g7','neural activity Heatmap');
end

%% ========================================================================
%  BOUCLE PRINCIPALE : j=1 -> Left paw / right DH , j=2 -> Right paw / right DH 
%  ========================================================================
for j = 1:2
    clearvars -except allneufqmeanall j files fileramp modality params
    dist_neufreqmean = cell(numel(files), 1);
    timedistall = [];
    distshock_normalized_selectedall = [];

    for p = 1:numel(files)

        %% 1. data loading and curve/interval synchronization (synccurvesall/int_select)
        %     force curve is exctrated from RAMalgo data recording 
        %     stacked in the "ramptrial" files, alors que la courbe de
        %     the temperature gradient (considered as a linear increase) is built here from "logline", a
        %     linear fonction line
        switch modality
            case 'pressure'
                load(fileramp{p}, 'synccurvesall', 'PWall','PWL','VFL','PWR','VFR', 'VFall','intstimVF');
                load(files{p}, 'x','xall','listegoodtrials','listegoodL','listegoodR');
                distshock = synccurvesall(:, 2);
                listegoodtrials = {listegoodL, listegoodR};
                Int_selec = [VFall, PWall]; 
                Int_selec_j     = Int_selec(listegoodtrials{j}, :);
                synccurvesall_j = synccurvesall;
                distshock_j     = distshock;

            case 'temperature'
                load(files{p}, 'x','xall','intThermalL','retthermalL','intThermalR','retthermalR');
                load(files{1}, 'logline');

                synccurvesallL = [];
                for i = 1:length(intThermalL(:,1))
                    t = logline(:,1) + intThermalL(i,1)';
                    mergetimecurve = [t logline(:,1)];
                    synccurvesallL = [synccurvesallL; mergetimecurve];
                end
                synccurvesallR = [];
                for i = 1:length(intThermalR(:,1))
                    t = logline(:,1) + intThermalR(i,1)';
                    mergetimecurve = [t logline(:,1)];
                    synccurvesallR = [synccurvesallR; mergetimecurve];
                end

                distshockL = synccurvesallL(:,2);
                distshockR = synccurvesallR(:,2);

                synccurvesall_cell = {synccurvesallL, synccurvesallR};
                distshock_cell     = {distshockL, distshockR}; 
                Int_selec_cell     = {intThermalL, intThermalR};

                synccurvesall_j = synccurvesall_cell{j};
                distshock_j     = distshock_cell{j};
                Int_selec_j     = Int_selec_cell{j};
        end

        %% 2. trial curve Normalization and trial concatenation 
        distshock_normalized_selected = [];
        distshock_raw_selected = [];
        binrangesall = [];
        synccurvesall_selected = [];
        for m = 1:size(Int_selec_j, 1)
            valid_indices = (synccurvesall_j(:,1) >= (Int_selec_j(m,1) + params.interval_offset) ...
                            & synccurvesall_j(:,1) <= Int_selec_j(m,2));
            distshock_selected = distshock_j(valid_indices);
            lower_end = min(distshock_selected);
            upper_end = max(distshock_selected);
            distshock_normalized = (distshock_selected - lower_end) / (upper_end - lower_end + eps);
            synccurves_selected = synccurvesall_j(valid_indices,2);
            binranges = synccurvesall_j(valid_indices,1);

            synccurvesall_selected  = [synccurvesall_selected; synccurves_selected];
            distshock_raw_selected  = [distshock_raw_selected; distshock_selected];
            distshock_normalized_selected = [distshock_normalized_selected; distshock_normalized];
            binrangesall = [binrangesall; binranges];
        end

        figure;
        subplot(2,1,1);
        plot(binrangesall, synccurvesall_selected);
        xlabel('time (sec)');
        ylabel(params.varname);
        subplot(2,1,2);
        plot(distshock_normalized_selected);
        xlabel('concatenated time (sec)');
        ylabel(['normalized ' params.varname]);

        %% 3. units Binning, normalization, interval selection
        if j == 1; x = xall.LDH; end  % to change for contralateral correlation
        if j == 2; x = xall.RDH; end

        maxselec = max(distshock_raw_selected);
        minselec = min(distshock_raw_selected);
        fieldNames = fieldnames(x);
        max_time = max(cellfun(@(field) max(x.(field)), fieldNames));
        Duration_Bin = 0.05;

        edges = 0:Duration_Bin:max_time;
        A = zeros(length(edges)-1, length(fieldNames));
        for f = 1:length(fieldNames)
            curr_neu = x.(fieldNames{f});
            if curr_neu <= 1
                binned_spike = zeros(size(edges)-1);
            else
                binned_spike = histcounts(curr_neu, edges);
            end
            A(:, f) = binned_spike;
        end

        relA = zeros(size(A));
        for k = 1:size(A, 2)
            curr_spikes = A(:, k);
            lower_end = min(curr_spikes);
            upper_end = max(curr_spikes);
            relA(:, k) = (curr_spikes - lower_end) / (upper_end - lower_end + eps);
        end
        relA = A; 

        Aall = [];
        for m = 1:size(Int_selec_j, 1)
            valid_indices = (edges(1:end-1) >= Int_selec_j(m,1) & edges(1:end-1) <= Int_selec_j(m,2));
            goodA = relA(valid_indices, :);
            Aall = [Aall; goodA];
        end
        relA = Aall;

        %% 4. pressure binning and occurence (timedistall), synchronized to neural activity (fqmeanall)
        binpressure = 0.05;
        edgespressure = 0:binpressure:1;

        indices_originaux = 1:length(distshock_normalized_selected);
        indices_reduits = linspace(1, length(distshock_normalized_selected), length(relA));
        distshock_normalized_selectedbinned = interp1(indices_originaux, distshock_normalized_selected, indices_reduits, 'linear');

        if length(relA) > length(distshock_normalized_selectedbinned)
            relA = relA(1:length(distshock_normalized_selectedbinned), :);
        end
        if length(relA) < length(distshock_normalized_selectedbinned)
            distshock_normalized_selectedbinned = distshock_normalized_selectedbinned(1:length(relA), :);
        end

        timedist = zeros(length(edgespressure)-1, 1);
        distfqmean = zeros(length(edgespressure)-1, size(relA, 2));
        for i = 1:(length(edgespressure)-1)
            idx = (distshock_normalized_selectedbinned >= edgespressure(i)) & (distshock_normalized_selectedbinned < edgespressure(i+1));
            dist = relA(idx, :);
            timedist(i) = size(dist, 1);
            distfqmean(i,:) = mean(dist, 1); %mean sum
        end
        timedistall = [timedistall timedist];

        distfqmean(isnan(distfqmean)) = 0;
        dist_neufreqmean{p} = distfqmean;
        distshock_normalized_selectedall{1,p} = distshock_normalized_selected;
    end

    %% 5. Concaténation de tous les fichiers pour ce côté (j)
    allneufqmean2 = cellfun(@(x) x, dist_neufreqmean, 'UniformOutput', false);
    allneufqmean2 = cat(2, allneufqmean2{:});
    allneufqmeanall{1,j} = allneufqmean2;
end

% Left paw (LDH) + right paw (RDH) concatenation
allneufqmean = [allneufqmeanall{1,1} allneufqmeanall{1,2}];

%% ========================================================================
            %  Analysis of all data and graphs
%  ========================================================================
x1 = edgespressure(1:end-1);   % harmonisé (longueur = size(allneufqmean,1) dans les 2 cas)
mean_rel_dist = mean(allneufqmean, 2);
sem_rel_dist = std(allneufqmean, 0, 2) / sqrt(size(allneufqmean, 2));

figure;

% Graph 1: cumulative stimulus occurence
subplot(3, 3, 1);
bar(x1, timedistall, 'b');
title(params.titles.g1);
axis([0 max(edgespressure(1:end-1)) 0 max(timedistall(1))]);

% Graph 2: mean neural activity (± SEM), all neurons
subplot(3, 3, 2);
shadedErrorBar(x1, smooth(mean_rel_dist,0.5,'loess'), smooth(sem_rel_dist,0.5,'loess'), ...
    'lineprops', {'Color', [0.2 0.7 0.8]}, 'transparent', 1, 'patchSaturation', 0.53);
xlabel(['normalized ' params.varname]);
ylabel('normalized neuronal activity');
title(params.titles.g2);
axis([0 max(edgespressure(1:end-1)) 0 params.axis2_max]);
legend('mean ± SEM', 'Location', 'best');

% Graph 3: Scatter plots with regression line for r>0.5
subplot(3, 3, 3);
hold on;
rall = [];
pall = [];
for i = 1:size(allneufqmean, 2)
    [r, p] = corr(x1', allneufqmean(:, i));
    rall = [rall, r];
    pall = [pall, p];
    if r > params.corr_threshold
        scatter(x1, allneufqmean(:, i), 'DisplayName', ['Colonne ', num2str(i)]);
        coeffs = polyfit(x1, allneufqmean(:, i), 1);
        xFit = linspace(min(x1+0.2), max(x1-0.2), 100);
        yFit = polyval(coeffs, xFit);
        plot(xFit, yFit, 'k-', 'LineWidth', 0.5, 'DisplayName', ['Pente: ', num2str(coeffs(1), 3)]);
        line(x1, allneufqmean(:, i));
    end
end
hold off;
xlabel(['normalized ' params.varname]);
ylabel('normalized neural activity');
title(params.titles.g3);
grid off;
rall(isnan(rall)) = 0;
pall(isnan(pall)) = 1;

% Graph 4: coefficient correlation histogram r
subplot(3, 3, 4);
histogram(rall, 16);
xlim([-1 1]);
ylim([0 20]);
xlabel('r');
ylabel('neuron count');
title(params.titles.g4);

% Graph 5: neural activity Heatmap
subplot(3, 3, 7:9);
heatmap(allneufqmean);
xlim([1 size(allneufqmean,2)]);
ylim([1 length(edgespressure)-1]);
caxis([0 0.2]);
xlabel('neuron #');
ylabel(['normalized ' params.varname]);
title(params.titles.g7);

% Graph 6: Cumulative activity of each neuron 
subplot(3, 3, 6);
corrneuron     = intersect(find(rall >  params.corr_threshold), find(pall < 0.05));
anticorrneuron = intersect(find(rall < -params.corr_threshold), find(pall < 0.05));
if params.plot_all_cumsum
    M_cumsum = cumsum(allneufqmean(:, :));
else
    M_cumsum = cumsum(allneufqmean(:, corrneuron));
end
hold on;
for i = 1:size(M_cumsum, 2)
    plot(x1, M_cumsum(:, i), 'DisplayName', ['Colonne ', num2str(i)]);
end
hold off;
axis(params.axis6);
xlabel(['normalized ' params.varname]);
ylabel('sum activity');
title(params.titles.g6);

% Graph 7: mean activity (± SEM) of selected neurons 
subplot(3, 3, 5);
mean_selec  = mean(allneufqmean(:, corrneuron), 2);
sem_selec   = std(allneufqmean(:, corrneuron), 0, 2) / sqrt(size(allneufqmean(:, corrneuron), 2));
mean_selec2 = mean(allneufqmean(:, anticorrneuron), 2);
sem_selec2  = std(allneufqmean(:, anticorrneuron), 0, 2) / sqrt(size(allneufqmean(:, corrneuron), 2));
shadedErrorBar(x1, smooth(mean_selec,0.5,'loess'), smooth(sem_selec,0.5,'loess'), ...
    'lineprops', {'Color', [0.2 0.7 0.8]}, 'transparent', 1, 'patchSaturation', 0.53);
hold on
shadedErrorBar(x1, smooth(mean_selec2,0.5,'loess'), smooth(sem_selec2,0.5,'loess'), ...
    'lineprops', {'Color', [0.2 0.2 0.8]}, 'transparent', 1, 'patchSaturation', 0.53);
xlabel(['normalized ' params.varname]);
ylabel('normalized neuronal activity');
title(params.titles.g5);
axis([0 1 0 params.axis7_max]);
legend('Mean ± SEM', 'Location', 'best');
hold off