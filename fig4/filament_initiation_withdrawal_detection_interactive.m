% --- CODE PRINCIPAL ---

%%% in this code, the Start time and withdrawal time is automatically
%%% detected from the rise and drop of the force curve, but can be reajusted
%%% manually if needed (e.g. artefact issues)

clear all;
close all;

fileramp = 'm31ramptrials.mat';
load(fileramp);

listetrialsgood = find(contains(orderfile, 'Mecha'));
dataorderedall = data(listeramptrialordertorec); % for 43 

PW = [];
VF = [];
quality = cell(length(listetrialsgood), 1);
shuttle_run_all = cell(length(listetrialsgood), 1);

% Déclare les variables globales pour les callbacks
global hLineGlobalStart;
global hLineGlobalThresh;

for g = 1:length(dataorderedall)  %(listetrialsgood)
    clear var varlist acceleration acc_window acc_speed inst_acc inst_speed time_limits

    Gtrial = dataorderedall{g, 1};
    change_time_bin = 1;
    varlist = [Gtrial.trial(:, 1)];
    time_limits = (1:change_time_bin:size(varlist, 1));

    for i = 1:size(time_limits, 2)-1
        var(i, :) = mean(varlist(time_limits(i)+1:time_limits(i+1), 1), 1);
    end

    for i = 2:size(var, 1)
        acceleration(i) = (var(i, :) - var(i-1, :)) * 1000;
    end

    acceleration = acceleration';
    acc_window = acceleration;
    acc_window_speed = var;

    % Traitement pour trouver start et thresh
    buff = 50;
    start_found = false;
    thresh_found = false;

    for j = (1+buff):size(acc_window, 1)
        inst_acc = acc_window;
        inst_speed = acc_window_speed;
        maxpress = max(inst_speed);
        meanbuff = mean(acc_window_speed(j-buff:j, :));
        inst_acc_thresh = min(inst_acc);

        if ~start_found && inst_speed(j) >= meanbuff*1.3
            start = time_limits(j-(buff-40));
            start_found = true;
        end

        if inst_speed(j-1) >= inst_speed(j)*1.4
            thresh = time_limits(j);
            thresh_found = true;
            break;
        end
    end

    % Affichage de la figure avec curseurs interactifs
    fig = figure;
    plot(time_limits(1:end-1), acc_window);
    hold on;
    plot(time_limits(1:end-1), acc_window_speed * 1000);
    plot([thresh, thresh], [min(acc_window), max(acc_window)], 'r--', 'LineWidth', 1.5);

    if ~exist('start', 'var')
        start = 0;
    end

    plot([start, start], [min(acc_window), max(acc_window)], 'b--', 'LineWidth', 1.5);
    xlabel('Temps (ms)');
    ylabel('Amplitude');
    title('Sélectionnez un temps avec les curseurs');

    % Ajout du curseur vertical pour "start"
    hSliderStart = uicontrol('Style', 'slider', ...
        'Min', min(time_limits), 'Max', max(time_limits), ...
        'Value', start, ...
        'Position', [20 50 200 20]);

    % Ajout du curseur vertical pour "thresh"
    hSliderThresh = uicontrol('Style', 'slider', ...
        'Min', min(time_limits), 'Max', max(time_limits), ...
        'Value', thresh, ...
        'Position', [20 20 200 20]);

    % Ligne verticale liée au curseur "start"
    hLineGlobalStart = plot([start, start], [min(acc_window), max(acc_window)], 'b-', 'LineWidth', 2);

    % Ligne verticale liée au curseur "thresh"
    hLineGlobalThresh = plot([thresh, thresh], [min(acc_window), max(acc_window)], 'g-', 'LineWidth', 2);

    % Assigne les callbacks
    set(hSliderStart, 'Callback', @updateCursorStart);
    set(hSliderThresh, 'Callback', @updateCursorThresh);

    % Ajoute un bouton de validation
    uicontrol('Style', 'pushbutton', 'String', 'Valider', ...
              'Position', [250 20 100 20], ...
              'Callback', 'uiresume(gcbf)');

    % Attend que l'utilisateur clique sur "Valider"
    uiwait(fig);

    % Récupère les valeurs des sliders après validation
    selected_time_start = hSliderStart.Value;
    selected_time_thresh = hSliderThresh.Value;

    % Boîte de dialogue pour confirmer la validation
    choice = questdlg('Valider ces temps ?', 'Validation', 'Oui', 'Non', 'Oui');

    if strcmp(choice, 'Oui')
        quality{g, 1} = 1;
    else
        quality{g, 1} = 0;
        selected_time_start = start; % Par défaut
        selected_time_thresh = thresh; % Par défaut
    end

    % Sauvegarde des temps sélectionnés
    start1 = intstimVF(1:end, 1);
    VF1 = selected_time_start/1000 + start1(g);
    PW1 = selected_time_thresh/1000 + start1(g);
    VF = [VF; VF1];
    PW = [PW; PW1];
    avoid_run = acc_window_speed; % À adapter
    shuttle_run_all{g, 1} = {avoid_run, time_limits, VF, PW, selected_time_start, selected_time_thresh};
    hold off;
end
%%

synccurvesall=[];

for i=1:length(dataorderedall)
    clear mergetimecurve t synccurve PW1 VF1

    if ~isempty(shuttle_run_all{i,1})
    t=(shuttle_run_all{i,1}{1,2}(1:end-1)/1000)+intstimVF(i,1)';
    synccurvet{i,1}=t;
    mergetimecurve=[t' shuttle_run_all{i,1}{1,1}];
    synccurve{i,1}=mergetimecurve;
    synccurvesall=[synccurvesall;mergetimecurve];
    
    else
        continue
    end
end
%
%%
VFall=shuttle_run_all{length(dataorderedall),1}{1,3};
PWall=shuttle_run_all{i,1}{1,4};
%save
%(fileramp,'shuttle_run_all','synccurvet','synccurve','synccurvesall','PWall','VFall','-append')


% --- FONCTIONS LOCALES ---
function updateCursorStart(~, ~)
    global hLineGlobalStart;
    % Récupère la valeur du slider "start"
    hSliderStart = gcbo;
    newTimeStart = hSliderStart.Value;
    % Met à jour la ligne bleue
    set(hLineGlobalStart, 'XData', [newTimeStart, newTimeStart]);
    % Met à jour le titre
    title(sprintf('Temps "start" sélectionné: %.2f ms', newTimeStart));
end

function updateCursorThresh(~, ~)
    global hLineGlobalThresh;
    % Récupère la valeur du slider "thresh"
    hSliderThresh = gcbo;
    newTimeThresh = hSliderThresh.Value;
    % Met à jour la ligne verte
    set(hLineGlobalThresh, 'XData', [newTimeThresh, newTimeThresh]);
    % Met à jour le titre
    title(sprintf('Temps "thresh" sélectionné: %.2f ms', newTimeThresh));
end

