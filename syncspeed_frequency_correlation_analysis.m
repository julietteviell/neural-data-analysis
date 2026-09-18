
clear all
close all
filename={'M93_spontaneous1preSNI'}; %,'M96_spontaneous2preSNI','M31_spontaneousconcatpreSNI'
    rall=[];
    pall=[];

for i= 1:length(filename)
      load(filename{i},'recstartTTL','recstopTTL','TTLstart','TTLstop','valeurs_fq','valeur_speed','rectime','exculsionint','valeurs_speed')
      if ~exist('valeur_speed')
          valeur_speed=valeurs_speed;
      end
       valeur_speed_sec=length(valeur_speed)/19.44;
    if i==3
        exclusionint=[1995.49118076616	2069.44749186667;...
            2090.71553539072	2144.69247633333;...
            2172.74542314509	2773.77672666667;...
            2801.15640417465	2828.398035;...
            2851.84549012405	3189.35718733333];
        startint= 1.972717835000000e+03;
        intrec=[[startint; exclusionint(1:end-2,2)] exclusionint(2:end,1) ];
        intreccum=[0;cumsum(intrec(:,2)-intrec(:,1))];

        j = 1; % indice pour rectime_fq et recvaleur_fq

        rectime_fq = cell(1, floor(size(valeurs_fq, 2)/2));
        recvaleur_fq = cell(1, floor(size(valeurs_fq, 2)/2));


        for i = 1:2:size(valeurs_fq,2)
            temp_rectime = [];
            temp_recvaleur = [];
            for k=1:length(intrec)
                clear valeurs_fqTTL
                valeurs_fqTTL = valeurs_fq(:,i) > intrec(k,1)&valeurs_fq(:,i) < intrec(k,2);
                temp_rectime = [temp_rectime; valeurs_fq(valeurs_fqTTL, i) - intrec(k, 1) + intreccum(k)];
                temp_recvaleur = [temp_recvaleur; valeurs_fq(valeurs_fqTTL, i + 1)];
            end
            rectime_fq{j} = temp_rectime;
            recvaleur_fq{j} = temp_recvaleur;
            j = j + 1;
        end


    else 

   
    close all
    % Paramètres

    % cutting neural data to video rec
    n_cols = floor(size(valeurs_fq,2)/2);
    rectime_fq = cell(1, n_cols); % On utilise des cellules pour stocker des vecteurs de tailles variables
    recvaleur_fq = cell(1, n_cols);

    j = 1; % indice pour rectime_fq et recvaleur_fq
    for i = 1:2:size(valeurs_fq,2)
        valeurs_fqTTL = valeurs_fq(:,i) > recstartTTL&valeurs_fq(:,i) < recstopTTL;
        rectime_fq{j} = valeurs_fq(valeurs_fqTTL, i)-recstartTTL;
        recvaleur_fq{j} = valeurs_fq(valeurs_fqTTL, i+1);
        j = j + 1;
    end
    end

    % cutting video data to TTL time

     if i==3
       recvaleur_speed= valeur_speed;
     else
    %videotime=162.672;
    % TTLstart=17.644;
    % TTLstop=162.672;
    timespeed=length(valeur_speed)*0.103; %0.103 = recalibration index
    t_speed=0.103:0.103:timespeed;
    valeurs_speedTTL = t_speed> TTLstart&t_speed < TTLstop;
    rectime_speed = t_speed(valeurs_speedTTL)-TTLstart;
    recvaleur_speed = valeur_speed(valeurs_speedTTL);
     end

    % Créer une grille temporelle commune (0.5s d'intervalle)
    t_min = min(rectime_speed);
    t_max = max(rectime_speed);
    t_common = 0:1:t_max;
    edges = [t_common, t_common(end)+0.5]; % Bins pour histcounts

    recvaleur_fq_binned = cell(1, length(recvaleur_fq));
    for j = 1:length(recvaleur_fq)
        % Trouver l'indice du bin pour chaque point
        [~, ~, bin_idx] = histcounts(rectime_fq{j}, edges);

        % Initialiser le vecteur binné
        binned_values = zeros(1, length(t_common));
        bin_counts = zeros(1, length(t_common));

        % Pour chaque point, ajouter sa valeur au bin correspondant
        for i = 1:length(bin_idx)
            if bin_idx(i) > 0 && bin_idx(i) <= length(t_common)
                binned_values(bin_idx(i)) = binned_values(bin_idx(i)) + recvaleur_fq{j}(i);
                bin_counts(bin_idx(i)) = bin_counts(bin_idx(i)) + 1;
            end
        end

        % Calculer la moyenne par bin (ou autre statistique)
        binned_values(bin_counts > 0) = binned_values(bin_counts > 0) ./ bin_counts(bin_counts > 0);
        % Si un bin est vide, on peut le laisser à 0 ou interpoler plus tard

        recvaleur_fq_binned{j} = binned_values;
    end

    % Initialiser les cellules pour les valeurs interpolées
    recvaleur_fq_interp = zeros(length(t_common),size(recvaleur_fq,2));

    for j = 1:length(recvaleur_fq)-1
        % Interpolation linéaire pour chaque canal
        %neu=interp1(rectime_fq{j}, recvaleur_fq{j}, t_common, 'next')';
        neu=recvaleur_fq_binned{j};
        neu(isnan(neu))=0;
        recvaleur_fq_interp(:,j) = neu';
    end

    [~, ~, bin_idx_speed] = histcounts(rectime_speed, edges);
    binned_speed = zeros(1, length(t_common));
    bin_counts_speed = zeros(1, length(t_common));
    for i = 1:length(bin_idx_speed)
        if bin_idx_speed(i) > 0 && bin_idx_speed(i) <= length(t_common)
            binned_speed(bin_idx_speed(i)) = binned_speed(bin_idx_speed(i)) + recvaleur_speed(i);
            bin_counts_speed(bin_idx_speed(i)) = bin_counts_speed(bin_idx_speed(i)) + 1;
        end
    end
    binned_speed(bin_counts_speed > 0) = binned_speed(bin_counts_speed > 0) ./ bin_counts_speed(bin_counts_speed > 0);

    valeurs_speed_continue=interp1(rectime_speed,recvaleur_speed,t_common,'next');
    valeurs_speed_continue(isnan(valeurs_speed_continue))=0;
    %'linear', 'nearest', 'next', 'previous', 'pchip', 'cubic', 'v5cubic', 'makima', or 'spline'


    %%
    % Tracé des deux courbesfigure;
    figure;
    plot(t_common,binned_speed,'b-', 'LineWidth', 2, 'DisplayName', 'velocity');
    hold on
    plot(t_common,recvaleur_fq_interp', 'LineWidth', 1, 'DisplayName', 'instant frequency');
    axis([0 160 0 200])
    xlabel('Temps (s)');
    ylabel('Fq (Hz)');
    title('correlation mobility neurone freq');
    %legend('show');
    grid off;

    clear r p  
    [r p]=corr(binned_speed',recvaleur_fq_interp);
    rall=[rall; r'];
    pall=[pall; p'];
    recvaleurall{i}=recvaleur_fq_interp;
    recspeedall{i}=binned_speed;
    end
   
FqImmobile=recvaleur_fq_interp(find(binned_speed<0.1),:);
meanFqImmobile=mean(FqImmobile,2);
FqMoving=recvaleur_fq_interp(find(binned_speed>0.1),:);
meanFqMoving=mean(FqMoving,2);
%
figure;
histogram(rall, 16);
xlim([-0.5 0.5]);
ylim([0 15]);
title('correlation histogram');
xlabel('r');
ylabel('neuron #');

correlated=rall>0.15&pall<0.05;

corrneuall=rall(correlated);
x1=t_common;

figure; 
mean_selec = mean(recvaleur_fq_interp, 2);
sem_selec = std(recvaleur_fq_interp, 0, 2)/sqrt(size(recvaleur_fq_interp, 2));
shadedErrorBar(x1, smooth(mean_selec,0.03,'loess'), smooth(sem_selec,0.03,'loess'), ...
    'lineprops', {'Color', [0.2 0.7 0.8]}, 'transparent', 1, 'patchSaturation', 0.53);
title('mean Frequency all neurons')
xlabel('Temps (s)');
ylabel('Fq (Hz)');

figure; 
plot(t_common,smooth(binned_speed,0.03,'loess'))
title('speed')
xlabel('Temps (s)');
ylabel('cm/s');





