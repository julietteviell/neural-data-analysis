%%
clear all
close all

%%

fileSCRAMAlgo={'m31'};%,'m40','m41','m43','m93','m96'}; %,
filename='processed_data_classification_script_50msbin';
clusterdata='clusterwaveSC';
files={fileSCRAMAlgo};
sessstring={'clus'};
bstring={'Left','right'};

%%
% Initialisation des structures et variables
fqrc = struct();
aucrc = struct();
hswrc = struct();
clusnumreal=zeros(1000,1);
%for o=1:1000 
for j = size(sessstring, 2)
    file = files{j};
    for k = 1:size(bstring, 2)
        % Réinitialisation des variables pour chaque itération
        fqall = [];
        aucall = [];
        hswall = [];
        waveall = [];
        biphasicall1=[];
        monophasicall1=[];
        hyper=[];
        % Boucle sur chaque fichier
        for p =1:numel(file)
            [auc, hsw, fq, wave, monophasicall, biphasicall, wenor1] = aucHswFqSCRamAlgo(file, p, k);
            fqall = [fqall, fq];
            aucall = [aucall, auc'];
            hswall = [hswall; hsw];
            waveall = [waveall; wave'];
            monophasicall1=[monophasicall1;monophasicall];
            biphasicall1=[biphasicall1;biphasicall];
            hyper=[hyper; wenor1];

        end

        % Stockage des résultats dans des structures
        fqrc.(bstring{k}) = fqall;
        aucrc.(bstring{k}) = aucall;
        hswrc.(bstring{k}) = hswall';
        biphasic.(bstring{k}) = biphasicall1;
        monophasic.(bstring{k}) = monophasicall1;
        waves.(bstring{k})=waveall;
        hyperall.(bstring{k})=hyper;
    end
        % Préparation des données pour le clustering
        datall = [
            [hswrc.(bstring{1})'; hswrc.(bstring{2})'], ...
            [aucrc.(bstring{1})';aucrc.(bstring{2})'] , ...
            [fqrc.(bstring{1})';  fqrc.(bstring{2})'], ...
            ];
          %  [biphasic.(bstring{1});biphasic.(bstring{2})],...
           % [monophasic.(bstring{1}); monophasic.(bstring{2})],...
       % ];
    
        datall(datall(:, 2) > 2000, :) = [];
        % datawave=[ waves.(bstring{1});  waves.(bstring{2})];
        % datahyper=[ hyperall.(bstring{1});  hyperall.(bstring{2})];
        % datamerge=[datahyper datall];



         % Détermination du nombre optimal de clusters par silhouette
        maxclust = 5; % Nombre maximal de clusters à tester
        silhouetteVals = zeros(1, maxclust-1);
        for nclust = 2:maxclust
            [idx1, ~] = kmeans(datall, nclust, 'dist', 'sqeuclidean');
            silhouetteVals(nclust-1) = mean(silhouette(datall, idx1));
        end
        [~, optimalClust1] = max(silhouetteVals);
        optimalClust1 = optimalClust1 + 1; % car on a commencé à 2 

   % K-means Clustering avec le nombre optimal de clusters
        [cidx1, cmeans] = kmeans(datall, optimalClust1, 'dist', 'sqeuclidean');
        

        
       
end
 % clusnumreal(o)=optimalClust1 ;
%end
% figure;
% histcounts(clusnumreal)

%%
% Visualisation 3D des clusters avec datall
cidx={cidx1};
optimalClust={optimalClust1};

for j=1
    figure;
    ptsymb = {'o', '^', '*','gd','.'};
    edgecolor = {[0.9660 0.640 0.1880],[0.9660 0.440 0.6880], [0.3290 0.9940 0.9250] [0.1290 0.3940 0.7250], [0.6290 0.2940 0.6250],[0.6290 0.2940 0.1250]};
    facecolor = {[0.9660 0.640 0.1880],[0.9660 0.440 0.2880], [0.3290 0.9940 0.9250] [0.1290 0.3940 0.7250], [0.6290 0.2940 0.6250],[0.6290 0.2940 0.1250],'k'};
    edgecolorcondi ={'w','k'};
    facecolorcondi ={'w','k'};
    hold on;

    for i = 1:optimalClust{j}
        % figure;
        plot3(datall(cidx{j} == i, 1), datall(cidx{j} == i, 2), datall(cidx{j} == i, 3), ...
            ptsymb{i}, 'MarkerEdgeColor', edgecolorcondi{2}, ...
            'Markersize', 8, 'MarkerFaceColor', facecolorcondi{1});
        legendStrings = arrayfun(@(x) sprintf('class%d (n=%d)', x, sum(cidx{j} == x)), 1:optimalClust{j}, 'UniformOutput', false);
        legend(legendStrings);

    end
     
    xlabel('Spike half-width (µs)', 'fontsize', 15, 'rotation', -17);
    ylabel('Area Under the Curve (mV²)', 'rotation', 18, 'fontsize', 15);
    zlabel('Firing rate (Hz)', 'fontsize', 15);
    xlim([-20 700]);
    ylim([-20 300]);
    zlim([0 15]);
    view(-115, 50);
end
%% 
load(filename,'LA','TLA','MLA','PNLAspe','MNLAspe','highNLA')
load(clusterdata)
condinew_class={LA,TLA,MLA,PNLAspe,MNLAspe,highNLA};
condistr={'PLA','TLA','MLA','PNLAspe','MNLAspe','highNLA'};
condi=condinew_class;
clearvars -except datall cidx condi condistr optimalClust
for j = 2
    figure;
    ptsymb = {'o', '^', '*', 'gd', '.'};
    edgecolor = {[0.9660 0.640 0.1880], [0.9660 0.440 0.6880], [0.3290 0.9940 0.9250], [0.1290 0.3940 0.7250], [0.6290 0.2940 0.6250], [0.6290 0.2940 0.1250]};
    facecolor = {[0.9660 0.640 0.1880], [0.9660 0.440 0.2880], [0.3290 0.9940 0.9250], [0.1290 0.3940 0.7250], [0.6290 0.2940 0.6250], [0.6290 0.2940 0.1250], 'k'};
    edgecolorcondi = {'w', 'k'};
    facecolorcondi = {'w', 'k'};
    hold on;

    % --- Boucle i : clusters ---
    h1 = gobjects(1, optimalClust{j});
    legendStrings = cell(1, optimalClust{j});
    for i = 1:optimalClust{j}
        h1(i) = plot3(datall(cidx{j} == i, 1), datall(cidx{j} == i, 2), datall(cidx{j} == i, 3), ...
            ptsymb{i}, 'MarkerEdgeColor', edgecolorcondi{2}, ...
            'Markersize', 8, 'MarkerFaceColor', facecolorcondi{1});
        legendStrings{i} = sprintf('class%d (n=%d)', i, sum(cidx{j} == i));
    end

    xlabel('Spike half-width (µs)', 'fontsize', 15, 'rotation', -17);
    ylabel('Area Under the Curve (mV²)', 'rotation', 18, 'fontsize', 15);
    zlabel('Firing rate (Hz)', 'fontsize', 15);
    xlim([-20 700]);
    ylim([-20 300]);
    zlim([0 15]);
    view(-115, 50);
 

    % Première légende (clusters)
    lgd1 = legend(h1, legendStrings, 'Location', 'northwest','AutoUpdate','off');

    % --- Boucle f : conditions ---
    h2 = gobjects(1, numel(condi));
    legendStrings1 = cell(1, numel(condi));
    for f = 1:numel(condi)
        h2(f) = plot3(datall(condi{f}, 1), datall(condi{f}, 2), datall(condi{f}, 3), ...
            ptsymb{5}, 'MarkerEdgeColor', edgecolor{f}, ...
            'Markersize', 15, 'MarkerFaceColor', facecolor{f});
        legendStrings1{f} = sprintf('%s (n=%d)', condistr{f}, numel(condi{f}));
    end


    xlabel('Spike half-width (µs)', 'fontsize', 15, 'rotation', -17);
    ylabel('Area Under the Curve (mV²)', 'rotation', 18, 'fontsize', 15);
    zlabel('Firing rate (Hz)', 'fontsize', 15);
    xlim([-20 700]);
    ylim([-20 300]);
    zlim([0 15]);
    view(-115, 50);
    grid on
    % Deuxième légende (conditions), sur un axe invisible superposé
    ax1 = gca;
    ax2 = axes('position', get(ax1, 'position'), 'visible', 'off');
    lgd2 = legend(ax2, h2, legendStrings1, 'Location', 'northeast');
    uistack(ax2, 'top');   % <-- force ax2 (et sa légende) au premier plan

    hold on;
end



