clear all
close all

fileSCRAMAlgo={'m31'};%'m40','m41','m43','m93','m96'
fileAnesAwake={'m31'}; %'m93','m96'
files={fileSCRAMAlgo}; % fileAnesAwake % for Fig 3 data  analysis
stringfile={'preSC'};
%%%make intervals and cell types %%%
for g=1:numel(files)
    clearvars -except files stringfile g  miceRVMSC  normiceRVMSC bstring strDH varstring
    file=files{1,g};
    strDH={'LDH','RDH'};
    mice2DH = struct();
    normice2DH = struct();
                          
    for m = 1:length(strDH)
        mice = cell(numel(file), 1);
        normice = cell(numel(file), 1);
        
  for p =1:numel(file)
                        clearvars -except files stringfile g p file miceRVMSC normiceRVMSC mice normice mice2DH normice2DH strDH m
    
            
            load(file{p},'stimanes','retthermalL','retthermalR','PWL','PWR','VFL','VFR','intThermalL','intThermalR','mergecurves_ind','mergecurves_ts,','retmechanicalL','retmechanicalR','pinchL1','pinchR1')
            if ~exist('PWL'), PWL=0; end ;if ~exist('PWR'), PWR=0; end; if ~exist('VFL'), VFL=0; end; if ~exist('VFR'), VFR=0; end

                var3dL={intThermalL(:,1),retthermalL,VFL,PWL}; %,stimanes
                var3dR={intThermalR(:,1),retthermalR,VFR,PWR};
                var3d.var3dL=var3dL;
                var3d.var3dR=var3dR;
                varstring={'heatstart','heatW','pressurestart','PressureW'};%,'stimanes'              
                bstring={'L','R'};
                bin=0.05;         %% time bin to change
                edges=(-15:bin:15);  %% egdes before and after 0
                limi=[-2.5 7];
                
                load(file{p},'x')
                [x]=crosco(file,p); %removes doubles based on cofiring probability
                %%
                names=fieldnames(x);
                L=struct();
                R=struct();
                for i = 1:length(fieldnames(x))
                    clear chaine deux chiffres
                    chaine=x.(names{i});
                    listeL={'01','02', '03', '04', '05', '06', '07', '08'};
                    listeR={'09', '10', '11' ,'12', '13' ,'14', '15' ,'16'};
                    deuxchiffres = regexp(names{i}, '\d{2}', 'match');
                    l=ismember(listeL,deuxchiffres)==1;
                    r=ismember(listeR,deuxchiffres)==1;
                    if any(l)      %(l==1))
                        L.(names{i})=chaine;
                    end
                    if  any(r)     %~isempty(find(r==1))
                        R.(names{i})=chaine;
                    end
                end
                if  p==2 || p==6 %for Ramalgo file  %p==3 for AnesAwake file
                    xall.LDH=R;
                    xall.RDH=L;
                else
                    
                    xall.LDH=L;
                    xall.RDH=R;
                end
                
                b=struct2cell(xall.(strDH{m}));%
                if isempty(b)
                    continue
                else
                    [total_counts_mean_var total_counts_mean_var_norm]= freqVarRAMalg(var3d,varstring,b,bstring,p,file,edges, bin);%creates spike counts per bins around each event 
                    mice{p,1}=total_counts_mean_var;
                    normice{p,1}=total_counts_mean_var_norm;

                end
                mice2DH.(strDH{m})=mice;
                normice2DH.(strDH{m})=normice;
        end
        miceSC.(stringfile{g})=mice2DH;
        normiceSC.(stringfile{g})=normice2DH;
    end
end

%% concatenates all mice neuronal binned data for each dorsal horn, on each paw, for each event, in normalized and raw mean frequency
clearvars -except edges bstring varstring files stringfile g limi normiceSC miceSC micevar var3d strDH
      
for m = 1:length(strDH)
    for k = 1:length(bstring)
        for i = 1:length(varstring)  
            normiceConcat = [];
            miceConcat = [];
            for idx = 1:length(normiceSC.preSC.(strDH{m}))
                if isempty(normiceSC.preSC.(strDH{m}){idx,1})
                    continue
                else
                normiceConcat = [normiceConcat, normiceSC.preSC.(strDH{m}){idx,1}.(bstring{k}).(varstring{i})];          
                miceConcat = [miceConcat, miceSC.preSC.(strDH{m}){idx,1}.(bstring{k}).(varstring{i})];  
                end
            end
            normicepreSC.(strDH{m}).(bstring{k}).(varstring{i}) = normiceConcat;
            micepreSC.(strDH{m}).(bstring{k}).(varstring{i}) = miceConcat;
        end
    end
end           
normiceall.SC=normicepreSC;
miceall.SC=micepreSC;

%% %%%SC neurons classification%%%q

for g=1:numel(files)
    for   m=1:numel(strDH)
        clearvars -except mice normice edges bstring varstring files stringfile g limi normiceSC miceSC micevar var3d indexall orderall...
            responseheatmapall zscoredall actcells_all inhcells_all miceall normiceall stringfile2 strDH m zscored_all_SC

        %%SCneurons classification
        mice = miceSC.(stringfile{g}).(strDH{m});
        estNonVide = cellfun(@(x) ~isempty(x) && ~(ischar(x) && all(isspace(x))), mice); mice=mice(estNonVide);
        normice = normiceSC.(stringfile{g}).(strDH{m});
        [zscored_all_varcr, basal] = zscoreVarSCRAMalgo(varstring,bstring,mice,g);
        [order]=heatMapZscoreSCRAMalgo(zscored_all_varcr,normiceall,bstring,varstring,stringfile,edges,limi,g); % stimuli and paw withdrawal response heatmaps ordered

        [actcells_allcr inhcells_allcr]= cellClassBootStrapRAMalgo(varstring,bstring,normice,mice,files,edges); % activation inhibition cell classification to each event 
        actcells_all.(stringfile{g}).(strDH{m})=actcells_allcr;
        inhcells_all.(stringfile{g}).(strDH{m})=inhcells_allcr;
        orderall.(['order' stringfile{g}])=order;
        zscored_all_SC.(strDH{m})=zscored_all_varcr;
        zscoredall.(['zscored' stringfile{g}])=zscored_all_SC;
        [responseheatmap percentagecombined indexLR]=heatmapResponsivenessMatriceOfEvents2(zscored_all_varcr, varstring, bstring, stringfile, actcells_allcr, inhcells_allcr, g)  %response panel for each event
        indexall.(['index' stringfile{g}])=indexLR;
        responseheatmapall.(['responseheatmap' strDH{m}])=responseheatmap;
    end
end


%% %%%% Analysis of cell response to paw and heat, for the Ipsi and Contralateral Paw
load ('processed_data_classification_script_50msbin')

names = varstring;
bin = 0.05;
edges = (-15:bin:15);
x1 = edges(1:end-1);
colors = {
    [0.8 0.3 0.8], [0.2 0.3 0.7], ... % Ipsi/Contra for zscore data
    [0.2 0.5 0.8], [0.2 0.8 0.8]       % Ipsi/Contra for normalized data
};
lineProps = {'transparent', 1, 'patchSaturation', 0.33};
smoothParam = 0.07;
smoothMethod = 'loess';   %% smoothing parameters

conditions =  {
    struct('All',sensorycells,'type', 'sensory', 'axisZscore', [-1 15], 'axisNorm', [0 0.5]),...
    struct('All',unresponsive,'type', 'UR', 'axisZscore', [-1 15], 'axisNorm', [0 0.5]),...
    struct('All',PNLAspe,'type', 'PNLA', 'axisZscore', [-1 15], 'axisNorm', [0 0.5]),... 
    struct('All',MNLAspe,'type', 'MNLA', 'axisZscore', [-1 15], 'axisNorm', [0 0.5]),...
    struct('All',TLA,'type', 'TLA', 'axisZscore', [-1 15], 'axisNorm', [0 0.5]),...
    struct('All',MLA,'type', 'MLA', 'axisZscore', [-1 15], 'axisNorm', [0 0.5]),...
    struct('All',PLA,'type', 'PLA', 'axisZscore', [-1 15], 'axisNorm', [0 0.5]),...
    };

% selected events (thermal withdrawal ; mechanical withdrawal)
selectedNames = [2 4];
nRows = length(selectedNames);  
nCols = length(conditions);     


figZscore = figure('Name', ['Ipsi contra DH/Paw responses - Zscore ' stringfile{1}], ...
    'NumberTitle', 'off', 'Position', [50 50 250*nCols 250*nRows]);
figNorm = figure('Name', ['Ipsi contra DH/Paw responses - Normalized ' stringfile{1}], ...
    'NumberTitle', 'off', 'Position', [50 50 250*nCols 250*nRows]);

% one iteration per category
for c = 1:length(conditions)
    cond = conditions{c};

    % one iteration per selected event
    for r = 1:nRows
        i = selectedNames(r);

        %%% bilateral PAW
        % Zscored data
        y1 = [zscored_all_SC.LDH.L.(names{i})'; zscored_all_SC.RDH.R.(names{i})']; %% ipsilateral paw stimulation
        y1 = y1(cond.All,:); 
        y3 = [zscored_all_SC.LDH.R.(names{i})'; zscored_all_SC.RDH.L.(names{i})'];  %% contralateral paw stimulation
        y3 = y3(cond.All, :);

        % Normalized data
        y2 = [normiceall.SC.LDH.L.(names{i})'; normiceall.SC.RDH.R.(names{i})']; %% ipsilateral paw stimulation
        y2 = y2(cond.All,:);
        y4 = [normiceall.SC.LDH.R.(names{i})'; normiceall.SC.RDH.L.(names{i})']; %% contralateral paw stimulation
        y4 = y4(cond.All,:);

        % mean and SEM computing
        [mean1, sem1] = computeMeanSEM(y1);
        [mean2, sem2] = computeMeanSEM(y2);
        [mean3, sem3] = computeMeanSEM(y3);
        [mean4, sem4] = computeMeanSEM(y4);

        % Index de subplot en grille (r,c) -> conditions en colonnes, names en lignes
        subplotIdx = (r-1)*nCols + c;

        %% --- Figure Zscore ---
        figure(figZscore);
        subplot(nRows, nCols, subplotIdx);
        shadedErrorBar(x1, smooth(mean1, smoothParam, smoothMethod), smooth(sem1, smoothParam, smoothMethod), ...
            'lineProps', {'Color', colors{1}}, lineProps{:});
        hold on;
        shadedErrorBar(x1, smooth(mean3, smoothParam, smoothMethod), smooth(sem3, smoothParam, smoothMethod), ...
            'lineProps', {'Color', colors{2}}, lineProps{:});
        axis([-2.5 10 -1 6]); %axis([-2 2 cond.axisZscore]);
        if r == 1
            title(cond.type, 'FontSize', 9); % titre de colonne en haut
        end
        if c == 1
            ylabel({names{i}, 'mean activity (zscore)'}, 'FontSize', 9); % nom de ligne à gauche
        end
        if r == 1 && c == nCols
            legend('Ipsi', 'Contra');
        end

        %% --- Figure Normalized ---
        figure(figNorm);
        subplot(nRows, nCols, subplotIdx);
        shadedErrorBar(x1, smooth(mean2, smoothParam, smoothMethod), smooth(sem2, smoothParam, smoothMethod), ...
            'lineProps', {'Color', colors{3}}, lineProps{:});
        hold on;
        shadedErrorBar(x1, smooth(mean4, smoothParam, smoothMethod), smooth(sem4, smoothParam, smoothMethod), ...
            'lineProps', {'Color', colors{4}}, lineProps{:});
        axis([-2.5 10 cond.axisNorm]);
        if r == 1
            title(cond.type, 'FontSize', 8);
        end
        if c == 1
            ylabel({names{i}, 'normalized activity'}, 'FontSize', 8);
        end
        if r == 1 && c == nCols
            legend('Ipsi', 'Contra');
        end
    end
end


%  mean and SEM computing function
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

