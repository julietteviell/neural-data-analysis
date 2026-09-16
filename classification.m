clear all
close all

fileSCRAMAlgo={'m31','m40','m41','m43','m93','m96'};%,,,'SST190preSNISC','SST243preSNISC'};%'preSNI96' ,,'preSNI96' 
fileAnesAwake={'m31','m93','m96'};
files={fileSCRAMAlgo}; % fileAnesAwake for Fig 3 data  analysis
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
                %[x]=crosco(file,p);
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
                    [total_counts_mean_var total_counts_mean_var_norm]= freqVarRAMalg(var3d,varstring,b,bstring,p,file,edges, bin); %var3dL,j
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

%%
clearvars -except edges bstring varstring files stringfile g limi normiceSC miceSC micevar var3d strDH
      
for m = 1:length(strDH)
    for k = 1:length(bstring)
        for i = 1:length(varstring)
            % Initialiser des matrices/cellules vides pour stocker les résultats
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

%% %%%RVM and SC neurons classification%%%q

for g=1:numel(files)
    for   m=1:numel(strDH)
        clearvars -except mice normice edges bstring varstring files stringfile g limi normiceSC miceSC micevar var3d indexall orderall...
            responseheatmapall zscoredall actcells_all inhcells_all miceall normiceall stringfile2 strDH m zscored_all_SC

        %%SCneurons classification
        mice = miceSC.(stringfile{g}).(strDH{m});
        estNonVide = cellfun(@(x) ~isempty(x) && ~(ischar(x) && all(isspace(x))), mice); mice=mice(estNonVide);
        normice = normiceSC.(stringfile{g}).(strDH{m});
        [zscored_all_varcr, basal] = zscoreVarSCRAMalgo(varstring,bstring,mice,g);
        [order]=heatMapNormRAMalgo(mice,bstring,varstring, stringfile,edges,limi,g); % 
        heatMapZscoreSCRAMalgo(zscored_all_varcr,normiceall,bstring,varstring,stringfile,edges,limi,g,order);

        [actcells_allcr inhcells_allcr]= cellClassBootStrapRAMalgo(varstring,bstring,normice,mice,files,edges);
        actcells_all.(stringfile{g}).(strDH{m})=actcells_allcr;
        inhcells_all.(stringfile{g}).(strDH{m})=inhcells_allcr;
        orderall.(['order' stringfile{g}])=order;
        zscored_all_SC.(strDH{m})=zscored_all_varcr;
        zscoredall.(['zscored' stringfile{g}])=zscored_all_SC;
        %[responseheatmap percentagecombined indexLR]=heatmapResponsivenessMatriceOfEvents2(zscored_all_varcr, varstring, bstring, stringfile, actcells_allcr, inhcells_allcr, g)
        %indexall.(['index' stringfile{g}])=indexLR;
        %responseheatmapall.(['responseheatmap' strDH{m}])=responseheatmap;      
    end
end



%%

%%%% Analysis of cell response to paw and heat, for the Ipsi and Contralateral DH/Paw
figure('Name', ['Ipsi contra DH/Paw responses ' stringfile{1}], 'NumberTitle', 'off');%, 'Position', [100 100 1200 800]);

names = varstring;
bin = 0.05;
edges = (-15:bin:15);
x1 = edges(1:end-1);

% Définition des couleurs et paramètres graphiques communs
colors = {
    [0.8 0.3 0.8], [0.2 0.3 0.7], ... % Ipsi/Contra pour zscore
    [0.2 0.5 0.8], [0.2 0.8 0.8]       % Ipsi/Contra pour normalized
};
lineProps = {'transparent', 1, 'patchSaturation', 0.33};
smoothParam = 0.07;
smoothMethod = 'loess';


 % conditions = {
 %            struct('LDH', [PressurePolyLDH PolymodalLDH], 'RDH', [PressurePolyRDH PolymodalRDH], 'type', 'Polycorr', 'axisZscore', [-1 5], 'axisNorm', [0 0.5]), ...
 %            struct('LDH', PolymodalLDH, 'RDH', PolymodalRDH, 'type', 'poly', 'axisZscore', [-1 5], 'axisNorm', [0 0.5]), ...
 %            struct('LDH', ThermalLDH, 'RDH', ThermalRDH, 'type', 'heat', 'axisZscore', [-1 5], 'axisNorm', [0 0.5]), ...
 %            struct('LDH', MechanicalLDH, 'RDH', MechanicalRDH, 'type', 'mecha', 'axisZscore', [-1 12], 'axisNorm', [0 0.5]), ...
 %            % struct('LDH', hcorripsiLDH, 'RDH', hcorripsiRDH, 'type', 'mecha', 'axisZscore', [-1 5], 'axisNorm', [0 0.5]), ...
 %            %  struct('LDH', mcorripsiLDH, 'RDH', mcorripsiRDH, 'type', 'mecha', 'axisZscore', [-1 12], 'axisNorm', [0 0.5]), ...          
 %            };

       conditions =  {
                struct('LDH',corrneuron,'type', 'corr', 'axisZscore', [-1 15], 'axisNorm', [0 0.5]),...         
                struct('LDH',PolyPoly,'type', 'TMLA', 'axisZscore', [-1 15], 'axisNorm', [0 0.5]),...
                struct('LDH',noncorrneuron,'type', 'nocorr', 'axisZscore', [-1 15], 'axisNorm', [0 0.5]),... 
                struct('LDH',Polynotemp,'type', 'TNLA', 'axisZscore', [-1 15], 'axisNorm', [0 0.5]),...
                struct('LDH',Polytemp,'type', 'TLA', 'axisZscore', [-1 15], 'axisNorm', [0 0.5]),...
                struct('LDH',Polynopressure,'type', 'MNLA', 'axisZscore', [-1 15], 'axisNorm', [0 0.5]),...
                struct('LDH',Polypressure,'type', 'MLA', 'axisZscore', [-1 15], 'axisNorm', [0 0.5]),...
                struct('LDH',Mechanical,'type', 'Mecha', 'axisZscore', [-1 15], 'axisNorm', [0 0.5]),...
                struct('LDH',Heatonly,'type', 'heat', 'axisZscore', [-1 15], 'axisNorm', [0 0.5]),...
                struct('LDH',Poly,'type', 'Poly', 'axisZscore', [-1 15], 'axisNorm', [0 0.5]),...   
                };

% conditions = {
%  struct('LDH', LDHPolycorr, 'RDH', RDHPolycorr, 'type', 'Polycorr', 'axisZscore', [-1 4], 'axisNorm', [0 0.5]), ...   
%  %struct('LDH', corrneuron, 'RDH', corrneuronR, 'type', 'Polycorr', 'axisZscore', [-1 12], 'axisNorm', [0 0.5]), ...
%  };
    % struct('LDH', LDHPolycorr, 'RDH', RDHPolycorr, 'type', 'Polycorr', 'axisZscore', [-1 12], 'axisNorm', [0 0.5]), ...
    % struct('LDH', LDHPolymodalpressure_corr, 'RDH', RDHPolymodalpressure_corr, 'type', 'poly_corr', 'axisZscore', [-1 12], 'axisNorm', [0 0.5]), ...
    % struct('LDH', LDHMecha_corr, 'RDH', RDHMecha_corr, 'type', 'mecha_corr', 'axisZscore', [-1 12], 'axisNorm', [0 0.5]), ...
    % struct('LDH', LDHPoly, 'RDH', RDHPoly, 'type', 'poly', 'axisZscore', [-1 12], 'axisNorm', [0 0.5]) ,...
    % struct('LDH', LDHHeatonly, 'RDH', RDHHeatonly, 'type', 'heat', 'axisZscore', [-1 12], 'axisNorm', [0 0.5]), ...
    % struct('LDH', LDHPressureonly, 'RDH', LDHPressureonly, 'type', 'pressure', 'axisZscore', [-1 12], 'axisNorm', [0 0.5]), ...
    % struct('LDH', LDHMecha, 'RDH', RDHMecha, 'type', 'mecha', 'axisZscore', [-1 12], 'axisNorm', [0 0.5]), ...
    % struct('LDH', LDHpolynoci, 'RDH', RDHpolynoci, 'type', 'polynoci', 'axisZscore', [-1 12], 'axisNorm', [0 0.5]), ...
    % struct('LDH', LDHmechanoci, 'RDH', RDHmechanoci, 'type', 'mechanoci', 'axisZscore', [-1 12], 'axisNorm', [0 0.5]), ...   

% 
% sensorytype1=setdiff(find(cidx2==1),type1other);
% sensorytype2=setdiff(find(cidx2==2),type2other);
% sensorytype3=setdiff(find(cidx2==3),type3other);

% Boucle sur chaque condition
for c = 1:length(conditions)
    cond = conditions{c};

    % Boucle sur chaque nom dans names
    for i = [2 4];%2%:length(names)
        %     %%% bilateral DH
        % % Zscored data
        % y1 = [zscored_all_SC.LDH.L.(names{i})(:, cond.LDH)'; zscored_all_SC.RDH.R.(names{i})(:, cond.RDH)'];
        % y3 = [zscored_all_SC.RDH.L.(names{i})(:, cond.RDH)'; zscored_all_SC.LDH.R.(names{i})(:, cond.LDH)'];
        % 
        % % Normalized data
        % y2 = [normiceall.SC.LDH.L.(names{i})(:, cond.LDH)'; normiceall.SC.RDH.R.(names{i})(:, cond.RDH)'];
        % y4 = [normiceall.SC.RDH.L.(names{i})(:, cond.RDH)'; normiceall.SC.LDH.R.(names{i})(:, cond.LDH)'];

            %%%   %%% bilateral PAW
                % Zscored data
        %y1 = [zscored_all_SC.LDH.L.(names{i})(:, cond.LDH)'; zscored_all_SC.RDH.R.(names{i})(:, cond.RDH)'];
        y1 = [zscored_all_SC.LDH.L.(names{i})'; zscored_all_SC.RDH.R.(names{i})'];
        % y1 =y1(cond.LDH,:); 
        y1 =y1(putact,:); 
        %y3 = [zscored_all_SC.LDH.R.(names{i})(:, cond.LDH)'; zscored_all_SC.RDH.L.(names{i})(:, cond.RDH)'];
        y3 = [zscored_all_SC.LDH.R.(names{i})'; zscored_all_SC.RDH.L.(names{i})'];
        %y3 =y3(cond.LDH,:); 
        y3 =y3(putact,:);
        % Normalized data
      %  y2 = [normiceall.SC.LDH.L.(names{i})(:, cond.LDH)'; normiceall.SC.RDH.R.(names{i})(:, cond.RDH)'];
      %  y4 = [normiceall.SC.LDH.R.(names{i})(:, cond.LDH)'; normiceall.SC.RDH.L.(names{i})(:, cond.RDH)'];

        % Calcul des moyennes et SEM
        [mean1, sem1] = computeMeanSEM(y1);
       % [mean2, sem2] = computeMeanSEM(y2);
        [mean3, sem3] = computeMeanSEM(y3);
       % [mean4, sem4] = computeMeanSEM(y4);

        % Subplot pour zscored data
      %  subplot(length(conditions), 1*length(names), (c-1)*1*length(names) + i);
        figure;
        shadedErrorBar(x1, smooth(mean1, smoothParam, smoothMethod), smooth(sem1, smoothParam, smoothMethod), 'lineProps', {'Color', colors{1}}, lineProps{:});
        hold on;
        shadedErrorBar(x1, smooth(mean3, smoothParam, smoothMethod), smooth(sem3, smoothParam, smoothMethod), 'lineProps', {'Color', colors{2}}, lineProps{:});
        ylabel('mean activity (zscore)', 'FontSize', 6);
        legend('Ipsi', 'Contra');
        axis([-2 2 -1 4]) %axis([-2 2 cond.axisZscore]);
        title([names{i} ' - ' cond.type],'Fontsize',6);

        % % Subplot pour normalized data
        % subplot(length(conditions), 2*length(names), (c-1)*2*length(names) + length(names) + i);
        % shadedErrorBar(x1, smooth(mean2, smoothParam, smoothMethod), smooth(sem2, smoothParam, smoothMethod), 'lineProps', {'Color', colors{3}}, lineProps{:});
        % hold on;
        % shadedErrorBar(x1, smooth(mean4, smoothParam, smoothMethod), smooth(sem4, smoothParam, smoothMethod), 'lineProps', {'Color', colors{4}}, lineProps{:});
        % ylabel('normalized activity', 'FontSize', 6);
        % %legend('Ipsi', 'Contra');
        % axis([-2.5 10 cond.axisNorm]);
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

%%

% varmean1=mean( zscored_all_varcr.(bstring{i}).(varstring{2})(1:100,:),1)';
% varmean2=mean( zscored_all_varcr.(bstring{i}).(varstring{2})(100:200,:),1)';
% varmean3=mean( zscored_all_varcr.(bstring{i}).(varstring{2})(200:300,:),1)';
% varmean4=mean(y1(:,200:300),2);
% varmean5=mean( zscored_all_varcr.(bstring{i}).(varstring{2})(400:500,:),1)';
% varmean6=mean( zscored_all_varcr.(bstring{i}).(varstring{2})(500:600,:),1)';
% ordervalues=[ varmean4 ];  % varmean1 varmean2 varmean3 varmean3 varmean5 varmean6
% [a1 index1]=sortrows(ordervalues,[-1 ]);  %-4 -5 -6 -2  -3 -4 
%  
%          meanActivation = mean(zscored_all_varcr.(bstring{i}).(varstring{j}),1);
%         [~, sortedIndices] = sort(meanActivation, 'descend');
%         neunorm_sorted = zscored_all_varcr.(bstring{i}).(varstring{j})(:,order);
%        %('Name','Heatmap','NumberTitle','off');
%         subplot(1,length(varstring),j)
% 
%         figure;
%         %data=zscored_all_varcr.(bstring{i}).(varstring{j})(:,index1)';
%         smoothdata=imgaussfilt(y1(index1,:),0.7); %order.(bstring{i}).(varstring{j})
%         %smoothdata=imgaussfilt(neunorm_sorted',0.7); %order.(bstring{i}).(varstring{j})
%         %h1=imagesc(data);
%         h1=imagesc(smoothdata)
%         set(h1,'XData',edges(1:end-1))
%         %caxis([-4 4])
%         caxis([-3 6])
%         xlim(limi)
%         %xlim([-2 23])
%         colorbar
%         title('selected cells')
%         ylabel('# neurons')
%         xlabel('time (sec)')
%          
%        % colormap('gray')
% beta = .05;
% brighten(beta)
% colormap(CustomColormap)