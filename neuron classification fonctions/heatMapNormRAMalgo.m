
function [order]=heatMapNormRAMalgo(mice,bstring,varstring, ~,edges,limi,~,varargin)

size(varargin)

Matneunorm=cell2mat(mice);

for k=1:size(bstring,1)
    figure('Name','Heatmap','NumberTitle','off');
    for j=1:size(varstring,2)
        neunorm=[];
        for i = 1:size(Matneunorm,1)
            neunorm=[neunorm Matneunorm(i).(bstring{k}).(varstring{j})];
        end
        neunorm_var.(bstring{k}).(varstring{j})=neunorm;
        meanActivation = mean(neunorm_var.(bstring{k}).(varstring{j}),1);
        [~, sortedIndices] = sort(meanActivation, 'descend');
        neunorm_sorted = neunorm_var.(bstring{k}).(varstring{1})(:,sortedIndices);

        if j == 1
            limi = [-6 10];
        elseif j == 2
            limi = [-2 7];
        elseif j == 3
            limi = [-5 10];
        end
        % varmean1=mean( neunorm_var.(bstring{k}).(varstring{j})(1:100,:),1)';
        % varmean2=mean( neunorm_var.(bstring{k}).(varstring{j})(100:200,:),1)';
        % varmean3=mean( neunorm_var.(bstring{k}).(varstring{j})(200:300,:),1)';
        % varmean4=mean( neunorm_var.(bstring{k}).(varstring{j})(300:400,:),1)';
        % varmean5=mean( neunorm_var.(bstring{k}).(varstring{j})(400:500,:),1)';
        % varmean6=mean( neunorm_var.(bstring{k}).(varstring{j})(500:600,:),1)';
        % ordervalues=[varmean4 varmean5 varmean6];  % varmean1 varmean2 varmean3
        % [a1 index1]=sortrows(ordervalues,[-1 -2 -3   ]);  %-4 -5 -6
        order=sortedIndices;

        subplot(1,length(varstring),j)
        h1=imagesc(neunorm_sorted');
        %h1=imagesc(neunorm_var.(bstring{k}).(varstring{1})(:,index1)');
        set(h1,'XData',edges(1:end-1))
        %caxis([0 20])
        xlim(limi)
        caxis([0.05 0.3])
        title(varstring{j})
        ylabel('# neurons')
        xlabel('time (sec)')
        colormap summer
        colorbar
    end

end


%suptitle([bstring{i} ' ' stringfile{g}])


end