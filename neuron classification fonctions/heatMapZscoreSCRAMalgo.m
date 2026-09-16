function [order]=heatMapZscoreSCRAMalgo(zscored_all_varcr,normiceall,bstring,varstring,stringfile,edges,limi,g,order,varargin)

size(varargin)

%varstring= fieldnames(zscored_all_varcr.(bstring{1}))'; 

load('Mycolormap')
for i=1:size(bstring,2)

    figure('Name',['zscore SC neurons' ' ' bstring{i} ],'NumberTitle','off');
    for j=1:size(varstring,2)


       if j == 2
            limi = [-5 10];
        elseif j == 3
            limi = [-2 10];
        elseif j == 4
            limi = [-5 10];
        end
 y1 = [zscored_all_SC.LDH.L.(varstring{1}) zscored_all_SC.RDH.R.(varstring{1})]'; 
 y2 = [zscored_all_SC.LDH.R.(varstring{1}) zscored_all_SC.RDH.L.(varstring{1})]'; 
%varmean1=mean( zscored_all_varcr.(bstring{i}).(varstring{1})(40:80,:),1)';
% varmean2=mean( zscored_all_varcr.(bstring{i}).(varstring{1})(100:200,:),1)';
% varmean3=mean( zscored_all_varcr.(bstring{i}).(varstring{1})(200:300,:),1)';
% varmean4=mean( zscored_all_varcr.(bstring{i}).(varstring{1})(300:400,:),1)';
varmean4=mean( y1(:,37:60),2);
%varmean5=mean( zscored_all_varcr.(bstring{i}).(varstring{2})(400:500,:),1)';
%varmean6=mean( zscored_all_varcr.(bstring{i}).(varstring{2})(500:600,:),1)';
ordervalues=[ varmean4 ];  % varmean1 varmean2 varmean3 varmean3 varmean5 varmean6

[a1 index1]=sortrows(ordervalues,[-1 ]);  %-4 -5 -6 -2  -3 -4 

        %meanActivation = mean(zscored_all_varcr.(bstring{i}).(varstring{j}),1);
        %[~, sortedIndices] = sort(meanActivation, 'descend');
        neunorm_sorted = y2(index1,:);
       %('Name','Heatmap','NumberTitle','off');
        %subplot(1,length(varstring),j)
        %figure;
        %data=zscored_all_varcr.(bstring{i}).(varstring{j})(:,index1)';
        %smoothdata=imgaussfilt(zscored_all_varcr.(bstring{i}).(varstring{j})(:,index1)',0.7); %order.(bstring{i}).(varstring{j})
        %smoothdata=imgaussfilt(zscored_all_varcr.(bstring{i}).(varstring{j})',0.7); %order.(bstring{i}).(varstring{j})
        smoothdata=imgaussfilt(neunorm_sorted,0.7); %order.(bstring{i}).(varstring{j})
        %h1=imagesc(data);
        h1=imagesc(smoothdata)
        set(h1,'XData',edges(1:end-1))
        caxis([-3 3])
        %caxis([0 0.8])
        %xlim(limi)
        xlim([-2 2])
        colorbar
        title(varstring{1})
        ylabel('# neurons')
        xlabel('time (sec)')
        
       colormap('gray')
beta = .05;
brighten(beta)
colormap(CustomColormap2)
colormap(mymap2)
% colormap winter
% colormap parula
% colormap(infernoColors)
% colormap(viridisColors)
% colormap(plasmaColors)


    end
    
end



end


