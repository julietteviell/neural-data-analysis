function [actcells_allcr inhcells_allcr]= cellClassBootStrapRVM(varstring,bstring,normice,mice,files,edges,bin,varargin)

for k=1:size(bstring,2)
    for j=1:size(varstring,2)
        clear active inhib
        %zero=find(edges>0&edges<0.01);
        bin=0.050;         %% time bin to change         %% egdes before and after 0
        edges=(-2:bin:2);
        limi=[-2 2];
        zero=find(edges==0);
        preevent=[];
        postevent=[];
        for l=1:size(mice,1)
            clear pre post
            if isempty(normice{l})
                continue
            else
            %zero=find(edges>0&edges<0.01);
            wind_selectpost=zero:zero+39; wind_selectpre=zero-39:zero; %wind_selectpost=zero:zero+100; wind_selectpre=zero-299:zero-200; ramalgo
            pre= normice{l}.(bstring{k}).(varstring{j})(wind_selectpre,:);
            post= normice{l}.(bstring{k}).(varstring{j})(wind_selectpost,:);
            preevent=[preevent pre];
            postevent=[postevent post];
            h2=mes(postevent,preevent,'U3','isDep',1,'nBoot',10000,'ROCtBoot',1,'doPlot',0,'missVal','listwise');%{'md'},'mdbysd','U1','U3','auroc','hedgesg','tailratio','glassdelta'});
            active=intersect((find(h2.t.p<0.05)),(find(h2.t.tstat>3.2)));
            inhib=intersect((find(h2.t.p<0.05)),(find(h2.t.tstat<-3.2)));
            actcells.(varstring{j})=active;
            inhcells.(varstring{j})=inhib;
            end
        end

        %         actcells.(varstring{j})=SC_neurons.ONcells;
        %         inhcells.(varstring{j})=SC_neurons.OFFcells;
    end
    actcells_allcr.(bstring{k})=actcells;
    inhcells_allcr.(bstring{k})=inhcells;

end
end
