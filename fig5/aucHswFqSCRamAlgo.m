function [auc hsw fq wave monophasicall biphasicall wenor1 varargout] = aucHswFqSCRamAlgo(file,p,k)  %auc hsw fq to_remove


clear SPK* wave waver wavec fq hsw auc fqspk auc1 liste Indexc Indexr use_spk count newliste truec
load(file{p},'waveL','waveR','xall')

    if exist('waveL','var')
        assert(size(waveL,2)==size(fieldnames(xall.LDH),1),('Lengths are not equal'));
    else
        waveL=[];
    end
    if exist('waveR','var')
        assert(size(waveR,2)==size(fieldnames(xall.RDH),1),('Lengths are not equal'));
    else
        waveR=[];
    end
    b={xall.LDH,xall.RDH};
    spikes=b{k};
    names= fieldnames(spikes);
    wa={waveL,waveR};
    wave=wa{k};
    
    if ~isempty(wave)
        Num_spikes_mice=length(names);
        
        for f=1:Num_spikes_mice, current_spike=spikes.(names{f});  edges=(1:1:2000); fqspk(f,:)= histcounts(current_spike,edges);  end
        fq=(mean(fqspk,2))';
        
        if size(wave,1)==60
            wave=wave(8:49,:);
        else
            wave=wave(1:42,:);
        end
        auc1 = zeros(size(wave,2),1);
        hsw = zeros(size(wave,2),1);
        xmax1=zeros(size(wave,2),1);
        xmin1=zeros(size(wave,2),1);
        monophasicall=zeros(size(wave,2),1);
        biphasicall=zeros(size(wave,2),1);
        wenor1=zeros(size(wave,2),1000);

        for c = 1:size(wave,2),
            clear wenor
            we = resample(wave(:,c),300,12);
            wf(:,c)=we';
            [xmin,idxmin]=min(we); [xmax,idxmax]=max(we(idxmin:end));
            %wf(:,c)=-wf(:,c)./min(we);
            m = max(we(find(we==min(we)):end));
            wenor=we./m;
            [m,idx]=min(we); we=we(idx:end); we=we(we>0); [mnor,idxnor]=min(wenor);  wenor=wenor(idxnor:end); wenor=wenor(wenor>0); 
            x_original = linspace(1, length(wenor), length(wenor))'; x_interpole = linspace(1, length(wenor), 1000)';
            wenor = interp1(x_original, wenor, x_interpole, 'linear');
            xmax1(c)=xmax;
            xmin1(c)=xmin;
            % 
            % figure;
            % plot(we)
            auc1(c) = trapz(we);
            auc = auc1;
            hsw(c) = idxmax;
            wenor1(c,:)=wenor';
        end
       monophasicall(find(xmax1<0.025))=1;
       biphasicall(find(xmax1>0.025))=1;
        % figure;
        % plot(wave)
        if auc == 0
            % nothing is returned
            return;
        else
        end
        if hsw == 1
            % nothing is returned
            return;
        else
        end
        if ~exist('fq')
            % nothing is returned
            return;
        else
        end
    else
        auc=[]; hsw=[]; fq=[]; wave=[]; monophasicall=[];  biphasicall=[]; wenor1=[];
    end
end