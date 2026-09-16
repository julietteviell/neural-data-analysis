function [total_counts_mean_var, total_counts_mean_var_norm] = freqVarRAMalg(var3d,varstring,b,bstring,p,file,edges,bin,varargin)
    results = struct();
    total_counts_mean_var = struct();
    total_counts_mean_var_norm = struct();

    for k = 1:size(fieldnames(var3d),1)
        clear varlist total_counts_mean total_counts_mean_norm
        Num_spikes_mice = length(b);

        if Num_spikes_mice == 0 || Num_spikes_mice == 1
            continue
        else
            for j = 1:size(var3d.(['var3d' bstring{k}]),2)
                varlist = var3d.(['var3d' bstring{k}]){j}(:,1);
                Yall = [];
                clear Counts_mean varnorm Varnorm varnorm1
                for f = 1:Num_spikes_mice
                    current_spike = b{f};
                    binrange = (0:bin:current_spike(end));
                    binned_spike = histcounts(current_spike,binrange);
                    maxvar = max(binned_spike);
                    minvar = min(binned_spike);
                    clear Binned varnorm var Allsess
                    for m = 1:length(varlist)
                        Rel = current_spike-varlist(m);
                        Binned(:,m) = histcounts(Rel,edges);
                    end
                    Counts_mean(:,f) = mean(Binned,2);
                    varnorm = sum(Binned,2);
                    elem1 = varnorm - min(varnorm);
                    elem2 = max(varnorm) - min(varnorm);
                    if elem2 == 0
                        elem2 = 0.001;
                        varnorm1 = elem1./elem2;
                    else
                        varnorm1 = elem1./elem2;
                    end
                    Yall = horzcat(Yall,varnorm1);
                    Yall(isnan(Yall)) = 0.001;
                end
                total_counts_mean.(varstring{j}) = Counts_mean;
                total_counts_mean_norm.(varstring{j}) = Yall;
            end
            total_counts_mean_var.(bstring{k}) = total_counts_mean;
            total_counts_mean_var_norm.(bstring{k}) = total_counts_mean_norm;
        end
    end
end
