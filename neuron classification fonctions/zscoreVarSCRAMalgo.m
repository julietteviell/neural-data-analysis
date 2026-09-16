function [zscored_all_varcr, basal] = zscoreVarSCRAMalgo(varstring,bstring,mice,N_type,g,varargin)

size(varargin);

%bstring={'L','R'};
for k=1:size(bstring,2)
    for j=1:size(varstring,2)
        zscored_all=[];
        for l=1:size(mice,1)
            if isfield(mice{l},bstring{k})
                var= mice{l}.(bstring{k}).(varstring{j});%(:,N_type{l,1}));
                varbasal=mice{l}.(bstring{k}).(varstring{j});%(:,N_type{l,1}));
            else
                continue
            end
            clear zscored
            for i=1:size(var,2)
                if j==1;                
                basal=varbasal(1:40,i) %bin 50ms %varbasal(100:150,i); %bin 100 %
                else
                basal=varbasal(1:200,i) %bin 50ms %varbasal(1:100,i); %bin 100 %
                end
                if basal == zeros(size(basal,1),1)
                    random_vector = randi([0, 1], 1, size(basal, 1));
                    basal = random_vector;%(random_vector * -0.5);
                    % Mapper les entiers 0 et 1 aux valeurs 0 et 0.5
                    mean_basal=mean(basal);
                    std_basal=std(basal);
                    zscored1=(var(:,i)-mean_basal)/(std_basal); %smooth(sem2,0.05,'loess'
                    random_vector1 = randi([0, 1], 1, size(zscored1, 1));
                    randiff=(random_vector1* -0.5)';
                    zscored(:, i) = zscored1-randiff;
                else               
                    %basal=varbasal(:,i);
                    mean_basal=mean(basal);
                    std_basal=std(basal);
                    %zscored(:,i)=smooth((var(:,i)-mean_basal)/std_basal,0.05,'loess'); %smooth(sem2,0.05,'loess'
                    zscored(:,i)=(var(:,i)-mean_basal)/(std_basal); %smooth(sem2,0.05,'loess'
                end
            end

            if ~exist('zscored')
                continue
            else
                zscored_all=[zscored_all zscored];
                % zscored_all(isnan(zscored_all))=0;
                % zscored_all(isinf(zscored_all))=0;
            end
            %             figure;
            %             plot(varbasal,'DisplayName','varbasal')
        end
        zscored_all_var.(varstring{j})=zscored_all;

    end
    zscored_all_varcr.(bstring{k})=zscored_all_var;
end
end
