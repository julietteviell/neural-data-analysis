
function [x]=crosco(file,p)

load(file{p},'*L','*R','Cha*','SPK*')

allVariables = whos('Cha*');
if isempty(allVariables)
    allVariables = whos('SPK*');
end

b = struct('name', {}, 'size', {}, 'bytes', {}, 'class', {}, 'global', {}, 'sparse', {}, 'complex', {}, 'nesting', {}, 'persistent', {});
x=struct();
% Parcourir les variables et filtrer celles avec les suffixes indésirables
for i = 1:length(allVariables)
    varName = allVariables(i).name;
    % Vérifier si la variable ne contient pas les suffixes à exclure
    if isempty(regexp(varName, '_wf$|_template$|_template_ts$|_wf_ts$'))
        b(end+1) = allVariables(i)';
        x.(varName)=eval(varName);
    end
end

b=b';

numberspikes=size(b,1);

cell_spk=cell(numberspikes,1);
count=0;

for nameLOOP=1:numberspikes;
    
 curr_var= b(nameLOOP,1).name
 
    if strcmp(curr_var(1:3),'Cha')
        count=count+1;
       
        a=eval(curr_var);
       cell_spk{count,1}=a;
    end
end

spikes=cell_spk;
binning=(0.001:0.001:4800); cellbin=repmat(binning,size(spikes,1),1);
cellbin1=mat2cell(cellbin,ones(1,size(cellbin,1)), size(cellbin,2));
res=cellfun(@ histc,spikes,cellbin1,'UniformOutput',false);
A=cell2mat(res');


m=1:size(A,2);
[CROSS,P] = corrcoef(A(:,m));

figure;
imagesc(CROSS)
xticks((1:1:length(CROSS)))
yticks((1:1:length(CROSS)))
caxis([0 0.5])
colorbar

[c,d]=find(CROSS>0.5);
DIAG=find(c==d);
c(DIAG)=[];
d(DIAG)=[];

if ~isempty(c)
    e=horzcat(c,d);
    % if isempty(c)

    A = e';
    [B index]=unique(A(2,:).','rows', 'stable');
    % if length(index)==length(e)
    %     listedel=A(1,1);
    % else
    e(index,:)=[];
    %      A = e';
    %     [B index]=unique(A(1,:).','rows', 'stable');
    %     e(index,:)=[];
    C=unique(sort(e(:,1),'ascend'));
    D=unique(sort(e(:,2),'ascend'));
    listedel=C';
    % end

    names = fieldnames(x);
    % listedel=unique(e(:,2),'rows')';
    for iname = listedel
        % x.(['del_', names{iname}]) = x.(names{iname});
        x = rmfield(x, names{iname});
    end
else
    x=x;
    % end
end
    save(file{p},'x','-append')
end

