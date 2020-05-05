%SCRIPT TO FIND THE MINIMAL CODEBOOK OF GENE EXPRESSION THAT ESTABLISHES
%NEURON IDENTITY GIVEN A SET OF GENE EXPRESSION FOR A GROUP OF CELLS

%Copyright Molly B. Reilly, Cyril Cros, Erdem Varol, Eviatar Yemini, and Oliver Hobert 2020   


clear all
clc
close all

%PARAMETERS
iter = 100;
set(0, 'DefaultTextInterpreter', 'none')
set(0, 'DefaultAxesTickLabelInterpreter', 'none')

%INPUT FILE:
filename = 'all conserved hbox 040420.xlsx'; % Input excel sheet should have cells (rows) x genes (columns) with headers and row-names
table=readtable(filename);
T=table2array(table(:,2:end));
tf_names=table.Properties.VariableNames(2:end);
tf_names=strrep(tf_names,'_','-');
tf_names=cellfun(@(x) ['{\it ' x '}'], tf_names, 'UniformOutput', false);
neurons=table2array(table(:,1));

% IDENTIFY UNIQUE NEURONS --- if the input table includes neurons that are
% indistinguishable by genetic expression, only the unique ones will be
% kept
[a,b,c]=unique(T,'rows','legacy');
T = T(b,:);
neurons=neurons(b);


% MINIMAL CODEBOOK SOLUTION
[minimal_codebook,indices]=optimal_codebook_solver(T);
Topt=T(:,indices); %retreiving the minimal set of genes necessary for identity
tf_names_opt=tf_names(indices); %retreiving the minimal set of genes necessary for identity


% display minimal codebook
figure('units','normalized','outerposition',[0 0 1 1])
h=heatmap(Topt);
h.Colormap=colormap('jet');
h.XDisplayLabels=tf_names_opt;
h.YDisplayLabels=neurons;
title(['Minimal codebook - ' num2str(length(indices)) ' genes'])
% display full codebook just for sanity check
figure('units','normalized','outerposition',[0 0 1 1])
h=heatmap(T);
h.Colormap=colormap('jet');
h.XDisplayLabels=tf_names;
h.YDisplayLabels=neurons;
title(['Full codebook - ' num2str(length(tf_names)) ' genes'])

% display codebook correlation to all Tf's
figure('units','normalized','outerposition',[0 0 1 1])
h=heatmap(corr(Topt,T));
h.Colormap=colormap('jet');
h.XDisplayLabels=tf_names;
h.YDisplayLabels=tf_names_opt;
title('Minimal codebook vs. all genes correlation')

% WRITE RESULTS TO FILES

writetable(cell2table(tf_names_opt'),'minimal_gene_set.csv');
Tmin(:,1)=array2table(neurons);
Tmin=[Tmin array2table(Topt)];
Tmin.Properties.VariableNames(1)={'neurons'};
Tmin.Properties.VariableNames(2:end)=matlab.lang.makeValidName(tf_names_opt);
writetable(Tmin,'minimal_codebook.xlsx');