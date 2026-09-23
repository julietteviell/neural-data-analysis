load('processed_data_classification_script_50msbin.mat','dataIPsicontrasanky2','dataIPsicontrasanky')
%% Option 1: use a table

data=dataIPsicontrasankyLANLA; %dataIPsicontrasankyAllcategories % fig S5

%% Customizable options
% Colormap: can be the name of matlab colormaps or a matrix of (N x 3).
%   Important: N must be the max number of categories in a layer 
%   multiplied by the number of layers. 
%   In the example of Option 1, N should be 4 * 4 = 16;
%   In the example of Option 2, N should be 5 * 3 = 15;

layer_labels={'Ipsi,''contra'};

options.color_map = 'parula';      
options.flow_transparency = 0.2;   % opacity of the flow paths
options.bar_width = 40;            % width of the category blocks
options.show_perc = true;          % show percentage over the blocks
options.text_color = [1 1 1];      % text color for the percentages
options.show_layer_labels = true;  % show layer names under the chart
options.show_cat_labels = true;   % show categories over the blocks.
options.show_legend = false;        % show legend with the category names. 
                                   % if the data is not a table, then the
                                   % categories are labeled as catX-layerY
plotSankeyFlowChart(data,options);