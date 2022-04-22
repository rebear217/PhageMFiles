close all
clc

M = data.infection.phi';
M = floor(60*M/max(max(M)));
M = M(1:end-1,:);
M = M(:,[1:92 , 96]);

bp = Bipartite(M);
bp.row_labels = {data.infection.bacteria_names{1:end-1}}';
bp.col_labels = {data.infection.phage_names{[1:92 , 96]}}';

bp

%%

set(0,'RecursionLimit',100)

bp.community.Detect();
bp.nestedness.Detect();

bp.plotter.use_type_interaction = true;

%%

figure(1)
bp.plotter.font_size = 8.0;
bp.plotter.use_isocline = true;
bp.plotter.use_module_format = false;
bp.plotter.use_type_interaction = false;

bp.community = LPBrim(bp.matrix);
bp.community.optimize_by_component = true;
% optimize by components
bp.plotter.PlotModularMatrix();

%%

set(1,'pos',[742         397        1396         948])
export_fig('./figures/modularityPhiMatrix.pdf');

%%

figure(2)
set(2,'pos',[742           1         763        1344])
bp.plotter.PlotModularGraph();
drawnow
export_fig('./figures/modularityPhiGraph.pdf');

%%

set(0,'RecursionLimit',150)
bp.statistics.DoCompleteAnalysis(1000, @NullModels.EQUIPROBABLE);
bp.printer.PrintStructureStatistics();
