function [] = setfig()
% 02/11/20
% inspired from: https://www.mathworks.com/matlabcentral/answers/183311-setting-default-interpreter-to-latex
% see onenote if link doesn't work

fnm = 'lmodern';
set(groot,'defaultAxesFontName',fnm,'defaultColorbarFontName',fnm,...
    'defaultGeoaxesFontName',fnm,'defaultGraphplotEdgeFontName',fnm,...
    'defaultGraphplotNodeFontName',fnm,'defaultLegendFontName',fnm,...
    'defaultPolaraxesFontName',fnm,'defaultTextFontName',fnm,...
    'defaultTextarrowshapeFontName',fnm,'defaultTextboxshapeFontName',fnm)

fsz = 28;
set(groot,'defaultAxesFontSize',fsz,'defaultColorbarFontSize',fsz,...
    'defaultGeoaxesFontSize',fsz,'defaultGraphplotEdgeFontSize',fsz,...
    'defaultGraphplotNodeFontSize',fsz,'defaultLegendFontSize',fsz,...
    'defaultPolaraxesFontSize',fsz,'defaultTextFontSize',fsz,...
    'defaultTextarrowshapeFontSize',fsz,'defaultTextboxshapeFontSize',fsz)

inp = 'latex';
set(groot,'defaultAxesTickLabelInterpreter',inp,'defaultColorbarTickLabelInterpreter',inp,...
    'defaultGraphplotInterpreter',inp,'defaultLegendInterpreter',inp,'defaultPolaraxesTickLabelInterpreter',inp,...
    'defaultTextInterpreter',inp,'defaultTextarrowshapeInterpreter',inp,'defaultTextboxshapeInterpreter',inp);

lw = 2;
set(groot,'defaultAxesLineWidth',lw,'defaultGraphplotLineWidth',lw,'defaultLineLineWidth',lw,...
    'defaultBarLineWidth',lw,'defaultGeoaxesLineWidth',lw,'defaultLegendLineWidth',lw,...
    'defaultLineshapeLineWidth',lw,'defaultStairLineWidth',lw);

% set(groot,'defaultAxesXGrid','on','defaultAxesXMinorGrid','on',...
%     'defaultAxesYGrid','on','defaultAxesYMinorGrid','on',...
%     'defaultAxesXMinorGridMode','manual','defaultAxesMinorGridAlpha',0.,...
%     'defaultAxesYMinorGridMode','manual');

set(groot,'defaultAxesXGrid','off','defaultAxesXMinorGrid','off',...
    'defaultAxesYGrid','off','defaultAxesYMinorGrid','off',...
    'defaultAxesXMinorGridMode','manual','defaultAxesMinorGridAlpha',0.,...
    'defaultAxesYMinorGridMode','manual');

mrksz = 15;
pos = [100,100,1000,600]; % [left bottom width height]
set(groot,'defaultFigurePosition',pos);
set(groot,'defaultLineMarkerSize',mrksz);

set(0, 'DefaultFigureRenderer', 'opengl');


end