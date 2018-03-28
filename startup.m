% This Matlab startup.m file sets up default values
% comments added on 11.04.2017
% ******* R. Cardot *******************************

%% DEFAULT VALUES FOR FIGURES %%
set(groot,'DefaultAxesFontName','Arial'); % Font Type
set(groot,'DefaultAxesFontSize',22); % Font Size
set(groot,'DefaultFigureColor','w'); % Background Color
HeightScreen=get(groot,'ScreenSize'); % Get the Screen size in pixels
HeightScreen=HeightScreen(4); % Get th height of the screen
set(groot,'DefaultFigurePosition',[0 35 HeightScreen-110 HeightScreen-110]); % Set up default figure size as a square of 'height of the screen' -20 pixels
clear HeightScreen
% set(groot,'TitleFontSizeMultiplier',1.2);
% set(groot,'LabelFontSizeMultiplier',1.15);