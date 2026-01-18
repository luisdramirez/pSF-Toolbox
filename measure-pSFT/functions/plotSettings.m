% plotSettings
%   Returns a structure of plot settings for the pSFT toolbox

function plot_settings = plotSettings()

%% Font settings

plot_settings.font_type = 'Helvetica';
plot_settings.axes_label_font_size = 14;
plot_settings.axes_tick_font_size = 13;
plot_settings.tick_length = 0.020;
plot_settings.line_width = 1;

%% Color settings

plot_settings.colors.red = [204 0 0]/255;
plot_settings.colors.green = [0 153 0]/255;
plot_settings.colors.blue = [0 76 152]/255;
plot_settings.colors.purple = [102 51 204]/255;
plot_settings.colors.yellow = [0.9, 0.8, 0];
plot_settings.colors.orange = [0.8, 0.5, 0.1];
plot_settings.colors.black = [0 0 0];
plot_settings.colors.white = [1 1 1];
plot_settings.colors.gray = plot_settings.colors.white/2;

end