function [] = printsetup(FigHandle, figure_width, figure_height, font_size, line_width, export_ppi, print_png, print_pdf, savename)


% PRINT PRIORITY METHOD
% Advantage: Lets you define precisely how large the figure should be when
% exported (in inches and ppi)
% Disadvantage: You cannot influence how large the figure looks like on
% screen. Small figures will look very small.

% WE ALSO NEED
screen_ppi = 72; % inherent property of Matlab. CHANGING THIS WILL MAKE FIGURE LOOK INCONSISTENT BETWEEN SCREEN AND EXPORT!

% DERIVED PROPERTIES (cannot be changed; for info only)
screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels

% SET UP FIGURE SIZE
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

% FONT
allaxes = findall(FigHandle, 'type', 'axes');
set(allaxes,'FontSize', font_size)
set(findall(allaxes,'type','text'),'FontSize', font_size)

% LINES
lines = findall(FigHandle, 'type', 'line');
set(lines, 'LineWidth', line_width);

% PRINT
if print_pdf
    print(gcf, '-dpdf', strcat(savename, '.pdf'));
end
if print_png
    print(gcf, '-dpng', strcat('-r',num2str(export_ppi)), strcat(savename, '.png'));
end

end