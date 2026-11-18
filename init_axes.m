% Function to initialise axes with latex axis labels and optional title
function ax = init_axes(ax,xlabeltext,ylabeltext,varargin)

cla(ax);
hold(ax, 'on');
ax.FontSize = 14;
ax.TickLabelInterpreter = 'latex';

xlabel(ax, xlabeltext, 'FontSize', 14, 'Interpreter', 'Latex')
ylabel(ax, ylabeltext, 'FontSize', 14, 'Interpreter', 'Latex')

if nargin >= 4
    xlims = varargin{1};
    xlim(ax, xlims)
end

if nargin >= 5
    ylims = varargin{2};
    ylim(ax, ylims)
end
    
end