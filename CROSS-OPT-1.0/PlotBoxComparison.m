% -----------------------------------------------------------------
%  PlotBoxComparison.m
% -----------------------------------------------------------------
%  programmers: Marcos Vinicius Issa
%               Anderson Pereira
%               Americo Cunha Jr
%
%  Originally programmed in: Sep 18, 2024
%           Last updated in: Sep 18, 2024
% -----------------------------------------------------------------
% This function plots a boxplot comparison of truss mass results 
% from different optimization methods.
%
% Input:
% x1       - results series 1
% x2       - results series 2
% x3       - results series 3
% graphobj - struct containing graph configuration parameters
%
% Output:
% fig      - handle to the created figure
% ----------------------------------------------------------------- 

function fig = PlotBoxComparison(x1,x2,x3,graphobj)

    % Check number of arguments
    if nargin < 4
        error('Too few inputs.');
    elseif nargin > 4
        error('Too many inputs.');
    end

    % Ensure all input vectors are row vectors
    if ~iscolumn(x1) 
        x1 = x1';
    end
    if ~iscolumn(x2)
        x2 = x2';
    end
    if ~iscolumn(x3)
        x3 = x3';
    end

    % Prepare the data for the boxplot
    x = [x1 x2 x3];

    % Create figure
    fig = figure('Name', graphobj.gname, 'NumberTitle', 'off');

    % Create the boxplot
    boxplot(x,'Notch','off',...
              'Labels',{graphobj.leg1,graphobj.leg2,graphobj.leg3});

    % Set font and box
    set(gcf,'color','white');
    set(gca,'position',[0.2 0.2 0.7 0.7]);
    set(gca,'Box','on');
    set(gca,'TickDir','out','TickLength',[.02 .02]);
    set(gca,'XMinorTick','off','YMinorTick','off');
    set(gca,'XGrid','off','YGrid','off');
    set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
    set(gca,'FontName','Helvetica');
    set(gca,'FontSize',18);
    box on;
    %grid on;

    % Set axis limits
    if ( strcmp(graphobj.xmin,'auto') || strcmp(graphobj.xmax,'auto') )
        xlim('auto');
    else
        xlim([graphobj.xmin graphobj.xmax]);
    end

    if ( strcmp(graphobj.ymin,'auto') || strcmp(graphobj.ymax,'auto') )
        ylim('auto');
    else
        ylim([graphobj.ymin graphobj.ymax]);
    end

    % Set labels
    xlabel(graphobj.xlab, 'FontSize', 20, 'FontName', 'Helvetica');
    ylabel(graphobj.ylab, 'FontSize', 20, 'FontName', 'Helvetica');

    % Set the title
    title(graphobj.gtitle, 'FontSize', 24, 'FontName', 'Helvetica');

    % Save the plot if required
    if strcmp(graphobj.print, 'yes')
        print('-depsc2', [graphobj.gname, '.eps']);
        print('-dpng', [graphobj.gname, '.png']);
    end

    % Close the figure if requested
    if strcmp(graphobj.close, 'yes')
        close(fig);
    end

end
% -----------------------------------------------------------------