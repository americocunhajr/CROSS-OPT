% -----------------------------------------------------------------
%  PlotTruss37.m
% -----------------------------------------------------------------
%  programmers: Marcos Vinicius Issa
%               Anderson Pereira
%               Americo Cunha Jr
%
%  Originally programmed in: Apr 04, 2024
%           Last updated in: Sep 13, 2024
% -----------------------------------------------------------------
%  Plot a 37 bars truss given the elements cross-section area.
% -----------------------------------------------------------------
function PlotTruss37(x,MyTruss,MyTitle)

    % truss structure parameters
    a       = MyTruss.a;
    ELEM    = MyTruss.ELEM;
    SYMBARS = MyTruss.SYMBARS;
    Nelem   = MyTruss.Nelem;

    % number of design variables
    NvarsA = MyTruss.NvarsA;
    NvarsH = MyTruss.NvarsH;

    % design variables vectors
    Area = zeros(Nelem,1);
    for i=1:NvarsA
        % area design variables
        Area(SYMBARS{i}) = x(i);
    end
    if NvarsH > 0
        % height design variables
        Height = x(NvarsA+1:end);

        % mesh points in global coordinates (Nnodes x 2)
        NODES = [    0 0;
                     a 0;
                     a Height(1);
                   2*a 0;
                   2*a Height(2);
                   3*a 0;
                   3*a Height(3);
                   4*a 0;
                   4*a Height(4);
                   5*a 0;
                   5*a Height(5);
                   6*a 0;
                   6*a Height(4);
                   7*a 0;
                   7*a Height(3);
                   8*a 0;
                   8*a Height(2);
                   9*a 0;
                   9*a Height(1);
                  10*a 0];
    else
        % mesh points in global coordinates (Nnodes x 2)
        NODES = MyTruss.NODES;
    end

    % custom color
    grayColor = [.7 .7 .7];

    % Open figure
    figure('DefaultAxesFontSize',10)
    clf
    hold on
    
    % Support 1
    xl1 = [-0.3 0.3]; 
    yl1 = [-0.1 -0.075];
    [X1,Y1] = hatch_coordinates(xl1, yl1, 0.025);
    plot(X1,Y1,'Color', grayColor,'linewidth',0.01);
    a11 = [0; -0.3; 0.3];
    a12 = [0; -0.075; -0.075];
    patch(a11,a12,grayColor,'EdgeColor','none');
    plot([a11(2) a11(3)],[a12(2) a12(3)],'k','linewidth',2.5);

    % Support 2
    xl1 = [9.7 10.3];
    yl1 = [-0.135 -0.11];
    [X1,Y1] = hatch_coordinates(xl1, yl1, 0.025);
    plot(X1,Y1,'Color', grayColor,'linewidth',0.01);
    a21 = [10; 9.7; 10.3];
    a22 = [0; -0.075; -0.075];
    patch(a21,a22,grayColor,'EdgeColor','none');
    plot([a21(2) a21(3)],[a22(2) a22(3)],'k','linewidth',2.5);
    plot([a21(2) a21(3)],[-0.11 -0.11]  ,'k','linewidth',2.5);
    f = [1 2];
    v = [9.8 -0.0925; 10.2 -0.0925];
    patch(        'Faces'          ,f,...
                  'Vertices'       ,v,...
                  'EdgeColor'      ,'none',...
                  'FaceColor'      ,'none', ...
                  'MarkerEdgeColor','k',...
                  'Marker'         ,'o',...
                  'MarkerFaceColor','white',...
                  'MarkerSize'     ,10,... 
                  'LineWidth'      ,0.01);
    
    % Custom color scheme
    mycolor = plasma(NvarsA-1);
    
    % Draw vertical and top horizontal bars
    for i=1:NvarsA-1
        idx1 = SYMBARS{i}(1);
        patch('Faces'    ,ELEM(idx1,:)   ,...
              'Vertices' ,NODES       ,...
              'EdgeColor',mycolor(i,:),...
              'FaceColor','none'      ,...
              'LineWidth',Area(i));

        idx2 = SYMBARS{i}(2);
        patch('Faces'    ,ELEM(idx2,:)   ,...
              'Vertices' ,NODES       ,...
              'EdgeColor',mycolor(i,:),...
              'FaceColor','none'      ,...
              'LineWidth',Area(i));
    end
    
    % Draw horizontal bars at the bottom
    for i = SYMBARS{NvarsA}
        patch('Faces'    ,ELEM(i,:),...
              'Vertices' ,NODES    ,...
              'EdgeColor',grayColor,...
              'FaceColor','none'   ,...
              'LineWidth',Area(i));
    end
    
    % Draw nodes
    for i=1:Nelem
        patch('Faces'          ,ELEM(i,:),...
              'Vertices'       ,NODES,...
              'EdgeColor'      ,'none',...
              'FaceColor'      ,'none', ...
              'MarkerEdgeColor','b',...
              'Marker'         ,'o',...
              'MarkerFaceColor','white',...
              'MarkerSize'     ,10,... 
              'LineWidth'      ,3);
    end
    
    % Draw added mass
    for i = 2:2:18
        plot([NODES(i,1)],[NODES(i,2)], ...
             'Marker'                 , 'o',...
             'MarkerEdgeColor'        ,[0.9290 0.6940 0.1250],...
             'LineWidth'              , 3.1 , ....
             'MarkerSize'             , 15);
    end
    
    set(gca,'xtick',[],'ytick',[])
    set(gca,'XColor', 'none','YColor','none')
    title(MyTitle,'FontSize',18)
    pause(1)
end
% ------------------------------------------------------------

% ------------------------------------------------------------
% This function returns the coordinates for plotting a hatch
% pattern. The angle of the lines can be adjusted by varying
% the ratio xstep/ystep. xlim and ylim are vectors with two 
% elements, where the first element needs to be smaller than 
% the second.
% ------------------------------------------------------------
function [X,Y] = hatch_coordinates(xlim,ylim,xstep,ystep,merge)

    % set default options
    if nargin < 3 ; xstep = 1     ; end
    if nargin < 4 ; ystep = xstep ; end
    if nargin < 5 ; merge = true  ; end

    % define base grid
    xpos = xlim(1):xstep:xlim(2) ; nx = numel(xpos) ;
    ypos = ylim(1):ystep:ylim(2) ; ny = numel(ypos) ;

    % Create the coordinates
    nanline = NaN*ones(1,nx+ny-3) ;
    X = [ [xpos(1)*ones(1,ny-2) xpos(1:end-1)]; ...
          [xpos(2:end) xpos(end)*ones(1,ny-2)]; ...
           nanline];
    Y = [ [ypos(end-1:-1:1) ylim(1)*ones(1,nx-2)]; ...
          [ypos(end)*ones(1,nx-1) ypos(end-1:-1:2)]; ...
           nanline];

    % merge if asked too
    if merge
        X = X(:) ;
        Y = Y(:) ;
    end
end
% ------------------------------------------------------------