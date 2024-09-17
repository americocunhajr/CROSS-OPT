% -----------------------------------------------------------------
%  Main_Truss10_MassMin_FrequencyBounds_CE.m
% -----------------------------------------------------------------
%  programmers: Marcos Vinicius Issa
%               Anderson Pereira
%               Americo Cunha Jr
%
%  Originally programmed in: Jun 25, 2024
%           Last updated in: Aug 29, 2024
% -----------------------------------------------------------------
%  ﻿Mass minimization with frequency constraints for a 10 bars
%  truss structure using the cross-entropy method.
% -----------------------------------------------------------------

clc; clear; close all;

disp(' ----------------------------------------- ')
disp(' Main_Truss10_MassMin_FrequencyBounds_CE.m ')
disp(' ----------------------------------------- ')

% random number generator (fix the seed for reproducibility)
%rng_stream = RandStream('mt19937ar','Seed',30081984);
%RandStream.setGlobalStream(rng_stream);

% convertion factors
Inch2Meter  = 0.0254;    % inch to meter  factor

% truss parameters
l1         = 360.0*Inch2Meter;  % 1st length (m)
l2         = 360.0*Inch2Meter;  % 2nd length (m)
h          = 360.0*Inch2Meter;  % height (m)
rho        = 2770.0;            % material density (kg/m^3)
E          = 69.8e9;            % elastic modulus (Pa)
AddedMass  = 454.0;             % added mass (kg)
omegaTresh = 2*pi*[7 15 20]';   % treshold frequencies (rad/s)

% sctruct with the truss parameters
MyTruss.l1         = l1;           % 1st length (m)
MyTruss.l2         = l2;           % 2nd length (m)
MyTruss.h          = h;            % height (m)
MyTruss.rho        = rho;          % material density (kg/m^3)
MyTruss.E          = E;            % elastic modulus (Pa)
MyTruss.AddedMass  = AddedMass;    % added mass (kg)
MyTruss.omegaTresh = omegaTresh;   % treshold frequencies (rad/s)
MyTruss.FixedDoFs  = [9 10 11 12]; % Diriclet boundary condition DoFs
MyTruss.NODES      = [(l1+l2),h;   % nodes coordinates
                      (l1+l2),0;
                       l1    ,h;
                       l1    ,0;
                        0    ,l1;
                        0    ,0];
MyTruss.ELEM       = [5,3;         % element nodes
                      3,1; 
                      6,4; 
                      4,2; 
                      4,3; 
                      2,1; 
                      5,4; 
                      3,6; 
                      3,2; 
                      1,4];
MyTruss.Nnodes     = size(MyTruss.NODES,1); % # of nodes
MyTruss.Nelem      = size(MyTruss.ELEM ,1); % # of elements
MyTruss.Ndofs      = 2*MyTruss.Nnodes;      % # of DoFs

% objective and constraint functions
fun     = @(x)TrussMass10           (x,MyTruss);
nonlcon = @(x)TrussFrequencyBounds10(x,MyTruss);

% number of variables
Nvars = 10;

% lower and upper bounds for design variables
lb =  0.1*ones(1,Nvars)*Inch2Meter^2;
ub = 20.0*ones(1,Nvars)*Inch2Meter^2;

% initial mean and standard deviation
xmean0 = 0.5*(ub+lb);
sigma0 = 5*(ub-lb);

% CE solver
tic
[Xopt,Fopt,ExitFlag,CEstr] = CEopt(fun,xmean0,sigma0,lb,ub,nonlcon);
toc

% check constraint violation
disp(' ')
disp(' Check if inequality constraints are <= 0')
disp(' ')
G = nonlcon(Xopt);
for i=1:length(G)
    fprintf('  G(%d) =  %5.2f\n',i,G(i));
end
disp(' ')
if sum(G > CEstr.TolCon) ~= 0
    disp(' One or more constraints violated :-( ')
else
    disp(' All constraints respected :-) ')
end

% display optimized mass
fprintf('\n  Original truss mass (kg): %5.2f',fun(xmean0));
fprintf('\n Optimized truss mass (kg): %5.2f',Fopt       );
disp(' ')

% plot the nominal truss structure
PlotTruss10(xmean0*400,MyTruss,'Non-optimal Truss Structure');

% plot the optimal truss structure
PlotTruss10(  Xopt*400,MyTruss,'Optimal Truss Structure');
