% -----------------------------------------------------------------
%  Main_Truss37_MassMin_FrequencyBounds_CE.m
% -----------------------------------------------------------------
%  programmers: Marcos Vinicius Issa
%               Anderson Pereira
%               Americo Cunha Jr
%
%  Originally programmed in: Apr 04, 2024
%           Last updated in: Sep 14, 2024
% -----------------------------------------------------------------
%  ﻿Mass minimization with frequency constraints for a 37 bars
%  truss structure using the cross-entropy method.
% -----------------------------------------------------------------

clc; clear; close all;

disp(' ----------------------------------------- ')
disp(' Main_Truss37_MassMin_FrequencyBounds_CE.m ')
disp(' ----------------------------------------- ')

% random number generator (fix the seed for reproducibility)
%rng_stream = RandStream('mt19937ar','Seed',30081984);
%RandStream.setGlobalStream(rng_stream);

% truss parameters
a          = 1.0;               % bar length (m)
h          = 1.0;               % truss height (m)
rho        = 7800.0;            % material density (kg/m^3)
E          =  210.0e9;          % elastic modulus (Pa)
AddedMass  = 10.0;              % added mass (kg)
omegaTresh = 2*pi*[20 40 60]';  % treshold frequencies (rad/s)

% sctruct with the truss parameters
MyTruss.a          = a;            % bar length (m)
MyTruss.h          = h;            % truss height (m)
MyTruss.rho        = rho;          % material density (kg/m^3)
MyTruss.E          = E;            % elastic modulus (Pa)
MyTruss.AddedMass  = AddedMass;    % added mass (kg)
MyTruss.omegaTresh = omegaTresh;   % treshold frequencies (rad/s)
MyTruss.FixedDoFs  = [1 2 40];     % Diriclet boundary condition DoFs
MyTruss.NODES      = [  0 0;       % nodes coordinates
                        a 0;
                        a h;
                      2*a 0;
                      2*a h;
                      3*a 0;
                      3*a h;
                      4*a 0;
                      4*a h;
                      5*a 0;
                      5*a h;
                      6*a 0;
                      6*a h;
                      7*a 0;
                      7*a h;
                      8*a 0;
                      8*a h;
                      9*a 0;
                      9*a h;
                     10*a 0];
MyTruss.ELEM       = [ 1  3;         % element nodes
                       2  3;
                       3  4;
                       3  5;
                       4  5;
                       5  6;
                       5  7;
                       6  7;
                       7  8;
                       7  9;
                       8  9;
                       9 10;
                       9 11;
                      10 11;
                      10 13;
                      11 13;
                      12 13;
                      12 15;
                      13 15;
                      14 15;
                      14 17;
                      15 17;
                      16 17;
                      16 19;
                      17 19;
                      18 19;
                      19 20;
                       1  2;
                       2  4;
                       4  6;
                       6  8;
                       8 10;
                      10 12;
                      12 14;
                      14 16;
                      16 18;
                      18 20];
MyTruss.SYMBARS    = {[ 1 27];       % symmetrical bars
                      [ 2 26];
                      [ 3 24];
                      [ 4 25];
                      [ 5 23];
                      [ 6 21];
                      [ 7 22];
                      [ 8 20];
                      [ 9 18];
                      [10 19];
                      [11 17];
                      [12 15];
                      [13 16];
                      [14 14];
                       28:37};
MyTruss.Nnodes     = size(MyTruss.NODES,1); % # of nodes
MyTruss.Nelem      = size(MyTruss.ELEM ,1); % # of elements
MyTruss.Ndofs      = 2*MyTruss.Nnodes;      % # of DoFs

% number of design variables of each kind
NvarsA         = size(MyTruss.SYMBARS,1);
NvarsH         = 5;
MyTruss.NvarsA = NvarsA;
MyTruss.NvarsH = NvarsH;

% objective and constraint functions
fun     = @(x)TrussMass37           (x,MyTruss);
nonlcon = @(x)TrussFrequencyBounds37(x,MyTruss);

% lower and upper bounds for design variables
lb_A =  100.0*ones(1,NvarsA)*1e-6;  % unit: m^2
ub_A = 1000.0*ones(1,NvarsA)*1e-6;  % unit: m^2
lb_H =    0.1*ones(1,NvarsH);       % unit: m
ub_H =   10.0*ones(1,NvarsH);       % unit: m
lb   = [lb_A lb_H];
ub   = [ub_A ub_H];

% initial mean and standard deviation
xmean0 = 0.5*(ub+lb);
sigma0 = 5*(ub-lb);

% CE solver
tic
[Xopt,Fopt,ExitFlag,CEstr] = CEopt(fun,[],[],lb,ub,nonlcon);
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
PlotTruss37([xmean0(1:NvarsA)*5e3 xmean0(NvarsA+1:end)],MyTruss,'Non-optimal Truss Structure');

% plot the optimal truss structure
PlotTruss37([  Xopt(1:NvarsA)*5e3   Xopt(NvarsA+1:end)],MyTruss,'Optimal Truss Structure');
