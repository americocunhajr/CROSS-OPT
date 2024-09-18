% -----------------------------------------------------------------
%  Main_Truss37_MassMin_FrequencyBounds_BoxPlot.m
% -----------------------------------------------------------------
%  programmers: Marcos Vinicius Issa
%               Anderson Pereira
%               Americo Cunha Jr
%
%  Originally programmed in: Jun 25, 2024
%           Last updated in: Sep 15, 2024
% -----------------------------------------------------------------
% ﻿ Mass minimization with frequency constraints for a 37 bars
%  truss structure. Box plot calculation for CE, GA and SQP.
% -----------------------------------------------------------------


clc; clear; close all;

disp(' ---------------------------------------------- ')
disp(' Main_Truss37_MassMin_FrequencyBounds_BoxPlot.m ')
disp(' ---------------------------------------------- ')

% random number generator (fix the seed for reproducibility)
rng_stream = RandStream('mt19937ar','Seed',30081984);
RandStream.setGlobalStream(rng_stream);

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
Nvars          = NvarsA + NvarsH;

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

% number of box plot samples
Nbp = 64;

% boxplot for CE
% ------------------------------------------------------------------
tic
disp(' ')
disp(' --------------- ')
disp(' CE Box Plot ... ')
disp(' --------------- ')
disp(' ')

        t = 0;
  Nrun_CE = 0;
  Xopt_CE = zeros(Nbp,Nvars);
  Fopt_CE = zeros(Nbp,1);
Fcount_CE = zeros(Nbp,1);

% optional parameters
CEstr.Verbose = 0;

while t < Nbp
        % CE solver
        [Xopt,Fopt,ExitFlag,CEstr] = CEopt(fun,[],[],lb,ub,nonlcon,CEstr);

        % update solver running counter
        Nrun_CE = Nrun_CE + 1;

        % check constraint violation
        if max(nonlcon(Xopt)) <= CEstr.TolCon
                         t = t + 1;
              Xopt_CE(t,:) = Xopt;
              Fopt_CE(t,:) = Fopt;
            Fcount_CE(t,:) = CEstr.Fcount;
        end
end
Xopt_CE_m = median(Xopt_CE)
toc
% ------------------------------------------------------------------


% boxplot for GA
% ------------------------------------------------------------------
tic
disp(' ')
disp(' --------------- ')
disp(' GA Box Plot ... ')
disp(' --------------- ')
disp(' ')

        t = 0;
  Nrun_GA = 0;
  Xopt_GA = zeros(Nbp,Nvars);
  Fopt_GA = zeros(Nbp,1);
Fcount_GA = zeros(Nbp,1);

% optional parameters
opt = optimoptions('ga','Display','off');

while t < Nbp
        % GA solver
        [Xopt,Fopt,ExitFlag,GAstr] = ga(fun,Nvars,[],[],[],[],lb,ub,nonlcon,opt);

        % update solver running counter
        Nrun_GA = Nrun_GA + 1;

        % check constraint violation
        if max(nonlcon(Xopt)) <= opt.TolCon
                         t = t + 1;
              Xopt_GA(t,:) = Xopt;
              Fopt_GA(t,:) = Fopt;
            Fcount_GA(t,:) = GAstr.funccount;
        end
end
Xopt_GA_m = median(Xopt_GA)
toc
% ------------------------------------------------------------------


% boxplot for SQP
% ------------------------------------------------------------------
tic
disp(' ')
disp(' ---------------- ')
disp(' SQP Box Plot ... ')
disp(' ---------------- ')
disp(' ')

         t = 0;
  Nrun_SQP = 0;
  Xopt_SQP = zeros(Nbp,Nvars);
  Fopt_SQP = zeros(Nbp,1);
Fcount_SQP = zeros(Nbp,1);

% optional parameters
opt = optimoptions('fmincon','Display','off','Algorithm','sqp');

while t < Nbp
        % random initial guess
        x0 = lb + (ub-lb).*rand(1,Nvars);

        % SPQ solver
        [Xopt,Fopt,ExitFlag,SQPstr] = fmincon(fun,x0,[],[],[],[],lb,ub,nonlcon,opt);

        % update solver running counter
        Nrun_SQP = Nrun_SQP + 1;

        % check constraint violation
        if max(nonlcon(Xopt)) <= opt.TolCon
                          t = t + 1;
              Xopt_SQP(t,:) = Xopt;
              Fopt_SQP(t,:) = Fopt;
            Fcount_SQP(t,:) = SQPstr.funcCount;
        end
end
Xopt_SPQ_m = median(Xopt_SQP)
toc
% ------------------------------------------------------------------


% save simulation results
% -----------------------------------------------------------
tic
disp(' ')
disp(' --- saving simulation results --- ');
disp(' ');
disp('    ... ');
disp(' ');

save('Truss37_MassMin_FrequencyBounds_BoxPlot.mat');

toc
% -----------------------------------------------------------


% post processing
% ------------------------------------------------------------------
disp(' ')
disp(' ------------------- ')
disp(' Post processing ... ')
disp(' ------------------- ')
disp(' ')

% ..........................................................
graphobj1.gname     = 'Truss37_MassMin_FrequencyBounds_BoxPlot_Fopt';
graphobj1.gtitle    = 'Objetive Function';
graphobj1.xmin      = 'auto';
graphobj1.xmax      = 'auto';
graphobj1.ymin      = 'auto';
graphobj1.ymax      = 'auto';
graphobj1.xlab      = [];
graphobj1.ylab      = 'Truss Mass (kg)';
graphobj1.leg1      = 'CE';
graphobj1.leg2      = 'GA';
graphobj1.leg3      = 'SQP';
graphobj1.print     = 'yes';
graphobj1.close     = 'no';
Fig1 = PlotBoxComparison(Fopt_CE,Fopt_GA,Fopt_SQP,graphobj1);
% ..........................................................

% ..........................................................
graphobj2.gname     = 'Truss37_MassMin_FrequencyBounds_BoxPlot_Fcount';
graphobj2.gtitle    = 'Computational Cost';
graphobj2.xmin      = 'auto';
graphobj2.xmax      = 'auto';
graphobj2.ymin      = 'auto';
graphobj2.ymax      = 'auto';
graphobj2.xlab      = [];
graphobj2.ylab      = 'Function Evaluations';
graphobj2.leg1      = 'CE';
graphobj2.leg2      = 'GA';
graphobj2.leg3      = 'SQP';
graphobj2.print     = 'yes';
graphobj2.close     = 'no';
Fig2 = PlotBoxComparison(Fcount_CE,Fcount_GA,Fcount_SQP,graphobj2);
% ..........................................................
% ------------------------------------------------------------------
