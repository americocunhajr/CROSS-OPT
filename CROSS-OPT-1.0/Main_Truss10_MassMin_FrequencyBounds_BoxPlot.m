% -----------------------------------------------------------------
%  Main_Truss10_MassMin_FrequencyBounds_BoxPlot.m
% -----------------------------------------------------------------
%  programmers: Marcos Vinicius Issa
%               Anderson Pereira
%               Americo Cunha Jr
%
%  Originally programmed in: Jun 25, 2024
%           Last updated in: Sep 15, 2024
% -----------------------------------------------------------------
% ﻿ Mass minimization with frequency constraints for a 10 bars
%  truss structure. Box plot calculation for CE, GA and SQP.
% -----------------------------------------------------------------


clc; clear; close all;

disp(' ---------------------------------------------- ')
disp(' Main_Truss10_MassMin_FrequencyBounds_BoxPlot.m ')
disp(' ---------------------------------------------- ')

% random number generator (fix the seed for reproducibility)
rng_stream = RandStream('mt19937ar','Seed',30081984);
RandStream.setGlobalStream(rng_stream);

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

save('Truss10_MassMin_FrequencyBounds_BoxPlot.mat');

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
graphobj1.gname     = 'Truss10_MassMin_FrequencyBounds_BoxPlot_Fopt';
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
graphobj2.gname     = 'Truss10_MassMin_FrequencyBounds_BoxPlot_Fcount';
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
