% -----------------------------------------------------------------
%  TrussFrequencyBounds37.m
% -----------------------------------------------------------------
%  programmers: Marcos Vinicius Issa
%               Anderson Pereira
%               Americo Cunha Jr
%
%  Originally programmed in: Apr 04, 2024
%           Last updated in: Sep 14, 2024
% -----------------------------------------------------------------
%  This function computes the frequency bounds for a truss
% -----------------------------------------------------------------
function [G,H] = TrussFrequencyBounds37(x,MyTruss)

    % truss structure parameters
    a          = MyTruss.a;
    rho        = MyTruss.rho;
    E          = MyTruss.E;
    AddedMass  = MyTruss.AddedMass;
    omegaTresh = MyTruss.omegaTresh;
    FixedDoFs  = MyTruss.FixedDoFs;
    ELEM       = MyTruss.ELEM; 
    SYMBARS    = MyTruss.SYMBARS;
    Nelem      = MyTruss.Nelem;
    Ndofs      = MyTruss.Ndofs;

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

    % preallocate memory for constraints
    Nconstr = length(omegaTresh);
    G       = zeros(1,Nconstr);  % inequality constraints
    H       = [];                %   equality constraint

    % preallocate memory for FEM matrices
    M = zeros(Ndofs,Ndofs);
    K = zeros(Ndofs,Ndofs);

    % assembly global matrices
    for e = 1:Nelem
        % distance between nodes
        dx = NODES(ELEM(e,2),1) - NODES(ELEM(e,1),1);
        dy = NODES(ELEM(e,2),2) - NODES(ELEM(e,1),2);
        l  = sqrt(dx^2+dy^2);
        % strain-displacement matrix
        c = dx/l; 
        s = dy/l;
        B = [-c -s c s];
        % element DoFs
        eDof = [2*ELEM(e,1)-1, 2*ELEM(e,1),...
                2*ELEM(e,2)-1, 2*ELEM(e,2)];
        % local stiffness matrix
        Ke = (E*Area(e)/l)*(B'*B);
        % local mass matrix
        Me = (rho*Area(e)*l/6)*[2 0 1 0;...
                                0 2 0 1;...
                                1 0 2 0;...
                                0 1 0 2];
        % update global matrices
        K(eDof,eDof) = K(eDof,eDof) + Ke;
        M(eDof,eDof) = M(eDof,eDof) + Me;
    end

    % update mass matrix with the added mass
    M = M + AddedMass*eye(Ndofs,Ndofs);
    
    % free DoFs coordinates
    FreeDoFs = setdiff(1:Ndofs,FixedDoFs);

    % solve the generalized eigenvalue problem
    omega2 = eig(K(FreeDoFs,FreeDoFs),M(FreeDoFs,FreeDoFs));
    
    % sort frequencies (rad/s)
    omega = sort(sqrt(omega2));
    
    % frequency constraints
    for j = 1:Nconstr
        G(j) = 1 - omega(j)/omegaTresh(j);
    end
end
% -----------------------------------------------------------------