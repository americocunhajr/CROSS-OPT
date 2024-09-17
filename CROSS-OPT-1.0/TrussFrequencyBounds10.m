% -----------------------------------------------------------------
%  TrussFrequencyBounds10.m
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
function [G,H] = TrussFrequencyBounds10(Area,MyTruss)

    % truss structure parameters
    rho        = MyTruss.rho;
    E          = MyTruss.E;
    AddedMass  = MyTruss.AddedMass;
    omegaTresh = MyTruss.omegaTresh;
    FixedDoFs  = MyTruss.FixedDoFs;
    NODES      = MyTruss.NODES;
    ELEM       = MyTruss.ELEM; 
    Nelem      = MyTruss.Nelem;
    Ndofs      = MyTruss.Ndofs;

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