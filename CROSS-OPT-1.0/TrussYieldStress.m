% -----------------------------------------------------------------
%  TrussYieldStress.m
% -----------------------------------------------------------------
%  programmers: Marcos Vinicius Issa
%               Anderson Pereira
%               Americo Cunha Jr
%
%  Originally programmed in: Apr 04, 2024
%           Last updated in: Aug 29, 2024
% -----------------------------------------------------------------
%  This function computes the yield stress constraints for a truss
% -----------------------------------------------------------------
function [G,H] = TrussYieldStress(A,MyTruss)

    % truss structure parameters
    E         = MyTruss.E;
    P         = MyTruss.P;
    SY1       = MyTruss.SY1;
    SY2       = MyTruss.SY2;
    WeakBar   = MyTruss.WeakBar;
    FixedDoFs = MyTruss.FixedDoFs;
    LoadDoFs  = MyTruss.LoadDoFs;
    NODES     = MyTruss.NODES;
    ELEM      = MyTruss.ELEM; 
    Nelem     = MyTruss.Nelem;
    Ndofs     = MyTruss.Ndofs;

    % preallocate memory for constraints
    G = zeros(1,Nelem);  % inequality constraints
    H = [];              %   equality constraint

    % preallocate memory for FEM objects
    K = zeros(Ndofs,Ndofs);
    u = zeros(Ndofs,1);
    f = zeros(Ndofs,1);
    
    % assembly the stiffness matrix
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
        Ke = (E*A(e)/l)*(B'*B);
        % update global stiffness matrix
        K(eDof,eDof) = K(eDof,eDof) + Ke;
    end

    % apply concentrated loads
    f(LoadDoFs) = - P;

    % free DoFs coordinates
    FreeDoFs = setdiff(1:Ndofs,FixedDoFs);
    
    % compute displacement vector
    u(FreeDoFs) = K(FreeDoFs,FreeDoFs)\f(FreeDoFs);
    
    % normal stress vector
    sigma = zeros(Nelem,1);
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
        % normal stress vector
        sigma(e) = (E/l)*B*u(eDof);
     end
    
    % yield stress inequality constraints
    for e = 1:Nelem
        SY = SY1;
        if e == WeakBar
            SY = SY2;
        end
        G(e) = abs(sigma(e))/SY - 1;
    end
end
% -----------------------------------------------------------------