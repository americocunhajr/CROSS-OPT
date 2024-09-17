% -----------------------------------------------------------------
%  TrussMass10.m
% -----------------------------------------------------------------
%  programmers: Marcos Vinicius Issa
%               Anderson Pereira
%               Americo Cunha Jr
%
%  Originally programmed in: Apr 04, 2024
%           Last updated in: Sep 14, 2024
% -----------------------------------------------------------------
%  This function computes the mass of a 37 bars truss structure.
% -----------------------------------------------------------------
function Mass = TrussMass37(x,MyTruss)

    % truss structure parameters
    a         = MyTruss.a;
    rho       = MyTruss.rho;
    ELEM      = MyTruss.ELEM;
    Nelem     = MyTruss.Nelem;
    SYMBARS   = MyTruss.SYMBARS;

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
        
%     A([1 27]) = x(6);
%     A([2 26]) = x(7);
%     A(26) = x(7);
%     A(3) = x(8);
%     A(24) = x(8);
%     A(4) = x(9);
%     A(25) = x(9);
%     A(5) = x(10);
%     A(23) = x(10);
%     A(6) = x(11);
%     A(21) = x(11);
%     A(7) = x(12);
%     A(22) = x(12);
%     A(8,:) = x(13);
%     A(20,:) = x(13);
%     A(9,:) = x(14);
%     A(18,:) = x(14);
%     A(10,:) = x(15);
%     A(19,:) = x(15);
%     A(11,:) = x(16);
%     A(17,:) = x(16);
%     A(12,:) = x(17);
%     A(15,:) = x(17);
%     A(13,:) = x(18);
%     A(16,:) = x(18);
%     A(14,:) = x(19);
% for i=28:37
%     A(i,:) = 4e-3;
% end

    % compute the mass
    Mass = 0.0;
    for e = 1:Nelem
        dx = NODES(ELEM(e,2),1) - NODES(ELEM(e,1),1);
        dy = NODES(ELEM(e,2),2) - NODES(ELEM(e,1),2);
        l  = sqrt(dx^2+dy^2);
        Mass  = Mass + rho*Area(e)*l;
    end
end
% -----------------------------------------------------------------