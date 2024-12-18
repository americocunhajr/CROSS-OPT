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
%  This function computes the mass of a 10 bars truss structure.
% -----------------------------------------------------------------
function Mass = TrussMass10(Area,MyTruss)

    % truss structure parameters
    rho       = MyTruss.rho;
    NODES     = MyTruss.NODES;
    ELEM      = MyTruss.ELEM;
    Nelem     = MyTruss.Nelem;
    
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