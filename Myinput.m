function [z,MW,w,Tc,Pc,Co,f,a,b,k] = Myinput(Ref1,Ref2,z1,z2)

global Pu Bi

switch Ref1
    case {'R290'},   id1 = 1;
    case {'R600a'},  id1 = 2;
    case {'R600'},   id1 = 3;
    case {'R1270'},  id1 = 4;
    case {'RC270'},  id1 = 5;

end

switch Ref2
    case {'R290'},   id2 = 1;
    case {'R600a'},  id2 = 2;
    case {'R600'},   id2 = 3;
    case {'R1270'},  id2 = 4;
    case {'RC270'},  id2 = 5;

end

MP = Bi(1+(id2-1)*12:12*id2,id1);

% Properties of mixture
z  = [z1 z2]; % Refrigerant Mass Fraction
MW = [Pu(id1,3) Pu(id2,3)]; % Molecular Mass [Kg/Kgmole]
w  = [Pu(id1,4) Pu(id2,4)]; % Acentric Factor
Tc = [Pu(id1,5) Pu(id2,5)]; % Critical Temperature [K]
Pc = [Pu(id1,6) Pu(id2,6)]; % Critical Pressure [bar]

Co = [Pu(id1,7) Pu(id1,8) Pu(id1,9);Pu(id2,7) Pu(id2,8) Pu(id2,9)]; % Antoine Coefficients
f  = [Pu(id1,10) Pu(id1,11);Pu(id2,10) Pu(id2,11)]; % f = rho_suc*(h_oc-h_ic)
a  = MP(1:4,1)'; % Regression Parameters for T_e
b  = MP(5:8,1)'; % Regression Parameters for y_rm
k  = [MP(9:10,1)';MP(11:12,1)']; % Interaction Parameter
end




