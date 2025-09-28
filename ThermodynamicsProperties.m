function [Te,vapor_rho,drho,h_lv,ye]=ThermodynamicsProperties(id1,id2,x,P)
%% Input
% id1 = Refrigerant 1
% id2 = Refrigerant 2
% x   = Mole Fraction of Refrigerant 1
% P   = Suction Pressure [bar]
%% Output
% Te        = Bubble Temperature ['c]
% vapor_rho = Vapor Density [kg/m^3]
% drho      = drho/dP [kg/m^3/bar]
% h_lv      = Enthalpy of Vaporization [J/kg]
% ye        = Mole Fraction of Vapor in Equilibrium with Liquid 
%%
% NOTE : P = [0.5 3] bar
% Only Use for suction Properties
%%
global H D dr T y
h  = H (1+(id2-1)*9:9*id2,id1); % Regression Parameter for h_lv
r  = D (1+(id2-1)*9:9*id2,id1); % Regression Parameter for Vapor density
d  = dr(1+(id2-1)*9:9*id2,id1); % Regression Parameter for drho/dP
t  = T (1+(id2-1)*9:9*id2,id1); % Regression Parameter for Bubble Temperature
Y  = y (1+(id2-1)*9:9*id2,id1); % Regression Parameter for vapor fraction
%%
Te        = t(1) + t(2)*x + t(3)*P + t(4)*x^2 +t(5)*x*P + t(6)*P^2 + t(7)*x^2*P +t(8)*x*P^2 + t(9)*P^3;
vapor_rho = r(1) + r(2)*x + r(3)*P + r(4)*x^2 +r(5)*x*P + r(6)*P^2 + r(7)*x^2*P +r(8)*x*P^2 + r(9)*P^3;
drho      = d(1) + d(2)*x + d(3)*P + d(4)*x^2 +d(5)*x*P + d(6)*P^2 + d(7)*x^2*P +d(8)*x*P^2 + d(9)*P^3;
h_lv      = (h(1) + h(2)*x + h(3)*P + h(4)*x^2 +h(5)*x*P + h(6)*P^2 + h(7)*x^2*P +h(8)*x*P^2 + h(9)*P^3)*1000;
ye        = Y(1) + Y(2)*x + Y(3)*P + Y(4)*x^2 +Y(5)*x*P + Y(6)*P^2 + Y(7)*x^2*P +Y(8)*x*P^2 + Y(9)*P^3;
end