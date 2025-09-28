clc
clear
global Pu  Bi buble H D dr T y
x = [1 2 0]; % x=[Refrigerant1 Refrigerant2 fractionRefrigerant1]
Pu    = xlsread('number','Pure');
Bi    = xlsread('number','Binary');
H     = xlsread('Regression Parameter','h_lv');
D     = xlsread('Regression Parameter','Vapor Density');
dr    = xlsread('Regression Parameter','drhodP');
T     = xlsread('Regression Parameter','Te');
y     = xlsread('Regression Parameter','y');
buble = xlsread('number','buble');
expP  = xlsread('Experimental-Data','SP06-Pow'); %Experimental Data for Power
expT  = xlsread('Experimental-Data','SP06.Temp'); % Experimental Data for Inside Temperature
%% Selection of refrigerant 1
if x(1)== 1,Ref1 = 'R290';  id1 =1; end
if x(1)== 2,Ref1 = 'R600a'; id1 =2; end 
if x(1)== 3,Ref1 = 'R600';  id1 =5; end
if x(1)== 4,Ref1 = 'R1270'; id1 =3; end 
if x(1)== 5,Ref1 = 'RC270'; id1 =4; end 

% Selection of refrigerant 2
if x(2)== 1,Ref2 = 'R290';  id2 =1; end
if x(2)== 2,Ref2 = 'R600a'; id2 =2; end
if x(2)== 3,Ref2 = 'R600';  id2 =5; end
if x(2)== 4,Ref2 = 'R1270'; id2 =3; end
if x(2)== 5,Ref2 = 'RC270'; id2 =4; end

%%
z(1) = x(3); 
z(2)=1-z(1);
[z,MW,w,Tc,Pc,Co,f,a,b,k] = Myinput(Ref1,Ref2,z(1),z(2));

%%
t_sim = 3; % Simulation time [hr]
dt    = .3; % Time Step [sec]
Nt    = ceil(t_sim*60*60/dt);

%%
T_goods    = zeros(1,Nt); % Goods Temperature ['c]
T_wall     = zeros(1,Nt); % Evaporator Wall Temperature ['c]
T_air      = zeros(1,Nt); % Cabinet Air Temperature ['c]
Psuc       = zeros(1,Nt); % Suction Pressure [bar]
M_ref      = zeros(1,Nt); % Refrigerant Liquid Mass in Evaporator [Kg]
t          = zeros(1,Nt); % Simulation Time [sec]
u_comp     = zeros(1,Nt); % Compressor ON/OFF [-]
Power      = zeros(1,Nt); % Compressor Input Power [w]
T_sh0      = zeros(1,Nt); % Compressor Input Temperature [K]
F_comp     = zeros(1,Nt); % Compressor Output Volumetric Flow Rate [m^3/s]
Te         = zeros(1,Nt); % Evaporator Saturated Temperature ['c]
m_r11      = zeros(1,Nt+1);     % Display_Case_1 & M_ref_1 & m_ref_1    [kg]
m_in       = zeros(1,Nt);
%%
T_max  = 9; % Maximum Temperature ['c]
T_min  = 4; % Minimum Temperature ['c]
T_amb  = 22.8;  % Ambiant Temperature ['c]


UA_goods_air = 0;    % Overall Heat Transfer Coefficent of Goods to Air [W/k]
UA_air_wall  = 3; % Overall Heat Transfer Coefficent of Air to Wall  [W/k]
UA_wall_ref  = 20;   % Overall Heat Transfer Coefficent of Wall to Refrigerant [W/k]
UA_amb_cab   = 2.94; % Overall Heat Transfer Coefficent of Ambiant to Refrigerator [W/k]


C_goods = 211200; % Specific Heat times Mass of Goods [J/k]
C_air   = 16*1083;  % Specific Heat times Mass of Air   [J/k]
C_wall  = 921.096*.5;  % Specific Heat times Mass of Wall  [J/k]

MP = buble(1+(id2-1)*12:12*id2,id1); % buble pressure  and liquid density paramiter needs for calculation in condenser
C_v = MP(7,1)*z(1)+MP(8,1);
M_ref_max = 0.053;          % Maximum Mass of Refrigerant
Mr_max    = M_ref_max*0.2; % Maximum Mass of Refrigerant in evaporator
Eta_vol   = 0.82;         % Volumetric Efficiency of Compressor
Eta_is    = 1;          % Isentropic Efficiency of Compressor
V_sl      = 10.2;     % displacement [cm^3/rev]
V_suc     = 7.4925e-4;      % suctin volume [m^3]
V_cond    = 8e-4; %condesser volume [m^3]
RPM       = 1800; % rotation per min
M_ref(1,1)   = .012;    %[kg/s]
m_r11(1,1)   = .012;
T_air(1,1)   = 4.8;  %['c]
T_goods(1,1) = 0;    %['c]
T_wall(1,1)  = -14;    %['c]
Psuc(1,1)    = .25;      %[bar]
u_comp(1,1)  = 1;
T_sh0(1)     = 273; 
%%
PI_swc = 0;
P_con_sat = (MP(1,1)*z(1)^2+MP(2,1)*z(1)+MP(3,1))/100000; % buble pressure of condenser
rho_in    = (MP(4,1)*z(1)^2+MP(5,1)*z(1)+MP(6,1)); % outlet density of condenser
%% Simulation
    for n=2:Nt
        t(n) = t(n-1)+dt;
        if M_ref(n-1)<0,  M_ref(n-1)    = 0; end
        if m_r11(n-1)<0,  m_r11(n-1)    = 0; end

        % Evaporators
        % Evaporators: Calculate mass fractions in liquid phase
        x_r = zeros(2,2);
        if M_ref(n-1)>0 && m_r11(n-1)>0
            x_r(1,1) = m_r11(n-1)/M_ref(n-1);
            if x_r(1,1)>1,x_r(1,1)=1;end
            x_r(1,2) = 1-x_r(1,1);
        else 
            x_r(1,1) = z(1);
            x_r(1,2) = z(2);
        end
        


        % Convert mass fractions to mole fractions in liquid phase of evaporators
        x_rm = zeros(2,2);
        x_rm(1,1) = (x_r(1,1)*MW(2))/(MW(1)-x_r(1,1)*MW(1)+x_r(1,1)*MW(2)); % Display_Case_1 & M_ref_1 & m_ref_1
        if x_rm(1,1)>1
             x_rm(1,1)=1;
        end
        if x_rm(1,1)<0
            x_rm(1,1)=0;
        end
        x_rm(1,2) = 1-x_rm(1,1);                                            % Display_Case_1 & M_ref_1 & m_ref_2

        % Calculate mole fractions in gas phase of evaporators
        [~,~,~,~,y_rm] = ThermodynamicsProperties(id1,id2,x_rm(1,1),Psuc(n-1));  % Display_Case_1 & M_ref_1 & m_ref_1     
        y_rm(1,1) = y_rm;
        if y_rm(1,1)>1,y_rm(1,1)=1;end
        if y_rm(1,1)<0,y_rm(1,1)=0;end
        y_rm(1,2) = 1-y_rm(1,1);                                              % Display_Case_1 & M_ref_1 & m_ref_2

        % Convert mole fractions to mass fractions in gas phase of evaporators
        y_r = zeros(2,2);
        y_r(1,1) = (y_rm(1,1)*MW(1))/(y_rm(1,1)*MW(1)+y_rm(1,2)*MW(2));        % Display_Case_1 & M_ref_1 & m_ref_1
        if y_r(1,1)>1, y_r(1,1)=1; end
        y_r(1,2) = 1-y_r(1,1) ;                                                % Display_Case_1 & M_ref_1 & m_ref_2
        [Te(n),~,~,h_lv,~]=ThermodynamicsProperties(id1,id2,y_rm(1,1),Psuc(n-1));
        % Cabinet Model
        Q_goods   = UA_goods_air*(T_goods(n-1)-T_air(n-1));
        Q_wall    = UA_air_wall*(T_air(n-1)-T_wall(n-1));
        Qe        = UA_wall_ref*(M_ref(n-1)/Mr_max)*(T_wall(n-1)-Te(n));
        Q_airload = UA_amb_cab*(T_amb-T_air(n-1));

        dT_goods  = -(Q_goods)/C_goods;
        dT_air    = (Q_goods+Q_airload-Q_wall)/C_air;
        dT_wall   = (Q_wall-Qe)/C_wall;

        T_goods(n) = T_goods(n-1)+dT_goods*dt;
        T_air(n)   = T_air(n-1)+dT_air*dt;
        T_wall(n)  = T_wall(n-1)+dT_wall*dt;
        % Compressor ON/OFF
        if T_air(n-1)>T_max
             u_comp(n) =1;
        end
        if T_air(n-1)<T_min
            u_comp(n) = 0.04;
        end
        if T_air(n-1)>=T_min && T_air(n-1)<=T_max
            u_comp(n) = u_comp(n-1);
        end
        % Compressor Model

        f_in      = (Qe/h_lv);          % Total mass flow inter suction manifold [kg/sec]
        f_in_suc1 = y_r(1,1)*(Qe/h_lv); % Outlet mass flow of component 1 from evaporators [kg/sec]
        f_in_suc2 = f_in - f_in_suc1;   % Outlet mass flow of component 2 from evaporators [kg/sec]

        T_sh0(n)  = Te(n)+10+273.15; % Superheat temperature of the suction manifold;

        [~,rho_suc,drho_suc,~,~] = ThermodynamicsProperties(id1,id2,z(1),Psuc(n-1));
        F_comp(n)   = (u_comp(n)*V_sl*RPM/1000)/60000;
        dP    = (f_in-F_comp(n)*rho_suc)/V_suc/drho_suc;
        Psuc(n)=Psuc(n-1)+dP*dt;

        % Regression model for Power calculation
        f_suc1    = f(1,1)*Psuc(n-1) + f(1,2); 
        f_suc2    = f(2,1)*Psuc(n-1) + f(2,2); 
        f_suc     = z(1)*f_suc1 + z(2)*f_suc2;  % f = rho_suc*(h_oc-h_ic)
        Power(n)  = F_comp(n)*f_suc/Eta_is; % [W]
 
         % Evaporator Model
        m_in(n)   = F_comp(n)*rho_suc;
        dM     = (m_in(n))-(f_in);
        dm_r11 = z(1)*dM;
        M_ref(n) = M_ref(n-1)+ dM*dt;
        m_r11(n) = m_r11(n-1)+ dm_r11*dt;
        if M_ref(n)<0,M_ref(n) = 0;end
        if M_ref(n)> Mr_max, M_ref(n)= Mr_max;end
        if m_r11(n)> z(1)*Mr_max, m_r11(n)= z(1)*Mr_max;end
        if m_r11(n)< 0, m_r11(n)= 0;end
        if u_comp(n)~=u_comp(n-1), PI_swc = PI_swc+1;end
    end
    
    %%
    set(0,'DefaultAxesColorOrder',[0 0 1; 1 0 0],'DefaultAxesLineStyleOrder','-|--|:')
set(gcf, 'PaperSize', [4 12]);
figure(1)
plot(t/60,Power,'linestyle','-','linew',1.5)
hold on
plot(expP(:,8)/60,expP(:,9),'linestyle','--','linew',1.5)
% title('Power Consumption with R290:R600a')
xlabel('Time [s]')
ylabel('Power [W]'), box on, set(figure(1),'color','w')
legend('Power_{sim}','Power_{exp}','Location','SE');
% legend('boxoff')
% axis([0 inf 0 65])
pbaspect([1 .4 1])
set(gca,'fontname','Times New Roman','FontSize',10)
figure(2)
plot(t/60,T_air,'linestyle','-','linew',1.5)
hold on
plot(expT(:,13)/60,expT(:,14),'linestyle','--','linew',1.5)
hold on
plot([0 500 t(end)],[T_min T_min T_min],'color',[0 0 0],'linestyle',':','linew',1)
hold on
plot([0 500 t(end)],[T_max T_max T_max],'color',[0 0 0],'linestyle',':','linew',1)
% title('Air Temperature with R290:R600a')
xlabel('Time [min]')
ylabel('Temperature [^oc]'), box on, set(figure(2),'color','w')
legend('Temp_{sim}','Temp_{exp}','T_{min,max}','Location','SE');
axis([0 t(end)/60 -2 10])
pbaspect([1 .4 1])
set(gca,'fontname','Times New Roman','FontSize',10)