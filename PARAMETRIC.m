
clc
clear
global Pu Bi buble H D dr T y

% Load data
Pu = xlsread('number.xlsx', 'Pure');
Bi = xlsread('number.xlsx', 'Binary');
buble = xlsread('number.xlsx', 'buble');
H = xlsread('Regression Parameter.xlsx', 'h_lv');
D = xlsread('Regression Parameter.xlsx', 'Vapor Density');
dr = xlsread('Regression Parameter.xlsx', 'drhodP');
T = xlsread('Regression Parameter.xlsx', 'Te');
y = xlsread('Regression Parameter.xlsx', 'y');

% List of refrigerants
refs = {'R290', 'R600a', 'R600', 'R1270', 'RC270'};

% Initialize results
results = struct('Ref1', {}, 'Ref2', {}, 'z1', {}, 'totalEnergy', {}, 'stdTair', {}, 'meanPsuc', {});

% Plot data
plotData = containers.Map();

idx = 1;

RPM = 3000;  % Typical compressor speed

for i = 1:length(refs)-1
    for j = i+1:length(refs)
        Ref1 = refs{i};
        Ref2 = refs{j};
        
        switch Ref1
            case 'R290', id1 = 1;
            case 'R600a', id1 = 2;
            case 'R600', id1 = 3;
            case 'R1270', id1 = 4;
            case 'RC270', id1 = 5;
        end
        switch Ref2
            case 'R290', id2 = 1;
            case 'R600a', id2 = 2;
            case 'R600', id2 = 3;
            case 'R1270', id2 = 4;
            case 'RC270', id2 = 5;
        end
        
        z_vals = 0:0.1:1;
        
        sorted = sort({Ref1, Ref2});
        pairKey = sprintf('%s-%s', sorted{:});
        
        if ~isKey(plotData, pairKey)
            plotData(pairKey) = struct('z1', [], 'totalEnergy', [], 'stdTair', [], 'meanPsuc', []);
        end
        
        for zz = 1:length(z_vals)
            z1 = z_vals(zz);
            z2 = 1 - z1;
            
            [z, MW, w, Tc, Pc, Co, f, a, b, k] = Myinput(Ref1, Ref2, z1, z2);
            
            MP = buble(1+(id2-1)*12:12*id2,id1);
            C_v = MP(7,1)*z1 + MP(8,1);  % Not used, but included
            
            t_sim = 3;
            dt = 0.3;
            Nt = ceil(t_sim*60*60/dt);
            
            T_goods = zeros(1,Nt);
            T_wall = zeros(1,Nt);
            T_air = zeros(1,Nt);
            Psuc = zeros(1,Nt);
            M_ref = zeros(1,Nt);
            t = zeros(1,Nt);
            u_comp = zeros(1,Nt);
            Power = zeros(1,Nt);
            T_sh0 = zeros(1,Nt);
            F_comp = zeros(1,Nt);
            Te = zeros(1,Nt);
            m_r11 = zeros(1,Nt+1);
            m_in = zeros(1,Nt);
            
            T_max = 9;
            T_min = 4;
            T_amb = 22.8;
            
            UA_goods_air = 0;
            UA_air_wall = 3;
            UA_wall_ref = 20;
            UA_amb_cab = 2.94;
            
            C_goods = 211200;
            C_air = 16*1083;
            C_wall = 921.096*0.5;
            
            M_ref_max = 0.053;
            Mr_max = M_ref_max*0.2;
            Eta_vol = 0.82;
            Eta_is = 1;
            V_sl = 10.2;  % cm^3/rev
            V_suc = 7.4925e-4;
            
            M_ref(1) = 0.012;
            T_air(1) = 6.5;
            T_goods(1) = 6.5;
            T_wall(1) = 4;
            Psuc(1) = 1;
            u_comp(1) = 1;
            m_r11(1) = z(1)*M_ref(1);
            
            PI_swc = 0;
            
            for n = 2:Nt
                t(n) = t(n-1) + dt;
                if M_ref(n-1) < 0, M_ref(n-1) = 0; end
                if m_r11(n-1) < 0, m_r11(n-1) = 0; end
                
                x_r = [0 0];
                if M_ref(n-1) > 0 && m_r11(n-1) > 0
                    x_r(1) = m_r11(n-1) / M_ref(n-1);
                    x_r(2) = 1 - x_r(1);
                else
                    x_r = z;
                end
                
                [Te(n), ~, ~, h_lv, ye] = ThermodynamicsProperties(id1, id2, x_r(1), Psuc(n-1));
                
                y_rm = [ye, 1 - ye];
                y_r = [0 0];
                y_r(1) = (y_rm(1) * MW(1)) / (y_rm(1)*MW(1) + y_rm(2)*MW(2));
                y_r(2) = 1 - y_r(1);
                
                Q_goods = UA_goods_air * (T_goods(n-1) - T_air(n-1));
                Q_wall = UA_air_wall * (T_air(n-1) - T_wall(n-1));
                Qe = UA_wall_ref * (M_ref(n-1)/Mr_max) * (T_wall(n-1) - Te(n));
                Q_airload = UA_amb_cab * (T_amb - T_air(n-1));
                
                dT_goods = -Q_goods / C_goods;
                dT_air = (Q_goods + Q_airload - Q_wall) / C_air;
                dT_wall = (Q_wall - Qe) / C_wall;
                
                T_goods(n) = T_goods(n-1) + dT_goods * dt;
                T_air(n) = T_air(n-1) + dT_air * dt;
                T_wall(n) = T_wall(n-1) + dT_wall * dt;
                
                if T_air(n-1) > T_max
                    u_comp(n) = 1;
                end
                if T_air(n-1) < T_min
                    u_comp(n) = 0.04;
                end
                if T_air(n-1) >= T_min && T_air(n-1) <= T_max
                    u_comp(n) = u_comp(n-1);
                end
                
                f_in = Qe / h_lv;
                f_in_suc1 = y_r(1) * f_in;
                f_in_suc2 = f_in - f_in_suc1;
                
                T_sh0(n) = Te(n) + 10 + 273.15;
                
                [~, rho_suc, drho_suc, ~, ~] = ThermodynamicsProperties(id1, id2, z(1), Psuc(n-1));
                
                F_comp(n) = (u_comp(n) * V_sl * RPM / 1000) / 60000;
                
                dP = (f_in - F_comp(n) * rho_suc) / V_suc / drho_suc;
                Psuc(n) = Psuc(n-1) + dP * dt;
                
                f_suc1 = f(1,1) * Psuc(n-1) + f(1,2);
                f_suc2 = f(2,1) * Psuc(n-1) + f(2,2);
                f_suc = z(1) * f_suc1 + z(2) * f_suc2;
                Power(n) = F_comp(n) * f_suc / Eta_is;
                
                m_in(n) = F_comp(n) * rho_suc;
                dM = m_in(n) - f_in;
                dm_r11 = z(1) * dM;
                M_ref(n) = M_ref(n-1) + dM * dt;
                m_r11(n) = m_r11(n-1) + dm_r11 * dt;
                if M_ref(n) < 0, M_ref(n) = 0; end
                if M_ref(n) > Mr_max, M_ref(n) = Mr_max; end
                if m_r11(n) > z(1) * Mr_max, m_r11(n) = z(1) * Mr_max; end
                if u_comp(n) ~= u_comp(n-1), PI_swc = PI_swc + 1; end
            end
            
            totalEnergy = sum(Power) * dt;  % Joules
            
            meanPsuc = mean(Psuc);
            
            results(idx).Ref1 = Ref1;
            results(idx).Ref2 = Ref2;
            results(idx).z1 = z1;
            results(idx).totalEnergy = totalEnergy;
            results(idx).stdTair = std(T_air);
            results(idx).meanPsuc = meanPsuc;
            
            pd = plotData(pairKey);
            if strcmp(Ref1, sorted{1})
                this_z = z1;
            else
                this_z = 1 - z1;
            end
            pd.z1 = [pd.z1 this_z];
            pd.totalEnergy = [pd.totalEnergy totalEnergy];
            pd.stdTair = [pd.stdTair std(T_air)];
            pd.meanPsuc = [pd.meanPsuc meanPsuc];
            plotData(pairKey) = pd;
            
            idx = idx + 1;
        end
    end
end

% Table
resultsTable = struct2table(results);
disp(resultsTable);
writetable(resultsTable, 'refrigerant_parametric_results.csv');
%%
% Plots
keys = plotData.keys;
for k = 1:length(keys)
    key = keys{k};
    pd = plotData(key);
    
    [pd.z1, sortIdx] = sort(pd.z1);
    pd.totalEnergy = pd.totalEnergy(sortIdx);
    pd.stdTair = pd.stdTair(sortIdx);
    pd.meanPsuc = pd.meanPsuc(sortIdx);
    
    [pd.z1, uniqueIdx] = unique(pd.z1);
    pd.totalEnergy = pd.totalEnergy(uniqueIdx);
    pd.stdTair = pd.stdTair(uniqueIdx);
    pd.meanPsuc = pd.meanPsuc(uniqueIdx);
    
    sorted = strsplit(key, '-');
    title_base = sprintf('%s (z1) and %s (1-z1)', sorted{1}, sorted{2});
    
    figure;
    subplot(3,1,1);
    plot(pd.z1, pd.totalEnergy / 3600 / 1000, 'o-');  % to kWh
    xlabel('Mass Fraction z1');
    ylabel('Total Energy Consumption over 3 hours (kWh)');
    title(['Energy Consumption vs Mass Fraction for ' title_base]);
    grid on;
    
    subplot(3,1,2);
    plot(pd.z1, pd.stdTair, 'o-');
    xlabel('Mass Fraction z1');
    ylabel('Std Deviation of Air Temperature (°C)');
    title(['Temperature Deviation vs Mass Fraction for ' title_base]);
    grid on;
    
    subplot(3,1,3);
    plot(pd.z1, pd.meanPsuc, 'o-');
    xlabel('Mass Fraction z1');
    ylabel('Mean Suction Pressure (bar)');
    title(['Suction Pressure vs Mass Fraction for ' title_base]);
    grid on;
end


