%{
    1  : example 1, p775, V_REF by resistor-MOSFET
%}
clc;

selector = 1;

%---------------------------------------------------------------------------------------
if (selector == 1)
    %a) a 1MEG resistor on a 10/2 NMOS
    resistor = 1e6; % ohms
    tc = 2e3; % temperature coefficient in ppm / C
    KPn = 120; % uA/V^2
    V_DD = 5;
    m_wid = 10;
    m_len = 2;
    
    % you have to know previous chapters to solve
    V_REF = 900e-3; % V_REF about =  V_GS
    I_D = (V_DD - V_REF) / resistor; % both must be 4.1000e-06  "4.1 uA"
    V_thresh = V_REF-1.17e-4;
    I_D = (KPn/2)*(m_wid/m_len)*(V_REF-V_thresh).^2
    
    part_T_Vthresh = -1; % mV/C
    part_T_Vref = part_T_Vthresh;
    
    % make a function out of it
end


%---------------------------------------------------------------------------------------
if (selector == 2)
end