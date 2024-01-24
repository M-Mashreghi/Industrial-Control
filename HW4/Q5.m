clear all
close all
clc


%%
t = 0:0.01:100;

s = tf('s')
G5 = (2*exp(-0.663*s))/((s+1)*(0.2*s+1)*(0.04*s+1)*(0.008*s+1))
[y,~] = step(G5,t);
margin(G5)
[Gm,Pm,Wcg,Wcp] = margin(G5)
nyquist(G5)
GG_estimate = 2*exp(-0.82*s)/(1.08*s+1)
step(G5)
fig = gca;
grid on;
figure;
g1 =2*exp(-0.82*s)/((2.3-0.82)*s+1) ;
[y1,~] = step(g1,t);
step(g1,G5)
legend("G_{max slope}","G")
figure
g2 = 2*exp(-0.82*s)/((1.95-0.82)*s+1) ;
[y2,~] = step(g2,t);
step(g2,G5)
legend("G_{one point}","G")
figure
g3 = 2*exp(-0.82*s)/(1.08*s+1) ;
[y3,~] = step(g3,t);
step(g3,G5)
legend("G_{two point}","G")
figure
g4 = 2*exp(-0.82*s)/(0.9*s+1) ;
[y4,~] = step(g4,t);
step(g4,G5)
legend("G_{ART}","G")

err_max_slope = mse(t.',y,y1)
err_one_point = mse(t.',y,y2)
err_two_point = mse(t.',y,y3)
err_ART = mse(t.',y,y4)

fig.XMinorGrid = 'on';
fig.YMinorGrid = 'on';
fig.XMinorTick = 'on';
fig.YMinorTick = 'on';
%%
k = 2;
tau_d = 0.82
tau = 1.08
alpha = tau_d/tau

% Initialize an empty list to store structures
resultsList = [];
dataTimeList0 = {};
dataTimeList1 = {};
dataTimeList2 = {};
dataTimeList3 = {};

%%
% ZN open loop
file_name = 'untitled.slx'
method_name ='ZN_open_loop'

kp = 1.2/(k*alpha)
Ti = 2 *tau_d
Td = 0.5 *tau_d

PID1 = sim(file_name);
Step_inf = stepinfo(PID1.Step_response.Data,PID1.Step_response.Time);
Step_inf.IE = PID1.IE.Data(end);
Step_inf.IAE = PID1.IAE.Data(end);
Step_inf.ISE = PID1.ISE.Data(end);
Step_inf.ITAE = PID1.ITAE.Data(end);
Step_inf.Sum_Abs_control_signal = PID1.Sum_Abs_control_signal.Data(end);
Step_inf.Max_control_signal = PID1.Max_control_signal.Data(end);
Step_inf.IE_D = PID1.IE_D.Data(end);
Step_inf.IAE_D = PID1.IAE_D.Data(end);
Step_inf.ISE_D = PID1.ISE_D.Data(end);
Step_inf.ITAE_D = PID1.ITAE_D.Data(end);
Step_inf.method_name = method_name;
Step_inf.kp = kp;
Step_inf.Ti = Ti;
Step_inf.Td = Td;
dataTimeList0{end+1} = struct('Data', PID1.Step_response.Data, 'Time', PID1.Step_response.Time , 'method_name',method_name);
dataTimeList1{end+1} = struct('Data', PID1.Step_response1.Data, 'Time', PID1.Step_response1.Time , 'method_name',method_name);
dataTimeList2{end+1} = struct('Data', PID1.Control_signal.Data, 'Time', PID1.Control_signal.Time , 'method_name',method_name);
dataTimeList3{end+1} = struct('Data', PID1.Control_signal1.Data, 'Time', PID1.Control_signal1.Time , 'method_name',method_name);

plotAndSaveStepResponse(PID1.Step_response.Data,PID1.Step_response.Time,method_name)
plotAndSaveStepResponsewithDis(PID1.Step_response1.Data,PID1.Step_response1.Time,method_name)
plotAndSavecontrolsignal(PID1.Control_signal.Data, PID1.Control_signal.Time, method_name);
plotAndSavecontrolsignalwithDis(PID1.Control_signal1.Data, PID1.Control_signal1.Time, method_name);

% Append the structure to the list
resultsList = [resultsList, Step_inf];
T = struct2table(Step_inf)
writetable(T, 'test.xlsx', 'Sheet', method_name, 'Range', 'A1');
 
 
%%
% CC
kp = (1/k)*(tau/tau_d)*((4/3)+(tau_d/(4*tau)))
Ti = tau_d *(32+6*(tau_d/tau))/(13+8*(tau_d/tau))
Td = tau_d*(4/(11+2*(tau_d/tau)))

file_name = 'untitled.slx'
method_name ='method_cc'


PID1 = sim(file_name);
Step_inf = stepinfo(PID1.Step_response.Data,PID1.Step_response.Time);
Step_inf.IE = PID1.IE.Data(end);
Step_inf.IAE = PID1.IAE.Data(end);
Step_inf.ISE = PID1.ISE.Data(end);
Step_inf.ITAE = PID1.ITAE.Data(end);
Step_inf.Sum_Abs_control_signal = PID1.Sum_Abs_control_signal.Data(end);
Step_inf.Max_control_signal = PID1.Max_control_signal.Data(end);
Step_inf.IE_D = PID1.IE_D.Data(end);
Step_inf.IAE_D = PID1.IAE_D.Data(end);
Step_inf.ISE_D = PID1.ISE_D.Data(end);
Step_inf.ITAE_D = PID1.ITAE_D.Data(end);
Step_inf.method_name = method_name;
Step_inf.kp = kp;
Step_inf.Ti = Ti;
Step_inf.Td = Td;

dataTimeList0{end+1} = struct('Data', PID1.Step_response.Data, 'Time', PID1.Step_response.Time , 'method_name',method_name);
dataTimeList1{end+1} = struct('Data', PID1.Step_response1.Data, 'Time', PID1.Step_response1.Time , 'method_name',method_name);
dataTimeList2{end+1} = struct('Data', PID1.Control_signal.Data, 'Time', PID1.Control_signal.Time , 'method_name',method_name);
dataTimeList3{end+1} = struct('Data', PID1.Control_signal1.Data, 'Time', PID1.Control_signal1.Time , 'method_name',method_name);

plotAndSaveStepResponse(PID1.Step_response.Data,PID1.Step_response.Time,method_name)
plotAndSaveStepResponsewithDis(PID1.Step_response1.Data,PID1.Step_response1.Time,method_name)
plotAndSavecontrolsignal(PID1.Control_signal.Data, PID1.Control_signal.Time, method_name);
plotAndSavecontrolsignalwithDis(PID1.Control_signal1.Data, PID1.Control_signal1.Time, method_name);


% Append the structure to the list
resultsList = [resultsList, Step_inf];
T = struct2table(Step_inf)
writetable(T, 'test.xlsx', 'Sheet', method_name, 'Range', 'A1'); 
 %%
 
% CHR

% Reference input tracking and 0% overshoot\
kp = (0.6/k)*(tau/tau_d)
Ti = tau
Td = 0.5*tau_d

file_name = 'untitled.slx'
method_name ='CHR ref track 0%OS'


PID1 = sim(file_name);
Step_inf = stepinfo(PID1.Step_response.Data,PID1.Step_response.Time);
Step_inf.IE = PID1.IE.Data(end);
Step_inf.IAE = PID1.IAE.Data(end);
Step_inf.ISE = PID1.ISE.Data(end);
Step_inf.ITAE = PID1.ITAE.Data(end);
Step_inf.Sum_Abs_control_signal = PID1.Sum_Abs_control_signal.Data(end);
Step_inf.Max_control_signal = PID1.Max_control_signal.Data(end);
Step_inf.IE_D = PID1.IE_D.Data(end);
Step_inf.IAE_D = PID1.IAE_D.Data(end);
Step_inf.ISE_D = PID1.ISE_D.Data(end);
Step_inf.ITAE_D = PID1.ITAE_D.Data(end);
Step_inf.method_name = method_name;
Step_inf.kp = kp;
Step_inf.Ti = Ti;
Step_inf.Td = Td;

dataTimeList0{end+1} = struct('Data', PID1.Step_response.Data, 'Time', PID1.Step_response.Time , 'method_name',method_name);
dataTimeList1{end+1} = struct('Data', PID1.Step_response1.Data, 'Time', PID1.Step_response1.Time , 'method_name',method_name);
dataTimeList2{end+1} = struct('Data', PID1.Control_signal.Data, 'Time', PID1.Control_signal.Time , 'method_name',method_name);
dataTimeList3{end+1} = struct('Data', PID1.Control_signal1.Data, 'Time', PID1.Control_signal1.Time , 'method_name',method_name);


plotAndSaveStepResponse(PID1.Step_response.Data,PID1.Step_response.Time,method_name)
plotAndSaveStepResponsewithDis(PID1.Step_response1.Data,PID1.Step_response1.Time,method_name)
plotAndSavecontrolsignal(PID1.Control_signal.Data, PID1.Control_signal.Time, method_name);
plotAndSavecontrolsignalwithDis(PID1.Control_signal1.Data, PID1.Control_signal1.Time, method_name);


% Append the structure to the list
resultsList = [resultsList, Step_inf];
T = struct2table(Step_inf)
writetable(T, 'test.xlsx', 'Sheet', method_name, 'Range', 'A1'); 
 
 %%
 
% CHR


% Reference input tracking and 20% overshoot

kp = (0.95/k)*(tau/tau_d)
Ti = 1.4*tau
Td = 0.47*tau_d


file_name = 'untitled.slx'
method_name ='CHR ref track 20%OS'


PID1 = sim(file_name);
Step_inf = stepinfo(PID1.Step_response.Data,PID1.Step_response.Time);
Step_inf.IE = PID1.IE.Data(end);
Step_inf.IAE = PID1.IAE.Data(end);
Step_inf.ISE = PID1.ISE.Data(end);
Step_inf.ITAE = PID1.ITAE.Data(end);
Step_inf.Sum_Abs_control_signal = PID1.Sum_Abs_control_signal.Data(end);
Step_inf.Max_control_signal = PID1.Max_control_signal.Data(end);
Step_inf.IE_D = PID1.IE_D.Data(end);
Step_inf.IAE_D = PID1.IAE_D.Data(end);
Step_inf.ISE_D = PID1.ISE_D.Data(end);
Step_inf.ITAE_D = PID1.ITAE_D.Data(end);
Step_inf.method_name = method_name;
Step_inf.kp = kp;
Step_inf.Ti = Ti;
Step_inf.Td = Td;

dataTimeList0{end+1} = struct('Data', PID1.Step_response.Data, 'Time', PID1.Step_response.Time , 'method_name',method_name);
dataTimeList1{end+1} = struct('Data', PID1.Step_response1.Data, 'Time', PID1.Step_response1.Time , 'method_name',method_name);
dataTimeList2{end+1} = struct('Data', PID1.Control_signal.Data, 'Time', PID1.Control_signal.Time , 'method_name',method_name);
dataTimeList3{end+1} = struct('Data', PID1.Control_signal1.Data, 'Time', PID1.Control_signal1.Time , 'method_name',method_name);


plotAndSaveStepResponse(PID1.Step_response.Data,PID1.Step_response.Time,method_name) 
plotAndSaveStepResponsewithDis(PID1.Step_response1.Data,PID1.Step_response1.Time,method_name)
plotAndSavecontrolsignal(PID1.Control_signal.Data, PID1.Control_signal.Time, method_name);
plotAndSavecontrolsignalwithDis(PID1.Control_signal1.Data, PID1.Control_signal1.Time, method_name);

% Append the structure to the list
resultsList = [resultsList, Step_inf];
T = struct2table(Step_inf)
writetable(T, 'test.xlsx', 'Sheet', method_name, 'Range', 'A1');
 
 
 %%
 
% CHR


% Removal of disturbance and 0% overshoot

kp = (0.95/k)*(tau/tau_d)
Ti = 2.4*tau
Td = 0.42*tau_d


file_name = 'untitled.slx'
method_name ='CHR rem dis 0%OS'


PID1 = sim(file_name);
Step_inf = stepinfo(PID1.Step_response.Data,PID1.Step_response.Time);
Step_inf.IE = PID1.IE.Data(end);
Step_inf.IAE = PID1.IAE.Data(end);
Step_inf.ISE = PID1.ISE.Data(end);
Step_inf.ITAE = PID1.ITAE.Data(end);
Step_inf.Sum_Abs_control_signal = PID1.Sum_Abs_control_signal.Data(end);
Step_inf.Max_control_signal = PID1.Max_control_signal.Data(end);
Step_inf.IE_D = PID1.IE_D.Data(end);
Step_inf.IAE_D = PID1.IAE_D.Data(end);
Step_inf.ISE_D = PID1.ISE_D.Data(end);
Step_inf.ITAE_D = PID1.ITAE_D.Data(end);
Step_inf.method_name = method_name;
Step_inf.kp = kp;
Step_inf.Ti = Ti;
Step_inf.Td = Td;
Step_inf.method_name = method_name;

dataTimeList0{end+1} = struct('Data', PID1.Step_response.Data, 'Time', PID1.Step_response.Time , 'method_name',method_name);
dataTimeList1{end+1} = struct('Data', PID1.Step_response1.Data, 'Time', PID1.Step_response1.Time , 'method_name',method_name);
dataTimeList2{end+1} = struct('Data', PID1.Control_signal.Data, 'Time', PID1.Control_signal.Time , 'method_name',method_name);
dataTimeList3{end+1} = struct('Data', PID1.Control_signal1.Data, 'Time', PID1.Control_signal1.Time , 'method_name',method_name);

plotAndSaveStepResponse(PID1.Step_response.Data,PID1.Step_response.Time,method_name) 
plotAndSaveStepResponsewithDis(PID1.Step_response1.Data,PID1.Step_response1.Time,method_name)
plotAndSavecontrolsignal(PID1.Control_signal.Data, PID1.Control_signal.Time, method_name);
plotAndSavecontrolsignalwithDis(PID1.Control_signal1.Data, PID1.Control_signal1.Time, method_name);

% Append the structure to the list
resultsList = [resultsList, Step_inf];
T = struct2table(Step_inf)
writetable(T, 'test.xlsx', 'Sheet', method_name, 'Range', 'A1');
 
 
 %%
 
% CHR

% Removal of disturbance and 20% overshoot

kp = (1.2/k)*(tau/tau_d)
Ti = 2*tau
Td = 0.42*tau_d


file_name = 'untitled.slx'
method_name ='CHR rem dis 20%OS'


PID1 = sim(file_name);
Step_inf = stepinfo(PID1.Step_response.Data,PID1.Step_response.Time);
Step_inf.IE = PID1.IE.Data(end);
Step_inf.IAE = PID1.IAE.Data(end);
Step_inf.ISE = PID1.ISE.Data(end);
Step_inf.ITAE = PID1.ITAE.Data(end);
Step_inf.Sum_Abs_control_signal = PID1.Sum_Abs_control_signal.Data(end);
Step_inf.Max_control_signal = PID1.Max_control_signal.Data(end);
Step_inf.IE_D = PID1.IE_D.Data(end);
Step_inf.IAE_D = PID1.IAE_D.Data(end);
Step_inf.ISE_D = PID1.ISE_D.Data(end);
Step_inf.ITAE_D = PID1.ITAE_D.Data(end);
Step_inf.method_name = method_name;
Step_inf.kp = kp;
Step_inf.Ti = Ti;
Step_inf.Td = Td;

dataTimeList0{end+1} = struct('Data', PID1.Step_response.Data, 'Time', PID1.Step_response.Time , 'method_name',method_name);
dataTimeList1{end+1} = struct('Data', PID1.Step_response1.Data, 'Time', PID1.Step_response1.Time , 'method_name',method_name);
dataTimeList2{end+1} = struct('Data', PID1.Control_signal.Data, 'Time', PID1.Control_signal.Time , 'method_name',method_name);
dataTimeList3{end+1} = struct('Data', PID1.Control_signal1.Data, 'Time', PID1.Control_signal1.Time , 'method_name',method_name);

plotAndSaveStepResponse(PID1.Step_response.Data,PID1.Step_response.Time,method_name) 
plotAndSaveStepResponsewithDis(PID1.Step_response1.Data,PID1.Step_response1.Time,method_name)
plotAndSavecontrolsignal(PID1.Control_signal.Data, PID1.Control_signal.Time, method_name);
plotAndSavecontrolsignalwithDis(PID1.Control_signal1.Data, PID1.Control_signal1.Time, method_name);

% Append the structure to the list
resultsList = [resultsList, Step_inf];
T = struct2table(Step_inf)
writetable(T, 'test.xlsx', 'Sheet', method_name, 'Range', 'A1');
 

 %%
 
% Damped Ocsillation

k_dmp = 0.78
T_dmp = 3.06
kp = 1.1*k_dmp
Ti = T_dmp/3.6
Td = T_dmp/9

file_name = 'untitled.slx'
method_name ='Damped Ocsillation'


PID1 = sim(file_name);
Step_inf = stepinfo(PID1.Step_response.Data,PID1.Step_response.Time);
Step_inf.IE = PID1.IE.Data(end);
Step_inf.IAE = PID1.IAE.Data(end);
Step_inf.ISE = PID1.ISE.Data(end);
Step_inf.ITAE = PID1.ITAE.Data(end);
Step_inf.Sum_Abs_control_signal = PID1.Sum_Abs_control_signal.Data(end);
Step_inf.Max_control_signal = PID1.Max_control_signal.Data(end);
Step_inf.IE_D = PID1.IE_D.Data(end);
Step_inf.IAE_D = PID1.IAE_D.Data(end);
Step_inf.ISE_D = PID1.ISE_D.Data(end);
Step_inf.ITAE_D = PID1.ITAE_D.Data(end);
Step_inf.method_name = method_name;
Step_inf.kp = kp;
Step_inf.Ti = Ti;
Step_inf.Td = Td;

dataTimeList0{end+1} = struct('Data', PID1.Step_response.Data, 'Time', PID1.Step_response.Time , 'method_name',method_name);
dataTimeList1{end+1} = struct('Data', PID1.Step_response1.Data, 'Time', PID1.Step_response1.Time , 'method_name',method_name);
dataTimeList2{end+1} = struct('Data', PID1.Control_signal.Data, 'Time', PID1.Control_signal.Time , 'method_name',method_name);
dataTimeList3{end+1} = struct('Data', PID1.Control_signal1.Data, 'Time', PID1.Control_signal1.Time , 'method_name',method_name);

plotAndSaveStepResponse(PID1.Step_response.Data,PID1.Step_response.Time,method_name) 
plotAndSaveStepResponsewithDis(PID1.Step_response1.Data,PID1.Step_response1.Time,method_name)
plotAndSavecontrolsignal(PID1.Control_signal.Data, PID1.Control_signal.Time, method_name);
plotAndSavecontrolsignalwithDis(PID1.Control_signal1.Data, PID1.Control_signal1.Time, method_name);

% Append the structure to the list
resultsList = [resultsList, Step_inf];
T = struct2table(Step_inf)
writetable(T, 'test.xlsx', 'Sheet', method_name, 'Range', 'A1');


 %%
 
% ZN close loop

ku = 1.376
Tu = 2.75
kp = 0.6*ku
Ti = Tu/2
Td = Tu/8


file_name = 'untitled.slx'
method_name = 'ZN_close_loop'


PID1 = sim(file_name);
Step_inf = stepinfo(PID1.Step_response.Data,PID1.Step_response.Time);
Step_inf.IE = PID1.IE.Data(end);
Step_inf.IAE = PID1.IAE.Data(end);
Step_inf.ISE = PID1.ISE.Data(end);
Step_inf.ITAE = PID1.ITAE.Data(end);
Step_inf.Sum_Abs_control_signal = PID1.Sum_Abs_control_signal.Data(end);
Step_inf.Max_control_signal = PID1.Max_control_signal.Data(end);
Step_inf.IE_D = PID1.IE_D.Data(end);
Step_inf.IAE_D = PID1.IAE_D.Data(end);
Step_inf.ISE_D = PID1.ISE_D.Data(end);
Step_inf.ITAE_D = PID1.ITAE_D.Data(end);
Step_inf.method_name = method_name;
Step_inf.kp = kp;
Step_inf.Ti = Ti;
Step_inf.Td = Td;

dataTimeList0{end+1} = struct('Data', PID1.Step_response.Data, 'Time', PID1.Step_response.Time , 'method_name',method_name);
dataTimeList1{end+1} = struct('Data', PID1.Step_response1.Data, 'Time', PID1.Step_response1.Time , 'method_name',method_name);
dataTimeList2{end+1} = struct('Data', PID1.Control_signal.Data, 'Time', PID1.Control_signal.Time , 'method_name',method_name);
dataTimeList3{end+1} = struct('Data', PID1.Control_signal1.Data, 'Time', PID1.Control_signal1.Time , 'method_name',method_name);

plotAndSaveStepResponse(PID1.Step_response.Data,PID1.Step_response.Time,method_name) 
plotAndSaveStepResponsewithDis(PID1.Step_response1.Data,PID1.Step_response1.Time,method_name)
plotAndSavecontrolsignal(PID1.Control_signal.Data, PID1.Control_signal.Time, method_name);
plotAndSavecontrolsignalwithDis(PID1.Control_signal1.Data, PID1.Control_signal1.Time, method_name);

% Append the structure to the list
resultsList = [resultsList, Step_inf];
T = struct2table(Step_inf)
writetable(T, 'test.xlsx', 'Sheet', method_name, 'Range', 'A1');


 %%
 
% Fertik ref track

alpha_f = tau_d/(tau+tau_d)
Tps = tau + tau_d
kp = (1/k) * (0.6)
Ti = Tps * (0.45)
Td = Tps * (0.36)


file_name = 'untitled.slx'
method_name = 'Fertik ref track'


PID1 = sim(file_name);
Step_inf = stepinfo(PID1.Step_response.Data,PID1.Step_response.Time);
Step_inf.IE = PID1.IE.Data(end);
Step_inf.IAE = PID1.IAE.Data(end);
Step_inf.ISE = PID1.ISE.Data(end);
Step_inf.ITAE = PID1.ITAE.Data(end);
Step_inf.Sum_Abs_control_signal = PID1.Sum_Abs_control_signal.Data(end);
Step_inf.Max_control_signal = PID1.Max_control_signal.Data(end);
Step_inf.IE_D = PID1.IE_D.Data(end);
Step_inf.IAE_D = PID1.IAE_D.Data(end);
Step_inf.ISE_D = PID1.ISE_D.Data(end);
Step_inf.ITAE_D = PID1.ITAE_D.Data(end);
Step_inf.method_name = method_name;
Step_inf.kp = kp;
Step_inf.Ti = Ti;
Step_inf.Td = Td;

dataTimeList0{end+1} = struct('Data', PID1.Step_response.Data, 'Time', PID1.Step_response.Time , 'method_name',method_name);
dataTimeList1{end+1} = struct('Data', PID1.Step_response1.Data, 'Time', PID1.Step_response1.Time , 'method_name',method_name);
dataTimeList2{end+1} = struct('Data', PID1.Control_signal.Data, 'Time', PID1.Control_signal.Time , 'method_name',method_name);
dataTimeList3{end+1} = struct('Data', PID1.Control_signal1.Data, 'Time', PID1.Control_signal1.Time , 'method_name',method_name);

plotAndSaveStepResponse(PID1.Step_response.Data,PID1.Step_response.Time,method_name) 
plotAndSaveStepResponsewithDis(PID1.Step_response1.Data,PID1.Step_response1.Time,method_name)
plotAndSavecontrolsignal(PID1.Control_signal.Data, PID1.Control_signal.Time, method_name);
plotAndSavecontrolsignalwithDis(PID1.Control_signal1.Data, PID1.Control_signal1.Time, method_name);


% Append the structure to the list
resultsList = [resultsList, Step_inf];
T = struct2table(Step_inf)
writetable(T, 'test.xlsx', 'Sheet', method_name, 'Range', 'A1');

 %%
 
% Fertik rem dis

alpha_f = tau_d/(tau+tau_d)
Tps = tau + tau_d
kp = (1/k) * (0.73)
Ti = Tps * (0.45)
Td = Tps * (0.36)


file_name = 'untitled.slx'
method_name = 'Fertik rem dis'


PID1 = sim(file_name);
Step_inf = stepinfo(PID1.Step_response.Data,PID1.Step_response.Time);
Step_inf.IE = PID1.IE.Data(end);
Step_inf.IAE = PID1.IAE.Data(end);
Step_inf.ISE = PID1.ISE.Data(end);
Step_inf.ITAE = PID1.ITAE.Data(end);
Step_inf.Sum_Abs_control_signal = PID1.Sum_Abs_control_signal.Data(end);
Step_inf.Max_control_signal = PID1.Max_control_signal.Data(end);
Step_inf.IE_D = PID1.IE_D.Data(end);
Step_inf.IAE_D = PID1.IAE_D.Data(end);
Step_inf.ISE_D = PID1.ISE_D.Data(end);
Step_inf.ITAE_D = PID1.ITAE_D.Data(end);
Step_inf.method_name = method_name;
Step_inf.kp = kp;
Step_inf.Ti = Ti;
Step_inf.Td = Td;

dataTimeList0{end+1} = struct('Data', PID1.Step_response.Data, 'Time', PID1.Step_response.Time , 'method_name',method_name);
dataTimeList1{end+1} = struct('Data', PID1.Step_response1.Data, 'Time', PID1.Step_response1.Time , 'method_name',method_name);
dataTimeList2{end+1} = struct('Data', PID1.Control_signal.Data, 'Time', PID1.Control_signal.Time , 'method_name',method_name);
dataTimeList3{end+1} = struct('Data', PID1.Control_signal1.Data, 'Time', PID1.Control_signal1.Time , 'method_name',method_name);


plotAndSaveStepResponse(PID1.Step_response.Data,PID1.Step_response.Time,method_name) 
plotAndSaveStepResponsewithDis(PID1.Step_response1.Data,PID1.Step_response1.Time,method_name)
plotAndSavecontrolsignal(PID1.Control_signal.Data, PID1.Control_signal.Time, method_name);
plotAndSavecontrolsignalwithDis(PID1.Control_signal1.Data, PID1.Control_signal1.Time, method_name);

% Append the structure to the list
resultsList = [resultsList, Step_inf];
T = struct2table(Step_inf)
writetable(T, 'test.xlsx', 'Sheet', method_name, 'Range', 'A1');





 %%
 
% Ciancone-Marline

alpha_f = tau_d/(tau+tau_d)
kp = (1/k) * (0.8)
Ti = (tau+tau_d) * (0.73)
Td = (tau+tau_d) * (0.6)


file_name = 'untitled.slx'
method_name = 'Ciancone-Marline ref track'


PID1 = sim(file_name);
Step_inf = stepinfo(PID1.Step_response.Data,PID1.Step_response.Time);
Step_inf.IE = PID1.IE.Data(end);
Step_inf.IAE = PID1.IAE.Data(end);
Step_inf.ISE = PID1.ISE.Data(end);
Step_inf.ITAE = PID1.ITAE.Data(end);
Step_inf.Sum_Abs_control_signal = PID1.Sum_Abs_control_signal.Data(end);
Step_inf.Max_control_signal = PID1.Max_control_signal.Data(end);
Step_inf.IE_D = PID1.IE_D.Data(end);
Step_inf.IAE_D = PID1.IAE_D.Data(end);
Step_inf.ISE_D = PID1.ISE_D.Data(end);
Step_inf.ITAE_D = PID1.ITAE_D.Data(end);
Step_inf.method_name = method_name;
Step_inf.kp = kp;
Step_inf.Ti = Ti;
Step_inf.Td = Td;

dataTimeList0{end+1} = struct('Data', PID1.Step_response.Data, 'Time', PID1.Step_response.Time , 'method_name',method_name);
dataTimeList1{end+1} = struct('Data', PID1.Step_response1.Data, 'Time', PID1.Step_response1.Time , 'method_name',method_name);
dataTimeList2{end+1} = struct('Data', PID1.Control_signal.Data, 'Time', PID1.Control_signal.Time , 'method_name',method_name);
dataTimeList3{end+1} = struct('Data', PID1.Control_signal1.Data, 'Time', PID1.Control_signal1.Time , 'method_name',method_name);

plotAndSaveStepResponse(PID1.Step_response.Data,PID1.Step_response.Time,method_name) 
plotAndSaveStepResponsewithDis(PID1.Step_response1.Data,PID1.Step_response1.Time,method_name)
plotAndSavecontrolsignal(PID1.Control_signal.Data, PID1.Control_signal.Time, method_name);
plotAndSavecontrolsignalwithDis(PID1.Control_signal1.Data, PID1.Control_signal1.Time, method_name);

% Append the structure to the list
resultsList = [resultsList, Step_inf];
T = struct2table(Step_inf)
writetable(T, 'test.xlsx', 'Sheet', method_name, 'Range', 'A1');





 %%
 
% Ciancone-Marline
alpha_f = tau_d/(tau+tau_d)
kp = (1/k) * (0.95)
Ti = (tau+tau_d) * (0.68)
Td = (tau+tau_d) * (0.4)


file_name = 'untitled.slx'
method_name = 'Ciancone-Marline rem dis'


PID1 = sim(file_name);
Step_inf = stepinfo(PID1.Step_response.Data,PID1.Step_response.Time);
Step_inf.IE = PID1.IE.Data(end);
Step_inf.IAE = PID1.IAE.Data(end);
Step_inf.ISE = PID1.ISE.Data(end);
Step_inf.ITAE = PID1.ITAE.Data(end);
Step_inf.Sum_Abs_control_signal = PID1.Sum_Abs_control_signal.Data(end);
Step_inf.Max_control_signal = PID1.Max_control_signal.Data(end);
Step_inf.IE_D = PID1.IE_D.Data(end);
Step_inf.IAE_D = PID1.IAE_D.Data(end);
Step_inf.ISE_D = PID1.ISE_D.Data(end);
Step_inf.ITAE_D = PID1.ITAE_D.Data(end);
Step_inf.method_name = method_name;
Step_inf.kp = kp;
Step_inf.Ti = Ti;
Step_inf.Td = Td;

dataTimeList0{end+1} = struct('Data', PID1.Step_response.Data, 'Time', PID1.Step_response.Time , 'method_name',method_name);
dataTimeList1{end+1} = struct('Data', PID1.Step_response1.Data, 'Time', PID1.Step_response1.Time , 'method_name',method_name);
dataTimeList2{end+1} = struct('Data', PID1.Control_signal.Data, 'Time', PID1.Control_signal.Time , 'method_name',method_name);
dataTimeList3{end+1} = struct('Data', PID1.Control_signal1.Data, 'Time', PID1.Control_signal1.Time , 'method_name',method_name);

plotAndSaveStepResponse(PID1.Step_response.Data,PID1.Step_response.Time,method_name) 
plotAndSaveStepResponsewithDis(PID1.Step_response1.Data,PID1.Step_response1.Time,method_name)
plotAndSavecontrolsignal(PID1.Control_signal.Data, PID1.Control_signal.Time, method_name);
plotAndSavecontrolsignalwithDis(PID1.Control_signal1.Data, PID1.Control_signal1.Time, method_name);

% Append the structure to the list
resultsList = [resultsList, Step_inf];
T = struct2table(Step_inf)
writetable(T, 'test.xlsx', 'Sheet', method_name, 'Range', 'A1');


 %%
 
% process response
ku = 1.376
Tiu = 0.65
Tdu = 0.95
kp = 0.5*ku
Ti = 3.3*Tiu
Td = 0.3*Tdu



file_name = 'untitled.slx'
method_name = 'process response'


PID1 = sim(file_name);
Step_inf = stepinfo(PID1.Step_response.Data,PID1.Step_response.Time);
Step_inf.IE = PID1.IE.Data(end);
Step_inf.IAE = PID1.IAE.Data(end);
Step_inf.ISE = PID1.ISE.Data(end);
Step_inf.ITAE = PID1.ITAE.Data(end);
Step_inf.Sum_Abs_control_signal = PID1.Sum_Abs_control_signal.Data(end);
Step_inf.Max_control_signal = PID1.Max_control_signal.Data(end);
Step_inf.IE_D = PID1.IE_D.Data(end);
Step_inf.IAE_D = PID1.IAE_D.Data(end);
Step_inf.ISE_D = PID1.ISE_D.Data(end);
Step_inf.ITAE_D = PID1.ITAE_D.Data(end);
Step_inf.method_name = method_name;
Step_inf.kp = kp;
Step_inf.Ti = Ti;
Step_inf.Td = Td;

dataTimeList0{end+1} = struct('Data', PID1.Step_response.Data, 'Time', PID1.Step_response.Time , 'method_name',method_name);
dataTimeList1{end+1} = struct('Data', PID1.Step_response1.Data, 'Time', PID1.Step_response1.Time , 'method_name',method_name);
dataTimeList2{end+1} = struct('Data', PID1.Control_signal.Data, 'Time', PID1.Control_signal.Time , 'method_name',method_name);
dataTimeList3{end+1} = struct('Data', PID1.Control_signal1.Data, 'Time', PID1.Control_signal1.Time , 'method_name',method_name);

plotAndSaveStepResponse(PID1.Step_response.Data,PID1.Step_response.Time,method_name) 
plotAndSaveStepResponsewithDis(PID1.Step_response1.Data,PID1.Step_response1.Time,method_name)
plotAndSavecontrolsignal(PID1.Control_signal.Data, PID1.Control_signal.Time, method_name);
plotAndSavecontrolsignalwithDis(PID1.Control_signal1.Data, PID1.Control_signal1.Time, method_name);

% Append the structure to the list
resultsList = [resultsList, Step_inf];
T = struct2table(Step_inf)
writetable(T, 'test.xlsx', 'Sheet', method_name, 'Range', 'A1');






% Convert the list of structures to a table
T = struct2table(resultsList);

% Rearrange columns: move last 4 columns to the beginning
T = [T(:, end-3:end), T(:, 1:end-4)];

% Write the table to an Excel file with the specified sheet name
writetable(T, 'test2.xlsx', 'Sheet', 'finall', 'Range', 'A1');


plotAndSaveMultipleControlSignals(dataTimeList0,'steprespone');
plotAndSaveMultipleControlSignals(dataTimeList1,'stepresponewithdis');
plotAndSaveMultipleControlSignals(dataTimeList2,'controlsignal');
plotAndSaveMultipleControlSignals(dataTimeList3,'controlsignalwithdis');

%%
function plotAndSaveStepResponse(Data, Time, method_name)
    % Plot the Step Response
    figure;
    plot(Time, Data);
    title(method_name);
    xlabel('Time');
    ylabel('Amplitude');
    grid on;

    % Save the plot as an image file (e.g., PNG)
    saveas(gcf, ['Step_Response_Plot_', method_name, '.png']);
end

 
function plotAndSaveStepResponsewithDis(Data, Time, method_name)
    % Plot the Step Response
    figure;
    plot(Time, Data);
    title(method_name);
    xlabel('Time');
    ylabel('Amplitude');
    grid on;

    % Save the plot as an image file (e.g., PNG)
    saveas(gcf, ['Step_Response with Dis_Plot_', method_name, '.png']);
 end

function plotAndSavecontrolsignal(Data, Time, method_name)
    % Plot the Step Response
    figure;
    plot(Time, Data);
    title(method_name);
    xlabel('Time');
    ylabel('signal control');
    grid on;

    % Save the plot as an image file (e.g., PNG)
    saveas(gcf, ['control_signal_Plot_', method_name, '.png']);
end

function plotAndSavecontrolsignalwithDis(Data, Time, method_name)
    % Plot the Step Response
    figure;
    plot(Time, Data);
    title(method_name);
    xlabel('Time');
    ylabel('signal control');
    grid on;

    % Save the plot as an image file (e.g., PNG)
    saveas(gcf, ['control_signal_with_Dis_Plot_', method_name, '.png']);
end

function plotAndSaveMultipleControlSignals(dataTimeList, group_name)
    % Plot all control signals together
    figure;
    
    for i = 1:length(dataTimeList)
        currentData = dataTimeList{i}.Data;
        currentTime = dataTimeList{i}.Time;
        method_name = dataTimeList{i}.method_name;

        % Plot each set of data
        plot(currentTime, currentData, 'DisplayName', method_name);
        hold on;
    end

    title(['Multiple Control Signals',group_name]);
    xlabel('Time');
    ylabel('Signal Control');
    legend('show');
    grid on;

    % Save the plot as an image file (e.g., PNG)
    saveas(gcf,['multiple_control_signals_plot', group_name, '.png']);
end

 