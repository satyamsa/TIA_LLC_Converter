%% TIA model 
function [fsw,IL_rms] = TIA_LLC_FB(Pload,Lr,Cr,Lm,Vin,Vout,n)

m = Lm/Lr;
fun = @(x)root2d(x,Vin,Pload,Lr,Cr,Vout,n,m);

x0 = [0,1.5338];

options = optimoptions('fsolve'); 
options.MaxIterations = 75000;
options.MaxFunctionEvaluations = 75000;
x = fsolve(fun,x0,options);

    A = x(1);
    B = x(2);
    fsw = 1/(4*B*sqrt(Lr*Cr));
    t1 = A*sqrt(Lr*Cr) ;


% above resonant frequency full bridge LLC converter 
% I1 = -n*Vout/(4*Lm*fsw);
% Ilm1 = I1;
% V1 = -(2*Lr*Vin*m - Lr*Vin*m*cos(A) - Lr*Vin*m*cos(A - 2*B) - Lr*Vout*m*n*cos(A) + Lr*Vout*m*n*cos(A - 2*B) + B*Vout*n*sin(A)*(Lr/Cr)^(1/2)*(Cr*Lr)^(1/2) + B*Vout*n*sin(A - 2*B)*(Lr/Cr)^(1/2)*(Cr*Lr)^(1/2))/(Lr*m*(cos(A) + cos(A - 2*B)));
% Ilm0 = n*Vout*(A/B -1)/(4*Lm*fsw);
% I0 = -(Lr*Vin*m*sin(2*A)*(Cr/Lr)^(1/2) + B*Vout*n*cos(A)*(Cr*Lr)^(1/2) + B*Vout*n*cos(A - 2*B)*(Cr*Lr)^(1/2) + Lr*Vout*m*n*(Cr/Lr)^(1/2)*(sin(3*A - 2*B)/2 - sin(A - 2*B)/2) + Lr*Vout*m*n*cos(2*B)*sin(A)*(Cr/Lr)^(1/2) - B*Vout*n*sin(2*B)*sin(A)*(Cr/Lr)^(1/2)*(Lr/Cr)^(1/2)*(Cr*Lr)^(1/2))/(Lr*m*cos(A)*(cos(A) + cos(A - 2*B)));
% V0 = (Lr*Vin*m*cos(A - 2*B) - Lr*Vin*m*cos(A) + Lr*Vout*m*n*cos(A) + Lr*Vout*m*n*cos(A - 2*B) - Lr*Vout*m*n*cos(2*A - 2*B) - Lr*Vout*m*n*cos(2*B) + B*Vout*n*sin(2*B)*(Lr/Cr)^(1/2)*(Cr*Lr)^(1/2))/(Lr*m*(cos(A) + cos(A - 2*B)));
% above resonant frequency full bridge LLC Converter kharan 

I1 = -n*Vout/(4*Lm*fsw);
Ilm1 = I1;
Ilm0 = n*Vout*(A/B -1)/(4*Lm*fsw);
V0= n*Vout + ((n*Vout*B*sin(2*B)/m) + (Vin*(cos(2*B-A)-cos(A))) -(2*n*Vout*cos(2*B-A)*cos(A)))/(2*cos(B)*cos(A-B));
V1 = Vin + (n*Vout*(cos(A)-cos(2*B-A))/(2*cos(B)*cos(A-B))) -((B*n*Vout/m)*tan(A-B))-(2*Vin/(2*cos(B)*cos(A-B)));
I0 = -1*((2*(Vin + (n*Vout*cos(2*B-A))))*sqrt(Cr/Lr)*sin(A)/(2*cos(B)*cos(A-B))+(n*Vout*cos(B)/(4*Lm*fsw*cos(A-B))));


% Pout_1 = Vin*Cr*(cos(A)-cos(2*B -A))/(cos(B)*cos(A-B));
% Pout_2 = -2*n*Vout*Cr*B*sin(B)*cos(B)/(m*cos(A-B)*cos(B));
% Pout_3 = n*Cr*Vout*((2*cos(A)*cos(2*B - A))-cos(A)-cos(2*B -A))/(cos(B)*cos(A-B));
% Pout_der = (Pout_1 + Pout_2 + Pout_3)*Vin*2*fsw;
% at resonance frequency 

%% Reconstrcuting the waveform of ILr,ILm,Vcr
time = 0:1e-10:(1/(2*fsw));
t_fl = 0:1e-10:(1/fsw);
[~,time_col] = size(time);
[~,t_fl_col] = size(t_fl);
Vcr = zeros(1,t_fl_col);
ILr = zeros(1,t_fl_col);
ILm = zeros(1,t_fl_col);

for i = 1:1:time_col 
    if( time(1,i) <= t1)
        Vcr(1,i) = Vin + n*Vout + (I0*sqrt(Lr/Cr)*sin(time(1,i)/sqrt(Lr*Cr))) + (V0 - Vin - n*Vout)*cos(time(1,i)/sqrt(Lr*Cr)); 
        ILr(1,i) = I0*cos(time(1,i)/sqrt(Lr*Cr)) - (V0 - Vin -n*Vout)*sqrt(Cr/Lr)*sin(time(1,i)/sqrt(Lr*Cr));
        ILm(1,i) = -(n*Vout*time(1,i)/Lm) + Ilm0;
    else
        Vcr(1,i) = Vin - n*Vout + (I1*sqrt(Lr/Cr)*sin((time(1,i)-t1)/sqrt(Lr*Cr))) + (V1 - Vin + n*Vout)*cos((time(1,i)-t1)/sqrt(Lr*Cr));
        ILr(1,i) = I1*cos((time(1,i)-t1)/sqrt(Lr*Cr)) - (V1 - Vin + n*Vout)*sqrt(Cr/Lr)*sin((time(1,i)-t1)/sqrt(Lr*Cr));
        ILm(1,i) = (n*Vout*(time(1,i)-t1)/Lm) + Ilm1;
    end
    
end

for i = (time_col):1:t_fl_col
    ILr(1,i) = -1*ILr(1,i-time_col+1);
    ILm(1,i) = -1*ILm(1,i-time_col+1);
    Vcr(1,i) = -1*Vcr(1,i-time_col+1);
end
figure
plot(t_fl,ILr);
hold on;
plot(t_fl,ILm);
IL_rms = rms(ILr);

figure
plot(Vcr,ILr)

%% Functions 
function F = root2d(x,Vin,Pload,Lr,Cr,Vout,n,m)

    F(1) = (16*Cr*(Vin^2)*sin(x(1)-x(2))*sin((x(1)/2)-x(2))*sin(x(1)/2))/(4*x(2)*sqrt(Lr*Cr)*((x(2)*(cos(2*x(2))+1)/m)+sin(2*x(2))))-Pload;
    F(2) = (-1/n)*(sin(x(1)-x(2)))/((x(2)*cos(x(2))/m)+sin(x(2))) - Vout/Vin ;
 
end


end

