function [t1_ini,t2_ini,fsw_ini] = TIA_LLC_FB_BR(Pload,Lr,Cr,Lm,Vin,Vout,n)
    m = Lm/Lr;
    fun = @(x)root2d(x,Vin,Pload,Lr,Cr,Lm,Vout,n,m);
    
    

    x0 = [3.33,0];
    options = optimoptions('fsolve'); 
    options.MaxIterations = 75000;
    options.MaxFunctionEvaluations = 75000;
    x = fsolve(fun,x0,options);
    F = root2d(x,Vin,Pload,Lr,Cr,Lm,Vout,n,m);
    A = x(1);
    T = x(2);
    t1 = A*sqrt(Lr*Cr) ;
    fsw = 1/(2*(T*sqrt((Lr+Lm)*Cr) + t1));
   

    %below resonance

    V0 = Vin - (n*Vout) +  ((A/m)*n*Vout*((cos(T)*sin(A)) + (sqrt(m+1)*sin(T)*cos(A)))/(((cos(A)-1)*(cos(T)-1)) - (sqrt(m+1)*sin(A)*sin(T)))) +  ((((2*Vin) - (n*Vout*(cos(T)+1)))*(cos(A)-1))/(((cos(A)-1)*(cos(T)-1)) - (sqrt(m+1)*sin(A)*sin(T)))); 
    V1 = (2*Vin) -(2*n*Vout) - V0 - (sqrt(Lr/Cr)*((n*Vout*t1*cot(A/2))/Lm));
    I0 = -((n*Vout*t1)/((1 - cos(A))*Lm)) - ((V0 - Vin +(n*Vout))*sqrt(Cr/Lr)*cot(A/2));
    I1 = (n*Vout*t1/Lm) + I0;
    Ilm1 = I1;
    Ilm0 = I0;
    
    %% Reconstrcuting the waveform of ILr,ILm,Vcr
    time = 0:1e-10:(1/(2*fsw));
    t_fl = 0:1e-10:(1/fsw);
    [~,time_col] = size(time);
    [~,t_fl_col] = size(t_fl);
    Vcr = zeros(1,t_fl_col);
    ILr = zeros(1,t_fl_col);
    ILm = zeros(1,t_fl_col);
    VLm = zeros(1,t_fl_col);
    t_req = 0;
    for i = 1:1:time_col 
        if( time(1,i) <= t1)
            Vcr(1,i) = Vin - n*Vout + (I0*sqrt(Lr/Cr)*sin(time(1,i)/sqrt(Lr*Cr))) + (V0 - Vin + n*Vout)*cos(time(1,i)/sqrt(Lr*Cr)); 
            ILr(1,i) = I0*cos(time(1,i)/sqrt(Lr*Cr)) - (V0 - Vin + n*Vout)*sqrt(Cr/Lr)*sin(time(1,i)/sqrt(Lr*Cr));
            ILm(1,i) = (n*Vout*time(1,i)/Lm) + Ilm0;
            VLm(1,i) = n*Vout;
            if ILr(1,i) < ILm(1,i) 
                t_req = time(1,i);
            end
        else
            Vcr(1,i) = Vin  + (I1*sqrt((Lr+Lm)/Cr)*sin((time(1,i)-t1)/sqrt((Lr+Lm)*Cr))) + (V1 - Vin)*cos((time(1,i)-t1)/sqrt((Lr+Lm)*Cr));
            ILr(1,i) = I1*cos((time(1,i)-t1)/sqrt((Lr+Lm)*Cr)) - (V1 - Vin)*sqrt(Cr/(Lr+Lm))*sin((time(1,i)-t1)/sqrt((Lr + Lm)*Cr));
            ILm(1,i) = ILr(1,i);
            VLm(1,i) = -1*Lm*(I1*sin((time(1,i)-t1)/sqrt((Lr+Lm)*Cr)) + (V1 - Vin)*sqrt(Cr/(Lr+Lm))*cos((time(1,i)-t1)/sqrt((Lr + Lm)*Cr)))/(sqrt((Lr+Lm)*Cr));
        end

    end

    for i = (time_col + 1):1:t_fl_col
        ILr(1,i) = -1*ILr(1,i-time_col);
        ILm(1,i) = -1*ILm(1,i-time_col);
        Vcr(1,i) = -1*Vcr(1,i-time_col);
        VLm(1,i) = -1*VLm(1,i-time_col);

    end
    figure 
    plot(t_fl,ILr);
    hold on;
    plot(t_fl,ILm);

    t1_ini = t1 - t_req;
    a = ((1/(2*fsw))-t_req)/(1/(2*fsw));
    t2_ini = a*(1/(2*fsw));
    fsw_ini = fsw;

    %% Functions 
    function F = root2d(x,Vin,Pload,Lr,Cr,Lm,Vout,n,m)

    rr = (sqrt(Lr/Cr));
    t1 = x(1)*sqrt(Lr*Cr) ;
    fsw = 1/(2*(x(2)*sqrt((Lr+Lm)*Cr) + t1));

    A11 = 1;
    A12 = 0;
    A13 = -cos(x(1));
    A14 = (sin(x(1))/rr);

    A21 = 0;
    A22 = 1;
    A23 = -sin(x(1))*rr;
    A24 = -cos(x(1));

    A31 = -cos(x(2));
    A32 = (sin(x(2))/(rr*sqrt(m+1)));
    A33 = -1;
    A34 = 0 ;

    A41 = -(sin(x(2))*rr*sqrt(m+1));
    A42 = -cos(x(2));
    A43 = 0;
    A44 = -1;


    MA = [ A11 A12 A13 A14 ;
           A21 A22 A23 A24 ;
           A31 A32 A33 A34 ;
           A41 A42 A43 A44 ;];

    B1 = (-n*Vout+Vin)*sin(x(1))/rr;
    B2 = (Vin-n*Vout)*(-cos(x(1))+1);
    B3 =  Vin*sin(x(2))/(rr*sqrt(m+1));
    B4 = Vin*(-cos(x(2))+1);



    MB = [B1;
          B2;
          B3;
          B4;];


    MR = (MA\MB);

    I1 = MR(1); V1 = MR(2); I0 = MR(3); V0 = MR(4);
    P1 = (-Vin+n*Vout+V0)*Cr*(-1+cos(x(1))) + I0*sqrt(Lr*Cr)*sin(x(1));
    P2 = I1*sqrt(Cr*(Lr+Lm))*sin(x(2)) + (Vin-V1)*Cr*(1-cos(x(2)));
    % 2*fsw*Vin*((Vin-n*Vout-V0)*Cr*(1-cos(x(1))) + I0*sqrt(Lr*Cr)*sin(x(1)) + I1*sqrt(Cr*(Lr+Lm))*sin(x(2)) + (Vin-V1)*Cr*(1-cos(x(2))))
    % Vin*(((2*m)/(x(1)*n))/ ((2*m/x(1)) + cot(x(1)/2) - (sqrt(m+1)*tan(x(2)/2))))
    F(1) = 2*fsw*Vin*(P1+P2) - Pload;
    F(2) = (((2*m)/(x(1)*n))/ ((2*m/x(1)) + cot(x(1)/2) - (sqrt(m+1)*tan(x(2)/2)))) - Vout/Vin ;
    end
    
end





