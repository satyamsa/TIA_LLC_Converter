function [fsw_req,ILr_rms,VLm_fswby2,Isec_rms] = TIA_LLC_FB_BR_3interval_FUNC(Pload,Lr,Cr,Lm,Vin,Vout,n,t1_ini,t2_ini,fsw_ini)

m = Lm/Lr;
fun = @(x)root2d(x,Vin,Pload,Lr,Cr,Lm,Vout,n,m);

    x0 = [t1_ini,t2_ini,fsw_ini]; %.85


%.75 .85
    options = optimoptions('fsolve'); 
    options.MaxIterations = 50000;
    options.MaxFunctionEvaluations = 10000;
    x = fsolve(fun,x0,options);
    [F,I1,V1,I2,V2,I0,V0] = root2d(x,Vin,Pload,Lr,Cr,Lm,Vout,n,m);
    % y = lsqnonlin(fun,x0,[4.70588 4.70588 0.85], [5 5 1],options)
    % [F1,~,~,~,~,~,~] = root2d(y,Vin,Pload,Lr,Cr,Lm,Vout,n,m)
    
    
    rr = (sqrt(Lr/Cr));
    t1 = x(1)*1e-6;
    t2 = x(2)*1e-6;
    fsw = x(3)*1e5;
    fsw_req = fsw;
    A = (t1/sqrt(Lr*Cr));
    T = (t2 - t1)/sqrt((Lr +Lm)*Cr);
    Ts = 1/fsw;
    B = ((Ts/2)-t2)/sqrt(Cr*(Lr+Lm));

    Ilm0 = I0;
    Ilm1 = I1; 
    %% Reconstrcuting the waveform of ILr,ILm,Vcr
    time = 0:1e-10:(1/(2*fsw));
    t_fl = 0:1e-10:(1/fsw);
    [~,time_col] = size(time);
    [~,t_fl_col] = size(t_fl);
    Vcr = zeros(1,t_fl_col);
    ILr = zeros(1,t_fl_col);
    ILm = zeros(1,t_fl_col);
    VLm = zeros(1,t_fl_col);
    for i = 1:1:time_col 
        if( time(1,i) <= t1)
            Vcr(1,i) = Vin - n*Vout + (I0*sqrt(Lr/Cr)*sin(time(1,i)/sqrt(Lr*Cr))) + (V0 - Vin + n*Vout)*cos(time(1,i)/sqrt(Lr*Cr)); 
            ILr(1,i) = I0*cos(time(1,i)/sqrt(Lr*Cr)) - (V0 - Vin + n*Vout)*sqrt(Cr/Lr)*sin(time(1,i)/sqrt(Lr*Cr));
            ILm(1,i) = (n*Vout*time(1,i)/Lm) + Ilm0;
            VLm(1,i) = n*Vout;
        end
        if(time(1,i)> t1 && time(1,i)<= t2)
            Vcr(1,i) = Vin  + (I1*sqrt((Lr+Lm)/Cr)*sin((time(1,i)-t1)/sqrt((Lr+Lm)*Cr))) + (V1 - Vin)*cos((time(1,i)-t1)/sqrt((Lr+Lm)*Cr));
            ILr(1,i) = I1*cos((time(1,i)-t1)/sqrt((Lr+Lm)*Cr)) - (V1 - Vin)*sqrt(Cr/(Lr+Lm))*sin((time(1,i)-t1)/sqrt((Lr + Lm)*Cr));
            ILm(1,i) = ILr(1,i);
            VLm(1,i) = -1*Lm*(I1*sin((time(1,i)-t1)/sqrt((Lr+Lm)*Cr)) + (V1 - Vin)*sqrt(Cr/(Lr+Lm))*cos((time(1,i)-t1)/sqrt((Lr + Lm)*Cr)))/(sqrt((Lr+Lm)*Cr));
        end
        if(time(1,i)> t2 && time(1,i)<=(1/(2*fsw)))
            Vcr(1,i) = -Vin  + (I2*sqrt((Lr+Lm)/Cr)*sin((time(1,i)-t2)/sqrt((Lr+Lm)*Cr))) + (V2 + Vin)*cos((time(1,i)-t2)/sqrt((Lr+Lm)*Cr));
            ILr(1,i) = I2*cos((time(1,i)-t2)/sqrt((Lr+Lm)*Cr)) - (V2 + Vin)*sqrt(Cr/(Lr+Lm))*sin((time(1,i)-t2)/sqrt((Lr + Lm)*Cr)); 
            ILm(1,i) = ILr(1,i);
            VLm(1,i) = -1*Lm*(I2*sin((time(1,i)-t2)/sqrt((Lr+Lm)*Cr)) + (V2 + Vin)*sqrt(Cr/(Lr+Lm))*cos((time(1,i)-t2)/sqrt((Lr + Lm)*Cr)))/(sqrt((Lr+Lm)*Cr));


        end

    end

    for i = (time_col + 1):1:t_fl_col
        ILr(1,i) = -1*ILr(1,i-time_col);
        ILm(1,i) = -1*ILm(1,i-time_col);
        Vcr(1,i) = -1*Vcr(1,i-time_col);
        VLm(1,i) = -1*VLm(1,i-time_col);
    end
    Isec = n*(ILr - ILm);
    Isec_rms = rms(Isec);
    figure
    plot(t_fl,ILr);
    hold on;
    plot(t_fl,ILm);
    figure
    plot(t_fl,Vcr)
    hold on 
    plot(t_fl,VLm);
    ILr_rms = rms(ILr);
    VLm_fswby2 = VLm(1,time_col);
    figure 
    plot(t_fl,Isec);
    figure
    plot(Vcr,ILr);
    
    %% Functions 
    function [F,I1,V1,I2,V2,I0,V0] = root2d(x,Vin,Pload,Lr,Cr,Lm,Vout,n,m)

    rr = (sqrt(Lr/Cr));
    t1 = x(1)*1e-6;
    t2 = x(2)*1e-6;
    fsw = x(3)*1e5;


    A = (t1/sqrt(Lr*Cr));
    T = (t2-t1)/sqrt((Lr +Lm)*Cr);
    Ts = 1/fsw;
    B = (Ts/2-t2)/sqrt(Cr*(Lr+Lm));


    A11 = 1;
    A12 = 0;
    A13 = 0;
    A14 = 0;
    A15 = -cos(A);
    A16 = (sin(A)/rr);

    A21 = 0;
    A22 = 1;
    A23 = 0;
    A24 = 0;
    A25 = -rr*sin(A);
    A26 = -cos(A);

    A31 = -cos(T);
    A32 = (sin(T)/(rr*sqrt(m+1)));
    A33 = 1 ;
    A34 = 0 ;
    A35 = 0 ;
    A36 = 0 ;

    A41 = -(sin(T)*rr*sqrt(m+1));
    A42 = -cos(T);
    A43 = 0;
    A44 = 1;
    A45 = 0;
    A46 = 0;

    A51 = 0;
    A52 = 0;
    A53 = -cos(B);
    A54 = sin(B)/(rr*sqrt(m+1));
    A55 = -1;
    A56 = 0;

    A61 = 0;
    A62 = 0;
    A63 = -sin(B)*rr*sqrt(m+1);
    A64 = -cos(B);
    A65 = 0;
    A66 = -1;


    MA = [ A11 A12 A13 A14 A15 A16;
           A21 A22 A23 A24 A25 A26;
           A31 A32 A33 A34 A35 A36;
           A41 A42 A43 A44 A45 A46;
           A51 A52 A53 A54 A55 A56;
           A61 A62 A63 A64 A65 A66;];

    B1 = -((n*Vout-Vin)*sin(A)/rr);
    B2 = (Vin-n*Vout)*(1 - cos(A));
    B3 =  Vin*sin(T)/(rr*sqrt(m+1));
    B4 =  Vin*(1-cos(T));
    B5 = -(Vin*sin(B)/(rr*sqrt(m+1)));
    B6 =  Vin*(cos(B)-1);

    MB = [B1;
          B2;
          B3;
          B4;
          B5;
          B6;];

    MR = (MA\MB);

    I1 = MR(1); V1 = MR(2); I2 = MR(3); V2 = MR(4); I0 = MR(5);V0 = MR(6);

    P1 = (-Vin+(n*Vout)+V0)*Cr*(cos(A)-1) + I0*sqrt(Lr*Cr)*sin(A);
    P2 = I1*sqrt(Cr*(Lr+Lm))*sin(T) + (V1-Vin)*Cr*(cos(T)-1);
    P3 = I2*sqrt(Cr*(Lr+Lm))*sin(B) + (V2+Vin)*Cr*(cos(B)-1);

    % 2*fsw*Vin*(P1+P2 - P3)
    % I1 - I0*cos(A) + (V0 - Vin + (n*Vout))*sin(A)*sqrt(Cr/Lr)
    % (Lm*((I1 - I0)/(t1 - 0))) - n*Vout
    F(1) = 2*fsw*Vin*(P1+P2-P3) - Pload;
    F(2) = I1 - I0*cos(A) + (V0 - Vin + (n*Vout))*sin(A)*sqrt(Cr/Lr);
    F(3) = (Lm*((I1 - I0)/(t1 - 0))) - n*Vout ;
    end
end