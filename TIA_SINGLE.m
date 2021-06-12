clear
clc

Pload = 4000;


Lr = 19.196e-06;
Cr = 72.403e-9;
Lm = 23.035e-06;
Vin = 250;
Vout = 15.5;
n = 30;
fres = 1/(2*pi*sqrt(Lr*Cr));
% [~,col] = size(Cr);
% fsw = zeros(1,col);
% ILr_rms = zeros(1,col);
% Isec_rms = zeros(1,col);

% for i = 1:1:col
%     
%     if n*Vout > Vin
%     
%         [t1_ini,t2_ini,fsw_ini]=TIA_LLC_FB_BR(Pload,Lr(1,i),Cr(1,i),Lm(1,i),Vin,Vout,n);
%         t1_ini = 1.0*t1_ini * 1e6;
%         t2_ini = 1.0*t2_ini * 1e6;
%         fsw_ini = fsw_ini/1e5;
%         [fsw(1,i),ILr_rms(1,i),~,Isec_rms(1,i)] = TIA_LLC_FB_BR_3interval_FUNC(Pload,Lr(1,i),Cr(1,i),Lm(1,i),Vin,Vout,n,t1_ini,t2_ini,fsw_ini);
%     %     per_chnge_req = (abs(VLm_fswby2) - (n*Vout))/(n*Vout);
%     %     t1_ini = (1+per_chnge_req)*t1_ini;
%     %     t2_ini = (1+per_chnge_req)*t2_ini;
%     %     [fsw,ILr_rms,VLm_fswby2] = TIA_LLC_FB_BR_3interval_FUNC(Pload,Lr,Cr,Lm,Vin,Vout,n,t1_ini,t2_ini,fsw_ini);
% 
%     end
% 
% 
% 
%     if n*Vout < Vin 
%         [fsw(1,i),ILr_rms(1,i)] = TIA_LLC_FB(Pload,Lr(1,i),Cr(1,i),Lm(1,i),Vin,Vout,n);
%     end
% 
%     if n*Vout == Vin            
%         fsw(1,i) = 1/(2*pi*sqrt(Lr(1,i)*Cr(1,i)));               
%     end
%     
% end


if n*Vout > Vin
    
    [t1_ini,t2_ini,fsw_ini]=TIA_LLC_FB_BR(Pload,Lr,Cr,Lm,Vin,Vout,n);
    t1_ini = 1.1*t1_ini * 1e6;
    t2_ini = 1.1*t2_ini * 1e6;
    fsw_ini = fsw_ini/1e5;
    [fsw,ILr_rms,VLm_fswby2,Isec_rms] = TIA_LLC_FB_BR_3interval_FUNC(Pload,Lr,Cr,Lm,Vin,Vout,n,t1_ini,t2_ini,fsw_ini);
    per_chnge_req = (abs(VLm_fswby2) - (n*Vout))/(n*Vout);
%     t1_ini = (1+per_chnge_req)*t1_ini;
%     t2_ini = (1+per_chnge_req)*t2_ini;
%     [fsw,ILr_rms,VLm_fswby2] = TIA_LLC_FB_BR_3interval_FUNC(Pload,Lr,Cr,Lm,Vin,Vout,n,t1_ini,t2_ini,fsw_ini);

end



if n*Vout < Vin 
    [fsw,ILr_rms] = TIA_LLC_FB(Pload,Lr,Cr,Lm,Vin,Vout,n);
end

if n*Vout == Vin            
    fsw = 1/(2*pi*sqrt(Lr*Cr));               
end

