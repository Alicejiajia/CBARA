function [CRB_dq] = cpu_CRB_dq(d_u,d_bs,Rcs,lamdac,kap,beam_gafc,Nr,sigma)
%计算位置估计的CRB
d=norm(d_u-d_bs);
beta = 7.3*10^(-12);
pathloss = Rcs*lamdac^2/((4*pi)^3*d^4); 
sens_chgain = pathloss * kap^4 * beam_gafc^4 / (Nr*sigma);
CRB_dq = beta / ((abs(sens_chgain))^2);

end

