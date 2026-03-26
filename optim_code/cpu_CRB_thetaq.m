function [CRB_theta] = cpu_CRB_thetaq(d_u,d_bs,Rcs,lamdac,kap,beam_gafc,Nr,sigma)
%计算角度估计的CRB
d=norm(d_u-d_bs);
beta = 7.7*10^(-18);
pathloss = Rcs*lamdac^2/((4*pi)^3*d^4);
sens_chgain = pathloss * kap^4 * beam_gafc^4 / (Nr*sigma);
CRB_theta = beta / ((abs(sens_chgain))^2);

end

