clear all;
close all;
clc;
%%
%初始化参数
global Nt %发射天线
global Nr  %接收天线
global Q  %感知目标数
global K  %基站个数
global I  %通信用户数
global N  %时隙个数
global fc  %信号频率 
global Ptotal  %总功率
global Btotal  %总带宽
global Lmin  %最少选择基站数
global Lmax  %最大选择基站数
global Rcs  %雷达截面
global D_T  %时隙间隔
Rcs =1;
Nt=32;
Nr=32;
K=3;
Q=3;
I=1 ;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
fc=3e9 ;%3GHz
Ptotal= 30  ;%W 
Btotal= 65; %MHz 
Lmin = 2;  
Lmax = 3;
lamdac = 3*10^8/fc ;%波长
cnoise_pwr = 10^(-17.5) ;%-145dBm
kap_s = sqrt(Nt*Nr) ; %感知阵列增益因子
kap_c = sqrt(Nt*Nr) ; %通信阵列增益因子
beam_gafc = 1 ; %波束增益因子
D_T = 0.5 ;%s
Ts = D_T;  
N=30;
latency_error = 1;  % 同步误差系数
eta = 0.5;  % 通感权衡尺度因子
dist = 120;  % 基站间距
scale = 10; %用以保证目标函数中通信和感知指标的数量级一致

[base,BS] = base_init(dist);
init_pos = {
    [-dist/2, dist/4*sqrt(3),0,-dist/4*sqrt(3)+15]; 
     [66,12,66,84];
     [-72,0,72,0];
    };

state = cell(1,Q);
for q = 1:Q
    state{q} = single_state_init(N+1, Ts, init_pos{q}(1), init_pos{q}(2), init_pos{q}(3), init_pos{q}(4));
end

Xi_qn = zeros(4,Q,N+1);
reXi_qn = zeros(4,Q,N+1);%状态矩阵
for n=1:N+1
    for q = 1:Q
        reXi_qn(:,q,n) = state{q}(n, :).';
    end
end

Sigma1 = 1.5 ;%目标噪声方差
Sigma2 = 2;
Sigma3 = 1 ;

sigma = [Sigma1,Sigma2,Sigma3];

loca_qn = zeros(Q,2,N+1);
loca_cn = zeros(I,2,N+1);
for nn = 1:N+1
    for i = 1:I
         loca_cn(i,:,nn) = state{Q-I+i}(nn,1:2) ;
    end
    loca_qn(:,:,nn) = reXi_qn(1:2,:,nn).';
end

I2 = diag([1,1]);
T_matrx_F = [1 Ts;0 1];
F_Xi = kron(T_matrx_F,I2); %状态转移矩阵
T_matrx_Fei = [1/3*Ts^3,1/2*Ts^2; 1/2*Ts^2,Ts];

Fei = zeros(4,4,Q);%噪声方差矩阵
for q = 1:Q
     Fei(:,:,q) = kron(T_matrx_Fei, sigma(q) * I2);
end


%***********************初始化J*************************************
Jt = ones(4,4,Q,N); 
J0 = ones(4,4,Q);
J0i = kron(I2,[3200 16;16 6]); 
for q =  1:Q
    Jt(:,:,q,1)= J0i;
end
for qq = 1:Q
    J0(:,:,qq)= J0i;   
end

bt = zeros(K,Q,N);
bt(:,:,1)=Btotal/(Q*K);
b00 = bt(:,:,1);

%***********************优化*************************************
[pt, bt, F, Jt,u, Rate_n] = ...
    optim_jointUPB(N, I, Q, K, Btotal, Ptotal, reXi_qn(:,:,1:N+1), J0, base, Rcs, lamdac, fc, kap_s, kap_c, beam_gafc, Nt, Nr,...
    sigma, cnoise_pwr, Fei, F_Xi, b00, Lmin, Lmax, loca_cn(:,:,1:N), eta, 0.1, latency_error,scale);





