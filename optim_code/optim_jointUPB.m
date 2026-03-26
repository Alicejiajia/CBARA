function [pn, bn, F, J_nN, u_n, rate, obj_val_history] = optim_jointUPB(N_slots, I, Q, K, Btotal, Ptotal, reXi_qn, J0, base, Rcs, lamdac, fc, kap_s, kap_c, beam_gafc, Nt, Nr, sigma, cnoise_pwr, Fei, F_Xi, b00, Lmin, Lmax, loca_cn, eta, alh, latency_error,scale)

% 初始化全时隙矩阵
pn = zeros(K,Q,N_slots);
pun = zeros(K,Q,N_slots);
bn = zeros(K,Q,N_slots);

J = zeros(4,4,Q,N_slots+1);
J(:,:,:,1) = J0;
F = zeros(1,N_slots);
F_uni = zeros(1,N_slots);

u_n = zeros(K,Q,N_slots);
II = [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
rate = zeros(1,N_slots);

obj_val_history = cell(1, N_slots); 

for n = 1:N_slots

  Xi_qn(:,:,n) = F_Xi*reXi_qn(:,:,n) ;
  loca_qn(:,:,n) = Xi_qn(1:2,:,n).';
  inv_Jq = zeros(4,4,Q);
  jnl = J(:,:,:,n);
  
  %***************************************************************************优化基站选择*************************************
  cvx_begin sdp quiet 
  variable D(4,4,Q) symmetric
  variable  pu(K,Q) 
  expression J_q(4,4,Q)

  sum_D = 0;
  for i =1:Q
      E = inv(Fei(:,:,i) + F_Xi*jnl(:,:,i)*F_Xi.');
      J_k = 0;
      for k = 1:K  
            [CRB_d_k] = cpu_CRB_dq(loca_qn(i,:,n),base(k,:),Rcs,lamdac,kap_s,beam_gafc,Nr,sigma(i));
            [CRB_theta_k] = cpu_CRB_thetaq(loca_qn(i,:,n),base(k,:),Rcs,lamdac,kap_s,beam_gafc,Nr,sigma(i));
            [CRB_v_k] = cpu_CRB_v(loca_qn(i,:,n),base(k,:),Rcs,lamdac,kap_s,beam_gafc,Nr,sigma(i));
            
            Ap_k = diag([b00(k,i)*10^6*(latency_error * CRB_d_k)^(-1), (latency_error * CRB_theta_k)^(-1), (latency_error * CRB_v_k)^(-1)]);
            Jacobihq_k = cpu_Jacobih(Xi_qn(:,i,n).',base(k,:));

             H_k = Jacobihq_k;
             V_k = H_k.'*Ap_k*H_k;
             J_k = J_k + pu(k,i)*V_k;
      end
      J_q(:,:,i) =  E + J_k;
      sum_D = sum_D + trace(D(:,:,i));
  end   

 [Rate_nu] = sum_rate(I,K,loca_cn(:,:,n),b00(:,3),pu(:,3),Nt,Nr,base,fc,kap_c,cnoise_pwr);

  fuc = eta*sum_D -(1-eta)*Rate_nu/scale;
  minimize (fuc) 

  subject to
  for q = 1:Q
     [D(:,:,q) II;II J_q(:,:,q)] == semidefinite(8);
  end

  for k = 1:K
      sum_p_k = 0;
      for q = 1:Q
          pu(k,q)>= 0;  
          pu(k,q)<=0.85*Ptotal;
          sum_p_k = sum_p_k + pu(k,q);    
      end
      sum_p_k <= Ptotal;
  end
  cvx_end
  
  pun(:,:,n) = pu;
  F_n = 0;
  for ii = 1:Q
      inv_Jq(:,:,ii) = inv(J_q(:,:,ii));
      F_n = F_n + trace(inv_Jq(:,:,ii));
  end
  F_uni(1,n) = F_n;
  
%*********************************************************************************优化u******************************************************************************
  u = zeros(K,Q);
  pu(:,:) = pu(:,:)/Ptotal;
  p_lim = pu;

  for qq = 1:Q
      for lim = 1:Lmin 
            [~, idx] = max(p_lim(:,qq));
            p_lim(idx,qq) = 0;
            u(idx,qq) = 1;
      end
      for lim = 1:(Lmax-Lmin)
          [~, idx] = max(p_lim(:,qq));
          if p_lim(idx,qq)>=alh
              u(idx,qq) = 1;
          end
          p_lim(idx,qq) = 0;
      end
  end

 %*******************************************再次优化p和b (交替优化算法 AO)******************************************************************************
 b_it = b00;
 epsilon = 1e-3;     
 max_iter = 10;      
 it = 0;         
 converged = false;  
 prev_obj_val = inf; 
 
 current_slot_history = zeros(1, max_iter);
 
 while ~converged && it < max_iter
     it = it + 1;
     
     % ==================== 优化 P ====================
     cvx_begin sdp quiet
     variable D(4,4,Q) symmetric
     variable  p(K,Q) 

     expression J_q(4,4,Q)
     jnl = J(:,:,:,n);
     F_n = 0;
     sum_D = 0;

     for i =1:Q
         E = inv(Fei(:,:,i) + F_Xi*jnl(:,:,i)*F_Xi.');
         J_k = 0;
         for k = 1:K  
                [CRB_d_k] = cpu_CRB_dq(loca_qn(i,:,n),base(k,:),Rcs,lamdac,kap_s,beam_gafc,Nr,sigma(i));
                [CRB_theta_k] = cpu_CRB_thetaq(loca_qn(i,:,n),base(k,:),Rcs,lamdac,kap_s,beam_gafc,Nr,sigma(i));
                [CRB_v_k] = cpu_CRB_v(loca_qn(i,:,n),base(k,:),Rcs,lamdac,kap_s,beam_gafc,Nr,sigma(i));
                Ap_k = diag([b_it(k,i)*10^6*(latency_error * CRB_d_k)^(-1), (latency_error * CRB_theta_k)^(-1), (latency_error * CRB_v_k)^(-1)]);
                Jacobihq_k = cpu_Jacobih(Xi_qn(:,i,n).',base(k,:));

                 H_k = Jacobihq_k;
                 V_k = H_k.'*Ap_k*H_k;
                 J_k = J_k + u(k,i)*p(k,i)*V_k;
         end
         J_q(:,:,i) =  E + J_k;
         sum_D = sum_D + trace(D(:,:,i));
     end    

     [Rate_np] = sum_rate_u(I,K,loca_cn(:,:,n),b_it(:,3),p(:,3),Nt,Nr,base,fc,kap_c,cnoise_pwr,u(:,3));

      fuc = eta*sum_D -(1-eta)*Rate_np/scale;
      minimize (fuc) 

      subject to
      for q = 1:Q
         [D(:,:,q) II;II J_q(:,:,q)] == semidefinite(8);
      end

      for k = 1:K
          sum_p_k = 0;
          for q = 1:Q
              if u(k,q) == 1
                   p(k,q)>=0.05*Ptotal;
                   p(k,q)<=0.85*Ptotal;
                   sum_p_k = sum_p_k + u(k,q)*p(k,q);    
              else
                  p(k,q) == 0;
              end
          end
          sum_p_k <= Ptotal;
      end
      cvx_end
      
      pn(:,:,n) = p;
      F_n = 0;
      for ii = 1:Q
          inv_Jq(:,:,ii) = inv(J_q(:,:,ii));
          F_n = F_n + trace(inv_Jq(:,:,ii));
      end
      F(1,n) = F_n;
      
    % ==================== 优化 B ====================
     cvx_begin   sdp quiet
     variable  b(K,Q)
     expression J_qb(4,4,Q)
     variable D_b(4,4,Q) symmetric
     variable x nonnegative
     
     sum_D_b = 0;
     
     for i =1:Q
         E = inv(Fei(:,:,i) + F_Xi*jnl(:,:,i)*F_Xi.');
         J_k = 0;
         for k = 1:K   
                [CRB_d_k] = cpu_CRB_dq(loca_qn(i,:,n),base(k,:),Rcs,lamdac,kap_s,beam_gafc,Nr,sigma(i));
                [CRB_theta_k] = cpu_CRB_thetaq(loca_qn(i,:,n),base(k,:),Rcs,lamdac,kap_s,beam_gafc,Nr,sigma(i));
                [CRB_v_k] = cpu_CRB_v(loca_qn(i,:,n),base(k,:),Rcs,lamdac,kap_s,beam_gafc,Nr,sigma(i));

                Ap_k = diag([b(k,i)*10^6*(latency_error * CRB_d_k)^(-1), (latency_error * CRB_theta_k)^(-1), (latency_error * CRB_v_k)^(-1)]);
                Jacobihq_k = cpu_Jacobih(Xi_qn(:,i,n).',base(k,:));
                 H_k = Jacobihq_k;
                 V_k = H_k.'*Ap_k*H_k;
                 J_k = J_k + u(k,i)*p(k,i)*V_k;
         end
         J_qb(:,:,i) =  E + J_k;
         sum_D_b = sum_D_b + trace(D_b(:,:,i));
     end

     [Rate_nb] = bsum_rate(I,K,loca_cn(:,:,n),b(:,3),p(:,3),Nt,Nr,base,fc,kap_c,cnoise_pwr,u(:,3));
     Rate_nb = Rate_nb/10^6;

     fuc = (1-eta)*Rate_nb/scale-eta*sum_D_b;
     maximize(fuc)
     
      subject to
       for q = 1:Q
         [D_b(:,:,q) II;II J_qb(:,:,q)] == semidefinite(8); %半正定约束
      end

      sum_b = 0;
      for k = 1:K
          for q = 1:Q
              if u(k,q) == 1
                  b(k,q)>= 0.05*Btotal; 
                  b(k,q)<= 0.85*Btotal;
                  sum_b = sum_b + u(k,q)*b(k,q) ;   
              else
                  b(k,q) == 0;
              end
          end
      end
      sum_b <= Btotal;
      cvx_end
      
      % ==================== 检查收敛性并记录 ====================
      current_obj_val = -cvx_optval;
      
      current_slot_history(it) = current_obj_val;
      
      if abs(current_obj_val - prev_obj_val) <= epsilon
          converged = true;
      end
      
      prev_obj_val = current_obj_val;
      b_it = b;
 end  

  obj_val_history{n} = current_slot_history(1:it);
 
  F_n = 0;
  for ii = 1:Q
      inv_Jq(:,:,ii) = inv(J_qb(:,:,ii));
      F_n = F_n + trace(inv_Jq(:,:,ii));
  end
  
  F(1,n) = F_n;
  bn(:,:,n)=b;
  u_n(:,:,n) = u;
  J(:,:,:,n+1)=inv_Jq;
  rate(1,n) = Rate_nb;
end

J_nN = J(:,:,:,2:N_slots+1);
end