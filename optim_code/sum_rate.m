function [Rate_n] = sum_rate(I,K,loca_cn,b,p,Nt,Nr,base,fc,kap,cnoise_pwr)

 Rate_n  = 0;
  for cc = 1:I
      d = loca_cn(cc,:);
      for k = 1:K
              d_bs = base(k,:);
              bk = b(k,cc)*10^6;
              pk = p(k,cc);
              at = steering_vector(Nt,d,d_bs); 
              ar = steering_vector(Nr,d,d_bs); 
              paloss_c = 32.4+20*log(norm(d-d_bs)*10^(-3))+20*log(fc*10^(-6));
              paloss_c = 10^(-paloss_c/20);
              sc = (abs(kap*paloss_c*(ar')*ar*(at')*at))^2/cnoise_pwr ;     
              Rate_n = Rate_n+b(k,cc)*log(1+pk*sc*inv_pos(bk))/log(2);

      end

  end
end

