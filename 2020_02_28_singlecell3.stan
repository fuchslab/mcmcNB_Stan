functions{
  vector rowSum(matrix Mat){
    vector[rows(Mat)] res =  Mat * rep_vector(1,cols(Mat));
    return res;
  }
  
  vector colSum(matrix Mat){
    vector[cols(Mat)] res =  (rep_row_vector(1,rows(Mat)) * Mat)';
    return res;
  }
  
  
  int[] compto0(int[] x){
    int res[num_elements(x)];
    for(i in 1:num_elements(x)){
      res[i] = (x[i] != 0) ;
    }
    return res;
  }
  
  vector reduceAlpha(int[] x, vector alpha){
    int x_new[num_elements(x)] = compto0(x);
    vector[sum(x_new)] alpha_red ;
    int j = 1;
    for (i in 1:num_elements(x)){
      if(x_new[i] == 1){
      alpha_red[j]= x[i] .* alpha[i];
      j += 1;
      }
    }
    return alpha_red;
  }
  
    vector reduceBeta(int[] x, vector beta){
    int x_new[num_elements(x)] = compto0(x);
    vector[sum(x_new)] beta_red ;
    int j = 1;
    for (i in 1:num_elements(x)){
      if(x_new[i] == 1){
      beta_red[j]=  beta[i];
      j += 1;
      }
    }
    return beta_red;
  }
  
  int[,] comb_summands(int n, int k); 

  int[,] comb_summands(int n, int k){
    int length_res = choose(n+k-1,n);
    int res[length_res,k];
    if(k==1){
      res[length_res] = {n};
    }else{
      int j = 1;  
      for(i in 0:n){
        int length_input_help = choose(n-i+k-2,n-i);  
        int help_recursion[length_input_help,(k-1)] ;
        help_recursion = comb_summands(n-i,k-1);
        for(i2 in 1:length_input_help){
          res[j+i2-1] = append_array({i},help_recursion[i2]);
        }
      j +=length_input_help;
      }
    }
  return res;
  }

 vector Q_param_calc(vector alpha, vector beta){
   int Summands = num_elements(alpha);
   vector[Summands] prob_param;
   vector[Summands] q_param;
   vector[Summands] Q_param;
   real prob_max;
   
   prob_param = beta ./(beta + 1);
   q_param = 1-prob_param;
   prob_max = max(prob_param);
   Q_param=((1 - prob_max) * prob_param)./(q_param * prob_max);
   return Q_param;
 } 
  
  
 real l_R_calc(vector alpha, vector Q_param){
   real l_R;
   
   l_R = sum(alpha .* log(Q_param)); // eigentlich modulo anstelle von mal
   return l_R;
 }
 
 vector lse_AE_vector_calc(vector alpha, vector Q_param, int infi_start, int infi_end){
   int Summands = num_elements(alpha);
   int infi_help = infi_end - infi_start+1;
   int i_counter_2;
   vector[infi_help] lse_AE_vector;
   matrix[infi_help, Summands] A_matrix;
   matrix[infi_help, Summands] E_matrix;
   matrix[infi_help,Summands] E_matrix_2;
   matrix[infi_help,Summands] AE_matrix;
    
   A_matrix = rep_matrix(to_row_vector(alpha), infi_help);
   E_matrix = rep_matrix(to_row_vector(1-Q_param), infi_help);
   
   i_counter_2 = infi_start;
   for(i_counter in 1:infi_help){
      
        for(i_counter_3 in 1:cols(E_matrix_2)){
          E_matrix_2[i_counter,i_counter_3] = pow(E_matrix[i_counter, i_counter_3],i_counter_2); //every element int row i is taken to power i
        }
        i_counter_2 += 1;
      }
        
      AE_matrix = A_matrix .* E_matrix_2;
      
    
      for(i_counter in 1:infi_help){
        lse_AE_vector[i_counter] = log(sum(AE_matrix[i_counter]));  //Eigene Funktion? aus alpha und Q-param
      }
      return lse_AE_vector;
   
 }
 
 vector l_delta_calc(vector alpha, vector beta, vector lse_AE_matrix, vector l_delta_before, int infi_start, int infi_end, int infi){
      vector[infi_end] l_delta;
      vector[infi_end] l_delta_rev;
      real max_beta;
      l_delta = l_delta_before[1:infi_end];
      
      if(infi_start < 3){
        l_delta_rev[infi_end] = 0;
        l_delta_rev[infi_end-1] = lse_AE_matrix[1];} 
      else{
        for(kk_counter in 1 : (infi_start-1)) l_delta_rev[(infi_end+1-kk_counter)] = l_delta_before[kk_counter];
      }
      
      max_beta = max(beta ./(beta + 1))/(1-max(beta ./(beta + 1)));  

     // delta_zero in the Furmann Paper is delta_1 here, pay attention every think is shifted, I hopw it is right  

      for(k_counter in max(2,(infi_start-1)):max(3,(infi_end-1))){
         l_delta_rev[infi_end-k_counter]= -log(k_counter) + log_sum_exp(lse_AE_matrix[1 : (k_counter)] + tail(l_delta_rev,k_counter));
      }
      for(kk_counter in 1 : infi_end){
         l_delta[(infi_end+1-kk_counter)] = l_delta_rev[kk_counter];
      }
      return l_delta;
 }

   
 real neg_binomial_sum_lpmf(int y, vector alpha, vector beta, int infi){

    int Summands = num_elements(alpha);
    int i_start = 1;
    
    vector[Summands] Q_param;
    vector[infi] l_delta;
    vector[infi] l_delta_help;
    vector[infi] lse_AE_vector;
    vector[infi] l_NB_infi;
    int infi_red[10];
    int i_end;
    int infi_start;
    int infi_end;
    int infi_final;
    
    //real prob_max;
    real l_R;
    real max_beta;
    real res;
    //matrix[infi,Summands] A_matrix;
    //matrix[infi,Summands] E_matrix;
    //matrix[infi,Summands] E_matrix_2;
    //matrix[infi,Summands] AE_matrix;
    
    max_beta = max(beta ./(beta + 1))/(1-max(beta ./(beta + 1)));  
      Q_param = Q_param_calc(alpha,beta);
      l_R = l_R_calc(alpha,Q_param); 
      
     // nicht alle bis infi berechnen, nur mit check wie bei altem: Immer nur bis dahin befüllen, ggf mehr! for mit if schleife? keine while!
      
      // 10,50,100,500,1000,5000,10000,50000,100000 so viele Lagen! -> viele Ram!
      // round(c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5,1))
      infi_red = {0,10,50,100,500,1000,5000,10000,50000,infi} ;
      //for (i in 1:9)if (infi_red[i]==0) infi_red[i]=1;
      
      //counter = min(counter + 3 * (l_delta[infi-k_counter] > l_delta[infi-k_counter+1]), infi-1)
      
      // sucher ersten eintrag in infi_red >1
      l_delta[1]=0;
      l_delta[2]=0;
      infi_end = 2;
      for (i_infi in 2:9) {
        
        if((infi_red[i_infi])<=infi){
         
          i_start = i_infi-1; 
          i_end = i_infi;
      
          infi_start = infi_red[i_start]+1;
          infi_end = infi_red[i_end];
          if(infi_start == 1 || l_delta[infi_start-1]>l_delta[infi_start-2|| l_delta[infi_start-1]>0]){
            lse_AE_vector[infi_start:infi_end] = lse_AE_vector_calc(alpha, Q_param, infi_start, infi_end);
            l_delta_help = l_delta;
            l_delta[1:infi_end] = l_delta_calc(alpha, beta, lse_AE_vector, l_delta_help, infi_start, infi_end, infi);
            infi_final = infi_end;
          }
        }
      }
      
      
      for(i2 in 1:infi_final){
            l_NB_infi[i2] = neg_binomial_lpmf(y| sum(alpha) + i2 - 1, max_beta);
          }

      res = l_R + log_sum_exp(l_delta[1:infi_final]+l_NB_infi[1:infi_final]) ;
      return(res);
  }






  real neg_binomial_sum_all_lpmf(int[] y, vector alpha, vector beta, vector theta, int infi, int N_cell, int Pop){
    real res = 0;
    int length_possi = choose(N_cell+Pop-1,N_cell);
    int possi[length_possi,Pop] = comb_summands(N_cell, Pop);
    for(y_count in 1:size(y)){          
      vector[length_possi] lps ;
      for(i_counter in 1:length_possi){     // go through all possible combinations N cells from Pop Populations
        int possi_help[Pop] = possi[i_counter] ;
        if(max(possi_help)==N_cell){ // if in the combination only one population is left, this is directly a NB with alpha = N*alpha and beta
          real alpha_help =  sum(alpha .* to_vector(possi_help));
          real beta_help = sum(beta .* to_vector(compto0(possi_help)));
          lps[i_counter] = multinomial_lpmf(possi_help | theta);
          lps[i_counter] += neg_binomial_lpmf(y[y_count] | alpha_help, beta_help);
        
        } else{
          int Pop_red = sum(compto0(possi_help));
          vector[Pop_red] alpha_help =  reduceAlpha(possi_help, alpha);
          vector[Pop_red] beta_help =  reduceBeta(possi_help, beta);
          lps[i_counter] = multinomial_lpmf(possi_help | theta);
          lps[i_counter] += neg_binomial_sum_lpmf(y[y_count] | alpha_help, beta_help, infi);
        }
    }
    res += log_sum_exp(lps);
  }
  return res;
}


 real neg_binomial2_lpmf(int[] y, real alpha, real beta){
   real res = neg_binomial_lpmf(y | alpha, beta);
   return res;
 }


 real multinomial_2_lpmf(int[] possi_help, vector theta){
   real res = multinomial_lpmf(possi_help | theta);
   return res;
 }
 

}



data {
  int<lower=1> N_cell; // # of cells per sample
  int<lower=1> Pop;        // Number of different cell populations
  int<lower=1> N_sample;     // # samples
  int<lower=0> y[N_sample]; // # of failures in trial n
  int<lower=1> infi;         //needed for Furman sum of NB
}

parameters{
  vector<lower=0>[Pop] alpha;  // shape = # of successes until stopping
  vector<lower=0>[Pop] beta; // inverse scale = p / (1-p)
  simplex[Pop] theta; // mixing proportions (sum up to 1)
//  real<lower=0> alpha;  // shape = # of successes until stopping
//  real<lower=0> beta; // inverse scale = p / (1-p)
//  real<lower=0> theta; // mixing proportions (sum up to 1)

 }


model{
  vector[Pop] alpha_help;
  vector[Pop] log_theta = log(theta);  // cache log calculation
  
  //for(p in 1:Pop) alpha[p] ~ lognormal(20, 100);
  // for(p in 1:Pop) beta[p] ~ lognormal(20, 100);

    if(N_cell==1 || Pop == 1){
      if(N_cell==1){
        alpha_help = alpha;
      } else if(Pop == 1){
        alpha_help =  N_cell * alpha;
      }
      for (n in 1:N_sample) {    
        vector[Pop] lps = log_theta;
        for (k in 1:Pop){
          lps[k] += neg_binomial_lpmf(y[n] | alpha_help[k], beta[k]);
        }
        target += log_sum_exp(lps);
      }
    } 
    else {

        y ~ neg_binomial_sum_all_lpmf(alpha, beta, theta, infi, N_cell, Pop);

    }

}
  



