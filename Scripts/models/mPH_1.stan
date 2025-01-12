data {
  int<lower=0> N; // Número de observaciones
  int<lower=0> N_valid; // Número de ensayos válidos
  int<lower=1> J; // Número de sujetos
  array[N] int<lower=1, upper=J> subj; // ID del sujeto por ensayo
  array[N] real trial; // ID del sujeto por ensayo
  array[N] int<lower=1, upper=2> Singleton; // Tipo de distractor (1 o 2)
  array[N] real reward; // Recompensa por ensayo
  array[N] int<lower=0, upper=1> valid; // Ensayos válidos
  array[N_valid] int<lower=0, upper=1> Omission; // 1 = miró distractor, 0 = miró objetivo
}

parameters {
  // Hiperparámetros
  real mu_intercept; 
  real mu_beta;
  real mu_beta2;
  real mu_lr;   
  real mu_gamma;        
  real alpha_0; // Nivel inicial de atención (de 0 a 1)
  
  vector<lower=0>[6] sigma; 

  // Parámetros no centrados
  matrix[6, J] z; 
}

transformed parameters {
  // Parámetros transformados
  matrix[6, J] beta_j;
  for (j in 1:J) {
    beta_j[1, j] = mu_intercept + sigma[1] * z[1, j]; // Intercepto
    beta_j[2, j] = inv_logit(mu_beta + sigma[2] * z[2, j]) * 10; // Pendiente
    beta_j[3, j] = mu_beta2 + sigma[3] * z[3, j]; // Intercepto
    beta_j[4, j] = inv_logit(mu_lr + sigma[4] * z[4, j]); // Tasa de aprendizaje [0, 1]
    beta_j[5, j] = inv_logit(alpha_0 + sigma[5] * z[5, j]); // Tasa de aprendizaje [0, 1]
    beta_j[6, j] = inv_logit(mu_gamma + sigma[6] * z[6, j]); // Tasa de aprendizaje [0, 1]
  }
}

model {
  // Priors para hiperparámetros
  mu_intercept ~ std_normal();
  mu_beta ~ std_normal();
  mu_beta2 ~ std_normal();
  mu_lr ~ std_normal();
  mu_gamma ~ std_normal();
  alpha_0 ~ std_normal();

  // Prior para el nivel inicial de atención

  // Priors para sigma
  sigma ~ std_normal();

  // Priors para parametrización no centrada
  to_vector(z) ~ std_normal();
  
  // Inicialización
  array[J] vector[2] v;     // Value
  array[J] vector[2] alpha; // Atención
  int counter = 1; // counter to index valid trials
  
  for (j in 1:J) {
    v[j] = rep_vector(0.0, 2);
    alpha[j] = rep_vector(beta_j[5, j], 2); // Usar el nivel inicial de atención
  }

  // Cálculo de latentes para cada ensayo
  for (i in 1:N) {
    int subject = subj[i];
    int s = Singleton[i];

    // Guardar alpha para calcular el likelihood
    if (valid[i] == 1) {
      Omission[counter] ~ bernoulli_logit(beta_j[1, subject] + beta_j[2, subject] * alpha[subject, s] + beta_j[3, subject] * trial[i]);
      counter += 1;
    }

    // PE y actualización de V
    real pe = reward[i] - v[subject, s];
    v[subject, s] += beta_j[4, subject] * alpha[subject, s] * pe;

    // Actualización de alpha
    alpha[subject, s] = beta_j[6, subject] * abs(pe) + (1 - beta_j[6, subject]) * alpha[subject, s]; 
  }
}

generated quantities {
  // Vector para guardar
  vector[N_valid] log_lik;
  matrix[N, 2] y_rep; // Posterior predictions

  { // local block
    array[J] vector[2] v;
    array[J] vector[2] alpha; 
    int counter = 1;

    for (j in 1:J) {
      v[j] = rep_vector(0.0, 2);
      alpha[j] = rep_vector(beta_j[5, j], 2); // Usar el nivel inicial de atención
    }

    for (i in 1:N) {
      int subject = subj[i];
      int s = Singleton[i];
      
      // Guardar alpha para los ensayos válidos
      if (valid[i] == 1) {
        log_lik[counter] = bernoulli_logit_lpmf(Omission[counter] | beta_j[1, subject] + beta_j[2, subject] * alpha[subject, s] + beta_j[3, subject] * trial[i]);
        counter += 1;
      }
      y_rep[i, 1] = bernoulli_logit_rng(beta_j[1, subject] + beta_j[2, subject] * alpha[subject, 1] + beta_j[3, subject] * trial[i]);
      y_rep[i, 2] = bernoulli_logit_rng(beta_j[1, subject] + beta_j[2, subject] * alpha[subject, 2] + beta_j[3, subject] * trial[i]);


      // PE y actualización de V
      real pe = reward[i] - v[subject, s];
      v[subject, s] += beta_j[4, subject] * alpha[subject, s] * pe;

      // Actualización de alpha
      alpha[subject, s] = beta_j[6, subject] * abs(pe) + (1 - beta_j[6, subject]) * alpha[subject, s]; 
    }
  }
}
