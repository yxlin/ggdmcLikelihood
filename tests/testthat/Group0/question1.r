CDM 8.3-14 (2025-07-13 14:03:01) 

Call:
CDM::gdina(data = data.frame(wide[, -1]), q.matrix = Q, rule = "RRUM")

Date of Analysis: 2025-10-01 15:41:25.720222 
Time difference of 0.991 secs
Computation Time: 0.991 

Generalized DINA Model 

Number of iterations = 267
Iteration with minimal deviance = 267 

Estimation method: ML
Optimizer: optim
Monotonicity constraints: FALSE
Number of items at boundary monotonicity constraint: NA

Parameter regularization: FALSE

Deviance = 46764  | Log likelihood = -23382 

Number of persons = 10000 
Number of groups = 1 
Number of items = 5 
Number of estimated parameters = 14 
Number of estimated item parameters = 11 
Number of estimated skill class parameters = 3 ( 4 latent skill classes)

AIC = 46792  | penalty = 28    | AIC = -2*LL + 2*p  
BIC = 46893  | penalty = 129    | BIC = -2*LL + log(n)*p   
CAIC = 46907  | penalty = 143    | CAIC = -2*LL + [log(n)+1]*p  (consistent AIC)   

-----------------------------------------------------------------
Used Q-matrix 

   Algebra Geometry
E1       1        0
E2       0        1
E3       1        1
E4       1        0
E5       0        1

-----------------------------------------------------------------

Item Parameter Estimates 

   link item itemno partype rule    est     se partype.attr
1   log   E1      1       0 ACDM -2.252 0.0327             
2   log   E1      1       1 ACDM  1.818 0.0365      Algebra
3   log   E2      2       0 ACDM -1.323 0.0188             
4   log   E2      2       1 ACDM  0.861 0.0249     Geometry
5   log   E3      3       0 ACDM -3.373 0.0625             
6   log   E3      3       1 ACDM  1.009 0.0837      Algebra
7   log   E3      3       2 ACDM  0.575 0.0875     Geometry
8   log   E4      4       0 ACDM -2.163 0.0311             
9   log   E4      4       1 ACDM  0.729 0.0505      Algebra
10  log   E5      5       0 ACDM -1.523 0.0214             
11  log   E5      5       1 ACDM  0.518 0.0353     Geometry

RMSD (RMSEA) Item Fit Statistics
   E1    E2    E3    E4    E5 
0.001 0.001 0.002 0.002 0.002 

Mean of RMSEA item fit: 0.001 

****
RRUM Parametrization
      pi r_Algebra r_Geometry
E1 0.648     0.162           
E2 0.630                0.423
E3 0.167     0.365      0.563
E4 0.238     0.482           
E5 0.366                0.596

#                        pi1   pi2   pi3 pi_item1 pi_item2 pi_item3 pi_item4
#                      0.236 0.315 0.449    0.227     0.34    0.302     0.38
# name_sorted_p_vector 0.100 0.200 0.300    0.100     0.20    0.300     0.40
#                      pi_item5   r11   r22   r31   r32   r41   r52
#                         0.427 0.263 0.556 0.246 0.448 0.379 0.575
# name_sorted_p_vector    0.500 0.100 0.400 0.300 0.600 0.400 0.660

# Skill Pattern Probabilities 

#     00     10     01     11 
# 0.6166 0.1627 0.1813 0.0394 
-----------------------------------------------------------------
Model Implied Conditional Item Probabilities 

   item rule        nessskill itemno skillcomb   prob
1    E1 ACDM          Algebra      1        A0 0.1052
2    E1 ACDM          Algebra      1        A1 0.6476
3    E2 ACDM         Geometry      2        A0 0.2664
4    E2 ACDM         Geometry      2        A1 0.6303
5    E3 ACDM Algebra-Geometry      3       A00 0.0343
6    E3 ACDM Algebra-Geometry      3       A10 0.0941
7    E3 ACDM Algebra-Geometry      3       A01 0.0610
8    E3 ACDM Algebra-Geometry      3       A11 0.1672
9    E4 ACDM          Algebra      4        A0 0.1150
10   E4 ACDM          Algebra      4        A1 0.2383
11   E5 ACDM         Geometry      5        A0 0.2181
12   E5 ACDM         Geometry      5        A1 0.3662
-----------------------------------------------------------------

Skill Probabilities 

         skill.prob0 skill.prob1
Algebra        0.798       0.202
Geometry       0.779       0.221
-----------------------------------------------------------------

Polychoric Correlations 

Group 1
         Algebra Geometry
Algebra    1.000   -0.064
Geometry  -0.064    1.000

 -----------------------------------------------------------------

