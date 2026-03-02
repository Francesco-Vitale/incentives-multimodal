# Sets
################################################################################################
set OD; # All ODs
set OD_O within OD; # All original ODs
set OD_V; # All added ODs (virtual ODs)
set OD_M; # All ODs that can be multimodal (to be split into virtual ODs)
set OD_C; # All ODs carrying cars
set OD_PT; # All ODs carrying PT
set Links;        
set Paths;       

# Parameters
################################################################################################
param t0{Links};            
param cap{Links};          
param K;                    
param L;                   
param delta{Links, Paths, OD} binary;  
param pw{Paths, OD};
param beta0; 
param beta1; 
param beta2; 
param beta3; 
param beta4; 
param beta5; 
param beta6;
param D{OD};               
param ttpt{OD};            
param B;                    
param N{OD};                  
param h{OD};
param w_w{OD} symbolic; # w in OD_M that w was originated from
param w_c{OD} symbolic; # w in (OD_C intersect OD_V) that w originated
param w_pt{OD} symbolic; # w in (OD_PT intersect OD_V) that w originated
param trans_time{OD};

# Variables
################################################################################################
var X{Links} >= 0;            
var t{Links} >= 0;             
var f{Paths, OD} >= 0;
var u{OD} >= 0;       
var qpt{OD} >= 0;                         
var Upt{OD};                  
var Uc{OD}; 
var Um{OD};                   
var WT{OD} >= 0;          
var phi{OD} >= 0;              
var Ipt{OD} >= 0, <= B;              
var qc{OD} >= 0;   
var ttc{OD} >= 0;            

# Objective
################################################################################################
minimize Total_Travel_Time:
    sum{w in OD_PT} (qpt[w]*ttpt[w]) + sum{a in Links} (X[a] * t[a]);
    #sum{w in OD_PT} (qpt[w]*ttpt[w]) + sum{w in OD_C} (qc[w]*ttc[w]);

# Constraints
###############################################################################################

subject to Mode_Choice1_1{w in OD_PT inter (OD_O diff OD_M)}:
    qpt[w] = D[w] * (1 / (1 + exp(Uc[w] - Upt[w])));

subject to Mode_Choice1_2{w in OD_PT inter OD_V : w_w[w] not in OD_PT}:
    qpt[w] = D[w_w[w]] * (1 / (1 + exp(Uc[w_w[w]] - Um[w_w[w]])));

subject to Mode_Choice1_3{w in OD_PT inter OD_V : w_w[w] in OD_PT}:
    qpt[w] = D[w_w[w]] * (1 / (1 + exp(Uc[w_w[w]] - Um[w_w[w]]) + exp(Upt[w_w[w]] - Um[w_w[w]])));

subject to Mode_Choice1_4{w in OD_PT inter OD_M}:
    qpt[w] = D[w_w[w]] * (1 / (1 + exp(Uc[w_w[w]] - Upt[w_w[w]]) + exp(Um[w_w[w]] - Upt[w_w[w]])));

subject to Mode_Choice1_5{w in OD: not w in OD_PT}:
    qpt[w] = 0;

###############################################################################################

subject to Multimodal_Utility{w in OD_M}:
    Um[w] = (beta0 + beta5) / 2 + beta6 * u[w_c[w]] + beta1 * ttpt[w_pt[w]] + beta2 * WT[w_pt[w]] + beta3 * phi[w_pt[w]] - beta4 * Ipt[w_pt[w]];
    #Um[w] = beta0 + beta1 * ttpt[w_pt[w]] + beta2 * WT[w_pt[w]] + beta3 * phi[w_pt[w]] - beta4 * Ipt[w_pt[w]];
    #Um[w] = beta0 + beta1 * ttpt[w] + beta2 * WT[w] + beta3 * phi[w] - beta4 * Ipt[w_pt[w]];

subject to Public_Transport_Utility{w in (OD_PT inter OD_O) union OD_M}:
    Upt[w] = beta0 + beta1 * ttpt[w] + beta2 * WT[w] + beta3 * phi[w] - beta4 * Ipt[w];

subject to Car_Utility{w in (OD_PT inter OD_O) union OD_M}:
    Uc[w] = beta5 + beta6 * u[w];

subject to Car_Travel_Time{w in OD}:
    ttc[w] = u[w];

###############################################################################################

subject to Budget_Constraint:
    sum{w in OD_PT} (qpt[w] * Ipt[w]) <= B;

###############################################################################################

subject to Waiting_Time{w in OD_PT}:
    WT[w] = trans_time[w] + (qpt[w]/N[w]) * (h[w] / 2);

subject to Phi_Definition{w in OD_PT}:
    phi[w] = exp(qpt[w] / N[w]) - 1;

###############################################################################################

subject to equillibrium_cond1{p in Paths, w in OD_C}:
    ( sum{l in Links} (delta[l, p, w] * t[l]) - u[w]*pw[p, w] ) >= 0;

subject to equillibrium_cond2{p in Paths, w in OD_C}:
    ( sum{l in Links} (delta[l, p, w] * t[l]) - u[w]*pw[p, w] )*f[p, w] = 0;

subject to Link_Flows{l in Links}:
    X[l] = sum{p in Paths, w in OD_C} delta[l, p, w]*f[p, w];

subject to Link_Travel_Time{l in Links}:
    t[l] = t0[l] * (1 + K * (X[l] / cap[l])^L);

################################################################################################

subject to Flow_Conservation1_1{w in OD_O diff OD_M}:
    D[w] = qpt[w] + qc[w];

subject to Flow_Conservation1_2{w in OD_C inter OD_V}:
    qpt[w_pt[w_w[w]]] = qc[w];

subject to Flow_Conservation1_3{w in OD_M}:
    D[w] = qc[w] + qpt[w_pt[w]] + qpt[w];

subject to Flow_Conservation1_4{w in OD_PT inter OD_V}:
    qc[w] = 0;

subject to Flow_Conservation2{w in OD_C}:
    qc[w] = sum{p in Paths} f[p, w]*pw[p, w];

################################################################################################