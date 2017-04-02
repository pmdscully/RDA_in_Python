# ==============================================================
# Written by Peter M. D. Scully (pmdscully@gmail.com) April 2015
# ==============================================================
# -- Based on: 
# Hilder, J. A., Owens, N. D., Hickey, P. J., Cairns, S. N., Kilgour, D. P., Timmis, J., & Tyrrell, A. 
# (2011). Parameter optimisation in the receptor density algorithm. 
# In Artificial Immune Systems (pp. 226-239). 
# Springer Berlin Heidelberg.
# 
# ==============================================================
from __future__ import division
import scipy.stats as stats
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as patches
import csv,sys

def write(filename, rows):
        with open(filename, 'wba') as f:
            writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            for row in rows:
                writer.writerow(row)

def region(si,pe,ne):
        #(1): 0 to B.
        if(si < (1-b)*B):
                pe = si/(1-b)
                ne = 0
                #print "Region 1"
        #(2): B
        elif (si > (1-b)*B) and (si < (1-b)*B + (g/(1-d))):  
                pe = B                  # APPROX EQUALS .. what does that mean?
                ne = si - (1-b)*B       # APPROX EQUALS .. what does that mean?
                #print "Region 2"
                
        #(3): B to l.
        elif(si > (1-b)*B + (g/(1-d))) and (si < (1-b)*l + (g/(1-d))):
                #pe = B                 # "Will not be able to hold p = B"... what does that mean?
                ne = g/(1-d)
                #print "Region 3"
                
        #(4): p > l..... surely, p > l is an anomaly, what does is "not of interest" mean?
        elif(si > (1-b)*l + (g/(1-d))):
                #pe = B                 # Calculated same as region 3.
                ne = g/(1-d)
                #print "Region 4"
        else:
                print "Ohhh..."

        #if pe > B: 
        #then negative_feedback is linearly generated, reduce p.
        return pe, ne

def heaviside_function(x):
        return 0.5 * (np.sign(x) +1.0)

def recurrence_equations(p,n,pe,ne,prev_p,prev_n,si):
        p       = (b*prev_p) + si - prev_n
        n       = (d*prev_n) + (g * heaviside_function(prev_p - B))
        return p,n

def classify(p):
        if p >= l:      # then classification = ct=1
                ct = 1
        else:           # else ct=0
                ct = 0
        return ct


def anomaly_state(i,ct, ct_ta,ct_ts,ct_te,p_color,ct_distB,p):
        
        if(ct == 1 and ct_ta == 0):    # Starting anomaly period
                #print "Anomaly start at time step "+str(i)+" = "+str(ct)
                ct_ta = 1
                ct_ts = i
                p_color = 'r'
                ct_distB = []           # start new signature set
                ct_distB += [p - B]       # - append to signature
        elif ct == 1 and ct_ta == 1:    # In anomaly period
                #print "+ "+str(i)
                p_color = 'r'
                ct_distB += [p - B]       # - append to signature
        elif ct == 0 and ct_ta == 1:    # Finished anomaly period
                ct_ta = 0
                ct_te = i
                # ------- NOTE --------
                # I understand that Hilder defines the signature as based upon the duration of the anomaly experience.
                avg_distB = get_signature(ct_distB,ct_ts,ct_te)
                store_signature(ct_te,ct_ts,avg_distB,ct_distB)
                
        return ct_ta,ct_ts,ct_te,p_color,ct_distB

def get_signature(ct_distB,ct_ts,ct_te):
        result = 0
        if len(ct_distB) > 0:
                time_ratio = 1/(ct_te - ct_ts) 
                accum = 0
                for dist in ct_distB:
                        accum += dist * heaviside_function(dist)
                result = time_ratio * accum
        return result

def store_signature(ct_te,ct_ts,avg_distB,ct_distB):
        global bootstrapping
        global signatures
        global CT_Total_Count
        global TP_count
        global FN_count
        print "Known Signatures: "+str(len(signatures))
        print "Store: "+str(ct_te-ct_ts)+",  \t"+str('%.4f' % avg_distB)+",  \t"+str('%.4f' % np.std(ct_distB))
        k = avg_distB
        v = (ct_te,ct_ts,avg_distB,ct_distB)
        if bootstrapping:
                signatures[k] = v
        else:   
                matched = False
                for B in signatures.values():
                        correlation             = correlation_of_single_sensor_signature(v, B)
                        match                   = is_matched( correlation , B )
                        if match:
                                matched = True
                if not matched:
                        print "Failed to find matching signature with error threshold: "+str(ct_error)
                        FN_count += 1
                else:
                        TP_count += 1
                CT_Total_Count += 1
        print "--------------------------------------------------------------------"

def is_matched( correlation , possible_match_sig ):
        global trained_min_avg
        global trained_max_avg
        global ct_error
        matched                 = False
        normalised_correl       = normalise(correlation, trained_min_avg, trained_max_avg)
        if abs(normalised_correl) < ct_error:
                (Bct_te,Bct_ts,Bavg_distB,Bct_distB) = possible_match_sig
                print "-----: "+str(Bct_te-Bct_ts)+",  \t"+str('%.4f' % Bavg_distB)+",  \t"+str('%.4f' % np.std(Bct_distB))+"  \tCorrelation:\t"+str('%.4f' % correlation)+"  \tCorrelation:\t"+str('%.4f' % correlation)+"  \t("+str('%.4f' % normalised_correl)+")"
                matched = True
        return matched

def normalise(x, min_x, max_x):
        if min_x == -1 or max_x == -1:
                print "Normalise failed as min / max value(s) not available.. Exiting."
                sys.exit(1)
        return (x-min_x)/(max_x-min_x)

def get_trained_max_mins():
        global signatures
        global trained_min_avg
        global trained_max_avg
        avg_values = []
        for B in signatures.values():
                (Bct_te,Bct_ts,Bavg_distB,Bct_distB) = B
                avg_values += [Bavg_distB]
        trained_min_avg = min(avg_values)
        trained_max_avg = max(avg_values)
        return trained_min_avg, trained_max_avg
                

def correlation_of_single_sensor_signature(A, B):
        # -------- NOTE --------------
        # The representation of a signature (A,B) will need to change when taking the inputs from multiple sensors.
        # Currently it only handles 1 input sensor.
        (Act_te,Act_ts,Aavg_distB,Act_distB) = A
        (Bct_te,Bct_ts,Bavg_distB,Bct_distB) = B
        Anum_receptors = 1
        Bnum_receptors = 1
        return (1/(Anum_receptors*Bnum_receptors)) * (Aavg_distB - Bavg_distB)

def correlation_of_multiple_sensor_signature(A, B):
        print "not implemented."

def reset():
        global signatures
        global bootstrapping
        global trained_min_avg
        global trained_max_avg
        signatures = {}
        bootstrapping = True
        trained_min_avg = -1
        trained_max_avg = -1

def main(M_train, M_test, init_p, init_n ):
        global bootstrapping
        global signatures
        global CT_Total_Count

        data = np.concatenate((M_train,M_test)).tolist()
        num_folds = 10
        subset_size = len(data)/num_folds
        for fold in range(num_folds):
                testing_this_round  = data[int(fold*subset_size):][:int(subset_size)]
                training_this_round = data[:int(fold*subset_size)] + data[int((fold+1)*subset_size):]
                M = np.concatenate((training_this_round,testing_this_round))
                # train using training_this_round
                # evaluate against testing_this_round
                # save accuracy
                
                reset()
                
                p  = init_p
                n  = init_n
                pe = p
                ne = n
                prev_p  = p
                prev_n  = n
                ct_ta   = 0     # Anomaly time is now? (0/1)
                ct_ts   = 0     # Anomaly time start step
                ct_te   = 0     # Anomaly time end step
                ct_distB = []   # Anomaly signature, distance of p from B

                fig, ax = plt.subplots(1)
                if WATCH_IN_REALTIME:
                        plt.show()
                        plt.ion()
                colors = cm.rainbow(np.linspace(0, 1, 5))
                
                marker_size = 3
                anom_Pro = plt.scatter(0, 0, color='r')
                norm_Pro = plt.scatter(0, 0, color=colors[2])
                
                
                print "============ Training: =============="
                for i in range(len(M)):
                        if i == len(training_this_round):
                                print "============ Testing: =============="
                                bootstrapping = False
                                get_trained_max_mins()
                        
                        # --- Input: ---
                        si      = M[i]
                        
                        # --- Evaluate: ---
                        pe,ne   = region(si,pe,ne)
                        p,n     = recurrence_equations(p,n,pe,ne,prev_p,prev_n,si)
                        prev_p  = p
                        prev_n  = n
                        ct      = classify(p)
                        p_color = colors[2]
                        ct_ta,ct_ts,ct_te,p_color,ct_distB = anomaly_state(i,ct, ct_ta,ct_ts,ct_te,p_color,ct_distB,p)
                        
                        
                        
                        # --- Reporting in Plot: ---
                        if ct == 1:
                                anom_Pro = plt.scatter(i, p, color=p_color, s=marker_size)
                        else:
                                norm_Pro = plt.scatter(i, p, color=p_color, s=marker_size)
                        Thres   = plt.scatter(i, l, color=colors[0], s=1) 
                        Neg     = plt.scatter(i, n, color=colors[1], s=marker_size)
                        Sensor  = plt.scatter(i, si, color=colors[3], s=marker_size)
                        leg     = plt.legend((anom_Pro, norm_Pro, Neg, Thres, Sensor),("Progress Anomaly", "Progress Normal", "Negative Feedback","Threshold","Sensor Value"),
                                           scatterpoints=1,
                                           loc='upper right',
                                           ncol=1,
                                           fontsize=8)
                        leg.get_frame().set_alpha(0.5)
                        
                        if WATCH_IN_REALTIME:
                                plt.draw()
                                time.sleep(0.00001)
                
                # place a text box in upper left in axes coords
                props = dict(boxstyle='round', facecolor='wheat', alpha=0.25)    # these are matplotlib.patch.Patch properties
                
                # +('$\mathrm{M3:} \mathrm{n}=%.2f$, $\mu=%.2f$, $\sigma=%.2f$, $\mathrm{Norm}$\n'%(T,mu1,sig1)) \
                # +('$\mathrm{M4:} \mathrm{n}=%.2f, \mathrm{Uniform}$\n'%(T)) \
                # +('$\mathrm{M5:} \mathrm{n}=%.2f$, $\mu=%.2f$, $\sigma=%.2f$, $\mathrm{Norm}$\n'%(T,mu2,sig2)) \
                # +('$\mathrm{M6:} \mathrm{n}=%.2f$, $\mu=%.2f$, $\sigma=%.2f$, $\mathrm{Norm}$\n'%(T,mu1,sig1)) \

                textstr = "Bootstrap Model: (0-"+str(len(training_this_round))+")\n" \
                         +('$\mathrm{M1:} \mathrm{n}=%.2f$, $\mu=%.2f$, $\sigma=%.2f$, $\mathrm{Norm}$\n'%(T,mu1,sig1)) \
                         +('$\mathrm{M2:} \mathrm{n}=%.2f$, $\mu=%.2f$, $\sigma=%.2f$, $\mathrm{Norm}$\n'%(T,mu2,sig2)) \
                         +("Test Model: ("+str(len(training_this_round))+"-"+str(len(training_this_round)+len(testing_this_round))+")\n") \
                         +('$\mathrm{M1:} \mathrm{n}=%.2f$, $\mu=%.2f$, $\sigma=%.2f$, $\mathrm{Norm}$\n'%(T,mu1,sig1)) \
                         +('$\mathrm{M2:} \mathrm{n}=%.2f$, $\mu=%.2f$, $\sigma=%.2f$, $\mathrm{Norm}$\n'%(T,mu2,sig2)) \
                         +("RDA Parameters:\n") \
                         +('$\mathrm{l}=%.2f$, $\mathrm{b}=%.2f$, $\mathrm{d}=%.2f$, $\mathrm{g}=%.2f$, $\mathrm{B}=%.2f$\n'%(l,b,d,g,B)) \
                         +("1-Class Classification:\n") \
                         +('$\mathrm{\epsilon}=%.4f,  $'%(ct_error))
                         
                textstr += "$TP=%.2f$, "%(TP_count/CT_Total_Count)+"$TN=%.2f$, "%(TN_count/CT_Total_Count)+"$FP=%.2f$, "%(FP_count/CT_Total_Count)+"$FN=%.2f$\n"%(FN_count/CT_Total_Count)
                textstr += "$Fold Samples: %.0f$, $Accumulated Tests: %.0f$"%(len(signatures),CT_Total_Count)
                plt.text(0.02, 0.98, textstr, fontsize=8, transform=ax.transAxes, verticalalignment='top', bbox=props)
                plt.xlabel('Time series inputs from model')
                plt.ylabel('Values')
                plt.title("Receptor Density Algorithm (RDA)", fontsize=10)
                
                plt.savefig("result_plot_fold="+str(fold)+"-of-"+str(num_folds)+".png", format='png')
        # find mean accuracy over all rounds


# =============================================================
# ---  Global Data:
# =============================================================

WATCH_IN_REALTIME = False
PROFILE = "HILDER_GA_OPTIMISED"

signatures = {}
bootstrapping = True
trained_min_avg = -1
trained_max_avg = -1
FP_count = 0
FN_count = 0
TP_count = 0
TN_count = 0
CT_Total_Count = 0

if __name__ == "__main__":
        if PROFILE == "HILDER":
        #        w = 500
        #        h = 1.9
        #        e = 0.0015
        #        #a = 4.67x10^-7
        #        D = 1
                l = 400         # 0 < B < l
                b = 0.84        # 0 < b,d < 1
                d = 0.98        # 0 < b,d < 1
                g = 0.35        # g > 0
                B = 1           # 0 < B < l
                initial_p = 0
                initial_n = 0
        
        elif PROFILE == "HILDER_GA_OPTIMISED":
        #        w = 500
        #        h = 1.9
        #        e = 0.0015
        #        #a = 4.67x10^-7
        #        D = 1
                l = 5           # 0 < B < l
                gb = 0          # mentioned, but not sure if it is g*b or gb.
                b = 0.84        # 0 < b,d < 1
                d = 0.98        # 0 < b,d < 1
                g = 0.35        # g > 0
                B = 1           # 0 < B < l
                initial_p = 0
                initial_n = 0
        elif PROFILE == "OWENS":
        #        w = 500
        #        h = 1.9
        #        e = 0.0015
        #        #a = 4.67x10^-7
        #        D = 1
                l = 0.1
                b = 0.00195     #0.84
                d = 0.000391    #0.98
                g = 4           #0.35
                B = 0.01        #1
                initial_p = 0
                initial_n = 0

        # -- Match Threshold:
        ct_error = 0.008

        # -- Model
        T = 200
        mu1 = 6
        mu2 = 120
        sig1 = 8.9
        sig2 = 20
        #seed = 282629734
        #np.random.seed(seed)
        M1 = stats.norm.rvs(mu1,sig1,T)
        M2 = stats.norm.rvs(mu2,sig2,T)
        M3 = stats.norm.rvs(mu1,sig1,T)
        M4 = stats.uniform.rvs(size=T)
        M5 = stats.norm.rvs(mu2,sig2,T)
        M6 = stats.norm.rvs(mu1,sig1,T)
        
        M_train = np.concatenate((M1,M2))
        M_train = np.concatenate((M_train,M1))
        M_train = np.concatenate((M_train,M2))
        M_test = np.concatenate((M1,M2))
        M_test = np.concatenate((M_test,M1))
        M_test = np.concatenate((M_test,M2))


        main(M_train, M_test, initial_p, initial_n)

