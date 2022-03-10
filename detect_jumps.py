import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from sklearn.metrics import f1_score

def main():
    data_csv = "/Users/alfredo/Desktop/data/1_inst_exp_data.csv"
    data_csv = "test_data.csv"

    event_csv = "/Users/alfredo/Desktop/data/1_inst_exp_events.csv"
    # pred_csv = "/Users/alfredo/Desktop/data/inst_exp_detected_500msavg.csv"
    # pred_csv = "/Users/alfredo/Desktop/data/GroEL events F=2000 jt=0 t=100 tavg=500.csv"
    # pred_csv = "/Users/alfredo/Desktop/data/GroEL events F=600 jt=0 t=50 tavg=500.csv"
    pred_csv = "/Users/alfredo/Desktop/data/GroEL events Fthresh=2000 tavg=200.csv"
    # pred_csv = "/Users/alfredo/Desktop/data/GroEL events Fthresh=2000 tavg=500.csv"

    
    # pred_csv = "/Users/alfredo/Desktop/data/inst_exp_detected.csv"
    # pred_100_csv = "/Users/alfredo/Desktop/data/1_inst_exp_detected_F=600_500ms_avg.csv"




    data = np.genfromtxt(data_csv, delimiter=',')
    # print(data[-1])
    pred_times = find_jumps_loc(data)

    events = np.genfromtxt(event_csv, delimiter=',')[0:1000]
    actual_times = events[:,0]

    # print(pred_times)

    # print(actual_times)

    pred_events = np.genfromtxt(pred_csv, delimiter=',')
    pred_times = pred_events[:,0]
    

    # tp, fp, fn, mse = evaluate_pred(pred_events, events, buffer = 0.001)
    # recall = (tp)/(tp+fn)
    # precision = (tp)/(tp+fp)
    # f1_score = 2*(precision * recall)/(precision + recall)
    # print(f1_score)
    # print(mse)
    plot_jumps(events, pred_events)
    # print(fp)


    

def find_jumps_loc(data, datatype = "synthetic"):
    '''
    Datatype: choose from 'noise', 'data', 'synthetic'
    '''
    n_samples, n_feat = data.shape
    n_modes = n_feat - 1

    if n_modes == 2:
        tsample = .00025        
        tmeas_detect = .10      
        tjump = .06             
        tjump_pre = .03         
        tmeas = .2              
    elif n_modes == 3:
        tsample = 0.02
        tmeas_detect = 0.5     
        tjump = 0.1            
        tjump_pre = 0.05
        tmeas = 2 
    
    tjump_post = tjump - tjump_pre

    Ndetect = int(np.floor(tmeas_detect/tsample)); 
    Njump = np.floor(tjump/tsample);
    Npre = np.floor(tjump_pre/tsample);
        

    # time, mode_1, mode_2 = data[:,0], data[:,1], data[:,2]
    tvect = data[:,0]
    fvect = data[:,1:n_modes+1].T
    print(fvect[:,0])

    if datatype == "data":
        if n_modes == 2:
            # note this logic for getting tvectuse may be wrong in porting data
            tvectuse = tvect > 10 and tvect < 92
            tvect = tvect[tvectuse]
            tvect = tvect-tvect[1]
            fvect = fvect[:,tvectuse]
            tstart = 10
            tlag_mode2_per_second = .013/830
            tlag_mode3_per_second = 0
            tlag_buffer = 100
        elif n_modes == 3:
            tstart = 0
            tlag_mode2_per_second = 0
            tlag_mode3_per_second = 0
            tlag_buffer = 0
    elif datatype == "noise" or datatype == "synthetic":
        tstart = 0
        tlag_mode2_per_second = 0
        tlag_mode3_per_second = 0
        tlag_buffer = 0

    ttot = len(tvect);
    tfin = int(ttot-Ndetect*3-Njump);
    print(tfin)
    Fstats = np.zeros(shape = tfin);
    Fstat_thresh_detect = 300;
    Fstat_thresh_meas = 600;

    # TODO: Confirm length of Fstats, and length of for loop
    for ti in range(1+tlag_buffer,tfin-tlag_buffer-Ndetect):
        print(ti/ttot*100); 

        xi_range = np.arange(ti,ti+Ndetect-1, dtype="int")
        yi_range = np.arange(ti+Ndetect+Njump,ti+2*Ndetect+Njump-1, dtype="int");
        
        tlag_m2 = (tvect[ti]+tstart)*tlag_mode2_per_second;
        Nlag_m2 = int(np.floor(tlag_m2/tsample))

        # print(fvect[1,xi_range])
        # print(fvect[2,xi_range+Nlag_m2].shape)


        
        fvect_xi = np.stack((fvect[0,xi_range], fvect[1,xi_range+Nlag_m2]));
        fvect_yi = np.stack((fvect[0,yi_range], fvect[1,yi_range+Nlag_m2]));

        
        if n_modes == 3:
            tlag_m3 = (tvect[ti]+tstart)*tlag_mode3_per_second;
            Nlag_m3 = int(np.floor(tlag_m3/tsample));
            fvect_xi = np.concatenate((fvect_xi, [fvect[2,xi_range+Nlag_m3]]));
            fvect_yi = np.concatenate((fvect_yi, [fvect[2,yi_range+Nlag_m3]]));
        
        xbar = np.mean(fvect_xi,axis = 1);
        ybar = np.mean(fvect_yi,axis = 1);
        
        sigma_x = np.cov(fvect_xi);
        sigma_y = np.cov(fvect_yi);
        sigma_pool = sigma_x/2 + sigma_y/2;
        t2 = Ndetect/2*np.matmul(np.matmul((xbar-ybar),np.linalg.inv(sigma_pool)),np.array((xbar-ybar)).T);
        p_dim = n_modes;
        Fstat = (2*Ndetect-p_dim-1)/(p_dim*(2*Ndetect-2))*t2;
        Fstats[int(ti+Ndetect)] = Fstat;
        
    val = Fstats.shape[0]
    # print(val)
    # print(tfin)
    print("Plotting Fstats")
    plt.scatter(tvect[0:val], Fstats)
    plt.show()
    
    jumps_detected = np.array([[0,0,0,0]]);
    ti=1+tlag_buffer;
    while ti < tfin-tlag_buffer:
        
        if Fstats[ti] < Fstat_thresh_detect:
            ti = ti+1;
        else:
            tii = 0;
            while tii < tfin-ti and Fstats[ti+tii] > Fstat_thresh_detect:
                tii = tii+1;
            
            Fstatmax = max(Fstats[ti:ti+tii]);
            
            peak_left_i = ti;
            while Fstats[peak_left_i] < Fstatmax/2 and peak_left_i < tfin:
                peak_left_i = peak_left_i + 1;
            
            peak_right_i = ti+tii;
            while Fstats[peak_right_i] < Fstatmax/2 and peak_right_i > 1:
                peak_right_i = peak_right_i - 1;
            
            
            peak_width_fwhm_range = peak_right_i - peak_left_i;
            ti_jump = ti+tii-1;
            # jumps_detected = [jumps_detected; ti_jump tvect[ti_jump] tvect[peak_width_fwhm_range] Fstatmax];
            jumps_detected = np.concatenate((jumps_detected, [[ti_jump, tvect[ti_jump], tvect[peak_width_fwhm_range], Fstatmax]]))
            ti = ti + tii + 1;
        
    print(jumps_detected.shape)
    jumps_measured = measure_jumps(jumps_detected, tmeas, tsample, tvect, fvect, Fstat_thresh_meas, tjump_post, Npre, Njump, tstart, tlag_mode2_per_second, tjump, tjump_pre, n_modes, tlag_mode3_per_second)
    print(jumps_measured.shape)
    print(jumps_measured)

    
    return(True)

def measure_jumps(jumps_detected, tmeas, tsample, tvect, fvect, Fstat_thresh_meas, tjump_post, Npre, Njump, tstart, tlag_mode2_per_second, tjump, tjump_pre, nmodes, tlag_mode3_per_second):
    Nmeas = np.floor(tmeas/tsample);
    jumps_measured = [[0,0,0,0,0,0]];
    rel_jump_ts_1 = [];
    rel_jump_ts_2 = [];
    rel_jump_ts_3 = [];
    Fstats_ts = [];

    # check indices on this range
    for ji in range(1, jumps_detected.shape[0]-2):
            
        # print(jumps_detected[ji-1,0])
        t_prev = tvect[int(jumps_detected[ji-1,0])];
        t_curr = tvect[int(jumps_detected[ji,0])];
        t_next = tvect[int(jumps_detected[ji+1,0])];
        
        Fstatmax = jumps_detected[ji,3];
        if Fstatmax < Fstat_thresh_meas:
            continue

        
        if (t_curr - tjump_pre < t_prev + tmeas + tjump_post) or (t_next - tjump_pre < t_curr + tmeas + tjump_post): 
            continue
        
        # % can optionally exclude jumps with too wide peak width
        peak_width_fwhm = jumps_detected[ji,2];
        if peak_width_fwhm > tjump*1.5:
            continue
        
        ti_jump = int(jumps_detected[ji,0]);
        ti = int(ti_jump-Npre-Nmeas);

        xi_range = np.arange(ti,ti+Nmeas-1, dtype="int")
        yi_range = np.arange(ti+Nmeas+Njump,ti+2*Nmeas+Njump-1, dtype="int");

        all_range = np.arange(ti,ti+2*Nmeas+Njump-1, dtype="int");
        
        tlag_m2 = (tvect[ti]+tstart)*tlag_mode2_per_second;
        Nlag_m2 = int(np.floor(tlag_m2/tsample));
        fvect_xi = np.stack((fvect[0,xi_range], fvect[1,xi_range+Nlag_m2]));
        fvect_yi = np.stack((fvect[0,yi_range], fvect[1,yi_range+Nlag_m2]));
        
        if nmodes == 3:
            tlag_m3 = (tvect[ti]+tstart)*tlag_mode3_per_second;
            Nlag_m3 = int(np.floor(tlag_m3/tsample));
            fvect_xi = np.concatenate((fvect_xi, [fvect[2,xi_range+Nlag_m3]]));
            fvect_yi = np.concatenate((fvect_yi, [fvect[2,yi_range+Nlag_m3]]));
        
        xbar = np.mean(fvect_xi,axis = 1);
        ybar = np.mean(fvect_yi,axis = 1);
        
        rel_jump = np.divide(ybar,xbar)-1;
        if nmodes ==2:
            jumps_measured = np.concatenate((jumps_measured, [[tvect[ti_jump], peak_width_fwhm, Fstatmax, rel_jump[0], rel_jump[1]]]))
        elif nmodes ==3:
            jumps_measured = np.concatenate((jumps_measured, [[tvect[ti_jump], peak_width_fwhm, Fstatmax, rel_jump[0], rel_jump[1], rel_jump[2]]]))

    return jumps_measured
        # jumps_measured = [jumps_measured; tvect(ti_jump) peak_width_fwhm Fstatmax rel_jump'];
        

def get_event_times(time, F_values):
    return True


def evaluate_pred(pred_events, events, buffer):
    pred_times = pred_events[:,0]
    actual_times = events[:,0]
    
    true_pos_count = 0
    false_pos_count = 0
    diffs = []
    for i, pred in enumerate(pred_times):
        pred_range = [pred-buffer,pred+buffer]
        is_real = False

        for j, time in enumerate(actual_times):
            if time > pred_range[0] and time < pred_range[1]:
                is_real = True
                diffs.append(time-pred)

        if is_real:
            true_pos_count += 1
            
        else:
            false_pos_count += 1

    # print(np.mean(diffs))
    # print(np.std(diffs))
    n_pos = len(actual_times)
    false_negaive_count = n_pos - true_pos_count
    mse = np.sum((np.array(diffs)**2))/len(diffs)
    return(true_pos_count, false_pos_count, false_negaive_count, mse)



def plot_jumps(data, pred_data):

    mass_prob = [.68,.28,.04]
    # mass_prob = [1,0,0]

    mass_scale = np.random.choice([1,2,3],size = data.shape[0], p = mass_prob)

    plt.title("Scatter plot of jumps")
    plt.ylabel('df2')
    plt.xlabel('df1')
    # plt.ylim(-3e-6, 1e-6)
    # plt.xlim(-3e-6, 1e-6)

    true = plt.scatter(data[:,1]*mass_scale, data[:,2]*mass_scale)
    pred = plt.scatter(pred_data[:,0], pred_data[:,1])
    plt.legend([true, pred],['COMSOL', "200ms Avg"])
    plt.rcParams["figure.figsize"] = [10, 10]
    plt.savefig("tavg200.jpg", dpi=600)
    plt.show()


def univariate_t_test(mode, thresh = 5000):
    '''
    An initial approach for jump detection using the threshold-t-value
    thresh is currently arbitraritly defined
    '''
    jump_len = 10
    df_window = 500
    total_window = 2*df_window+5
    t_values = []


    for i in range(0, len(mode)-total_window, 1):
        first_range = mode[i:i+df_window]
        second_range = mode[i+df_window+jump_len:i+total_window]
        t_val = stats.ttest_ind(first_range, second_range).statistic
        t_values.append(t_val)

    t_values = np.array(t_values)
    pred_jumps_ind = np.nonzero(t_values > thresh)[0]
    return pred_jumps_ind.tolist()

if __name__ == "__main__":
    main()
