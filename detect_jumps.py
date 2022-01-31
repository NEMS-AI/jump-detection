import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

def main():
    # data_csv = "/Users/alfredo/Desktop/data/1_inst_exp_data.csv"
    event_csv = "/Users/alfredo/Desktop/data/1_inst_exp_events.csv"
    # pred_csv = "/Users/alfredo/Desktop/data/inst_exp_detected_500msavg.csv"
    # pred_csv = "/Users/alfredo/Desktop/data/GroEL events F=2000 jt=0 t=100 tavg=500.csv"
    # pred_csv = "/Users/alfredo/Desktop/data/GroEL events F=600 jt=0 t=50 tavg=500.csv"
    pred_csv = "/Users/alfredo/Desktop/data/GroEL events Fthresh=2000 tavg=200.csv"
    # pred_csv = "/Users/alfredo/Desktop/data/GroEL events Fthresh=2000 tavg=500.csv"


    # pred_csv = "/Users/alfredo/Desktop/data/inst_exp_detected.csv"
    # pred_100_csv = "/Users/alfredo/Desktop/data/1_inst_exp_detected_F=600_500ms_avg.csv"





    # data = np.genfromtxt(data_csv, delimiter=',')
    # pred_times = find_jumps_loc(data)

    events = np.genfromtxt(event_csv, delimiter=',')
    actual_times = events[:,0]

    # print(pred_times)

    # print(actual_times)

    pred_events = np.genfromtxt(pred_csv, delimiter=',')
    pred_times = pred_events[:,0]

    # tp, fp = evaluate_pred(pred_times, actual_times)
    plot_jumps(events, pred_events)
    # print(fp)

def multivariate_t_test(data, thresh = 1500, plot = True):
    '''
    An initial approach for jump detection using the threshold-t-value
    thresh is currently arbitraritly defined
    '''
    data = data.T
    jump_len = 10
    df_window = 500
    total_window = 2*df_window+jump_len
    F_values = []


    for i in range(0, data.shape[0]-total_window, 1):
        # print(str(i/(data.shape[0]-total_window)*100) + "%" )
        first_range = data[i:i+df_window,:]
        second_range = data[i+df_window+jump_len:i+total_window,:]
        F_val = TwoSampleT2Test(first_range, second_range)
        F_values.append(F_val)

    F_values = np.array(F_values)
    if plot:
        
        x = np.arange(len(F_values))

        plt.title("Scatter plot of jumps")
        plt.scatter(y = F_values, x = x)
        plt.rcParams["figure.figsize"] = [10, 10]
        # plt.savefig("plt500ms.jpg", dpi=600)
        plt.show()

    pred_jumps_ind = np.argwhere(F_values > thresh).flatten()
    return pred_jumps_ind.tolist()

def TwoSampleT2Test(X, Y):
    '''
    From https://www.r-bloggers.com/2020/10/hotellings-t2-in-julia-python-and-r/
    '''
    nx, p = X.shape
    ny, _ = Y.shape
    n_samples = nx+ny
    x_bar = np.mean(X, axis=0) 
    y_bar = np.mean(Y, axis=0)
    delta = x_bar - y_bar
    sigma_x = np.cov(X, rowvar=False)
    sigma_y = np.cov(Y, rowvar=False)
    sigma_pool = sigma_x/2 + sigma_y/2;
    t2 = n_samples/2*(delta.T@np.linalg.inv(sigma_pool))@(delta)
    statistic = t2 * (nx+ny-p-1)/(p*(nx+ny-2))
    F = (2*n_samples-p-1)/(p*(2*n_samples-2))*t2;  
    return F
    

def find_jumps_loc(data):
    time, mode_1, mode_2 = data[:,0], data[:,1], data[:,2]

    mode_1_jumps_loc = multivariate_t_test(np.array([mode_1, mode_2]))
    mode_2_jumps_loc = []

    all_jump_loc = list(set(mode_1_jumps_loc+mode_2_jumps_loc))
    pred_times = time[all_jump_loc]
    
    return(pred_times)

def get_event_times(time, F_values):
    return True


def evaluate_pred(pred_times, actual_times):
    true_pos_count = 0
    false_pos_count = 0
    diffs = []
    for pred in pred_times:
        pred_range = [pred-0.01,pred+0.01]
        is_real = False

        for time in actual_times:
            if time > pred_range[0] and time < pred_range[1]:
                is_real = True
                diffs.append(time-pred)

        if is_real:
            true_pos_count += 1
            
        else:
            false_pos_count += 1

    print(np.mean(diffs))
    print(np.std(diffs))

    return(true_pos_count, false_pos_count)

def plot_jumps(data, pred_data):
    plt.title("Scatter plot of jumps")
    plt.ylabel('df2')
    plt.xlabel('df1')
    plt.ylim(-3e-6, 1e-6)
    plt.xlim(-3e-6, 1e-6)

    true = plt.scatter(data[:,1], data[:,2])
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

