import numpy as np
from scipy import stats

def main():
    data_csv = "/Users/alfredo/Desktop/data/inst_exp_data_short.csv"
    event_csv = "/Users/alfredo/Desktop/data/inst_exp_events_short.csv"

    data = np.genfromtxt(data_csv, delimiter=',')
    pred_times = find_jumps_loc(data)

    events = np.genfromtxt(event_csv, delimiter=',')
    actual_times = events[:,0]

    print(pred_times)

    print(actual_times)


    tp, fp = evaluate_pred(pred_times, actual_times)
    print(fp)

def multivariate_t_test(data, thresh = 25):
    '''
    An initial approach for jump detection using the threshold-t-value
    thresh is currently arbitraritly defined
    '''
    jump_len = 10
    df_window = 500
    total_window = 2*df_window+5
    t_values = []


    for i in range(0, data.shape[0]-total_window, 1):
        first_range = data[i:i+df_window,:]
        second_range = data[i+df_window+jump_len:i+total_window,:]
        t_val, _ = TwoSampleT2Test(first_range, second_range)
        t_values.append(t_val)

    t_values = np.array(t_values)
    pred_jumps_ind = np.nonzero(t_values > thresh)[0]
    return pred_jumps_ind.tolist()

def TwoSampleT2Test(X, Y):
    '''
    From https://www.r-bloggers.com/2020/10/hotellings-t2-in-julia-python-and-r/
    '''
    nx, p = X.shape
    ny, _ = Y.shape
    delta = np.mean(X, axis=0) - np.mean(Y, axis=0)
    Sx = np.cov(X, rowvar=False)
    Sy = np.cov(Y, rowvar=False)
    S_pooled = ((nx-1)*Sx + (ny-1)*Sy)/(nx+ny-2)
    t_squared = (nx*ny)/(nx+ny) * np.matmul(np.matmul(delta.transpose(), np.linalg.inv(S_pooled)), delta)
    statistic = t_squared * (nx+ny-p-1)/(p*(nx+ny-2))
    F = f(p, nx+ny-p-1)
    p_value = 1 - F.cdf(statistic)
    print(f"Test statistic: {statistic}\nDegrees of freedom: {p} and {nx+ny-p-1}\np-value: {p_value}")
    return statistic, p_value
    

def find_jumps_loc(data):
    time, mode_1, mode_2 = data[:,0], data[:,1], data[:,1]
    mode_1_jumps_loc = univariate_t_test(mode_1)
    mode_2_jumps_loc = univariate_t_test(mode_2)
    all_jump_loc = list(set(mode_1_jumps_loc+mode_2_jumps_loc))
    pred_times = time[all_jump_loc]
    return(pred_times)

def univariate_t_test(mode, thresh = 40):
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

def evaluate_pred(pred_times, actual_times):
    true_pos_count = 0
    false_pos_count = 0
    for pred in pred_times:
        pred_range = [pred-0.01,pred+0.01]
        is_real = False

        for time in actual_times:
            if time > pred_range[0] and time < pred_range[1]:
                is_real = True

        if is_real:
            true_pos_count += 1
        else:
            false_pos_count += 1


    return(true_pos_count, false_pos_count)


if __name__ == "__main__":
    main()