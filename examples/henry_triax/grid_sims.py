import os
import numpy as np
import time
import multiprocessing as mp

# This file currently runs 50 simulations in the triaxial case.
# To run 50 simulations in the oblate case
#   1. rename dir_path to note difference
#   2. Change "triax" variable in dydt_rb() in problem.c to 0.0
#   3. To plot, make sure to change triax to True

dir_path = "./2body_data_triax"
# dir_path = "./2body_data_obl"

def run_sim(trial_num, n_trials=50):

    # max_omega = 4. # 2 because otherwise obliquity is excited # (1+(np.pi/2/np.arctan(1/Q_tide)))
    # omega_to_n = max_omega*np.random.default_rng().uniform()
    omegas = [1.6,3.1]
    omega_to_n = omegas[trial_num % len(omegas)]
    thetas = np.linspace(0.05,179.95,int(n_trials / len(omegas)))
    theta = thetas[trial_num//2]
    tf = 1e7 # in orbital periods

    # make output directory and file
    os.makedirs(dir_path, exist_ok=True)

    ### RUN SIMULATION ###
    outfile = str(dir_path)+"/trial_"+str(int(trial_num))
    if os.path.exists(outfile):
        print(outfile, 'already exists, skipping')
        return
    command = "./rebound "+str(omega_to_n)+" "+str(theta)+" "+str(tf)+" "+outfile
    print(command)
    os.system(command)

# main function
if __name__ == '__main__':
    # to change params, see keyword args in run_sim()
    n_trials = 50 # if I change this, make sure to also change in keyword args in run_sim()
    start = time.time()
    with mp.Pool(9) as pool:
        pool.map(run_sim, range(n_trials))

    tot_time = time.time() - start
    hrs = tot_time // 3600
    mins = (tot_time % 3600) // 60
    secs = int((tot_time % 3600) % 60)
    print(f"Total Runtime: {hrs} hours {mins} minutes {secs} seconds", flush=True)

