from loica import *
from genephlo import *

net, A_in, B_in, C_in, ref = construct_multi_bit_adder(N=1, deg_rate=1)

# Optimize switching concentrations for operators
optimize(net, d=100, alpha0=500)

# ---
# Test network
# Get inputs from job id
import sys
if len(sys.argv)>1:
    seed = int(sys.argv[1])
    np.random.seed(seed)

ref_period = 120
interval = 0.24
scale = 100
errors = 0
f = open(f'result_full_adder_{seed}.csv', 'wt')
f.write('seed,IN1,IN2,IN3,SUM,CARRY,result,error\n')



for IN1 in range(2):
    for IN2 in range(2):
        for IN3 in range(2):
            input_phases = np.array([(1 - IN1), (1 - IN2), (1 - IN3)])
            sample, prof = setup_adder_sample(net, 1, ref, A_in, B_in, C_in, IN1, IN2, IN3, ref_period, scale=100)
            meas = run_assay(sample, ref_period, interval, stochastic=False)
            phases, lags, corrs = compute_phases(net, meas, prof, interval, ref_period)
            result = (phases[2]<90)*1 + (phases[1]<90)*2
            sum_phase =  phases[2]
            carry_phase = phases[1]
            error = (IN1 + IN2 + IN3) - result
            print(f'IN1={IN1}, IN2={IN2}, IN3={IN3}, result={result}, error={error!=0}')
            f.write(f'{seed},{input_phases[0]*180}, {input_phases[1]*180}, {input_phases[2]*180}, {sum_phase}, {carry_phase}, {result},{error!=0}\n')

f.close()

# ---

