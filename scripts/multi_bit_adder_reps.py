from loica import *
from genephlo import *

net, A_in, B_in, C_in, ref = construct_multi_bit_adder(N=4, deg_rate=1)

# Optimize switching concentrations for operators
optimize(net, d=100, alpha0=500)

# ---
# Test network
# Get inputs from job id
import sys
AB = int(sys.argv[1]) - 1
A = (AB & 0b11110000) >> 4
B = AB & 0b00001111
C = 0

ref_period = 120
interval = 0.24
errors = 0
f = open(f'result_multi_bit_adder_{A}_{B}.csv', 'wt')
f.write('rep,A,B,sum,result,error\n')

sample, prof = setup_adder_sample(net, 4, ref, A_in, B_in, C_in, A, B, C, ref_period, scale=100)

for rep in range(1):
    meas = run_assay(sample, ref_period, interval, stochastic=True)
    phases, lags, corrs = compute_phases(net, meas, prof, interval, ref_period)
    result = (phases[2]<90)*1 + (phases[4]<90)*2 + (phases[6]<90)*4 + (phases[8]<90)*8 + (phases[7]<90)*16
    error = A + B - result
    print(f'A={A}, B={B}, A+B={A+B}, result={result}, error={error!=0}')
    f.write(f'{rep},{A},{B},{A+B},{result},{error!=0}\n')

f.close()

# ---
