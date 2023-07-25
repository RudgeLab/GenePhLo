from loica import *
from genephlo import *

def setup_sample(ref_period, phase, scale):
    def growth_rate(t):
        return 0 # gompertz_growth_rate(t, 0.01, 1, 1, 0.5)

    def biomass(t):
        return 1 # gompertz(t, 0.01, 1, 1, 0.5)

    metab = SimulatedMetabolism('metab', biomass, growth_rate)

    sample = Sample(genetic_network=net, 
                        metabolism=metab)

    def prof_input(period, phase=0):
        def f(t):
            p = period
            return ((t-phase*0.5*p/np.pi)%p)/p < 1/2
        return f

    # Input
    sample.set_supplement(ahl, concentration=scale, profile=prof_input(ref_period, phase))

    return sample, prof_input(ref_period, 0)

net = GeneticNetwork()

deg_rate = 1

repressors = [Regulator(name=f'Rep{i}', degradation_rate=deg_rate) for i in range(1)]

reporters = [Reporter(name=f'SFP{i}', color='red', degradation_rate=deg_rate) for i in range(2)]

ahl = Supplement(name='IN')
rec = Receiver(input=ahl, output=[repressors[0], reporters[0]], alpha=[0,1], K=0.5, n=2)
                        
not_ = Hill1(input=repressors[0], output=reporters[1], alpha=[1,0], K=0.5, n=4)

net.add_operator([rec, not_])
net.add_reporter(reporters)
net.add_regulator(repressors)
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
f = open(f'result_not_{seed}.csv', 'wt')
f.write('seed,IN,OUT,result,error\n')

for IN in range(2):
    in_phase = (1 - IN) * 180
    sample, prof = setup_sample(ref_period, (1-IN)*np.pi, scale)
    meas = run_assay(sample, ref_period, interval, stochastic=True)
    phases, lags, corrs = compute_phases(net, meas, prof, interval, ref_period)
    out_phase = phases[1]
    result = (out_phase<90)*1
    error = (1 - IN) - result
    print(f'IN={IN}, result={result}, error={error!=0}')
    f.write(f'{seed},{in_phase},{out_phase},{result},{error!=0}\n')

f.close()

# ---

