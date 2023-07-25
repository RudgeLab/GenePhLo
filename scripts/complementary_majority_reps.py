from loica import *
from genephlo import *

def setup_sample(ref_period, input_phases, scale):
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
    for phase,ahl in zip(input_phases,ahls):
        sample.set_supplement(ahl, concentration=scale, profile=prof_input(ref_period, phase))

    return sample, prof_input(ref_period, 0)

scale = 100

net = GeneticNetwork()

deg_rate = 1

activators = [Regulator(name=f'Act{i}', degradation_rate=deg_rate) for i in range(3)]
repressors = [Regulator(name=f'Act{i}', degradation_rate=deg_rate) for i in range(1)]

reporters = [Reporter(name=f'SFP{i}', color='red', degradation_rate=deg_rate) for i in range(5)]

ahls = [Supplement(name=f'IN{i}') for i in range(3)]
recs = [Receiver(input=ahls[i], output=[activators[i], reporters[i]], alpha=[0,1], K=1, n=2) for i in range(3)]
                        
sum_ = Sum(input=activators, output=[reporters[3], repressors[0]], alpha=[[0,1],[0,1],[0,1]], K=[1,1,1], n=[2,2,2])

not_ = Hill1(input=repressors[0], output=reporters[4], alpha=[1,0], K=1, n=4)

net.add_operator(recs + [sum_, not_])
net.add_reporter(reporters)
net.add_regulator(activators)
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
f = open(f'result_complementary_majority_{seed}.csv', 'wt')
f.write('seed,IN1,IN2,IN3,OUT,result,error\n')

for IN1 in range(2):
    for IN2 in range(2):
        for IN3 in range(2):
            input_phases = np.array([(1 - IN1), (1 - IN2), (1 - IN3)])
            sample, prof = setup_sample(ref_period, input_phases*np.pi, scale)
            meas = run_assay(sample, ref_period, interval, stochastic=True)
            phases, lags, corrs = compute_phases(net, meas, prof, interval, ref_period)
            out_phase = phases[4]
            result = (out_phase<90)*1
            error = ((IN1 + IN2 + IN3) > 1)*1 - result
            print(f'IN1={IN1}, IN2={IN2}, IN3={IN3}, result={result}, error={error!=0}')
            f.write(f'{seed},{input_phases[0]*180}, {input_phases[1]*180}, {input_phases[2]*180}, {out_phase},{result},{error!=0}\n')

f.close()

# ---




