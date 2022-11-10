from loica import *

def add_adder(net, A, B, C, scale, gate_mismatch=0):
    # Place a full adder in GeneticNetwork net with inputs A,B,C
    # A, B = input supplements
    # C = input activator
    deg_rate = 1

    activators = [Regulator(name='Act', degradation_rate=deg_rate) for i in range(6)]

    reporters = [Reporter(name=f'SFP{len(net.reporters)+i}', color='red', degradation_rate=deg_rate) for i in range(2)]

    repressors = [Regulator(name='Rep', degradation_rate=deg_rate) for i in range(5)]

    # Input A
    #ahl1 = Supplement(name='A')
    rec1 = Receiver(input=A, output=[activators[0], repressors[0]], alpha=[0.0125*scale,1.25*scale], K=0.5*scale, n=2)
    # Input B
    #ahl2 = Supplement(name='B')
    rec2 = Receiver(input=B, output=[activators[1]], alpha=[0.0125*scale,1.25*scale], K=0.5*scale, n=2)
    # Input C
    #ahl3 = Supplement(name='C')
    #rec3 = Receiver(input=ahl3, output=[reporters[2], activators[2]], alpha=[0,1.25], K=0.5, n=2)

    # NOT(A)
    not1 = Hill1(input=repressors[0], output=[activators[2]], alpha=[1*scale,0.01*scale], K=np.random.lognormal(sigma=gate_mismatch)*0.5*scale, n=4)

    # MAJ(NOT(A), B, C)
    K = np.random.lognormal(sigma=gate_mismatch, size=(3,))*0.5*scale
    sum1 = Sum(input=[activators[2], activators[1], C], output=[repressors[1]], alpha=[[0.01*scale,1*scale],[0.01*scale,1*scale],[0.01*scale,1*scale]], K=K, n=[2,2,2])

    # MAJ(A, B, C)
    K = np.random.lognormal(sigma=gate_mismatch, size=(3,))*0.5*scale
    sum2 = Sum(input=[activators[0], activators[1], C], output=[repressors[2]], alpha=[[0.01*scale,1*scale],[0.01*scale,1*scale],[0.01*scale,1*scale]], K=K, n=[2,2,2])

    # NOT(MAJ(A, B, C))
    K = np.random.lognormal(sigma=gate_mismatch)*1.25*scale
    not2 = Hill1(input=repressors[2], output=[repressors[3]], alpha=[1*scale,0.01*scale], K=1.25*scale, n=4)

    # NOT(NOT(MAJ(A, B, C))) = CARRY
    K = np.random.lognormal(sigma=gate_mismatch)*0.5*scale
    not3 = Hill1(input=repressors[3], output=[activators[3], reporters[0]], alpha=[1*scale,0.01*scale], K=K, n=4)

    # NOT(MAJ(NOT(A), B, C))
    K = np.random.lognormal(sigma=gate_mismatch)*1.25*scale
    not4 = Hill1(input=repressors[1], output=[activators[4]], alpha=[1*scale,0.01*scale], K=K, n=4)

    # MAJ(NOT(A), NOT(MAJ(NOT(A), B, C)), CARRY)
    K = np.random.lognormal(sigma=gate_mismatch, size=(3,))*0.5*scale
    sum3 = Sum(input=[activators[2], activators[4], activators[3]], output=repressors[4], alpha=[[0.01*scale,1*scale],[0.01*scale,1*scale],[0.01*scale,1*scale]], K=K, n=[2,2,2])

    # NOT(MAJ(NOT(A), NOT(MAJ(NOT(A), B, C)), CARRY)) = SUM
    K = np.random.lognormal(sigma=gate_mismatch)*1.25*scale
    not5 = Hill1(input=repressors[4], output=reporters[1], alpha=[1*scale,0.01*scale], K=K, n=4)

    net.add_operator([rec1, rec2, not1, sum1, sum2, not2, not3, not4, sum3, not5])
    net.add_reporter(reporters)
    net.add_regulator(activators)
    net.add_regulator(repressors) 

    # Return carry signal
    return activators[3]



def growth_rate(t):
    return 0 # gompertz_growth_rate(t, 0.01, 1, 1, 0.5)

def biomass(t):
    return 1 # gompertz(t, 0.01, 1, 1, 0.5)

metab = SimulatedMetabolism('Test', biomass, growth_rate)

def setup_sample(A, B, ref_period):
    sample = Sample(genetic_network=net, 
                        metabolism=metab)

    '''
    def prof_input(f, p):
        def prof(t):
            #return 0.5 + 0.5*np.cos(2 * np.pi * t * f + p)
            period = 1 / f
            tau = p / 2 / np.pi * period
            if (t+tau)%period < period/2:
                return 1
            else:
                return 0
        return prof
    '''
    def prof_input(period, phase=0):
        def f(t):
            p = period
            return ((t-phase*0.5*p/np.pi)%p)/p < 1/2
        return f

    scale = 100
    
    # A
    sample.set_supplement(A1, concentration=scale, profile=prof_input(ref_period, (1 - A&1) * np.pi))
    sample.set_supplement(A2, concentration=scale, profile=prof_input(ref_period, (1 - (A&2)/2) * np.pi))
    sample.set_supplement(A3, concentration=scale, profile=prof_input(ref_period, (1 - (A&4)/4) * np.pi))
    sample.set_supplement(A4, concentration=scale, profile=prof_input(ref_period, (1 - (A&8)/8) * np.pi))
    # B
    sample.set_supplement(B1, concentration=scale, profile=prof_input(ref_period, (1 - B&1) * np.pi))
    sample.set_supplement(B2, concentration=scale, profile=prof_input(ref_period, (1 - (B&2)/2) * np.pi))
    sample.set_supplement(B3, concentration=scale, profile=prof_input(ref_period, (1 - (B&4)/4) * np.pi))
    sample.set_supplement(B4, concentration=scale, profile=prof_input(ref_period, (1 - (B&8)/8) * np.pi))
    # ZERO
    sample.set_supplement(zero, concentration=scale, profile=prof_input(ref_period, np.pi))
    # Reference
    sample.set_supplement(ref, concentration=scale, profile=prof_input(ref_period, 0))

    # Sum
    sum_ = A + B
    #print(sum_)
    #print("{:0b}".format(sum_))
    
    return sample


def run_assay(sample, interval):
    assay = Assay([sample], 
                  n_measurements=int(ref_period*2/interval), 
                  interval=interval,
                  name='Loica phase logic',
                  description='Simulated phase logic circuit generated by loica'
                 )
    assay.run(stochastic=True)
    return assay.measurements

from scipy.signal import correlate, correlation_lags

def compute_phases(meas, interval, ref_period, plot=False):
    phases = np.zeros((len(net.reporters),))
    for i in range(len(net.reporters)):
        ref_sig = meas[meas.Signal=='SFP0']
        m1 = ref_sig[ref_sig.Time>ref_period].sort_values('Time').Measurement.values
        sig = meas[meas.Signal==f'SFP{i}']
        m2 = sig[sig.Time>ref_period].sort_values('Time').Measurement.values
        c = correlate(m1 - m1.mean(), m2 - m2.mean())
        lags = correlation_lags(len(m1), len(m2))
        phase_diff = lags[np.argmax(c)] * interval / ref_period * 360
        phase_diff = phase_diff % 360 
        if plot:
            plt.figure(figsize=(6,2))
            #plt.plot(lags[int(len(lags)/2):], c[int(len(c)/2):])
            plt.plot(lags * interval / ref_period * 360, c)
            plt.plot(phase_diff, c[np.argmax(c)], 'r*')
            plt.xticks([-360, -180, 0, 180, 360])
            plt.yticks([])
        
        phase_diff = abs(phase_diff)
        if phase_diff > 180:
            phase_diff = 360 - phase_diff
        phases[i] = phase_diff

    return phases


# ---
# Construct multi-bit adder network
scale = 100
mismatch = 0

net = GeneticNetwork()

A1 = Supplement(name='A1')
A2 = Supplement(name='A2')
A3 = Supplement(name='A3')
A4 = Supplement(name='A4')

B1 = Supplement(name='B1')
B2 = Supplement(name='B2')
B3 = Supplement(name='B3')
B4 = Supplement(name='B4')

zero = Supplement(name='ZERO')
C = Regulator(name='C', degradation_rate=1)
rec_zero = Receiver(input=zero, output=C, alpha=[0.0125*scale,1.25*scale], K=0.5*scale, n=2)
net.add_regulator(C)
net.add_operator(rec_zero)

ref = Supplement(name='REF')
rep_ref = Reporter(name='SFP0', degradation_rate=1)
rec_ref = Receiver(input=ref, output=rep_ref, alpha=[0.0125*scale,1.25*scale], K=0.5*scale, n=2)
net.add_reporter(rep_ref)
net.add_operator(rec_ref)

carry1 = add_adder(net, A1, B1, C, scale=scale, gate_mismatch=mismatch)
carry2 = add_adder(net, A2, B2, carry1, scale=scale, gate_mismatch=mismatch)
carry3 = add_adder(net, A3, B3, carry2, scale=scale, gate_mismatch=mismatch)
carry4 = add_adder(net, A4, B4, carry3, scale=scale, gate_mismatch=mismatch)
carries = [carry1, carry2, carry3, carry4]

#carry_reps = [Reporter(name=f'SFP{i}', degradation_rate=1) for i in range(len(net.reporters), len(net.reporters)+3)]
#carry_buffs = [Hill1(input=carries[i], output=carry_reps[i], alpha=[0,1], K=0.5, n=2) for i in range(3)]
#net.add_reporter(carry_reps)
#net.add_operator(carry_buffs)

# ---
# Test network
# Get inputs from job id
import sys
AB = int(sys.argv[1]) - 1
A = (AB & 0b11110000) >> 4
B = AB & 0b00001111


ref_period = 120
interval = 0.24
errors = 0
f = open(f'result_multi_bit_adder_{A}_{B}.csv', 'wt')
f.write('A,B,sum,result,error\n')

sample = setup_sample(A, B, ref_period)
meas = run_assay(sample, interval)
phases = compute_phases(meas, interval, ref_period)
result = (phases[2]<90)*1 + (phases[4]<90)*2 + (phases[6]<90)*4 + (phases[8]<90)*8 + (phases[7]<90)*16
error = A + B - result
print(f'A={A}, B={B}, A+B={A+B}, result={result}, error={error!=0}')
f.write(f'{A},{B},{A+B},{result},{error!=0}\n')

f.close()

# ---
