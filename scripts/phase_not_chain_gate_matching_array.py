import numpy as np
from loica import *


def build_network(N=1, gate_mismatch=1, scale=100, seed=0):
    net = GeneticNetwork()

    deg_rate = 1

    np.random.seed(seed)

    repressors = [Regulator(name=f'Rep{i}', degradation_rate=deg_rate) for i in range(N)]

    reporters = [Reporter(name=f'SFP{i}', color='red', degradation_rate=deg_rate) for i in range(2)]

    ahl = Supplement(name='IN')
    rec = Receiver(input=ahl, output=[repressors[0], reporters[0]], alpha=[0.0125*scale,1.25*scale], K=0.5*scale, n=2)

    nots = [Hill1(input=repressors[i], output=repressors[i+1], alpha=[1*scale,0.01*scale], K=np.random.lognormal(sigma=gate_mismatch)*0.1*scale, n=4) for i in range(N-1)]

    not_out = Hill1(input=repressors[N-1], output=reporters[1], alpha=[1*scale,0.01*scale], K=np.random.lognormal(sigma=gate_mismatch)*0.1*scale, n=4)

    net.add_operator(nots)
    net.add_operator(not_out)
    net.add_operator(rec)
    net.add_reporter(reporters)
    net.add_regulator(repressors)
    
    return net, ahl

def find_downstream_ops(op, ops):
    outputs = op.output
    if type(outputs)!=list:
        outputs = [outputs]
    downstream_ops = []
    for out in outputs:
        for op2 in ops:
            op_inputs = op2.input
            if type(op_inputs)!=list:
                op_inputs = [op_inputs]
            if out in op_inputs:
                downstream_ops.append((out,op,op2))
    return downstream_ops

def phi1(x, delta, n):
    return (1/delta + x**(n/2)) / (1 + x**(n/2))

def phi2(x, delta, n):
    return (1/delta + x**(-n/2)) / (1 + x**(-n/2))

def optimize_Ks(net, gate_mismatch=0):
    d = 100
    alpha0 = 500
    regs = net.regulators
    ops = net.operators
    recs = [op for op in ops if type(op)==Receiver]
    # Set up a list in each operator to store upstream gate output dynamic ranges
    for op in ops:
        if type(op)==Hill1:
            op.din = -1
            op.dout = -1
        elif type(op)==Sum:
            op.din = np.zeros((3,)) - 1
            op.dout = -1

    # Start from receivers
    downstream_ops = []
    for rec in recs:
        # Set optimal parameters assuming large input dynamic range
        rec.K = 1
        rec.Kout = np.sqrt(alpha0*alpha0/d)
        rec.din = d
        rec.dout = d
        rec.alpha = [alpha0/d, alpha0]
        # Find downstream operators
        downstream_ops += find_downstream_ops(rec, ops)
        
    while len(downstream_ops)!=0:
        new_downstream_ops = []
        # Set optimal parameters for downstream ops 
        for reg,upop,dop in downstream_ops:
            new_downstream_ops += find_downstream_ops(dop, ops)
            if type(dop)==Hill1:
                if upop.dout!=-1 and dop.dout==-1:
                    dop.din = upop.dout
                    pon = phi1(dop.din, d, dop.n)
                    poff = phi2(dop.din, d, dop.n)
                    alpha = 1 / np.sqrt(pon * poff)
                    dop.Kout = alpha0 / alpha
                    dop.K = upop.Kout * np.random.lognormal(sigma=gate_mismatch)
                    dop.alpha = [alpha0, alpha0/d] # [alpha, alpha/d]
                    dop.dout = pon/poff
                    print('NOT', dop.din, dop.dout, pon, poff, dop.Kout)
            elif type(dop)==Sum:
                idx = np.where(np.array(dop.input)==reg)
                idx = idx[0][0]
                if upop.dout!=-1:
                    dop.din[idx] = upop.dout
                    dop.alpha[idx] = [alpha0/d, alpha0] # [alpha/d, alpha]
                    dop.K[idx] = upop.Kout * np.random.lognormal(sigma=gate_mismatch)
                if not np.any(dop.din==-1) and dop.dout==-1:
                    p1A = phi1(dop.din[0], d, dop.n[idx])
                    p1B = phi1(dop.din[1], d, dop.n[idx])
                    p1C = phi1(dop.din[2], d, dop.n[idx])
                    p1s = np.sort([p1A, p1B, p1C])
                    p2A = phi2(dop.din[0], d, dop.n[idx])
                    p2B = phi2(dop.din[1], d, dop.n[idx])
                    p2C = phi2(dop.din[2], d, dop.n[idx])
                    p2s = np.sort([p2A, p2B, p2C])
                    
                    pon = p1s[:2].sum() + p2s[-1]
                    poff = p1s[-1] + p2s[:2].sum()
                    alpha = 1 / np.sqrt( pon * poff )
                    dop.Kout = alpha0 / alpha
                    dop.dout = pon / poff
                    print('SUM', dop.din, dop.dout, pon, poff, dop.Kout)
        downstream_ops = new_downstream_ops

def growth_rate(t):
    return 0 # gompertz_growth_rate(t, 0.01, 1, 1, 0.5)

def biomass(t):
    return 1 # gompertz(t, 0.01, 1, 1, 0.5)

metab = SimulatedMetabolism('Test', biomass, growth_rate)

def prof_input(period, phase=0):
    def f(t):
        p = period
        return ((t-phase*0.5*p/np.pi)%p)/p < 1/2
    return f
    
def define_sample(net, ahl, phase_based=True, ref_period=200, scale=100, logic=0):
    sample = Sample(genetic_network=net, 
                        metabolism=metab)

    if phase_based:
        prof = prof_input(ref_period, (1-logic)*np.pi)
        conc = scale
    else:
        prof = None
        conc = scale*logic

    # Input
    sample.set_supplement(ahl, concentration=conc, profile=prof)
    return sample

from scipy.signal import correlate, correlation_lags

def compute_phases(meas, interval, ref_period, plot=False):
    phases = np.zeros((len(net.reporters),))
    for i in range(len(net.reporters)):
        ref_sig = meas[meas.Signal=='SFP0']
        m1 = ref_sig[ref_sig.Time>ref_period].sort_values('Time').Measurement.values
        sig = meas[meas.Signal==f'SFP{i}']
        m2 = sig[sig.Time>ref_period].sort_values('Time').Measurement.values
        c = correlate(m1 - m1.mean(), m2 - m2.mean())
        c = c[int(len(c)/2):]
        lags = correlation_lags(len(m1), len(m2))
        lags = lags[int(len(lags)/2):]
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

import sys
if len(sys.argv)>=2:
    seeds = np.array(sys.argv[1:]).astype(int)
else:
    seeds = [0]

interval = 0.24
ref_period = 120
scale = 100

err0_phase_list = []
err1_phase_list = []
err0_level_mean_list = []
err1_level_mean_list = []
err0_level_list = []
err1_level_list = []
mismatch = np.log(2)/2
Ns = [1, 3, 5, 7, 9] #range(1, 6)
for N in Ns:
    errs0_phase = 0
    errs1_phase = 0
    errs0_level_mean = 0
    errs1_level_mean = 0
    errs0_level = 0
    errs1_level = 0
    for seed in seeds:
        net, ahl = build_network(N=N, gate_mismatch=mismatch, scale=scale, seed=seed)
        optimize_Ks(net, mismatch)
        
        # Phase-based
        sample0 = define_sample(net, ahl, phase_based=True, ref_period=ref_period, scale=scale, logic=0)
        sample1 = define_sample(net, ahl, phase_based=True, ref_period=ref_period, scale=scale, logic=1)
        assay = Assay([sample0, sample1], 
                      n_measurements=int(2*ref_period/interval), 
                      interval=interval,
                      name='Loica phase logic',
                      description='Simulated phase logic circuit generated by loica'
                     )
        assay.run(stochastic=True)
        meas = assay.measurements
        phases0 = compute_phases(meas[meas.Sample==0], interval, ref_period, plot=False)
        phases1 = compute_phases(meas[meas.Sample==1], interval, ref_period, plot=False)
        print(phases0, phases1)
        if phases0[1]<90:
            errs0_phase += 1
        if phases1[1]<90:
            errs1_phase += 1
            
        # Level-based
        sample0 = define_sample(net, ahl, phase_based=False, scale=scale, logic=0)
        sample1 = define_sample(net, ahl, phase_based=False, scale=scale, logic=1)
        assay = Assay([sample0, sample1], 
                      n_measurements=int(2*ref_period/interval), 
                      interval=interval,
                      name='Loica phase logic',
                      description='Simulated phase logic circuit generated by loica'
                     )
        assay.run(stochastic=True)
        meas = assay.measurements
        level0 = meas[meas.Signal=='SFP1'][meas.Time>ref_period][meas.Sample==0].Measurement
        level1 = meas[meas.Signal=='SFP1'][meas.Time>ref_period][meas.Sample==1].Measurement

        #print(level0, level1)
        if N%2==0:
            if level0.mean()>=scale/2:
                errs0_level_mean += 1
                print(f"Error mean(level0) = {level0.mean()}")
            if level1.mean()<scale/2:
                errs1_level_mean += 1            
                print(f"Error mean(level1) = {level1.mean()}")
            if level0.max()>=scale/2:
                errs0_level += 1
                print(f"Error max(level0) = {level0.max()}")
            if level1.min()<scale/2:
                errs1_level += 1            
                print(f"Error min(level1) = {level1.min()}")
        else:
            if level0.mean()<scale/2:
                errs0_level_mean += 1
                print(f"Error mean(level0) = {level0.mean()}")
            if level1.mean()>=scale/2:
                errs1_level_mean += 1
                print(f"Error mean(level1) = {level1.mean()}")
            if level0.min()<scale/2:
                errs0_level += 1
                print(f"Error min(level0) = {level0.min()}")
            if level1.max()>=scale/2:
                errs1_level += 1            
                print(f"Error max(level1) = {level1.max()}")
            
    err0_phase_list.append(errs0_phase)        
    err1_phase_list.append(errs1_phase)
    err0_level_mean_list.append(errs0_level_mean)        
    err1_level_mean_list.append(errs1_level_mean)
    err0_level_list.append(errs0_level)        
    err1_level_list.append(errs1_level)
    
fout = open(f'not_chain_seeds_{seeds}.csv', 'wt')
fout.write('N,seed,err0_phase,err1_phase,err0_level,err1_level,err0_level_mean,err1_level_mean\n')
for i in range(len(Ns)):
    fout.write(f'{Ns[i]}, {seeds}, {err0_phase_list[i]}, {err1_phase_list[i]}, {err0_level_list[i]}, {err1_level_list[i]}, {err0_level_mean_list[i]}, {err1_level_mean_list[i]}\n')
fout.close()
