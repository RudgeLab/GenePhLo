from loica import *
import numpy as np

def add_adder(net, A, B, C, deg_rate):
    # Place a full adder in GeneticNetwork net with inputs A,B,C
    # A, B = input supplements
    # C = input activator
    # Parameters of operators are place holders, which should be optimized using optimize()

    activators = [Regulator(name='Act', degradation_rate=deg_rate) for i in range(6)]

    reporters = [Reporter(name=f'SFP{len(net.reporters)+i}', color='red', degradation_rate=deg_rate) for i in range(2)]

    repressors = [Regulator(name='Rep', degradation_rate=deg_rate) for i in range(5)]

    # Input A
    #ahl1 = Supplement(name='A')
    rec1 = Receiver(input=A, output=[activators[0], repressors[0]], alpha=[0.0125,1.25], K=0.5, n=2)
    # Input B
    #ahl2 = Supplement(name='B')
    rec2 = Receiver(input=B, output=[activators[1]], alpha=[0.0125,1.25], K=0.5, n=2)
    # Input C
    #ahl3 = Supplement(name='C')
    #rec3 = Receiver(input=ahl3, output=[reporters[2], activators[2]], alpha=[0,1.25], K=0.5, n=2)

    # NOT(A)
    not1 = Hill1(input=repressors[0], output=[activators[2]], alpha=[1,0.01], K=1, n=4)

    # MAJ(NOT(A), B, C)
    K = [1,1,1] 
    sum1 = Sum(input=[activators[2], activators[1], C], output=[repressors[1]], alpha=[[0.01,1],[0.01,1],[0.01,1]], K=K, n=[2,2,2])

    # MAJ(A, B, C)
    K = [1,1,1] 
    sum2 = Sum(input=[activators[0], activators[1], C], output=[repressors[2]], alpha=[[0.01,1],[0.01,1],[0.01,1]], K=K, n=[2,2,2])

    # NOT(MAJ(A, B, C))
    K = 1
    not2 = Hill1(input=repressors[2], output=[repressors[3]], alpha=[1,0.01], K=1.25, n=4)

    # NOT(NOT(MAJ(A, B, C))) = CARRY
    K = 1
    not3 = Hill1(input=repressors[3], output=[activators[3], reporters[0]], alpha=[1,0.01], K=K, n=4)

    # NOT(MAJ(NOT(A), B, C))
    K = 1
    not4 = Hill1(input=repressors[1], output=[activators[4]], alpha=[1,0.01], K=K, n=4)

    # MAJ(NOT(A), NOT(MAJ(NOT(A), B, C)), CARRY)
    K = [1,1,1]
    sum3 = Sum(input=[activators[2], activators[4], activators[3]], output=repressors[4], alpha=[[0.01,1],[0.01,1],[0.01,1]], K=K, n=[2,2,2])

    # NOT(MAJ(NOT(A), NOT(MAJ(NOT(A), B, C)), CARRY)) = SUM
    K = 1
    not5 = Hill1(input=repressors[4], output=reporters[1], alpha=[1,0.01], K=K, n=4)

    net.add_operator([rec1, rec2, not1, sum1, sum2, not2, not3, not4, sum3, not5])
    net.add_reporter(reporters)
    net.add_regulator(activators)
    net.add_regulator(repressors) 

    # Return carry signal
    return activators[3]


def construct_multi_bit_adder(N, deg_rate):
    # ---
    # Construct multi-bit adder network
    # N = number of bits to add
    net = GeneticNetwork()

    A = [Supplement(name=f'A{i+1}') for i in range(N)]

    B = [Supplement(name=f'B{i+1}') for i in range(N)]

    C = Supplement(name='C')
    carry = Regulator(name='C', degradation_rate=deg_rate)
    rec_C = Receiver(input=C, output=carry, alpha=[0.01,1], K=1, n=2)
    net.add_regulator(carry)
    net.add_operator(rec_C)

    ref = Supplement(name='REF')
    rep_ref = Reporter(name='SFP0', degradation_rate=deg_rate)
    rec_ref = Receiver(input=ref, output=rep_ref, alpha=[0.01,1], K=1, n=2)
    net.add_reporter(rep_ref)
    net.add_operator(rec_ref)

    for i in range(N):
        carry = add_adder(net, A[i], B[i], carry, deg_rate=deg_rate)

    return net, A, B, C, ref

def setup_adder_sample(net, N, ref, A_in, B_in, C_in, A, B, C, ref_period, scale=100):
    def growth_rate(t):
        return 0 # gompertz_growth_rate(t, 0.01, 1, 1, 0.5)

    def biomass(t):
        return 1 # gompertz(t, 0.01, 1, 1, 0.5)

    metab = SimulatedMetabolism('Test', biomass, growth_rate)
    sample = Sample(genetic_network=net, 
                        metabolism=metab)

    def prof_input(period, phase=0):
        def f(t):
            p = period
            return ((t-phase*0.5*p/np.pi)%p)/p < 1/2
        return f

    for i in range(N):
        val = 2**i
        # A
        sample.set_supplement(A_in[i], concentration=scale, profile=prof_input(ref_period, (1 - (A & val)/val) * np.pi))
        # B
        sample.set_supplement(B_in[i], concentration=scale, profile=prof_input(ref_period, (1 - (B & val)/val) * np.pi))

    # C
    sample.set_supplement(C_in, concentration=scale, profile=prof_input(ref_period, (1 - C&1) * np.pi))
    # Reference
    sample.set_supplement(ref, concentration=scale, profile=prof_input(ref_period, 0))
    
    return sample, prof_input(ref_period, 0)


