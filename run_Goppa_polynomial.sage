import time

load("utils.sage")

McEliece0 = {'label': 0, 'n': 41,  'k': 11, 'w': 5, 'm': 6, 'f': x^6 + x + 1}
McEliece1 = {'label': 1, 'n': 3488, 'k': 2720, 'w': 64, 'm': 12, 'f':x^12+x^3+1}
McEliece2 = {'label': 2,'n': 4608, 'k': 3360, 'w': 96, 'm': 13, 'f':x^13+x^4+x^3+x+1}
McEliece3 = {'label': 3,'n': 6960, 'k': 5413, 'w': 119,'m': 13, 'f':x^13+x^4+x^3+x+1}
McEliece4 = {'label': 4,'n': 8192, 'k': 6528, 'w': 128,'m': 13, 'f':x^13+x^4+x^3+x+1}


param_set = McEliece4

R = ZZ['x']
p = 2
kappa = param_set['m']
n_max = p**kappa - 1
F = GF(p)
Fx = PolynomialRing(F, 'x')
defining_poly = param_set['f']
print('defining_poly:', defining_poly)
ff.<a> = FiniteField(p**kappa, modulus = defining_poly)
R = PolynomialRing(ff, 'x')
x = R.gen()

n = param_set['n']
k = param_set['k']
t = param_set['w']
sizeI = n-k+1
print('sizeI:', sizeI)

N_instances = 20
Ntrials_per_instance = 100
for i_instance in range(N_instances):
    print('running instance #', i_instance, 'on paramset ', param_set['label'])
    L_recovered = 0
    HInv_failures = 0
    L_failures = 0

    H, L, g = gen_instance(param_set,ff,F)
    Lset = set(L)
    time_total= 0
    for i_trial in range(Ntrials_per_instance):
        I = genI(sizeI, n)
        Lleaked =  genLabr(L, I)
        start = time.time()
        """
        if not someL==0:
            if someL==L:
                #print('L recovered!')
                L_recovered+=1
            else: print('wrong someL')

        """
        Hpriv_candidate = recover_Hpriv_complete(H, g, Lleaked, n-k, F, ff, kappa, I)
        if not Hpriv_candidate==0:
            completeL = recover_all_alpha(Lleaked, g, Hpriv_candidate, ff, kappa, k-1)
            if (not completeL==set()):
                L_recovered+=1
                end = time.time()
                time_total+=end-start
            else:
                print('recover_all_alpha failed')
                L_failures+=1
        else:
            print('Hpriv_candidate was not found')
            HInv_failures+=1
        """
        completeL = recoverL(H,I,g,Lleaked,F,R)
        end = time.time()
        if not completeL==0:
            print('L recovered!')
            L_recovered+=1

        Hpriv_candidate = recover_Hpriv_complete(H, g, Lleaked, n-k, F, ff, kappa, I)
        if not Hpriv_candidate==0:
            assert(Hpriv_candidate==Hpriv_real), 'Hpriv_candidate!=Hpriv_real'
            completeL = recover_all_alpha(Lleaked, g, Hpriv_candidate, ff, kappa)
            if not completeL==0:
                L_recovered+=1
            else: print('recover_all_alpha failed')
        else: print('Hpriv_candidate was not found')
        """

        if i_trial%10==0: print(i_trial, ' L_recovered:', L_recovered)
    print('# L recovered:', L_recovered, '; in averege time:', (time_total/L_recovered), 'HInv_failures:', HInv_failures, 'L_failures:', L_failures)
