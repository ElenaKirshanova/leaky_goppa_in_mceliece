import time

load("utils.sage")

McEliece0 = {'label': 0, 'n': 41,  'k': 11, 'w': 5, 'm': 6, 'f': x^6 + x + 1}
McEliece1 = {'label': 1, 'n': 3488, 'k': 2720, 'w': 64, 'm': 12, 'f':x^12+x^3+1}
McEliece2 = {'label': 2,'n': 4608, 'k': 3360, 'w': 96, 'm': 13, 'f':x^13+x^4+x^3+x+1}
McEliece3 = {'label': 3,'n': 6960, 'k': 5413, 'w': 119,'m': 13, 'f':x^13+x^4+x^3+x+1}
McEliece4 = {'label': 4,'n': 8192, 'k': 6528, 'w': 128,'m': 13, 'f':x^13+x^4+x^3+x+1}


param_set = McEliece3

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

N_instances = 10
Ntrials_per_instance = 50
for i_instance in range(N_instances):
    print('running instance #', i_instance, 'on paramset ', param_set['label'])
    L_recovered = 0
    HInv_failures = 0
    L_failures = 0
    false_positive = 0

    H, L, g = gen_instance(param_set,ff,F)
    time_total= 0
    for i_trial in range(Ntrials_per_instance):

        while True:
            I = genI(sizeI, n)
            Lleaked =  genLabr(L, I)
            H_abr = H[[i for i in range(H.nrows())],[j for j in range(H.ncols()) if j in I]]
            if rank(H_abr)==H.nrows():
                break

        #I = genI(sizeI, n)
        #Lleaked =  genLabr(L, I)

        start = time.time()
        Habr, Gabr = genHG_abr(H, I, F)
        kabr = Gabr.rank()
        assert(kabr>0)

        products = create_products(Lleaked, R.gen())
        #print('obtaining candidates for g')
        g_candidates = get_g(Gabr, kabr, products, t, R)

        if len(g_candidates)>1:
            print('False positive g is found')
            #print('correct g:', g)
            false_positive+=1

        for g_cand in g_candidates:
            #if len(g_candidates)>1: print('is g_cand==g?:', g==g_cand)

            Hpriv_candidate = recover_Hpriv_complete(H, g_cand, Lleaked, F, ff, kappa, I)
            if not Hpriv_candidate==0:
                completeL = recover_all_alpha(Lleaked, g_cand, Hpriv_candidate, ff, kappa, k-1)
                if (not completeL==set()) and len(completeL)==k-1:
                    L_recovered+=1
                    #if not g==g_cand:
                    #    print('L remaining:', sorted(set(L).difference(set(Lleaked))))
                    #    print('L returned:', sorted(completeL))
                    end = time.time()
                    time_total+=end-start
                else:
                    print('recover_all_alpha failed, |g_candidates| = ', len(g_candidates), 'completeL:', completeL)
                    L_failures+=1
            else:
                #print('Hpriv_candidate was not found')
                HInv_failures+=1

        if i_trial%10==0: print(i_trial, ' L_recovered:', L_recovered)
    print('# keys recoverd:', L_recovered, '; in averege time:', (time_total/L_recovered), 'HInv_failures:', HInv_failures, 'L_failures:', L_failures, 'false positive g:', false_positive)
