import time

load("utils.sage")

McEliece0 = {'label': 0, 'n': 41,  'k': 11, 'w': 5, 'm': 6, 'f': x^6 + x + 1}
McEliece1 = {'label': 1, 'n': 3488, 'k': 2720, 'w': 64, 'm': 12, 'f':x^12+x^3+1}
McEliece2 = {'label': 2,'n': 4608, 'k': 3360, 'w': 96, 'm': 13, 'f':x^13+x^4+x^3+x+1}
McEliece3 = {'label': 3,'n': 6960, 'k': 5413, 'w': 119,'m': 13, 'f':x^13+x^4+x^3+x+1}
McEliece4 = {'label': 4,'n': 8192, 'k': 6528, 'w': 128,'m': 13, 'f':x^13+x^4+x^3+x+1}


param_set = McEliece1

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
sizeI = n-k+2

N_instances = 1
Ntrials_per_instance = 1
for i_instance in range(N_instances):
    print('running instance #', i_instance, 'on paramset ', param_set['label'])
    truly_equiv_g = 0
    g_recovered = 0
    false_positive = 0
    H, L, g = gen_instance(param_set,ff,F)
    #print('g = ', g)

    for i_trial in range(Ntrials_per_instance):
        start = time.time()
        I = genI(sizeI, n)
        Lleaked =  genLabr(L, I)

        #print('generating Habr, Gabr')
        Habr, Gabr = genHG_abr(H, I, F)
        kabr = Gabr.rank()
        assert(kabr>0)

        #print('generating products...')
        products = create_products(Lleaked, R.gen())
        #print('obtaining candidates for g')
        g_candidates = get_g(Gabr, kabr, products, t, R)
        #g_candidates = get_gSet_plain_syndrome(Gabr, t, Lleaked, R.gen())
        was_g_in = g in g_candidates
        end = time.time()
        print('# candidates:', len(g_candidates), '; is correcnt g found?', was_g_in, '; in time:', end-start)
        if was_g_in==True:
            g_recovered+=1
        if len(g_candidates)>1:
            false_positive+=1
            print('False positive is found')
            """
            print('computing G from H...')
            G = getGfromH(H)
            print('creating the product tree...')
            products_all = create_products(L, R.gen())

            for g_cand in g_candidates:
                if g_cand == g: continue
                print('g_cand:', g_cand)
                test_g_full = test_g(g_cand, G, products_all, R)
                test_g_abr  = test_g(g_cand, Gabr, products, R)
                if (test_g_abr == 1 and test_g_full==0):
                    false_positive+=1
                    print('false positive g detected', 'rank Gabr:', Gabr.rank())
                elif (test_g_abr == 1 and test_g_full==1 and not g_cand==g):
                    truly_equiv_g+=1
                    print('equivalent g detected', 'rank Gabr:', Gabr.rank())
                else:
                    print('g = ', g, 'g_cand = ', g)
                    print('UNPREDICTED EVENT: test_g_abr:', test_g_abr, 'test_g_full:', test_g_full)
            assert(False)
            """
    print('g_recovered:', g_recovered, 'false_positive:', false_positive, 'truly_equiv_g:', truly_equiv_g)
