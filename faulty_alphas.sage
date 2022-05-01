import time

load("utils.sage")

McEliece0 = {'label': 0, 'n': 41,  'k': 5, 'w': 6, 'm': 6, 'f': x^6 + x + 1} #toy example
McEliece0 = {'label': 0, 'n': 41,  'k': 11, 'w': 5, 'm': 6, 'f': x^6 + x + 1}
McEliece1 = {'label': 1, 'n': 3488, 'k': 2720, 'w': 64, 'm': 12, 'f':x^12+x^3+1}
McEliece2 = {'label': 2,'n': 4608, 'k': 3360, 'w': 96, 'm': 13, 'f':x^13+x^4+x^3+x+1}
McEliece3 = {'label': 3,'n': 6960, 'k': 5413, 'w': 119,'m': 13, 'f':x^13+x^4+x^3+x+1}
McEliece4 = {'label': 4,'n': 8192, 'k': 6528, 'w': 128,'m': 13, 'f':x^13+x^4+x^3+x+1}


def spoil_L(L,p,ff):
	"""
		changes p elements of L to random field elements
	"""
	spoiled_positions = set()
	n = len(L)
	Lspoiled = copy(L)
	for i in range(p):
		while True:
			pos_tmp = ZZ.random_element(0, n)
			if not pos_tmp in spoiled_positions:
				new_element = ff.random_element()
				if not Lspoiled[pos_tmp] == new_element:
					Lspoiled[pos_tmp] = new_element
					break
	return Lspoiled

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


nspoiled = 1

N_instances = 100
L_recovered = 0
HInv_failures = 0
L_failures = 0
L_failures_set_differ = 0
false_positive = 0
empty_g_set = 0
empty_L_set = 0
wrong_g_found = 0
not_enough_alphas = 0
wrong_g_and_empty_L = 0
for i_instance in range(N_instances):
	#print('running instance #', i_instance, 'on paramset ', param_set['label'], 'with nspoiled = ', nspoiled)
	H, L, g = gen_instance(param_set,ff,F)
	Lset = set(L)

	while True:
		I = genI(n-k+1, n)
		Lleaked =  genLabr(L, I)
		H_abr = H[[i for i in range(H.nrows())],[j for j in range(H.ncols()) if j in I]]
		if rank(H_abr)==H.nrows():
			break


	Lspoiled = spoil_L(Lleaked,nspoiled,ff)
	Habr, Gabr = genHG_abr(H, I, F)
	kabr = Gabr.rank()
	assert(kabr>0)

	products = create_products(Lspoiled, R.gen())
	#print('obtaining candidates for g')
	g_candidates = get_g(Gabr, kabr, products, t, R)

	was_g_in = g in g_candidates
	if (not was_g_in) and len(g_candidates)>0: wrong_g_found+=1

	if len(g_candidates)==0:
		#print('no g_candidates found')
		empty_g_set+=1

	if len(g_candidates)>1:
		#print('False positive g is found')
		#print('correct g:', g)
		false_positive+=1

	for g_cand in g_candidates:
		#if len(g_candidates)>1: print('is g_cand==g?:', g==g_cand)

		Hpriv_candidate = recover_Hpriv_complete(H, g_cand, Lspoiled, F, ff, kappa, I)
		if not Hpriv_candidate==0:
			completeL = recover_all_alpha(Lspoiled, g_cand, Hpriv_candidate, ff, kappa, k-1)
			if (not completeL==set()):
				if (len(completeL)==k-1):
					if completeL==set(L).difference(Lleaked):
						L_recovered+=1
					else:
						L_failures_set_differ+=1
				else:
					not_enough_alphas+=1
			else:
				empty_L_set+=1
				if not was_g_in: wrong_g_and_empty_L+=1
		else:
			HInv_failures+=1

	if i_instance%10==0: print(i_instance, 'L_recovered:', L_recovered, 'wrong_g_found:', wrong_g_found, 'empty_g_set:', empty_g_set, 'false_positive:', false_positive, 'L_failures_set_differ:', L_failures_set_differ, 'not_enough_alphas:', not_enough_alphas, 'empty_L_set:', empty_L_set, 'wrong_g_and_empty_L:', wrong_g_and_empty_L)

print('Total stats: keys recoverd:', L_recovered, 'wrong_g_found:', wrong_g_found, 'empty_g_set:', empty_g_set, 'false_positive:', false_positive, 'L_failures_set_differ:', L_failures_set_differ, 'not_enough_alphas:', not_enough_alphas, 'empty_L_set:', empty_L_set, 'wrong_g_and_empty_L:', wrong_g_and_empty_L)
