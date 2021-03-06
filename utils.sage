"""
	helper function for the project
	Decoding McEliece with a Hint -- Secret Goppa Key Parts Reveal Everything
	https://eprint.iacr.org/2022/525
"""


def constructL(ff, n):
	"""
		construct a set of size of the finite field ff
	"""
	L = set()
	while len(L)<n:
		while True:
			tmp = ff.random_element()
			if not tmp in L:
				L.add(tmp)
				break
	return list(L)

def gen_parity(L, g):
	"""
		generates the parity check matrix of GoppaCode(L, g)
		over the extension field
	"""
	n = len(L)
	glist = list(g)
	t  = g.degree()
	G = matrix(parent(glist[0]),t)
	for i in range(t):
		for k in range(i+1):
			G[i,k] = glist[t-k]
	V = matrix.vandermonde(L).transpose()

	V1 = V[0:t]
	Gdiag = matrix.diagonal([1/g(L[i]) for i in range(n)])
	#print(G)
	return V1*Gdiag

def gen_parity_full(L, g, kappa, n, F, ff, echelonForm=True):
	"""
		generates the parity check matrix of GoppaCode(L, g)
		over the base (binary) field
	"""
	Hbar = gen_parity(L,g)
	V, from_V, to_V = ff.vector_space(F, map=True)
	H = matrix(F, kappa*Hbar.nrows(), n)
	for i in range(Hbar.nrows()):
		for j in range(n):
			tmp_vec = to_V(Hbar[i,j])
			for k in range(kappa):
				H[kappa*i+k, j] = tmp_vec[k]
	if echelonForm: H = H.echelon_form()
	return H

def gen_error(nerrs, n):
	"""
		generate a size-n binary vector with nerrs-many 1's
	"""
	error_vec = [0]*n
	error_pos = set()
	for i in range(nerrs):
		while True:
			error_pos_tmp = ZZ.random_element(0, n)
			if not error_pos_tmp in error_pos:
				error_pos.add(error_pos_tmp)
				break
		error_vec[error_pos_tmp] = 1
	return vector(GF(2),error_vec)

def getGfromH(H):
	"""
		generator matrix from parity check matrix
	"""
	G = matrix(H.transpose().kernel().basis())
	return G

def gen_message(k):
	"""
		generate a random binary vector of length k
	"""
	m = [0]*k
	for i in range(k):
		m[i] = ZZ.random_element(0,2)
	return vector(GF(2),m)

def compute_syndrome2(y,products):
	"""
		compute the syndrome of y using precomputed products
	"""
	n = len(y)
	s = 0
	for i in range(n):
		s +=y[i]*products[i]
	return s

def compute_syndrome_plain(y,L,x):
	"""
		compute the syndrome of y unsing L and inderminant x
	"""
	n = len(y)
	assert(n==len(L))
	s = 0
	for i in range(n):
		if y[i]==1:
			tmp_prod = 1
			for j in range(n):
				if not j==i: tmp_prod*=(x-L[j])
			s+=tmp_prod
	return s

def genLabr(L, I):
	"""
		returns the subset of L indexed by I
	"""
	T = len(I)
	Labr = [0]*(T)
	counter = 0
	for j in range(n):
		if j in I:
			Labr[counter] = L[j]
			counter+=1
	return Labr

def genHG_abr(H, I, F):
	"""
		generate the submatrix of H taking only columns indexed by I
		generate the corresponding abridged generator matrix
	"""
	T = len(I)
	n = H.ncols()
	Habr = matrix(F, H.nrows(), T)
	counter = 0
	for j in range(n):
		if j in I:
			for i in range(H.nrows()):
				Habr[i,counter] = H[i,j]
			counter+=1
	Gabr = matrix(Habr.transpose().kernel().basis())
	#kabr = Gabr.rank()
	return Habr, Gabr

def create_product_all(L, x):
	"""
		generate the product tree of elements from L
		similar to Fast Multipoint Evaluation
		(see von zu Gathen "Modern Computer Algerbra")
		and return prod_i (x-alpha_i)
	"""
    n = len(L)
    depth_max = ceil(log(n,2))
    #print('depth_max:', depth_max)
    current_depth = 0
    path = [None]*(depth_max)
    ind = 0
    while current_depth<depth_max-1:
        if ind<n-1:
            tmp = (x - L[ind])*(x-L[ind+1])
            ind+=2
            #i = 0
        elif ind == n-1:
            tmp = (x - L[ind])
            ind+=1
            #i = 0
        else:
            tmp = 1
            #if current_depth>0: path[current_depth-1] = 1
        i = 0
        #print('ind:', ind)
        while not path[i] == None:
            tmp *= path[i]
            i+=1
        path[i] = tmp
        current_depth = i
        #print('i:', i)
        if i>0:
            for j in range(i): path[j] = None

        #print(i, path, ind, current_depth)
    return(path[depth_max-1])

def create_products(L, x):
	"""
		generate product of the form prod_{i\neq j} (x-alpha_i) for all j
		using the  prod_i (x-alpha_i)
	"""
	all_prod = create_product_all(L, x)
	prods = [all_prod/(x-L[i]) for i in range(len(L))]
	return prods

def genI(T, n):
	'''
	generate I --  set of indicies of unknown elements in L of size T
	'''
	I = set()
	for i in range(T):
		while True:
			pos_tmp = ZZ.random_element(0, n)
			if not pos_tmp in I:
				I.add(pos_tmp)
				break
	return I

def test_g(g_candidate, G, prods, R):
	"""
	for each codeword represented by a row of G,
	compute its syndrome and check if g_candidate is divisible by it.
	This eliminates faulty Goppa polynomial g_candidate
	"""
	k = G.nrows()
	zerov = vector([0]*G.ncols())
	Rquo = R.quotient(g_candidate, 'x1')
	for i in range(k):
		c = G[i]
		assert(not c==zerov)
		syndrome = compute_syndrome2(c, prods)
		if not Rquo(syndrome) == 0:
			return 0
	return 1


def get_gSet(G, k, prods, t):
	"""
		returns a set of possible Goppa polynomial
		computed as factors of syndromes of G rows
		We do not check for multiplicity (it should be 2)
		of the factors, we check only for the degree. Any factor of degree >1
		appears with multiplicity 2
	"""
	assert(G.nrows()>0)
	c1 = G[0]
	syndrome = compute_syndrome2(c1, prods)
	g_out_factors = list(syndrome.factor())
	print('g_out_factors:', g_out_factors)
	Gset = set()
	for i, gpoly in enumerate(g_out_factors):
		if gpoly[0].degree() == t:
			Gset.add(gpoly[0])
			print(gpoly[0])
	return Gset


def get_g(G, k, prods, t, R):
	"""
		returns a set of Goppa polynomials that passed through test_g
	"""
	Gset = get_gSet(G, k, prods, t)
	good_g = set()
	for g_candidate in Gset:
		if test_g(g_candidate, G, prods, R):
			good_g.add(g_candidate)
	return good_g


def check_for_zero_on_S(vec, S):
	"""
		check if vec has 0-coordinates on indices S
	"""
	res = sum([ZZ(vec[i]) for i in range(len(vec)) if i in S])
	return not bool(res)

def positions_in_set(setA, setB):
	"""
		setB is a subset of setA
		returns positions of elements of setB in setA
	"""
	I = set()
	for i in range(len(setA)):
		if list(setA)[i] in setB:
			I.add(i)
	return I

def gen_instance(param_set,ff,F):
	"""
		generate a McEliece instance
		Returns public H, Goppa point L, and Goppa polynomial g
	"""

	n = param_set['n']
	k = param_set['k']
	t = param_set['w']

	g = PolynomialRing(ff, 'x').irreducible_element(t)

	#print('constructing L...')
	L = constructL(ff, n)
	#print('constructing H...')
	H = gen_parity_full(L, g, kappa, n, F, ff)


	return H, L, g

def recoverL(H, I, g, Labr,F,R):
	"""
	recover the complete set of Goppa points
	recieving on input public parity-check matrix H,
	Goppa polynomial g, some Goppa points Labr indexed by I,
	the finite field F and the polynomial ring R = F[x]


	The algorithm first finds an invertible submatrix in H with columns indexed by Labr
	(trying trials-many times)

	It then computes unkonwn Goppa points using Equation (5) from
	https://eprint.iacr.org/2022/525
	"""

	S.<xbar> = R.quo(g)
	n = H.ncols()
	k = n-H.nrows()
	Hpub_abr_ = H[[i for i in range(H.nrows())],[j for j in range(H.ncols()) if j in I]]
	assert(Hpub_abr_.ncols()==len(Labr))
	assert(len(I)>=H.nrows())
	trials = 200
	counter = 0
	while counter<trials:
		indicator = gen_error(H.nrows(),len(Labr))
		Hpub_abr = Hpub_abr_[[i for i in range(Hpub_abr_.nrows())],[j for j in range(len(Labr)) if indicator[j]==1]]
		assert(Hpub_abr.nrows()==Hpub_abr.ncols())
		if Hpub_abr.is_invertible():
			break
		counter+=1
	if counter==trials:
		print('could not find invertible Hpub_abr')
		return 0


	nSet = set([i for i in range(n)])
	while len(I)<n:
		exit_flag = False
		notI = nSet.difference(I)
		for istar in notI:
			if exit_flag: break

			Iaug = sorted(set(I).union({istar}))
			Hpub_abr = H[[i for i in range(n-k)],[j for j in range(n) if j in Iaug]]
			ind = list(positions_in_set(Iaug, {istar}))[0]
			zerov = vector([0]*Hpub_abr.ncols())

			for vec in Hpub_abr.transpose().kernel():

				# finding a good codeword
				if vec[ind]==1:



					c = vector(F, [0]*n)
					for i in range(len(vec)):
						if vec[i]==1:
							c[Iaug[i]] = 1
					assert(H*c == vector(F, [0]*(n-k))), 'non zero syndrome'

					syndrome = 0
					counter = 0
					for i in range(len(c)):
						if (not i == istar) and (i in I):
							syndrome+=c[i]/(x-Labr[counter])
							counter+=1

					# finding alpha unsing Eq.(5)
					smodg = S(syndrome.subs({x:xbar}))
					if not smodg==0:
						#if (1/smodg).lift().degree()>1:
							#print('smodg:', smodg)
							#print((1/smodg).lift())
						assert((1/smodg).lift().degree()==1), 'wrong degree in 1/smodg'
						element = (1/smodg).lift().roots()[0][0]
					else:
						element = 0
						#print('found 0 element in L')
					#print('found L[', istar, '] = ', element)
					Labr = Labr[:ind]+[element]+Labr[ind:]
					exit_flag = True
					I = set(Iaug)

					break
		if not exit_flag:
			print('have run through all notI, terminate')
			return 0
			break
	return Labr

def recover_Hpriv_complete(H, g, Labr, F, ff, kappa, I):
	"""
		computing a privite parity-check matrix
		by (1) creating an invertible submatrix of the public H and
		(2) generating a privite submatrix from the subset of Goppa points Labr
		and g. We then find an invertible transformation (U) between them.

		Assumes that len(Labr)==nrows+1. Uncomment the commented out lines for
		the function to work with any Labr>=nrows+1
	"""
	#assert(dim<len(Labr))
	nrows = H.nrows()
	assert(len(Labr)>nrows)
	Hpub_abr_ = H[[i for i in range(nrows)],[j for j in range(H.ncols()) if j in I]]
	#assert(Hpub_abr_.ncols()==len(Labr))

	assert(len(Labr)==nrows+1)
	for indicator in range(len(Labr)):
		#print(ind,[j for j in range(len(Labr)) if not j==ind])
		Hpub_abr = Hpub_abr_[[i for i in range(nrows)],[j for j in range(len(Labr)) if not j==indicator]]
		assert(Hpub_abr.nrows()==Hpub_abr.ncols())
		if Hpub_abr.is_invertible():
			break
	if indicator==len(Labr)-1:
		print('could not find invertible Hpub_abr')
		return 0


	"""
	trials = 500
	counter = 0
	while counter<trials:
		indicator = gen_error(nrows, len(Labr))
		#print(indicator)
		Hpub_abr = Hpub_abr_[[i for i in range(nrows)],[j for j in range(len(Labr)) if indicator[j]==1]]
		#print(Hpub_abr.nrows(), Hpub_abr.ncols(), rank(Hpub_abr))
		if Hpub_abr.is_invertible():
			break
		counter+=1
	if counter==trials:
		print('could not find invertible Hpub_abr')
		return 0
	"""

	Hpriv_abr = gen_parity_full(Labr, g, kappa, len(Labr), F, ff, echelonForm=False)
	Hpriv_abr = Hpriv_abr[[i for i in range(nrows)],[j for j in range(len(Labr)) if not j==indicator]]
	#print('Hpriv_abr is constructed')
	U = Hpriv_abr*Hpub_abr.inverse()
	#print('U is constructed')

	Hpriv_candidate = U*H
	return Hpriv_candidate

def recover_all_alpha(Lleaked, g, Hpriv_candidate, ff, kappa, target_size):
	"""
	Using a candindate for the privite parity-check matrix found in recover_Hpriv_complete,
	find all Goppa points by brute-forcing over the field ff
	"""
	V, from_V, to_V = ff.vector_space(F, map=True)
	LeakedSet = set(Lleaked)
	nr = Hpriv_candidate.nrows() / kappa
	L_remaining = set()
	column_set = Hpriv_candidate.columns()
	exit_flag = False
	for el in ff:
		if exit_flag: break
		if not el in LeakedSet:
			inv = 1/g(el)
			vector_candidate = vector([inv*el^i for i in range(nr)])
			vector_candidate_bin_ = [to_V(ff(vector_candidate[i])) for i in range(nr)]
			tmp = []
			for x in vector_candidate_bin_: tmp+=list(x)
			tmp = (vector(tmp))
			try:
				ind = column_set.index(tmp)
				L_remaining.add(el)
				if len(L_remaining)==target_size: exit_flag=True
			except ValueError:
				continue
				#print(el, 'is not in L')
				#return L_remaining

	return L_remaining


def recoverL_beyond_bound(H, I, g, Labr, R):
	"""
	same as recoverL()
	"""
	S.<xbar> = R.quo(g)
	n = H.ncols()
	k = n-H.nrows()
	# find linear dependency between some columns indexed by I and *one* vector indexed by notI
	while len(I)<n:
		#print('|I| = ', len(I))
		notI = set([i for i in range(n)]).difference(I)
		exit_flag = False
		for istar in notI:
			if exit_flag: break
			Iaug = sorted(set(I).union({istar}))
			Htmp = H[[i for i in range(n-k)],[j for j in range(n) if j in Iaug]]
			ind = list(positions_in_set(Iaug, {istar}))[0]

			if Htmp.rank()<len(I):
				# find a non-zero kernel in Htmp
				for vec in Htmp.transpose().kernel():
					# find the one with 1 on ind
					if vec[ind]==1:
						c = vector(F, [0]*n)
						counter = 0
						for i in range(len(vec)):
							if vec[i]==1:
								c[Iaug[i]] = 1
						assert(H*c == vector(F, [0]*(n-k)))

						syndrome = 0
						counter = 0
						for i in range(len(c)):
							if (not i == istar) and (i in I):
								syndrome+=c[i]/(x-Labr[counter])
								counter+=1

						smodg = S(syndrome.subs({x:xbar}))
						assert((1/smodg).lift().degree()==1)
						element = (1/smodg).lift().roots()[0][0]
						#print('found L[', istar, '] = ', element)
						Labr = Labr[:ind]+[element]+Labr[ind:]
						exit_flag = True
						I = set(Iaug)
						break
				#print('didnt find good vec with', istar, 'trying another one')
		if not exit_flag:
			print('have run through all notI, terminate')
			return 0
			break
	return Labr
