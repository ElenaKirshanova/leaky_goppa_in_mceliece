import time
import copy
import os.path
from os.path import isfile, join
import sys

load("utils.sage")

McEliece0 = {'label': 0, 'n': 41,  'k': 11, 'w': 5, 'm': 6, 'f': x^6 + x + 1}
McEliece01 = {'label': 1, 'n': 488, 'k': 120, 'w': 16, 'm': 12, 'f':x^12+x^3+1}
McEliece1 = {'label': 1, 'n': 3488, 'k': 2720, 'w': 64, 'm': 12, 'f':x^12+x^3+1}
McEliece2 = {'label': 2,'n': 4608, 'k': 3360, 'w': 96, 'm': 13, 'f':x^13+x^4+x^3+x+1}
McEliece3 = {'label': 3,'n': 6960, 'k': 5413, 'w': 119,'m': 13, 'f':x^13+x^4+x^3+x+1}
McEliece4 = {'label': 4,'n': 8192, 'k': 6528, 'w': 128,'m': 13, 'f':x^13+x^4+x^3+x+1}

param_set = McEliece1


def generate_system(f, u, ff, Rg, Rg_squared, t):
    """
    generate the system of the form [f | x*f | ... | x^(u-1)*f| -Id]
    """
    A = matrix(ff, 2*t, 2*u-1)
    for j in range(u):
        tmp_ = list(Rg_squared(f*x^j))
        for i in range(len(tmp_)):
            A[i,j] = tmp_[i]
    for i in range(u, 2*u-1): A[i-u, i] = -1
    return A

def generate_c2(H,I,t,m,F,c1):
    """
    generate  a codeword that has only one '1' coordintes on positions supp(c_1) \setminus I
    and has at most t ones on positions [1..n]\setminus I
    """
    #print(H)
    n = H.ncols()
    Ntrials = 700
    trial = 0
    full_rank_failures = 0

    supp_c1 = []
    for i in range(n):
        if (not i in I) and (c1[i]==1): supp_c1.append(i)
    Habr = matrix(F, H.nrows(), t*m+2) #+const to hope for full row-rank
    counter = 0
    for j in range(n):
        if (j in I):
            for i in range(H.nrows()):
                Habr[i,counter] = H[i,j]
            counter+=1
    save_counter = counter

    while trial < Ntrials:
        new_columns = []
        counter = save_counter

        while counter<t*m+2: #+the same const to hope for full row-rank
            new_column = ZZ.random_element(0, n)
            if (not new_column in I) and (not new_column in new_columns) and (not new_column in supp_c1):
                for i in range(H.nrows()):
                    Habr[i,counter] = H[i,new_column]
                counter+=1
                new_columns.append(new_column)
        if Habr.rank()<t*m:
            trial+=1
            full_rank_failures+=1
            continue

        #for each i s.t. c1[i]==1, try to express c1[i] as a lin. combination of Habr and this lin. combination has <t+2 1's on notI
        # such linear combination will form c2
        for i in supp_c1:

            assert(not i in new_columns)

            #ith column
            columnH = [0]*(t*m)
            for j in range(t*m): columnH[j] = H[j, i]
            #print(Habr.ncols(), Habr.nrows(), t*m, Habr.rank())

            res = Habr.solve_right(columnH)
            assert(len(res)==len(I)+len(new_columns))

            ctr = 0
            codeword = [0]*n
            codeword[i] = 1
            wt_res = 1
            for j in I:  # cols in Habr are positioned differently to cols in H, a loop over range(n) is wrong
                codeword[j] = res[ctr]
                ctr+=1
            for j in new_columns:
                codeword[j] = res[ctr]
                ctr+=1
                if codeword[j]==1: wt_res +=1

            assert(H*vector(F,codeword) == vector(F,[0]*H.nrows()))

            if wt_res < t+1:
                return codeword, i
        trial+=1
    print('Warning: ran out of trials in generate_c2')
    print('number of full_rank_failures:', full_rank_failures)
    return [],[]


def generate_codeword(H,I,t,m,F):
    """
    generate a codeword that has at most t+1 1's on [1..n]\setminus I positions
    """
    n = H.ncols()

    Ntrials = 300
    trial = 0
    while trial < Ntrials:

        notI = []#copy.deepcopy(i_star)
        counter = len(notI)


        while counter<(t*m+1)-len(I):
            tmp_pos = ZZ.random_element(0, n)
            if (not tmp_pos in I) and (not tmp_pos in notI):
                notI.append(tmp_pos)
                counter += 1

        Habr = matrix(F, H.nrows(), t*m+1)

        counter = 0
        for j in range(n):
            if (j in I) or (j in notI):
                for i in range(H.nrows()):
                    Habr[i,counter] = H[i,j]
                counter+=1

        Habr_kernel = Habr.transpose().kernel().basis()

        for v in Habr_kernel:

            counter = 0
            codeword = [0]*n
            for i in range(n):
                if (i in I) or (i in notI):
                    codeword[i] = v[counter]
                    counter+=1

            #if (not i_star==[]) and (codeword[i_star[0]]==0): continue

            n_ones = sum([int(codeword[i]) for i in notI ])
            #if (not i_star==[]):
            #    assert(codeword[i_star[0]]==1)

            assert(H*vector(F,codeword) == vector(F,[0]*H.nrows()))


            if n_ones < t+2:  #condition assures that we will be able to find the unknown alphas
                return codeword

        trial+=1
    print('Warning: ran out of trials in generate_codeword()')
    return []



def obtain_alphas(codeword, I, notI, Lleaked, ff, Rg, Rg_squared):
    """
        obtain the *set* A_codewords  -- alpha_i's for i \in notI
    """

    known_part = 0
    counter = 0
    for i in I:
        if codeword[i]==1:
            known_part+=1/Rg_squared(x-Lleaked[counter])
        counter+=1

    A = generate_system(known_part, len(notI), ff, Rg, Rg_squared, t)
    target = vector(ff, list(Rg_squared(known_part*x^(len(notI)) +(len(notI)%2)*x^(len(notI)-1))) )

    res = A.solve_right(target)
    res_ = res[:len(notI)]



    #root_fining_start = time.time()
    res_roots = R(list(res_)+[1]).roots(ff)
    #root_fining_end = time.time()

    alphas = [res_roots[i][0] for i in range(len(res_roots))]

    return alphas


p = 2
F = GF(p)
kappa = param_set['m']
defining_poly = param_set['f']
print('defining_poly:', defining_poly, 'param_set:', param_set['label'])
ff.<a> = FiniteField(p**kappa, modulus = defining_poly)
R = PolynomialRing(ff, 'x')
x = R.gen()


H, L, g = gen_instance(param_set,ff,F)
Rg = R.quotient(g, 'x')
Rg_squared = R.quotient(g*g, 'x')

#print('g:', g, 'g^2:', g*g)

n = param_set['n']
t = g.degree()
m = param_set['m']
sizeI =(t*m)+1-2*t
print('sizeI:', sizeI)
print('t:', t)

I = genI(sizeI, n)
I = list(I)
I.sort() # otherwise known_part is computed wrongly


Lleaked =  genLabr(L, I)

def obtain_common_alpha(c1, notI1, c2, notI2, I, Lleaked, ff, Rg, Rg_squared):
    """
        A_c1 \cap A_c2
    """

    alpha1 = obtain_alphas(c1, I, notI1, Lleaked, ff, Rg, Rg_squared)
    alpha2 = obtain_alphas(c2, I, notI2, Lleaked, ff, Rg, Rg_squared)


    if alpha1==[] or alpha2==[]:
        return -1


    intersection = list(set(alpha1) & set(alpha2))
    return  intersection

def notI(I, c):
    """
        return supp(c_1)\setminus I
    """
    notI = []
    for i in range(len(c)):
        if (not i in I) and c[i] == 1: notI.append(i)

    return notI

def update_Lleaked(I, i_star, Lleaked, new_alpha):
    """
    insert i_star into Lleaked into the correct positions
    """
    i = 0
    while i_star>I[i] and i<len(I):
        i+=1
    n = len(Lleaked)
    pos_insert = i

    Lleaked_tmp = copy.deepcopy(Lleaked)

    Lleaked.append(0)
    Lleaked[pos_insert] = new_alpha

    for i in range(pos_insert, n):
        Lleaked[i+1] = Lleaked_tmp[i]
    return Lleaked


start = time.time()
new_alphas = []
new_positions = []
sizeI_ = sizeI
time_c1 = 0
time_c2 = 0
time_newalpha = 0
time_per_alpha = 0
ctr_system_failure = 0
while sizeI_ < t*m+1:
    start_c1 = time.time()
    c1 = generate_codeword(H,I,t,m,F)
    finish_c1 = time.time()
    time_c1 +=finish_c1 - start_c1


    notI1 = notI(I, c1)
    start_c2 = time.time()
    c2, i_star = generate_c2(H,I,t,m,F,c1)
    finish_c2 = time.time()
    time_c2 +=finish_c2 - start_c2


    if i_star in new_positions:
        continue
    notI2 = notI(I, c2)

    start_new_alpha = time.time()
    new_alpha = obtain_common_alpha(c1, notI1, c2, notI2, I, Lleaked, ff, Rg, Rg_squared)
    if new_alpha == -1:
        ctr_system_failure+=1
        continue
    finish_new_alpha = time.time()
    time_newalpha+=finish_new_alpha-start_new_alpha

    if (not len(new_alpha)==1) or (len(new_alpha)==0):
        continue
    assert(L[i_star]==new_alpha[0])
    new_alphas.append(new_alpha[0])
    new_positions.append(i_star)
    sizeI_ += 1

    if sizeI_ < t*m - 10:
        Lleaked = update_Lleaked(I, i_star, Lleaked, new_alpha[0])
        I.append(i_star)
        I.sort()
    print(sizeI_, L[i_star], i_star, new_alpha)
    finish_one_alpha = time.time()
    time_per_alpha+= finish_one_alpha - start_c1
end = time.time()
print('obtained ', sizeI_-sizeI, ' new alphas in time', end-start)
print('average c1:', time_c1/(sizeI_-sizeI), 'average c2:', time_c2/(sizeI_-sizeI), 'time_newalpha:',time_newalpha/(sizeI_-sizeI), 'time per alpha:', time_per_alpha/(sizeI_-sizeI) )
print('ctr_system_failure:', ctr_system_failure)

"""
# checking the implementation above by multiplying A with the known solution
c1, _ = generate_codeword(H,I,t,m,F)
notI1 = notI(I, c1)
print('len(notI1):', len(notI1))
print([L[i] for i in notI1])
alpha1 = obtain_alphas(c1, I, notI1, Lleaked, ff, Rg, Rg_squared, L)
print('alpha1:', alpha1)
#print(sum([c1[i]/Rg_squared(x-L[i]) for i in range(n)]))
known = sum([c1[i]/Rg_squared(x-L[i]) for i in I])
print('known:', sum([c1[i]/Rg_squared(x-L[i]) for i in I]))
print('unknown:', sum([c1[i]/Rg_squared(x-L[i]) for i in notI1]))
print('lhs:', prod([x-L[i] for i in notI1])*known)
rhs = 0
for i in notI1:
    tmp = 1
    for j in notI1:
        if not j==i: tmp*=x-L[j]
    rhs+=tmp
print('rhs:', rhs)

p_poly = 1
for i in notI1:
    if c1[i]==1 :p_poly*=R(x - L[i])
print('p_poly:', p_poly)

p_prime = 0
for i in notI1:
    tmp = 1
    for j in notI1:
        if (not j==i) and (c1[i]==1):
            tmp *= R(x-L[j])
    p_prime +=R(tmp)
print('p_prime:', p_prime)
#print('p_poly*f:', Rg(p_poly*known))

solution = list(p_poly[:p_poly.degree()])
p_prime_l = list(p_prime)
if (not p_prime == 1):
    solution+=[p_prime_l[0]]
    for i in range(1, p_prime.degree()+1):
        if i<p_poly.degree()-1:
            solution+=[p_prime_l[i]]
    #if p_prime.degree()%2==0: solution+=[0]
    solution+=[0]*(p_poly.degree()-p_prime.degree()-2)

print('solution:', solution)
A = generate_system(known, len(notI1), ff, Rg, Rg_squared, t)
print(A.nrows(), A.ncols(), len(solution))
target2 = A*vector(ff,solution)
print(target2)
print(A.solve_right(target2))
print(A.right_kernel())
print(A.right_kernel().dimension())
"""
