pnfm=.999 #Probability of finding Noise-Free Matrix after A attempts
lenK=4096 #Corresponds to |K|, the number of key guesses
Tp=100    #Number of traces for sampling

W=20      #Window size
tau=1/16  #Noise rate

I=30      #Number of iterations for trying finding false positives / missing solutions

RF=RealField(10000) #To ensure correct probabilistic parameters computation for values very close to zero

#Preparing probabilistic parameters for every window sizes < W
mARRAY=[]
cARRAY=[]
AARRAY=[]
for w in range(1,W+1):
    alpha = RF(1)/(2**w)
    beta = (1/lenK)*(((RF(1)-tau)/2)**w)

    m = ((sqrt((3/2)*log(1/alpha))+sqrt(log(1/beta)))/((1/2)-tau))**2
    mARRAY.append(ceil(m))

    c = tau*m+sqrt(3*((1/2)-tau)*log(1/alpha)*m)
    cARRAY.append(round(c))

    A=ceil(log(RF(1)-pnfm)/(log(1-(RF(1)-tau)^w)))
    AARRAY.append(A)

def WBLPN(Mpool, Spool, Mv, Sv, mostProbableKey):
    #Ensuring linearly independent columns, for potential invertible sampled matrices
    linIndColumns = Mpool.pivots()

    #Adapting probabilistic parameters to the potentially reduced LPN dimension
    Wprime=len(linIndColumns)
    m=mARRAY[Wprime-1]
    c=cARRAY[Wprime-1]
    A=AARRAY[Wprime-1]

    numberOfAttempts=0
    while numberOfAttempts < A:
        #Sampling Ms from the pool
        sampledRows=Combinations(range(Tp), Wprime).random_element()
        Ms=Mpool[sampledRows,linIndColumns]

        if Ms.is_invertible():
            MsInv=Ms.inverse()
            numberOfAttempts+=1

            #We stranspose E since Sage computes hamming_weight faster for rows than columns
            E=(((Mv[range(m),linIndColumns]*MsInv)*Spool[sampledRows,range(lenK)])+Sv[range(m),range(lenK)]).transpose()

            for i in range(lenK):
                if E.row(i).hamming_weight() < c:
                    if i not in mostProbableKey:
                        mostProbableKey.append(i)
    return(mostProbableKey)

print("Testing the ability of WBLPN of not returning false-positives:\n")
falsePositives=0
for i in range(I):
    #Generating random traces
    Mpool=MatrixSpace(GF(2), Tp, W).random_element()
    Spool=MatrixSpace(GF(2), Tp, lenK).random_element()
    m=mARRAY[W-1]
    Mv=MatrixSpace(GF(2), m, W).random_element()
    Sv=MatrixSpace(GF(2), m, lenK).random_element()

    mostProbableKey=[]
    mostProbableKey=WBLPN(Mpool, Spool, Mv, Sv, mostProbableKey)
    if len(mostProbableKey)>0:
        falsePositives+=1
    print('\033[1A', end='\x1b[2K')
    print("%.2f%% done" % ((i+1)*100/I))
print('\033[1A', end='\x1b[2K')
print("Total false positives for %d attempts: %d\n" % (I, falsePositives))


print("Testing the ability of WBLPN of not missing a solution contained in a window:\n")
missedSolutions=0
for i in range(I):
    #Generating random traces
    Mpool=MatrixSpace(GF(2), Tp, W).random_element()
    Spool=MatrixSpace(GF(2), Tp, lenK).random_element()
    m=mARRAY[W-1]
    Mv=MatrixSpace(GF(2), m, W).random_element()
    Sv=MatrixSpace(GF(2), m, lenK).random_element()

    #Creating noisy solution: M.column(0)+M.column(1)+errorvector=S.column(0)
    #Pooled Gauss should be able to return the key guess "2" corresponding to S.column(2) each time
    noisyBits=0
    temp=[]
    for j in range(Tp):
        noise=uniform(0,1)
        if noise < tau:
            temp.append(Mpool[j][0]+Mpool[j][1]+1)
            noisyBits+=1
        else:
            temp.append(Mpool[j][0]+Mpool[j][1])
    Spool[:,2]=vector(GF(2),temp)
    temp=[]
    for j in range(m):
        noise=uniform(0,1)
        if noise < tau:
            temp.append(Mv[j][0]+Mv[j][1]+1)
        else:
            temp.append(Mv[j][0]+Mv[j][1])
    Sv[:,2]=vector(GF(2),temp)

    mostProbableKey=[]
    mostProbableKey=WBLPN(Mpool, Spool, Mv, Sv, mostProbableKey)

    if (len(mostProbableKey)==0) or (2 not in mostProbableKey):
        missedSolutions+=1
        print('\033[1A', end='\x1b[2K')
        print("Missed solution for tau = %.2f%%, and for in-practice tau in this generated Mpool = %.2f%%" % (tau*100, noisyBits*100/Tp))
        if (noisyBits>tau*Tp):
            print("To avoid having more in-practice noise in Mpool than tau*Tp, augment Tp value\n")
    print('\033[1A', end='\x1b[2K')
    print("%.2f%% done" % ((i+1)*100/I))
print('\033[1A', end='\x1b[2K')
print("Total missed solutions for %d attempts: %d" % (I, missedSolutions))
