from time import time

pnfm=.999        #Probability of finding Noise-Free Matrix after A attempts
lenK=4096        #Corresponds to |K|, the number of key guesses
Tp=100           #Number of traces for sampling

W=20             #Window size
O=2              #Order
tauOrd=1/(2**7)  #Noise rate of order O

I=30             #Number of iterations for trying finding false positives / missing solutions

RF=RealField(10000) #To ensure correct probabilistic parameters computation for values very close to zero

#Preparing probabilistic parameters for every window sizes < binomial(W,<=Order)
Wextended=0
for o in range(1,O+1):
    Wextended+=binomial(W,o)
mARRAY=[]
cARRAY=[]
AARRAY=[]
for w in range(1,Wextended+1):
    alpha = RF(1)/(2**w)
    beta = (1/lenK)*(((RF(1)-tauOrd)/2)**w)

    m = ((sqrt((3/2)*log(1/alpha))+sqrt(log(1/beta)))/((1/2)-tauOrd))**2
    mARRAY.append(ceil(m))

    c = tauOrd*m+sqrt(3*((1/2)-tauOrd)*log(1/alpha)*m)
    cARRAY.append(round(c))

    A=ceil(log(RF(1)-pnfm)/(log(1-(RF(1)-tauOrd)^w)))
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

def HOWBLPN(M, Spool, Sv, mostProbableKey):
    T=Tp+mARRAY[Wextended-1]
    #Ensuring linearly independent columns, for potential invertible sampled matrices
    linIndColumns = M.pivots()

    #Extending the array
    MextendedLIST=M[range(T), linIndColumns].columns()
    for combSize in range(2,O+1):
        for comb in Combinations(MextendedLIST[0:len(linIndColumns)], combSize):
            ANDedCombVect=[1]*T
            for vect in comb:
                for traceNumber in range(T):
                    ANDedCombVect[traceNumber]*=vect[traceNumber]
            MextendedLIST.append(ANDedCombVect)
    Mextended=MatrixSpace(GF(2), Wextended, T)(MextendedLIST).transpose()

    #Preparing WBLPN inputs
    Mpool=Mextended[range(Tp), range(Wextended)]
    Mv=Mextended[range(Tp,Tp+mARRAY[Wextended-1]), range(Wextended)]

    return(WBLPN(Mpool, Spool, Mv, Sv, mostProbableKey))

#We show below that the for the example given in the paper (P=x1+x2+x3+x4+x5x6+x7x8+x9x10x11x12x13x14x15),
#HOWBLPN of order 2 takes about 3.13s to be performed (mistake in the paper: we said 20s due to a previous error)
print("Measuring average time of HOWBLPN of order %d, %d-order noise rate=%.5f, window size %d for %d iterations\n" % (O,O,tauOrd,W,I))
average=0.0
for i in range(I):
    #Generating random traces
    Spool=MatrixSpace(GF(2), Tp, lenK).random_element()
    m=mARRAY[Wextended-1]
    M=MatrixSpace(GF(2), m+Tp, W).random_element()
    Sv=MatrixSpace(GF(2), m, lenK).random_element()

    mostProbableKey=[]
    t1 = time()
    mostProbableKey=HOWBLPN(M, Spool, Sv, mostProbableKey)
    t2 = time()

    average+=(t2-t1)
    print('\033[1A', end='\x1b[2K')
    print("%.2f%% done" % ((i+1)*100/I))

print('\033[1A', end='\x1b[2K')
print("Average for %d iterations: %.5fs" % (I,average/I))
