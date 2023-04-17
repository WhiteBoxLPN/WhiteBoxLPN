from time import time
import pickle #To register time samples

I=20      #Number of iterations of measuring time per window

pnfm=.999 #Probability of finding Noise-Free Matrix after A attempts
lenK=4096 #Corresponds to |K|, the number of key guesses
Tp=100    #Number of traces for sampling

RF=RealField(10000) #To ensure correct probabilistic parameters computation

def WBLPN(Mpool, Spool, Mv, Sv, mostProbableKey):
    #Ensuring linearly independent columns, for potential invertible sampled matrices
    linIndColumns = Mpool.pivots()

    #Adapting probabilistic parameters to the potentially reduced LPN dimension
    Wprime=len(linIndColumns)
    m=mARRAY[Wprime-2]
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


#Measuring time for tau in {1/4, 1/8, 1/16, 1/32}
for tauPower in range(2,6):

    #Preparation of the parameters that are computed before performing the sliding window method:
    print("tau=1/%d" % (2**tauPower))
    tau=RF(1)/(2**tauPower)
    mARRAY=[]
    cARRAY=[]
    AARRAY=[]
    for w in range(1,51):
        alpha = RF(1)/(2**w)
        beta = (1/lenK)*(((RF(1)-tau)/2)**w)

        m = ((sqrt((3/2)*log(1/alpha))+sqrt(log(1/beta)))/((1/2)-tau))**2
        mARRAY.append(ceil(m))

        c = tau*m+sqrt(3*((1/2)-tau)*log(1/alpha)*m)
        cARRAY.append(round(c))

        A=ceil(log(RF(1)-pnfm)/(log(1-(RF(1)-tau)^w)))
        AARRAY.append(A)

    W=5
    fp=0 #counter of false-positives
    average=0.0
    #We stop measuring time if we exceed 2min of computation
    while (average<120) and (W<51):
        average=0.0
        listOfResults=[]

        #sampling 20 times:
        for _ in range(I):

            #Generation of random inputs:
            Mpool=MatrixSpace(GF(2), Tp, W).random_element()
            Spool=MatrixSpace(GF(2), Tp, lenK).random_element()

            m=mARRAY[W-1]
            Mv=MatrixSpace(GF(2), m, W).random_element()
            Sv=MatrixSpace(GF(2), m, lenK).random_element()

            #Measuring time of PooledGauss:
            mostProbableKey=[]
            t1 = time()
            mostProbableKey=WBLPN(Mpool, Spool, Mv, Sv, mostProbableKey)
            t2 = time()

            listOfResults.append(t2-t1)

            average+=(t2-t1)
            fp+=len(mostProbableKey)

        #Saving our results
        with open("../timeResults/WBLPN/WBLPN_tau2pow%d/WBLPN_I%d_tau2pow%d_W%d.pkl" % (tauPower, 20, tauPower, W), "wb") as f:
            pickle.dump(listOfResults, f)

        #Displaying execution progress, since the whole execution takes hours
        average=average/I
        print("W=%d: " % W, end='')
        print(average)
        W+=1
    print("\nEncountered false-positives: %d\n" % fp)
