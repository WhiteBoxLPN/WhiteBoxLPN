from time import time
import pickle #To register time samples

I=20      #Number of iterations of measuring time per window

lenK=4096 #Corresponds to |K|, the number of key guesses
T=100     #Number of traces, the more traces, the more likely the best key candidate is
          #is the correct one. T impact linearly the time complexity of HODCA

def correlation(V1, V2):
    l=len(V1)
    matches=0
    mismatches=0
    for i in range(l):
        if V1[i] == V2[i]:
            matches+=1
        else:
            mismatches+=1
    return((matches-mismatches)/l)

for O in range(2,6): #O represents the Order of HODCA
    print("Order=%d" % O)
    W=5
    average=0.0
    #We stop measuring time if we exceed 2min of computation
    while (average<120) and (W<51):
        average=0.0
        listOfResults = []
        for _ in range(I):
            #Generation of a random Window:
            extendedWindow=[]
            for node in range(W):
                vectTemp=[]
                for traceNumber in range(T):
                    vectTemp.append(randrange(2))
                extendedWindow.append(vectTemp)

            #Generation of the random selection vectors:
            listOfSelVectors=[]
            for keyByteGuess in range(256):
                vectPos=[]
                for bytePosition in range(16):
                    vectTrace=[]
                    for traceNumber in range(T):
                        vectTrace.append(randrange(2))
                    vectPos.append(vectTrace)
                listOfSelVectors.append(vectPos)

            #mostProbableKey is constituted of 16 pair of values: the first element of the
            #   pair represents the key guess in {0,...,255}, the second its corresponding
            #   best correlation found so far. In practice it can be interesting to keep
            #   the top 10 best candidates.
            mostProbableKey = [[0, 0.4] for __ in range(16)]

            extendedWindowcp = extendedWindow.copy()
            #Actual HODCA
            t1 = time()
            #Expanding the window with XOR combinations
            for combSize in range(2,O+1):
                for comb in Combinations(extendedWindowcp[0:W], combSize):
                    XORedCombVect=[0]*T
                    for vect in comb:
                        for traceNumber in range(T):
                            XORedCombVect[traceNumber]^=vect[traceNumber]
                    extendedWindowcp.append(XORedCombVect)

            #Computing correlation of each vector of the extended window with each of the selection vectors
            for vect in extendedWindowcp:
                for keyByteGuess in range(256):
                    for bytePosition in range(16):
                        c=abs(correlation(listOfSelVectors[keyByteGuess][bytePosition], vect))
                        if c > mostProbableKey[bytePosition][1]:
                            mostProbableKey[bytePosition]=[keyByteGuess, c]
            t2 = time()

            listOfResults.append(t2-t1)
            average+=t2-t1

        with open("../timeResults/HODCA/HODCA_Order%d/HODCA_I20_O%d_W%d.pkl" % (O, O, W), "wb") as f:
            pickle.dump(listOfResults, f)

        #Displaying execution progress, since the whole execution takes hours
        average=average/I
        print("W=%d: " % W, end='')
        print(average)
        W+=1
    print()
