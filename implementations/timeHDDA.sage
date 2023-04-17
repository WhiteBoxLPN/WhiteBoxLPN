from time import time
import pickle #To register time samples

I=20      #Number of iterations of measuring time per window

lenK=4096 #Corresponds to |K|, the number of key guesses

for O in range(2,6):
    print("Order %d:" % O)
    W=5
    average=0.0
    fp=0 #counter of false-positives

    while (average<120) and (W<51):
        average=0.0
        listOfResults = []

        #T=binomial(W,<=O)+30
        #Adding 30 ensures that for only 2^-30 of the time we will have a false-positive
        T=30
        for o in range(1,O+1):
            T+=binomial(W,o)
        VS=VectorSpace(GF(2), T)

        #We stop measuring time if we exceed 2min of computation
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
            for keyByteGuess in range(lenK):
                listOfSelVectors.append(VS.random_element())

            extendedWindowcp = extendedWindow.copy()

            #Actual HDDA
            t1 = time()
            #Expanding the window with AND combinations
            for combSize in range(2,O+1):
                for comb in Combinations(extendedWindowcp[0:W], combSize):
                    ANDedCombVect=[0]*T
                    for vect in comb:
                        for traceNumber in range(T):
                            ANDedCombVect[traceNumber]&=vect[traceNumber]
                    extendedWindowcp.append(ANDedCombVect)

            #Solving the linear system by computing the parity check matrix, and
            #   verify whether a selection vector S verifies PCM*S==0 or not. If
            #   so, the key guess corresponding to S might be a correct key guess
            PCM=Matrix(GF(2), extendedWindowcp).right_kernel().basis_matrix()
            nrows=PCM.nrows()
            PCM=PCM.rows()
            for vectPos in range(lenK):
                correct=True
                row = 0
                while correct and (row<nrows):
                    col = 0
                    while col < T:
                        if PCM[row][col]*listOfSelVectors[vectPos][col]:
                            correct=False
                        col+=1
                    row+=1
                if correct:
                    #We should register the key byte guess corresponding to the
                    #   selection vector. But it does not have an impact on time
                    #   measurng. However, we are sure that we found a false-
                    #   positive.
                    fp+=1

            t2 = time()
            listOfResults.append(t2-t1)
            average+=t2-t1

        with open("../timeResults/HDDA/HDDA_Order%d/HDDA_O%d_W%d.pkl" % (O, O, W), "wb") as f:
            pickle.dump(listOfResults, f)

        #Displaying execution progress, since the whole execution takes hours
        average=average/I
        print("W=%d: " % W, end='')
        print(average)
        W+=1
    print("\nEncountered false-positives: %d\n" % fp)
