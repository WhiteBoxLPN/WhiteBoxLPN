import os
import io
import pickle #For reading results

donothing=0
#time measurments of WBLPN
print("WBLPN:\n")
for tauPower in range(2,6):
    print("\ntau=1/2^%d:\n" % tauPower)
    print("Average:")
    for W in range(5, 51):
        try:
            with open("../ourTimeResults/WBLPN/WBLPN_tau2pow%d/WBLPN_I20_tau2pow%d_W%d.pkl" % (tauPower, tauPower, W), "rb") as f:
                ListOfTimes = pickle.load(f)
            print(sum(ListOfTimes)/20)
        except IOError:
            donothing+=1
    print("\nMin:")
    for W in range(5, 51):
        try:
            with open("../ourTimeResults/WBLPN/WBLPN_tau2pow%d/WBLPN_I20_tau2pow%d_W%d.pkl" % (tauPower, tauPower, W), "rb") as f:
                ListOfTimes = pickle.load(f)
            print((sum(ListOfTimes)/20)-min(ListOfTimes))
        except IOError:
            donothing+=1
    print("\nMax:")
    for W in range(5, 51):
        try:
            with open("../ourTimeResults/WBLPN/WBLPN_tau2pow%d/WBLPN_I20_tau2pow%d_W%d.pkl" % (tauPower, tauPower, W), "rb") as f:
                ListOfTimes = pickle.load(f)
            print(max(ListOfTimes)-(sum(ListOfTimes)/20))
        except IOError:
            donothing+=1

#time measurments of HODCA
print("\n\nHODCA:\n")
for Order in range(2,6):
    print("\nOrder=%d:\n" % Order)
    print("Average:")
    for W in range(5, 51):
        try:
            with open("../ourTimeResults/HODCA/HODCA_Order%d/HODCA_I20_O%d_W%d.pkl" % (Order, Order, W), "rb") as f:
                ListOfTimes = pickle.load(f)
            print(sum(ListOfTimes)/20)
        except IOError:
            donothing+=1
    print("\nMin:")
    for W in range(5, 51):
        try:
            with open("../ourTimeResults/HODCA/HODCA_Order%d/HODCA_I20_O%d_W%d.pkl" % (Order, Order, W), "rb") as f:
                ListOfTimes = pickle.load(f)
            print((sum(ListOfTimes)/20)-min(ListOfTimes))
        except IOError:
            donothing+=1
    print("\nMax:")
    for W in range(5, 51):
        try:
            with open("../ourTimeResults/HODCA/HODCA_Order%d/HODCA_I20_O%d_W%d.pkl" % (Order, Order, W), "rb") as f:
                ListOfTimes = pickle.load(f)
            print(max(ListOfTimes)-(sum(ListOfTimes)/20))
        except IOError:
            donothing+=1

#time measurments of HDDA
print("\n\nHDDA:\n")
for Order in range(2,6):
    print("\nDegree=%d:\n" % Order)
    print("Average:")
    for W in range(5, 51):
        try:
            with open("../ourTimeResults/HDDA/HDDA_Order%d/HDDA_O%d_W%d.pkl" % (Order, Order, W), "rb") as f:
                ListOfTimes = pickle.load(f)
            print(sum(ListOfTimes)/20)
        except IOError:
            donothing+=1
    print("\nMin:")
    for W in range(5, 51):
        try:
            with open("../ourTimeResults/HDDA/HDDA_Order%d/HDDA_O%d_W%d.pkl" % (Order, Order, W), "rb") as f:
                ListOfTimes = pickle.load(f)
            print((sum(ListOfTimes)/20)-min(ListOfTimes))
        except IOError:
            donothing+=1
    print("\nMax:")
    for W in range(5, 51):
        try:
            with open("../ourTimeResults/HDDA/HDDA_Order%d/HDDA_O%d_W%d.pkl" % (Order, Order, W), "rb") as f:
                ListOfTimes = pickle.load(f)
            print(max(ListOfTimes)-(sum(ListOfTimes)/20))
        except IOError:
            donothing+=1
