import pandas as pd
import numpy as np
import os
from scipy.optimize import minimize_scalar

count = 1

#Error Handling, handles Subloading_tij.exe outputs that are NAN
def format_exp(x):
    if isinstance(x, float):
        return x
    if x == "-Infinity" or "Infinity":
        return 10**300.5
    tmpp = x.rfind('+', 1)
    tmpm = x.rfind('-', 1)
    if tmpp == -1:
        x = x[:tmpm-1] + "e" + x[tmpm:]
    if tmpm == -1:
        x = x[:tmpp-1] + "e" + x[tmpp:]
    return x

def tijmodelfn(e0, itr):

#Input of final stress state for layered condition in which soil layers are placed and compacted in succesion
    stress_states = np.array([[7.332556785, 0.332556785],
                     [7.997670356, 0.997670356],
                     [8.662783927, 1.662783927],
                     [9.327897497, 2.327897],
                     [9.993011068, 2.993011068],
                     [10.65812464, 3.658124638],
                     [11.32323821, 4.323238209],
                     [11.98835178, 4.98835178],
                     [12.65346535, 5.65346535],
                     [13.31857892, 6.318578921]])

    inputstress = np.array([stress_states[itr,0], stress_states[itr,1]])

# Parameters from optimization
    N = 0.3576
    e_max = 0.551
    pl = 1
    a = 103.32
    aAF = a
    aIC = a
    Qw0 = 0
    bw = 0
    Pa = 14.7
    powerIC = 2
    wgt_IC = 1
    pw = 1
    bxk = 0
    Kappa = 0.0016
    w_hpSoft = 1
    Poisson = 0.2
    Lamda = 0.0207
    Rcs = 5.085
    Beta = 1.052
    axk = 1

#Parameter Array
    tparams = [Lamda, Kappa, Rcs, N, Poisson, Beta, e_max, e0.item(), aAF, aIC, pl, axk, Qw0, bw, pw, bxk, Pa, powerIC, wgt_IC, w_hpSoft]
    tparams2 = np.array([0.003, 1.e-7, 0.0, 1.e-7, 0.0, 0.0, 1.e10, 1.0, 0.0])
    tparams3 = np.array([2.0, 5.0, 2.0])

#Running model for given iteration, model parameters, and stress state, Input for parameter text file
    with open('Para.txt', 'w') as para:
        arglist = ["Model       1", "iVersion    3", "iCyclic&iCreep      0", "iAssoc      3", "Elastic     N"]
        desc = ["[1: Subloading tij,    2: Cam clay]", "[0: original,    3: alternative, 4: alternative4]",
                "[0: no-creep,        1: creep,       10: cyclic]",
                "[0: (AF+IC) for original, 3:  (AF+IC) for alternative,   1: (AF) for both iVersions]",
                "[H: non-linear HooKean,    N: no-tension elastic]"]
        tempdf = pd.DataFrame(desc, index=arglist, columns=[''])
        tempdfAsString = tempdf.to_string(header=False, index=True)
        para.write(tempdfAsString)
        para.write(
            "\nLamda\tKappa\tRcs	N\tPoisson\tBeta\te_max\te0\taAF*\taIC*\tpl\taxk\tQw0\tbw*\tpw\tbxk\tPa\tpowerIC\twgt_IC\tw_hpSoft\n")
        for param in tparams:
            para.write("%s" % param)
            para.write("\t")
        para.write(
            "\n\nLamda_alfa\tVoiddot_ref\tVoiddot_min\tVoiddot_ini\tbR_faster\tbR_slower\tVoiddot_max\tw_Vcrp\tw_Psi\n")
        for param in tparams2:
            para.write("%s" % param)
            para.write("\t")
        para.write("\n\na_cyc\tkm_cyc\tkn_cyc\n")
        for param in tparams3:
            para.write("%s" % param)
            para.write("\t")

    para.close()

#Input for stress path text file
    with open('Path.txt', 'w') as path:
        path.write("NoPath= 2\n")
        path.write("InitialStresses= 0.332556785,\t0.332556785,\t0.332556785,\t0.0,\t0.0,\t0.0\n")
        path.write("Consolidation\tOneDimensional\t%s\t1000\t10\n" % inputstress[0])
        path.write("Consolidation\tOneDimensional\t%s\t1000\t10\n" % inputstress[1])

    path.close()

#Run Subloadingtij.exe from given file path (same path as Path and Para files)
    os.system(r'"C:\Users\Miles Drazkowski\Documents\Miles Drazkowski\Thesis\Ko Calculation\SubloadingTij_a.exe"')

#Pull resulting void ratio from results
    colspecs = [(231, 241)]
    names = ["e"]
    temp = pd.read_fwf('Ana.cal', encoding='shift_jis', colspecs=colspecs, names=names, skiprows=10)

#Error Handling
    temp.fillna(value=10 ** 101.5, inplace=True)
    if temp['e'].dtypes == 'object':
        print(temp['e'])
        try:
            temp['e'] = temp['e'].astype('float64')
        except:
            temp['e'] = temp['e'].map(format_exp)
        else:
            temp['e'] = temp['e'].astype('float64')

#Writing results to array
    SS = temp[['e']].to_numpy(dtype='float64')
    SStemp = SS[200]
#Determine if e0 (0.245) is reached to determine where to calculate K0
    if round(SStemp[0],4) == 0.245:
        #Reading horizontal and vertical stresses from calculated values
        colspecs = [(12, 25), (25, 38)]
        names = ["sx", "sy"]
        temp = pd.read_fwf('Ana.cal', encoding='shift_jis', colspecs=colspecs, names=names, skiprows=10)
        #Error handling
        if temp['sx'].dtypes == 'object':
            print(temp['sx'])
            try:
                temp['sx'] = temp['e'].astype('float64')
            except:
                temp['sx'] = temp['e'].map(format_exp)
            else:
                temp['sx'] = temp['e'].astype('float64')
        if temp['sy'].dtypes == 'object':
            print(temp['sy'])
            try:
                temp['sy'] = temp['e'].astype('float64')
            except:
                temp['sy'] = temp['e'].map(format_exp)
            else:
                temp['sy'] = temp['e'].astype('float64')
        sxtemp = temp[['sx']].to_numpy(dtype='float64')
        sytemp = temp[['sy']].to_numpy(dtype='float64')

        #Calculate K0
        K0 = sytemp[200]/sxtemp[200]
        #Write K0 to txt file
        with open('k0optvals.txt', 'a') as kopt:
            kopt.write(str(K0))
            kopt.write("\t")
            kopt.write(str(itr))
            kopt.write("\n")
        kopt.close()

    #Calculate residual
    res = abs(0.245 - SS[200])

    os.remove("Ana.cal")

    return res

bnds = (0.245, 0.551)
x0 = 0.3341

#For each lift calculate the precompaction void ratio in order to generate the final void ratio by minimizing the
#residual using SciPy
for i in range(10):
    k0paramopt = minimize_scalar(tijmodelfn, bounds=bnds, args=(i,), method='bounded')
    with open('k0e0opt.txt', 'a') as opt:
        resultx=str(k0paramopt.x)
        opt.write(resultx)
        opt.write("\n")
        resultsuc=str(k0paramopt.message)
        opt.write(resultsuc)
        opt.write("\n")
    opt.close()
