##This program performs a parameter optimization for Subloading tij geomaterials model based on triaxial testing data (axial strain, deviator stress, volumetric strain, and a known confining stress) 
##and isotropic consolidation resutls (void ratio and mean pressure). 
import pandas as pd
import numpy as np
import os
from lmfit import Model
import corner
import matplotlib.pyplot as plt

# Function for error handling for "Subloading_tij.exe" that cannot be intrepreted by the reader
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

#Defining the Subloading tij model that recieves parameter inputs and outputs calculated Triaxial data and Isotropic Consolidation data
def tijmodelfn(astrains, Lamda, Rcs, Beta, a, N, Kappa, axk):

    #Defining non iteration parameters/constants these are standard Subloading tij parameters
    e0 = 0.3341
    e_max = 0.551
    pl = 1
    aAF = a
    aIC = a
    Qw0 = 0
    bw = 0
    Pa = 14.7
    powerIC = 2
    wgt_IC = 1
    pw = 1
    bxk = 0
    w_hpSoft = 1
    Poisson = 0.2
    axk = int(round(axk))

  #Writing parameter values to array to be written to the input files
    tparams = [Lamda, Kappa, Rcs, N, Poisson, Beta, e_max, e0, aAF, aIC, pl, axk, Qw0, bw, pw, bxk, Pa, powerIC, wgt_IC, w_hpSoft]
    tparams2 = np.array([0.003, 1.e-7, 0.0, 1.e-7, 0.0, 0.0, 1.e10, 1.0, 0.0])
    tparams3 = np.array([2.0, 5.0, 2.0])

  #Writing to/creating the Subloading_tij.exe parameter input file
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

  #Writing to/creating the Subloading_tij.exe stress path input file for first confining stress 5.0 psi 
  #Will output 300 points with a final strain of 15%
    with open('Path.txt', 'w') as path:
        path.write("NoPath= 1\n")
        path.write("InitialStresses= 5.0\t5.0\t5.0\t0.0\t0.0\t0.0\n")
        path.write("Triaxial\tDrain\tsrConstant\t0.15\t300\t1")

    path.close()

    #Running the Subloading_tij.exe program
    os.system(r'"C:\Users\Miles Drazkowski\Documents\Miles Drazkowski\Thesis\Optimization Paramaters\Without 20psi\with consol\SubloadingTij_a.exe"')

    #Reading from the Subloading_tij.exe output file "Ana.cal" 
    colspecs = [(205, 217), (129, 141)]
    names = ["εv", "q"]
    temp = pd.read_fwf('Ana.cal', encoding='shift_jis', colspecs=colspecs, names=names, skiprows=10)
    temp.fillna(value=10 ** 101.5, inplace=True)

  #Error handling for non-standard/undefined Ana.cal outputs
    if temp['q'].dtypes == 'object':
        temp['q'] = temp['q'].map(format_exp)

    if temp['εv'].dtypes == 'object':
        print(temp['εv'])
        temp['εv'] = temp['εv'].map(format_exp)

    #Dimensioning a results array and writing to it.
    SS = np.empty([1013, 2])

    SS[0:301, :] = temp[['εv', 'q']].to_numpy(dtype='float64')

  #Deleting results allowing Subloading_tij.exe to be run again. 
    os.remove("Ana.cal")

    #Defining starting void ratio for second confining stress
    tparams[7] = 0.3331

    #Rewriting parameter file using updataed parameter values, this might be optimized by replacing with a text search and replace.
    with open('Para.txt', 'w') as para:
        arglist = ["Model       1", "iVersion    3", "iCyclic&iCreep      0", "iAssoc      3", "Elastic     N"]
        desc = ["[1: Subloading tij,    2: Cam clay]", "[0: original,    3: alternative, 4: alternative4]",
                "[0: no-creep,        1: creep,       10: cyclic]",
                "[0: (AF+IC) for original, 3:  (AF+IC) for alternative,   1: (AF) for both iVersions]",
                "[H: non-linear HooKean,    N: no-tension elastic]"]
        tempdf = pd.DataFrame(desc, index=arglist, columns=[''])
        tempdfAsString = tempdf.to_string(header=False, index=True)
        para.write(tempdfAsString)
        para.write("\nLamda\tKappa\tRcs	N\tPoisson\tBeta\te_max\te0\taAF*\taIC*\tpl\taxk\tQw0\tbw*\tpw\tbxk\tPa\tpowerIC\twgt_IC\tw_hpSoft\n")
        for param in tparams:
            para.write("%s" % param)
            para.write("\t")
        para.write("\n\nLamda_alfa\tVoiddot_ref\tVoiddot_min\tVoiddot_ini\tbR_faster\tbR_slower\tVoiddot_max\tw_Vcrp\tw_Psi\n")
        for param in tparams2:
            para.write("%s" % param)
            para.write("\t")
        para.write("\n\na_cyc\tkm_cyc\tkn_cyc\n")
        for param in tparams3:
            para.write("%s" % param)
            para.write("\t")

    para.close()

   #Rewriting stress path file, confining stress 10.0 psi 
    with open('Path.txt', 'w') as path:
        path.write("NoPath= 1\n")
        path.write("InitialStresses= 10.0\t10.0\t10.0\t0.0\t0.0\t0.0\n")
        path.write("Triaxial\tDrain\tsrConstant\t0.15\t300\t1")

    path.close()

    #Running Subloading_tij.exe
    os.system(r'"C:\Users\Miles Drazkowski\Documents\Miles Drazkowski\Thesis\Optimization Paramaters\Without 20psi\with consol\SubloadingTij_a.exe"')

    #Reading results and error handling
    temp = pd.read_fwf('Ana.cal', encoding='shift_jis', colspecs=colspecs, names=names, skiprows=10)
    temp.fillna(value=10 ** 101.5, inplace=True)

    if temp['q'].dtypes == 'object':
        temp['q'] = temp['q'].map(format_exp)

    if temp['εv'].dtypes == 'object':
        print(temp['εv'])
        temp['εv'] = temp['εv'].map(format_exp)

    #Writing to Results array
    SS[301:602, :] = temp[['εv', 'q']].to_numpy(dtype='float64')

    os.remove("Ana.cal")

    #Initial void ratio for next stress path
    tparams[7] = 0.3233

    #Rewriting parameter file using updataed parameter values
    with open('Para.txt', 'w') as para:
        arglist = ["Model       1", "iVersion    3", "iCyclic&iCreep      0", "iAssoc      3", "Elastic     N"]
        desc = ["[1: Subloading tij,    2: Cam clay]", "[0: original,    3: alternative, 4: alternative4]",
                "[0: no-creep,        1: creep,       10: cyclic]",
                "[0: (AF+IC) for original, 3:  (AF+IC) for alternative,   1: (AF) for both iVersions]",
                "[H: non-linear HooKean,    N: no-tension elastic]"]
        tempdf = pd.DataFrame(desc, index=arglist, columns=[''])
        tempdfAsString = tempdf.to_string(header=False, index=True)
        para.write(tempdfAsString)
        para.write("\nLamda\tKappa\tRcs	N\tPoisson\tBeta\te_max\te0\taAF*\taIC*\tpl\taxk\tQw0\tbw*\tpw\tbxk\tPa\tpowerIC\twgt_IC\tw_hpSoft\n")
        for param in tparams:
            para.write("%s" % param)
            para.write("\t")
        para.write("\n\nLamda_alfa\tVoiddot_ref\tVoiddot_min\tVoiddot_ini\tbR_faster\tbR_slower\tVoiddot_max\tw_Vcrp\tw_Psi\n")
        for param in tparams2:
            para.write("%s" % param)
            para.write("\t")
        para.write("\n\na_cyc\tkm_cyc\tkn_cyc\n")
        for param in tparams3:
            para.write("%s" % param)
            para.write("\t")

    para.close()

    #Rewriting stress path file, confining stress 30.0 psi 
    with open('Path.txt', 'w') as path:
        path.write("NoPath= 1\n")
        path.write("InitialStresses= 30.0\t30.0\t30.0\t0.0\t0.0\t0.0\n")
        path.write("Triaxial\tDrain\tsrConstant\t0.15\t300\t1")

    path.close()

     #Running Subloading_tij.exe
    os.system(r'"C:\Users\Miles Drazkowski\Documents\Miles Drazkowski\Thesis\Optimization Paramaters\Without 20psi\with consol\SubloadingTij_a.exe"')

    #Reading results and error handling
    temp = pd.read_fwf('Ana.cal', encoding='shift_jis', colspecs=colspecs, names=names, skiprows=10)
    temp.fillna(value=10 ** 101.5, inplace=True)

    if temp['q'].dtypes == 'object':
        temp['q'] = temp['q'].map(format_exp)

    if temp['εv'].dtypes == 'object':
        print(temp['εv'])
        temp['εv'] = temp['εv'].map(format_exp)

    #Writing to Results array  
    SS[602:903, :] = temp[['εv', 'q']].to_numpy(dtype='float64')

    os.remove("Ana.cal")

    #Initial void ratio for Isotropic consolidation. 
    tparams[7] = 0.3341
  
    #Rewriting parameter file using updataed parameter values
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

    #Writing Stress Path file for isotropic consolidation with intial consolidation of 5.0 psi, loading to 50.0 psi and unloading to 5.0 psi
    with open('Path.txt', 'w') as path:
        path.write("NoPath= 2\n")
        path.write("InitialStresses= 5.0\t5.0\t5.0\t0.0\t0.0\t0.0\n")
        path.write("Consolidation\tOneDimensional\t50\t540\t10\n")
        path.write("Consolidation\tOneDimensional\t5\t540\t10")

    path.close()

    #Running subloading_tij.exe
    os.system(r'"C:\Users\Miles Drazkowski\Documents\Miles Drazkowski\Thesis\Optimization Paramaters\Without 20psi\with consol\SubloadingTij_a.exe"')

    #Reading results and error handling
    colspecs = [(231, 241), (154, 167)]
    names = ["e", "s1"]
    temp = pd.read_fwf('Ana.cal', encoding='shift_jis', colspecs=colspecs, names=names, skiprows=10)
    temp.fillna(value=10 ** 101.5, inplace=True)

    if temp['s1'].dtypes == 'object':
        print(temp['s1'])
        try:
            temp['s1'] = temp['s1'].astype('float64')
        except:
            temp['s1'] = temp['s1'].map(format_exp)
        else:
            temp['s1'] = temp['s1'].astype('float64')

    if temp['e'].dtypes == 'object':
        print(temp['e'])
        try:
            temp['e'] = temp['e'].astype('float64')
        except:
            temp['e'] = temp['e'].map(format_exp)
        else:
            temp['e'] = temp['e'].astype('float64')

    #Writing to results array
    SS[904:1013, :] = temp[['e', 's1']].to_numpy(dtype='float64')

    os.remove("Ana.cal")

    SS = SS.astype('float64')

    #Converting from percent strain to decimal strain
    SS[0:903, 0] = SS[0:903, 0] / 100

    #Flattening results into 1D array
    SS = SS.flatten('F')

    return SS

#Reading triaxial data and loading in data arrays
data = np.loadtxt("Triaxial Test Data and Consol.20.txt", dtype='float64')
xdata = data[:, 0]
xdata = np.resize(xdata, 2026)
ydata = data[:, 1:3]
ydata = ydata.flatten('F')
#Loading weights for error calculation (optional)
fweights = np.loadtxt('weightsforfitting.4124.txt', dtype='float64')
fweights=fweights.flatten('F')

#Defining parameter array and Model class
tparams = ['Lamda', 'Rcs', 'a', 'Beta', 'N', 'Kappa', 'axk']
tijmodel = Model(tijmodelfn, param_names=tparams)

print(tijmodel.param_names)

#Define Parameter class and add to Model
paramin = tijmodel.make_params(Lamda=dict(value=0.02073268, max=0.025, min=0.019),
                               Rcs=dict(value=5.08521251, max=5.5, min=4.5),
                               Beta=dict(value=1.05154567, max=3, min=1), a=dict(value=103.324742, max=1000, min=10),
                               N=dict(value=0.35760210, max=0.551, min=0.3323), Kappa=dict(value=0.00160315, max=0.0019, min=0.0013),
                               axk=dict(value=1, max=50, min=1))

#Optimization algorithm keyword arguements and optimization 
emcee_kws = dict(steps=4500, nwalkers=500, burn=1000, is_weighted=False)
#difev_kws = dict(atol=0, tol=0, maxiter=10000000)
paramopt = tijmodel.fit(data=ydata, astrains=xdata, nan_policy='omit', params=paramin, method='emcee', fit_kws=emcee_kws)
#paramopt = tijmodel.fit(data=ydata, astrains=xdata, weights=fweights, nan_policy='omit', params=paramin, method='differential_evolution', fit_kws=difev_kws)


#Writing Optimization results to file and printing results as well as confidence intervals
with open('Optparatest20D.4-15-24.txt', 'w') as opt:
    opt.write(paramopt.fit_report())
opt.close()
print(paramopt.fit_report())
print(paramopt.ci_report())

#Plotting emcee emcee acceptance fraction
plt.plot(paramopt.acceptance_fraction, 'o')
plt.title('Test20.0')
plt.xlabel('walker')
plt.ylabel('acceptance fraction')
plt.savefig('AcceptanceFractest20.update.3.png')

#Plotting posterior distributions
para_plot = corner.corner(paramopt.flatchain, labels=paramopt.var_names, truths=list(paramopt.params.valuesdict().values()))
para_plot.savefig('paramplotvoltest20.update.3.png')

plt.show()
ara_plot.show()
