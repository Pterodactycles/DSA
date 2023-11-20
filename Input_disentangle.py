# shift-and-add & grid disentangling, by Tomer Shenar, with contributions from Matthias Fabry & Julia Bodensteiner
# 21.11.2022, V1.0; feel free to contact at T.Shenar@uva.nl or tomer.shenar@gmail.com for questions/inquires
# Input file

# INTRODUCE ROUND WHEN SAVING FILE NAMES!!!!


import numpy as np
import json
import argparse


###############################
##### STAR AND DATA INFO ######
###############################

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--system", help="Name of the system")
parser.add_argument("-rv", "--find_rv", help="Whether or not to find the radial velocity, boolean")
parser.add_argument("-o", "--order", help="Helps pick predefined values from dictionaries below.")
args = parser.parse_args()

### Name of object (just for file names)
StarName = args.system or input("System: ")
find_rv = bool(args.find_rv)
Order = f"order_{args.order}"


# Path to data (folder where all spectra are stored):
# ObsPath = f'/home/users/cmb255/DBiSS/data/{StarName}/ordered/{Order}/'
ObsPath = '/home/users/cmb255/DBiSS/4950_4965/'

# Output_dir = f'/home/users/cmb255/DeBiSS/data/{StarName}/disentangled'
Output_dir = "/home/users/cmb255/DBiSS/disentangled"

### JSON PATH
json_file = f'/home/users/cmb255/DBiSS/data/orbital_parameters.json'

### Type of data format. There are two options:
### OPTION 1: ObsFormat = 'TXT' 
### assumes that the observations are in ascii format, each file containing 2-column tables of wave & normalised flux. 
### In addition, the observation directory MUST contain a file called 'ObsDat.txt', which has the following format:
###   MJD          obsname
###   xxxx          PATH1
###   yyy           PATH2
### The paths should be either absolute or relative to the directory in which the script is stored.
### OPTION 2: ObsFormat = 'FITS' 
### The script will look for ALL fits files in the given directory. 
### The script will attempt to retrieve the dates from the fits headers using a user-specified header keyword
### IMPORTANT NOTES:
### 1. It doesn't matter if the dates are "MJD", "JD", "HJD", etc -- important is that the T0 provided by the user matches this!
### 2. For the "fits" option, I include a few built-in functions to read e.g. HERMES, X-SHOOTER, FEROS spectra.... 
### User should feel free to update the reading of the files!

ObsFormat = 'TXT'

#Only important if ObsFormat='FITS'
MJDHeader = 'MJD-OBS'



###############################
##### SYSTEM PROPERTIES #######
###############################


#Number of components; possible up to four components
CompNum = 2

# Orbital parameters 
### P, T0, ecc, omega, and gamma cannot be derived with current version and are assumed by the user
### K1, K2 can be explored via chi2 if required by the user, but initial guesses should be given.
### Important: omega is defined via vr(1) = Gamma1 + K1*(cos(nu) + ecc * cos(omega) )
### If results don't make sense, your best bet is to set omega --> omega + pi
Orbital_Params = {
####### MUST BE FILLED BELOW ALWAYS (inner binary) ######   
    'Period': 18,
    'T0': 0.,
    'ecc': 0.3,
    'omega': 60., #IN DEGREES!
    'Gamma': 0.,
    'K1': 87.,
    'K2': 135.  ,
####### Only for triples / quadruples: refers to outer orbit of companion/binary AROUND the inner binary ######       
####### IF outer period very long/negligible, set PeriodOut = 1E10 to neglect long-term motion     
    'PeriodOut': 200.,
    'T0Out' : 34.,
    'eccOut' : 0.3,
    'omegaOut' : 160.,
    'KOut' : 29.,
####### Only for Quadruples: refers to orbit of 2nd binary in the system  ######           
    'Period_2' : 0.,
    'T0_2' : 0.,
    'ecc_2' : 0.,
    'omega_2' : 0.,  
    'K3' : 0. ,
    'K4' : 0.}

### SUCH DATA IS PLACED IN A JSON FILE AND RETRIEVED BASED ON THE SYSTEM NAME, FOR CONVENIENCE
with open(json_file) as file:
    data = json.load(file)
Orbital_Params = data[f'{StarName}']
print(Orbital_Params)


# MUST INTRODUCE COMPNUM TO AVOID FUNNY RATIOS!
# Vector of light ratios, [l2, l3, l4], i.e. flux_i / sum(flux). Assumed constant throughout range.
lguessVec = [0.5, 0., 0.] 

if CompNum==2:    
    lguess1 = 1.-lguessVec[0]
elif CompNum==3:
    lguess1 = 1. - lguessVec[0] - lguessVec[1]
elif CompNum==4:
    lguess1 = 1. - np.sum(np.array(lguessVec))



#Where to measure S2N, only important for  defining continuum and weighting of spectra when co-adding (not critical)
continuumdictionary = {
    "order_11" : {
        "S2Nblue" : 653.7,
        "S2Nred" : 654.2
    },
    "order_21" : {
        "S2Nblue" : 582,
        "S2Nred" : 584
    },
    "order_28" : {
        "S2Nblue" : 544.9,
        "S2Nred" : 545.4
    },
    "order_29" : {
        "S2Nblue" : 543.7,
        "S2Nred" : 544.3
    },
    "order_30" : {
        "S2Nblue" : 535.5,
        "S2Nred" : 536.1
    },
    "order_34" : {
        "S2Nblue" : 514,
        "S2Nred" : 516
    },
    "order_39" : {
        "S2Nblue" : 495.5,
        "S2Nred" : 496
    },
    "order_51" : {
        "S2Nblue" : 448,
        "S2Nred" : 450
    }
}

S2Nblue = continuumdictionary[Order]["S2Nblue"]
S2Nred = continuumdictionary[Order]["S2Nred"]



################################
##### Disentangling options ####
################################


# Run grid disentangling? 
# If TRUE: will conduct grid disentangling and derive Ks
# If FALSE: will only peform separation using input K1,K2 
GridDis = True

# Define grid search (only important if GridDis = True). 
# For setting K1, K2, K3, K4 search arrays: Karr = np.arange(IniFacK*K, FinFacK*K, Dense)
# Current version only works for first two columns (K1, K2)
# If DenseKArr[i] = 1, then the search is "1D", i.e. K is fixed to the value specified by the user.
DenseKArr = [1, 1, 12, 1]
IniFacKArr = [0.1, 0.1, 0.1, .3]
FinFacKArr = [2., 2., 2., 2.]



# Number of iterations
### IMPORTANT NOTES:
### 1. Ideally convergence could be determined via a condition on EPS (See above). However, a suitable condition could not yet be developed
### --> User needs to judge when results are "sufficiently converged", by either comparing the results for different itr numbers, or using
###     options below.
### 2. itrnumlim is the number of iterations per K1,K2 pair; NumItrFinal is the number of iterations for the final separation, 
###    after K1, K2 have been derived / set. Often, itrnumlim < NumItrFinal, for speed, and since individual iterations occur on individual lines.
### 3. See documentation for tips and insights about number of iterations.

itrnumlim =30
NumItrFinal = 1000


# If StrictNegA = True, enforce disentangled spectra to be below continuum except for prespecified regions (given in array).
# Below continuum = ForceNegSigma "sigmas" below continuum. 
# HIGHLY RECOMMENDED for OB-type stars -- otherwise, output often exhibits cosmetic "wings" and continuum offsets.
# For WR stars typically "False" is better.
# "Positive regions" should be regions with expected emission lines etc. For WR stars, 

ForceNegSigma = 2.

StrictNegA = True

#Only relevant if StrictNegA=True
PosLimCondA = np.array([ 
            [3968., 3969.]    
              ])

# Same as StrictNegA for secondary
StrictNegB = True

PosLimCondB = np.array([ 
            [3968., 3969.]
              ])

# Same as StrictNegA for tertiary
StrictNegC = True

PosLimCondC = np.array([ 
            [3968., 3969.]
              ])

# Same as StrictNegA for fourth companion
StrictNegD = True

PosLimCondD = np.array([ 
            [3968., 3969.]
              ])



# Define regions where the solution is allowed to be above continuum (where emission is expected)
PosLimCondNeb = np.array([ 
            [3968., 3971.],
            [4025., 4027.],            
            [4100.5, 4103],
            [4143., 4145],            
            [4339., 4342],
            #[4340., 4345],            
            [4387., 4391.5],       
            #[4385., 4393],             
            #[4335., 4345.],            
            #[4465., 4482.]       
            [4470., 4473.]              
            #[4135., 4150.],            
            #[4228., 4238.],
            #[4330., 4355.], 
              #[4381., 4396.],              
              #[4465., 4485.],  
              #[4840., 4880.], 
              #[4916., 4930.],  [5010., 5024.],
              #[5160., 5180.], [5190., 5210.], [5225., 5240.], [5270., 5290.],
              #[5310, 5325.], [5255., 5370.], [5520., 5540.], [6140., 6160.], [6230., 6260.], [6310., 6330], 
              #[6340., 6390.], [6410., 6460.], [6510, 6520.], [6553., 6574.], [7505, 7525], [7766., 7786.]
              ])



# Plot fits between disentangled spectra, their sum, and the observations at RV extremes.
# Highly recommended for sanity checks and presentation in papers.
# The plot is shown for the K1, K2 pair most closely matching (Velo_plot_usrK1_ext, Velo_plot_usrK2_ext, ...) given by the user.
# Recommended: True
PLOTEXTREMES = True
Velo_plot_usrK1_ext =  Orbital_Params['K1']
Velo_plot_usrK2_ext=  Orbital_Params['K2']
Velo_plot_usrK3_ext=  Orbital_Params['KOut']
Velo_plot_usrK4_ext=  Orbital_Params['K4']

# line width and figsize for "Extreme plots"
linewidExt = 7
ExtremesFigSize = (10, 10)




# Plot fits between disentangled spectra, their sum, and each epoch of observation. 
# Useful to examine all data; rather tedious but important for critical systems (e.g., black holes!)
# The plot is shown for the K1, K2 pair most closely matching (Velo_plot_usrK1, Velo_plot_usrK2) given by the user.
# Recommended: False
PLOTFITS = False
Velo_plot_usrK1 = Orbital_Params['K1']
Velo_plot_usrK2 = Orbital_Params['K2']
Velo_plot_usrK3 = Orbital_Params['KOut']
Velo_plot_usrK4 = Orbital_Params['K4']

# Plot convergence plot
# If True, will produce converge plot, i.e. EPS vs. itr for each run.
# EPS = SUM(DisSpec[i+1] - DisSpec[i]), where the maximum on all components is taken. 
# Recommended: False
PLOTCONV = False


# Plot disentangled spectra after each "N_Iteration_Plot" iterations; helpful for judging convergence.
# Recommended: False
PLOTITR = False
N_Iteration_Plot = 100


# Type of interpolation in interp1d (see python doc for options);
# 'linear' can lead to artificial increase of S/N due to interpolation
# 'cubic' performs better, but is slower. 
InterKind='cubic'

# Region for fitting parabola of chi2 in index steps from minimum
ParbSize=2



################################
##### Disentangling lines ######
################################

# User chooses in which line/lines the K1,K2 search should occur.
# All that is required is:
# 1. Ranges = [ [l1, l2], [l3, l4], ...]
# 2. Rangestr = 'xxx' -- used below to pick range, but either way, needs to be specified for file-saving purposes.
# For convenience, typical lines (for massive stars) are provided below. 
# USERS: feel free to edit wavelength regions below!!!
# IMPORTANT: the final ranges used/plotted are NOT identical to those provided by the user: the script reduces them to ensure that edge issues are avoided.
# The reduction depends on K1, K2; the user should judge (e.g., using "PLOTEXTREMES") that the lines are well covered and reach continuum at both edges.
# Ideally, disentangled region should be line-dominated (to enhance the signal on chi2), but certainly reach continuum at the edges.

#Rangestr = 'Hdelta'
#Rangestr = 'Hgamma'
#Rangestr = 'Hbeta'
#Rangestr = 'Halpha'
#Rangestr = 'Balmer'
#Rangestr = 'Balmer_noHalpha'
#Rangestr = 'HeI'
# Rangestr = 'HeI4472'
#Rangestr = 'HeI4122'
#Rangestr = 'HeI4009'
#Rangestr = 'HeI4026'
#Rangestr = 'HeI4144'
#Rangestr = 'HeII4200'
#Rangestr = 'HeI4388'
#Rangestr = 'HeII4546'
#Rangestr = 'HeI5878'
#Rangestr = '4120Region'
#Rangestr = '4020Region'
#Rangestr = 'IronEmission'
#Rangestr = 'OIII'
#Rangestr = 'OIII8446'
#Rangestr = 'HI8367'
#Rangestr = 'Fe4584'
#Rangestr = 'Fe5168'     
#Rangestr = 'Fe5192'   
#Rangestr = 'Fe5234'   
#Rangestr = 'Fe5275'  
#Rangestr = 'Fe5316'
##Rangestr = 'Fe5362'
#Rangestr = 'AllHeI'
#Rangestr = 'AllHeII'
#Rangestr = 'AllHe'
#Rangestr = 'Indiv'
Rangestr = 'specific'

##### Define ranges corresponding too the strings above.... CHANGE IF NEEDED
waveranges = {
    "RangeHa" : [655.3, 657.0],
    "RangeHb" : [484.0, 487.7],
    "RangeHg" : [431.0, 437.0],
    "RangeHd" : [407.0, 414.0],
    "RangeHeI5878" : [586.9,588.1],
    "RangeHeI4472" : [445.7,448.9],
    "RangeHeI4144" : [412.0,417.0],
    "RangeOIIIJulia" : [776.0,778.5],
    "RangeFe4584" : [458.0,458.8],
    "RangeFe4584" : [458.0,458.6],
    "RangeFe5168" : [516.2,517.4],
    "RangeFe5192" : [519.0,520.5],
    "RangeFe5234" : [523.0,524.0],
    "RangeFe5275" : [526.8,528.2],
    "RangeFe5316" : [531.0,532.2],
    "RangeFe5362" : [535.8,536.7],
    "RangeOIII8446" : [843.8,845.5],
    "RangeHI8367" : [836.7,850.0],
    "RangeHeI4122" : [411.5,412.7],
    "RangeHeI4009" : [400.3,401.8],
    "RangeHeI4026" : [400.0,405.0],
    "RangeHeI4388" : [436.5,441.0],
    "RangeHeII4545" : [451.5,456.5],
    "RangeHeII4200" : [418.5,421.5],
    "order_11" : [654,658],
    "order_21" : [584,588],
    "order_28" : [543,550],
    "order_29" : [538,543],
    "order_30" : [535,539],
    "order_34" : [516,519.5],
    "order_39" : [495,496],
    "order_51" : [445,450],
}

# Define "Ranges" list based on user's choices from above.
Rangelists = {
    "specific" :[ waveranges[Order]],
    "Hdelta" : [waveranges["RangeHd"]],
    "Hgamma" : [waveranges["RangeHg"]],
    "Hbeta" : [waveranges["RangeHb"]],
    "Halpha" : [waveranges["RangeHa"]],
    "Balmer" : [
        waveranges["RangeHd"],
        waveranges["RangeHg"],
        waveranges["RangeHb"],
        waveranges["RangeHa"]
    ],
    "Balmer_noHalpha" : [
        waveranges["RangeHd"],
        waveranges["RangeHg"],
        waveranges["RangeHb"]
    ],
    "HeI" : [
        waveranges["RangeHeI4009"],
        waveranges["RangeHeI4026"],
        waveranges["RangeHeI4122"],
        waveranges["RangeHeI4144"],
        waveranges["RangeHeI4388"],
        waveranges["RangeHeI4472"]
    ],
    "HeI4472" : [waveranges["RangeHeI4472"]],
    "HeI5878" : [waveranges["RangeHeI5878"]],
    "HeI4144" : [waveranges["RangeHeI4144"]],
    "HeII4200" :  [waveranges["RangeHeII4200"]],
    "4120Region" : [waveranges["RangeHeI4144"]],
    "IronEmission" : [
        waveranges["RangeFe5168"],
        waveranges["RangeFe5192"],
        waveranges["RangeFe5275"],
        waveranges["RangeFe5316"]
    ],
    "Fe4584" : [waveranges["RangeFe4584"]],
    "Fe5168" : [waveranges["RangeFe5168"]],
    "Fe5192" : [waveranges["RangeFe5192"]],
    "Fe5234" : [waveranges["RangeFe5234"]],
    "Fe5275" : [waveranges["RangeFe5275"]],
    "Fe5316" : [waveranges["RangeFe5316"]],
    "Fe5362" : [waveranges["RangeFe5362"]],
    "OIII" : [waveranges["RangeOIIIJulia"]],
    "OIII8446" : [waveranges["RangeOIII8446"]],
    "HI8367" : [waveranges["RangeHI8367"]],
    "HeI4122" : [waveranges["RangeHeI4122"]],
    "HeI4009" : [waveranges["RangeHeI4009"]],
    "HeI4026" : [waveranges["RangeHeI4026"]],
    "HeI4388" : [waveranges["RangeHeI4388"]],
    "HeII4546" : [waveranges["RangeHeII4545"]],
    "AllHe" : [
        waveranges["RangeHeI4026"],
        waveranges["RangeHeII4200"],
        waveranges["RangeHeI4472"],
        waveranges["RangeHeII4545"]
    ],
    "AllHeI" : [
        waveranges["RangeHeI4026"],
        waveranges["RangeHeI4144"],
        waveranges["RangeHeI4388"],
        waveranges["RangeHeI4472"],
        waveranges["RangeHeII4545"]
    ],
    "AllHeII" : [
        waveranges["RangeHeII4200"],
        waveranges["RangeHeII4545"]
    ],
    "Indiv" : [
        waveranges["RangeHeII4200"], 
        waveranges["RangeHeI4472"]
    ]
}

Ranges = Rangelists[Rangestr]


################################
##### Fancy options   ##########
################################     


#Clean cosmics?
CleanCos = False


#Renormalise spectra at pre-specified points. 
Renormalise = False
NormPoints = [3961., 4006., 4016., 4038., 4088., 4116., 4129., 4138., 4154., 4195., 4210., 4328., 4362., 4386., 4400., 4462., 4490., 4494., 4530., 4557., 4560]

# Nebular line handling?
NebLines = False


####################
