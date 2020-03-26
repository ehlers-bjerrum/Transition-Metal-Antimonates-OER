#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 12:28:52 2019

@author: Anders Bjerrum
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from numpy.linalg import inv
from numpy.linalg import cholesky, det, lstsq
from scipy.optimize import minimize
from ase.db import connect
from TBModel import TBModel
from time import sleep
import sys
from calcAdsorptionEnergy import calcAdsorptionEnergy
from loadDatabases import loadDatabases
import constants as c

### Variable encoding
# The ordering is more important than the actual values
r = 0.5
varO = 0
varOH = 1 + r*1j
varH2O = 2 + 2*r*1j
varVac = 3


def GPpredictor(voltage,state=False,conf=False):
  ## Predicts possible groundstate coverages for the facet under study.
  # INPUT: Applied bias vs RHE.
  # OPTIONAL INPUT: conf, input a possible surface configuration to get the predicted energy
  # OPTIONAL INPUT: state, input a possible state to get the predicted energy
  # OUTPUT: (1): A list of candidates for the groundstate in descending order of energy. (2): The predicted energies.
  # (3) The uncertainty in the prediction
  
  ########## SETUP ##########

  # Only used for referencing the formation energies. The label refers to a row in the references database references.db
  material = 'Mn'

  ## Connect to ASE database. 
  # Each entry should have a meta data value flagged OCCUPANCY. In our case this is formatted as
  # eg. 'OH3A,OH3B,OH5,O7,O9A,O9B' where the letters indicate the occupancy, and the numbers
  # indicate the site. The flag should be called as row.data.OCCUPANCY
  
  db1 = connect('/home/cat/s144109/dataProcessing/MnSb2O6/databases/GPTraining2.db')
  dbs = [db1]


  ######## PREDICTIONS ############
  
  ## Call the function to train the Gaussian Processes (GP) model on the training data
  model = trainModel(dbs,material='Mn')
  
  ## Generate possible states. A state is configuration (coverage)  encoded using the variable encoding.
  X = getPossibleStates()

  ## Predict the energy of a given surface configuration (coverage)
  if conf is not False:
    return predictEnergy(getState(conf.split(',')),model,voltage)
  
  ## Predict the energy of a given state
  if state is not False:
    return predictEnergy(state,model,voltage)
  

  ## Get training states
  XD = getDoneCalc(dbs)
  ## Remove the training states from the trial states X
  XR = removeDoneCalc(X,XD)

  ## Find all states predicted to have an energy within 2 std. deviations from the hull
  [candidates,dG_T,S] = predictCandidates(model,XR,dbs,material,voltage)
  ## Sort the states according to their formation energy
  [[candidates,S],dG_T] = sortCandidates([candidates,S],dG_T)
  

  return [candidates,dG_T,S]

  

def getBestActivity(dbs,material,voltage,model):
  # Predicts the maximum activity of the training data states. We define the activity as
  # ACTIVITY = - (FORMATION_ENERGY + ACTIVATION_ENERGY)
  # A high activity implies a low formation energy and low activation energy for OO coupling
  [OS,I,R] = loadDatabases()
  [conf_full,dG_Ad_full,Y_full,X_full] = getTrainingData(dbs,R,material)
  
   
  # For each coverage, get the most active site (Occupied by Oxide) 
  activity = []
  activityError = []
  site = []
  for state in X_full:
    [A,Asig] = predictActivity(state,model,voltage)
    minIndex = A.index(min(A))
    activity.append(A[minIndex])
    activityError.append(Asig[minIndex])
    site.append(minIndex)
  

  # Sort the coverages according to their most active site
  [[VASPState,VASPError,site],VASPActivity] = sortCandidates([conf_full,activityError,site],activity)
  #print('Most active surfaces')
  #for j in range(0,20):
  #  print(VASPState[j])
  #  print(VASPActivity[j])
  #  print(VASPError[j])
  #  print(site[j])

  # Return the highest activity of the training set
  bestActivity = VASPActivity[0]
  error = VASPError[0]

  return [bestActivity,error]  

def predictActivity(state,model,voltage):
  # Get the activity of all sites for a given state 
  getBarrier = lambda v,d: -0.58*v-0.48*d+2.41 # Barrier calculation for OOH formation see C. Dickens and J. NÃrskov, J. Phys. Chem. C, volume 123, number 31, pages 18960-18977, 2019., DOI: 10/1021/acs.jpcc.9b03830

  [G,sig] = predictEnergy(state,model,voltage) # Predict the energy of the state
  A = []
  Asig = []
  for j in range(0,len(state)):
    tempState = list(state)
      
    # Calculate the activity if the site has an oxide (a possible active site for OOH formation)  
    if tempState[j] == varO:
      tempState[j] = varOH
      [GTemp,sigTemp] = predictEnergy(tempState,model,voltage)
      desc = G-GTemp + voltage
      barrier = getBarrier(voltage,desc)
      activity = G + barrier
      error = np.sqrt(0.48**2*(sigTemp**2 + sig**2)+sigTemp**2)

      A.append(activity)
      Asig.append(error)
    
  return [A,Asig]

def predictCandidates(model,XR,dbs,material,voltage):
  ## Selects states that should be submitted for a VASP calculation
  # Compares predicted formation energies with the most stable state found so far, and selects a state for calculation if the predicted formation energy of that state is within 2 standard deviations from the most stable state.
  # Compares predicted activites with the best activity of the training data, and selects states for VASP calculation if the predicted activity is within 1 standard deviation from the most active site
  
  ## Gets the most stable state from the training set
  minEnergyVASP = getMinimumEnergy(dbs,material,voltage)
  ## Gets the most active state from the training set
  [BestActivityVASP,ActivityErrorVASP] = getBestActivity(dbs,material,voltage,model)

  candidates = []
  dG_T = []
  S = []
  for cand in XR:
    [dG,sig] = predictEnergy(cand,model,voltage)
    [A,Asig] = predictActivity(cand,model,voltage)
    if (dG-minEnergyVASP) < sig*1.95: # Selects stable states
      candidates.append(getConf(cand))
      dG_T.append(dG)
      S.append(sig)
      continue
    for j in range(0,len(A)): # Selects active states
      if (A[j]-BestActivityVASP) < np.sqrt(Asig[j]**2 + ActivityErrorVASP**2):
        candidates.append(getConf(cand))
        dG_T.append(dG)
        S.append(sig)
        break

  return [candidates,dG_T,S]

def getConf(state):
  ## Gets the configuration of a given state
  conf = []
  sites = ['3A','3B','5','7','9A','9B']
  for j in range(0,len(state)):
    if state[j] == varOH:
      conf.append('OH'+sites[j])
    elif state[j] == varO:
      conf.append('O'+sites[j])
    elif state[j] == varVac:
      conf.append('V'+sites[j])
    elif state[j] == varH2O:
      conf.append('H2O'+sites[j])
  return conf


def getMinimumEnergy(dbs,material,voltage):
  ## Gets the  minimum formation energy from training states

  ## Read in the References database R
  [OS,I,R] = loadDatabases()
  ## Get the training data from the training SQL-lite database
  [conf_full,dG_Ad_full,Y_full,X_full] = getTrainingData(dbs,R,material)  
  
  ## Subtract the electrostatic potential energy from the 0V vs RHE formation energy
  energyAtRHE_VASP = []
  for j in range(0,len(dG_Ad_full)):
    U = 0
    for site in X_full[j]:
      if site==varO:
        U = U - 2*voltage
      elif site==varOH:
        U = U - voltage
    energyAtRHE_VASP.append(dG_Ad_full[j] + U)

  minEnergyVASP = min(energyAtRHE_VASP)
  
  return minEnergyVASP

def sortCandidates(toSort,ref):
  ## Sorts a list of elements toSort based on another list ref

  if toSort:
    sort = []
    for temp in toSort:
      temp = [x for _,x in sorted(zip(ref,temp))]
      sort.append(temp)
    ref.sort()
    return [sort,ref]
  else:
      return [False,False,False]



def kernel(X1, X2, l, sigma_f):
## Kernel definition for GP
'''
 # Isotropic squared exponential kernel. Computes 
  a covariance matrix from points in X1 and X2.
  Args:
    X1: Array of m points (m x d).
    X2: Array of n points (n x d).
  Returns:
    Covariance matrix (m x n).
  '''
  sqdist = np.sum(X1*np.conj(X1), 1).reshape(-1, 1) + np.sum(X2*np.conj(X2), 1) - np.dot(np.conj(X1), X2.T) - np.dot(X1, np.conj(X2).T)
  
  sqdist = np.real(sqdist) 
  
  return sigma_f**2 * np.exp(-0.5 / l**2 * sqdist)



def posterior_predictive(X_s, X_train, Y_train, l, sigma_f, noise):
## Makes a posterior prediction using the training data and the kernel
'''  
  Computes the suffifient statistics of the GP posterior predictive distribution 
  from m training data X_train and Y_train and n new inputs X_s.  
  Args:
    X_s: New input locations (n x d).
    X_train: Training locations (m x d).
    Y_train: Training targets (m x 1).
    l: Kernel length parameter.
    sigma_f: Kernel vertical variation parameter.
    sigma_y: Noise parameter.
  Returns:
  Posterior mean vector (n x d) and covariance matrix (n x n).
  '''
  K = kernel(X_train, X_train, l, sigma_f) + np.eye(len(Y_train))*noise**2 #+ 1e-2 * np.eye(len(Y_train))
  K_s = kernel(X_train, X_s, l, sigma_f)
  K_ss = kernel(X_s, X_s, l, sigma_f) + 1e-8 * np.eye(len(X_s))

  K_inv = inv(K)

  # Equation (4)
  mu_s = K_s.T.dot(K_inv).dot(Y_train)
  # Equation (5)
  cov_s = K_ss - K_s.T.dot(K_inv).dot(K_s)
  return mu_s, cov_s

def nll_fn(X_train, Y_train):
  ## Computes the negative log-likelihood of the model
  '''
  Returns a function that computes the negative marginal log-
  likelihood for training data X_train and Y_train and given 
  noise level. 
  Args:
    X_train: training locations (m x d).
    Y_train: training targets (m x 1).
    noise: known noise level of Y_train.
    naive: if True use a naive implementation of Eq. (7), if 
      False use a numerically more stable implementation.   
  Returns:
    Minimization objective.
  '''

  def nll_stable(theta):
    # Numerically more stable implementation of Eq. (7) as described
    # in http://www.gaussianprocess.org/gpml/chapters/RW2.pdf, Section
    # 2.2, Algorithm 2.1.
    K = kernel(X_train, X_train, l=theta[0], sigma_f=theta[1]) + theta[2]**2*np.eye(len(Y_train))
    
    L = cholesky(K)
    return np.sum(np.log(np.diagonal(L))) + 0.5 * Y_train.T.dot(lstsq(L.T, lstsq(L, Y_train,rcond=None)[0],rcond=None)[0]) + 0.5 * len(X_train) * np.log(2*np.pi)
    
  return nll_stable



def trainModel(dbs,material='Mn',leaveOneOutState=False):
  ## Trains the GP model on the training data (Conditions the multivariate normal distribution, see http://www.gaussianprocess.org/gpml/chapters/RW2.pdf) 


  # Load the reference databases
  [OS,I,R] = loadDatabases()
  # Get the possible states
  X = getPossibleStates()
  # Get the adsorption energies and states of the training data
  [conf_train,dG_Ad_train,Y_train,X_train] = getTrainingData(dbs,R,material)

  # Enforce the symmetries of the adsorption sites
  [X_train,Y_train] = enforceSymmetry(X_train,Y_train)
  
  # Calculate statistical parameters for standardizing the data
  Y_mean = np.mean(Y_train)
  Y_std = np.std(Y_train)
  Y_train = (Y_train - Y_mean)/Y_std

  # Minimize 
  res = minimize(nll_fn(np.asarray(X_train), np.asarray(Y_train)), x0=[1, 1, 1],bounds=((1e-1, None),(1e-1, None), (1e-3,None)),method='L-BFGS-B')
  dG_ML, cov_s = posterior_predictive(np.asarray(X), np.asarray(X_train), np.asarray(Y_train), *res.x)
  
  dG_ML = dG_ML*Y_std + Y_mean
  Y_train = Y_train*Y_std + Y_mean

  print(*res.x)
  error = np.sqrt(np.diag(cov_s))*Y_std
  for i in range(0,len(X)):
    if X[i] in X_train:
      j = X_train.index(X[i])
      dG_ML[i] = Y_train[j]

  return [X,dG_ML,error]

def enforceSymmetry(X,Y):
  ## For the sites under consideration, we have the symmetry that swapping 3A and 3B AND swapping 9A and 9B, yields the same coverage (state). This function ensures that both symmetrically equivalent states are sampled, when one of them has been.

  for j in range(0,len(X)):
    state = X[j]
    tempState = list(state)
    tempState[0] = state[1]
    tempState[1] = state[0]
    tempState[4] = state[5]
    tempState[5] = state[4]

    if tempState not in X:
      X.append(list(tempState))
      Y.append(Y[j])
    else:
      index = X.index(tempState)
      Ynew = min([Y[index],Y[j]])
      Y[j] = Ynew
      Y[index] = Ynew
    
  return [X,Y]

def getState(occ):
  ## Gets the encoded state corresponding to a given configuration
  state = [varO,varO,varO,varO,varO,varO]
  sites = ['3A','3B','5','7','9A','9B']
  for site in occ:
    if 'OH' in site:
      index = sites.index(site.replace('OH',''))
      state[index] = varOH
    elif 'V' in site:
      index = sites.index(site.replace('V',''))
      state[index] = varVac
    elif 'H2O' in site:
      index = sites.index(site.replace('H2O',''))
      state[index] = varH2O
  return state

def getPossibleStates():
  ## Generates all possible states
  bin = [varO,varOH]
  bin2 = [varO,varOH,varH2O,varVac]
  Xconf = [[a,b,c,d,e,f] for a in bin for b in bin for c in bin2 for d in bin for e in bin for f in bin]
  
  return Xconf

def removeDoneCalc(X,XD):
  ## Removes states in the training set from the list of states investigated by the GP routine
  XR = []
  for state in X:
    status = 1
    for stateD in XD:
      if state == stateD:
        status = 0
        break
    if status:
      XR.append(state)
  return XR

def getDoneCalc(dbs):
  ## Gets the states of the training data

  XD = [] # Contains all the surface configurations contained in the training data
  # Loops over databases containing training data
  for db in dbs:
    for row in db.select():
      if row.relaxed:
        conf = row.data.OCCUPANCY.replace('\n','').split(',') #Retrives the configuration from the occupancy flag stored in the database metadata
        XD.append(getState(conf))
  
  for j in range(0,len(XD)):
    state = XD[j]
    tempState = list(state)
    tempState[0] = state[1]
    tempState[1] = state[0]
    tempState[4] = state[5]
    tempState[5] = state[4]

    if tempState not in XD:
      XD.append(list(tempState))
  
  return XD


def getTrainingData(dbs,R,material):
  # Generates training data. Finds adsorption energies of surfaces calculated from VASP.
  
  # Variable definitions:
  X_full = [] # Stores the states of all the training data
  Y_full = [] # Stores the VASP formation energy for all training data
  dG_Ad_full = [] # Same as Y_full
  conf_full = [] # Contains all the configurations contained in the training data

  # Loops over databases containing training data
  for db in dbs:
    for row in db.select():
        if row.relaxed and row.stable:
          #Retrives the configuration from the occupancy flag stored in the database metadata
          conf = row.data.OCCUPANCY.replace('\n','').split(',') 
          # Retrive reference slab       
          encut = row.calculator_parameters['encut']
          refID = R['Name'].index(material+'-'+str(encut))
          refConf = R['Occupancy'][refID].split(',')
          refEnergy = R['Energy'][refID]
         
          dG_Ad = calcAdsorptionEnergy(conf,row.energy,refConf,refEnergy) # Calculates the adsorption energy
          conf_full.append(conf)
          dG_Ad_full.append(dG_Ad)
          Y_full.append(dG_Ad) 
          X_full.append(getState(conf))


  return [conf_full,dG_Ad_full,Y_full,X_full]


def predictEnergy(state,model,v):
  # Predicts the energy of a given state using the 
  X = np.asarray(state)
  
  
  index = model[0].index(state)
  U = 0
  for site in state:
    if site==varO:
      U = U - 2*v
    elif site==varOH:
      U = U - v
  
  return [model[1][index] + U, model[2][index]]



