#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 12:28:52 2019

@author: Anders
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
r = 0.5
varO = 0
varOH = 1 + r*1j
varH2O = 2 + 2*r*1j
varVac = 3


def GPpredictor(voltage,state=False,conf=False):
  ########## SETUP ##########
  material = 'Mn'

  ## Connect to databases
  db1 = connect('/home/cat/s144109/dataProcessing/MnSb2O6/databases/GPTraining2.db')
  dbs = [db1]


  ######## PREDICTIONS ############
  
  model = trainModel(dbs,material='Mn')
  
  # Generate possible configurations that have not already been calculated
  X = getPossibleStates()

  if conf is not False:
    return predictEnergy(getState(conf.split(',')),model,voltage)
  
  if state is not False:
    return predictEnergy(state,model,voltage)
  
  XD = getDoneCalc(dbs)
  XR = removeDoneCalc(X,XD)

  # Minimum energy of training set
  [candidates,dG_T,S] = predictCandidates(model,XR,dbs,material,voltage)
  [[candidates,S],dG_T] = sortCandidates([candidates,S],dG_T)
  
  return [candidates,dG_T,S]

  

def getBestActivity(dbs,material,voltage,model):
  # Predict maximum from training data
  [OS,I,R] = loadDatabases()
  [conf_full,dG_Ad_full,Y_full,X_full] = getTrainingData(dbs,R,material)
  
    
  activity = []
  activityError = []
  site = []
  for state in X_full:
    [A,Asig] = predictActivity(state,model,voltage)
    minIndex = A.index(min(A))
    activity.append(A[minIndex])
    activityError.append(Asig[minIndex])
    site.append(minIndex)
  

  [[VASPState,VASPError,site],VASPActivity] = sortCandidates([conf_full,activityError,site],activity)
  print('Most active surfaces')
  for j in range(0,20):
    print(VASPState[j])
    print(VASPActivity[j])
    print(VASPError[j])
    print(site[j])

  bestActivity = VASPActivity[0]
  error = VASPError[0]

  return [bestActivity,error]  


def predictDone(state,dbs):
  XD = getDoneCalc(dbs)
  Y = []
  for j in range(0,len(state)):
    tempState1 = list(state)
    tempState2 = list(state)
    tempState1[j] = varO
    tempState2[j] = varOH

    if tempState1 in XD and tempState2 in XD:
      Y.append(1)
    else:
      Y.append(0)
  
  return Y


def predictUnstable(state,dbs):
  XD = getUnstableCalc(dbs)
  Y = []
  for j in range(0,len(state)):
    tempState1 = list(state)
    tempState2 = list(state)
    tempState1[j] = varO
    tempState2[j] = varOH

    if tempState1 in XD or tempState2 in XD:
      Y.append(1)
    else:
      Y.append(0)

  return Y



def predictActivity(state,model,voltage):
  # get the activity of the site 
  getBarrier = lambda v,d: -0.58*v-0.48*d+2.41

  [G,sig] = predictEnergy(state,model,voltage)
  A = []
  Asig = []
  for j in range(0,len(state)):
    tempState = list(state)
      
    if tempState[j] == varOH:
      A.append(0)
      Asig.append(0.1)
      
      #tempState[j] = varO
      #[GTemp,sigTemp] = predictEnergy(tempState,model,voltage)   
      #desc = GTemp-G + voltage
      #barrier = getBarrier(voltage,desc)
      #activity = GTemp + barrier
      #error = np.sqrt(0.48**2*(sigTemp**2 + sig**2)+sigTemp**2)
      
      #A.append(activity)
      #Asig.append(error)
    
    elif tempState[j] == varO:
      tempState[j] = varOH
      [GTemp,sigTemp] = predictEnergy(tempState,model,voltage)
      desc = G-GTemp + voltage
      barrier = getBarrier(voltage,desc)
      activity = G + barrier
      error = np.sqrt(0.48**2*(sigTemp**2 + sig**2)+sigTemp**2)

      A.append(activity)
      Asig.append(error)
    
    elif tempState[j] == varH2O or tempState[j] == varVac:
      A.append(0)
      Asig.append(0.1)

      #tempState2 = list(state)
      #tempState[j] = varO
      #[GTemp1,sigTemp1] = predictEnergy(tempState,model,voltage)
      #tempState2[j] = varOH
      #[GTemp2,sigTemp2] = predictEnergy(tempState2,model,voltage)
      
      #desc = GTemp1 - GTemp2 + voltage
      #barrier = getBarrier(voltage,desc)
      #activity = GTemp1 + barrier
      #error = np.sqrt(0.48**2*(sigTemp1**2 + sigTemp2**2)+sigTemp1**2)

      #A.append(activity)
      #Asig.append(error)
    
  return [A,Asig]

def predictCandidates(model,XR,dbs,material,voltage):

  minEnergyVASP = getMinimumEnergy(dbs,material,voltage)  
  [BestActivityVASP,ActivityErrorVASP] = getBestActivity(dbs,material,voltage,model)

  candidates = []
  dG_T = []
  S = []
  for cand in XR:
    [dG,sig] = predictEnergy(cand,model,voltage)
    [A,Asig] = predictActivity(cand,model,voltage)
    if (dG-minEnergyVASP) < sig*1.95:
      candidates.append(getConf(cand))
      dG_T.append(dG)
      S.append(sig)
      continue
    for j in range(0,len(A)):
      if (A[j]-BestActivityVASP) < np.sqrt(Asig[j]**2 + ActivityErrorVASP**2):
        candidates.append(getConf(cand))
        dG_T.append(dG)
        S.append(sig)
        break

  return [candidates,dG_T,S]

def getConf(state):
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

  # Get minimum energy from training data
  [OS,I,R] = loadDatabases()
  [conf_full,dG_Ad_full,Y_full,X_full] = getTrainingData(dbs,R,material)  
  
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
  
  [[VASPC],VASPE] = sortCandidates([conf_full],energyAtRHE_VASP)
  print('Thermodynamic hull')
  for j in range(0,20):
    print(VASPC[j])
    print(VASPE[j])

  return minEnergyVASP

def sortCandidates(toSort,ref):
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

def nll_fn(X_train, Y_train, naive=False):
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
  def nll_naive(theta):
    # Naive implementation of Eq. (7). Works well for the examples 
    # in this article but is numerically less stable compared to 
    # the implementation in nll_stable below.
    K = kernel(X_train, X_train, l=theta[0], sigma_f=theta[1]) + noise**2*np.eye(len(Y_train))
    return 0.5 * np.log(det(K)) + 0.5 * Y_train.T.dot(inv(K).dot(Y_train)) + 0.5 * len(X_train) * np.log(2*np.pi)

  def nll_stable(theta):
    # Numerically more stable implementation of Eq. (7) as described
    # in http://www.gaussianprocess.org/gpml/chapters/RW2.pdf, Section
    # 2.2, Algorithm 2.1.
    K = kernel(X_train, X_train, l=theta[0], sigma_f=theta[1]) + theta[2]**2*np.eye(len(Y_train))
    
    L = cholesky(K)
    return np.sum(np.log(np.diagonal(L))) + 0.5 * Y_train.T.dot(lstsq(L.T, lstsq(L, Y_train,rcond=None)[0],rcond=None)[0]) + 0.5 * len(X_train) * np.log(2*np.pi)
    
  if naive:
    return nll_naive
  else:
    return nll_stable



def trainModel(dbs,material='Mn',leaveOneOutState=False):
  
  [OS,I,R] = loadDatabases()
  X = getPossibleStates()
  [conf_train,dG_Ad_train,Y_train,X_train] = getTrainingData(dbs,R,material)

  for i in range(0,len(X_train)):
    s1 = X_train[i]
    X_train_temp = list(X_train)
    X_train_temp.pop(i)
    for j in range(0,len(X_train_temp)):
      s2 = X_train_temp[j]
      if s1==s2:
        print(s1)
        print(i)
        print(s2)
        print(j)

        print('duplicate')

  [X_train,Y_train] = enforceSymmetry(X_train,Y_train)
  
  Y_mean = np.mean(Y_train)
  Y_std = np.std(Y_train)
  Y_train = (Y_train - Y_mean)/Y_std

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
  bin = [varO,varOH]
  bin2 = [varO,varOH,varH2O,varVac]
  Xconf = [[a,b,c,d,e,f] for a in bin for b in bin for c in bin2 for d in bin for e in bin for f in bin]
  
  return Xconf

def removeDoneCalc(X,XD):
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
  XD = [] # Contains all the surface configurations contained in the training data
  # Loops over databases containing training data
  for db in dbs:
    for row in db.select():
      if row.relaxed and row.stable:
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

def getUnstableCalc(dbs):
  XD = [] # Contains all the surface configurations contained in the training data
  # Loops over databases containing training data
  for db in dbs:
    for row in db.select():
      if row.stable is False:
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
  # Generates training data. Finds adsorption energies of surfaces calculated from VASP and subtracts the result of a tight binding model, using on-site energies and and pair interaction energies (optional)
  
  # Variable definitions:
  X_full = [] # Stores the mean field configuration of all the training data
  Y_full = [] # Stores the VASP subtracted TB model energy for all training data
  dG_Ad_full = [] # Contains the formation energy of all surfaces
  conf_full = [] # Contains all the surface configurations contained in the training data

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
          Y_full.append(dG_Ad) # The training data should capture higher order interactions, so we remove the tight binding model
          X_full.append(getState(conf))


  return [conf_full,dG_Ad_full,Y_full,X_full]


def predictEnergy(state,model,v):
  X = np.asarray(state)
  
  
  index = model[0].index(state)
  U = 0
  for site in state:
    if site==varO:
      U = U - 2*v
    elif site==varOH:
      U = U - v
  
  return [model[1][index] + U, model[2][index]]



