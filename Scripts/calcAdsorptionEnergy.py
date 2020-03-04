from ase.db import connect
import sys
import constants as c
from loadDatabases import loadDatabases

def getFormationEnergy(energy,conf):

  H2O = c.G_H2O + c.SV_H2O - c.G_H2O_GAS
  OOH = c.G_OOH_3H + c.SV_OOH - 2*c.G_H2O_GAS #- 2*1.23
  OO = c.G_OO_4H + c.SV_OO - 2*c.G_H2O_GAS # - 2*1.23
  OH = c.G_OH_H + c.SV_OH - c.G_H2O_GAS
  O = c.G_O_H2 + c.SV_O - c.G_H2O_GAS
  H = c.G_H_H + c.SV_H - c.G_H2_GAS
  

  for site in conf:
    if 'H2O' in site:
      energy = energy + H2O
    elif 'OOH' in site:
      energy = energy + OOH
    elif 'OH' in site:
      energy = energy + OH
    elif 'OO' in site:
      energy = energy + OO
    elif 'O' in site:
      energy = energy + O
    elif 'H' in site:
      energy = energy + H
  return energy


def calcAdsorptionEnergy(conf,energy,refConf,refEnergy):

  ### Calculate and return adsorption energies
  energy = getFormationEnergy(energy,conf)
  refEnergy = getFormationEnergy(refEnergy,refConf)


  return energy - refEnergy
