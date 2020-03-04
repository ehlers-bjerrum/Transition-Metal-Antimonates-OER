
#!/usr/bin/python

import os,sys,subprocess
from ase.io import read
from ase.calculators.vasp import Vasp
import ase.calculators.vasp as vasp_calculator
import numpy as np
from ase.optimize import BFGS
from pymatgen import MPRester
from ase.db import connect
from os.path import exists
import os
import subprocess
from ase.calculators.vasp import Vasp2

def assignMagMom(atoms):
  ### Asign magnetic moments ###
  sign = 1 # Alternating sign
  for i in range(len(atoms)):
    if atoms[i].symbol in elements:
      for j in range(len(elements)):
        if atoms[i].symbol == elements[j]:
          if MO == 'AFM':
            if atoms[i].symbol in magMetals:
              atoms[i].magmom = sign*magmom[j]
              sign = -sign
          else:
            atoms[i].magmom = magmom[j]
    else:
      atoms[i].magmom = 0
  return atoms


### BASIC SETUP
ntasks=40
DATABASENAME = DBtemp
ISIF = ISIFtemp
ENCUT = ENCUTtemp
KPTS = KPTStemp
MO = MOtemp # Magnetic order
ID = IDtemp

magMetals = ['Ni','Mn'] # Magnetic metals
elements = ['Ni','Mn','Sb','O','H']
magmom = [5,8,0.6,0.6,0.6]


### Database setup
db = connect(DATABASENAME)
name=0
if exists('config'):
  f = open("config", "r")
  for line in f:
    flag=line.split(':')
    if flag[0]=='NAME':
      name=flag[1].replace('\n','')
      f.close()
      break
    f.close()
if name==0:
  name=os.getcwd()


### Calculation parameters ###

calc = Vasp2(
  istart = 0,
  ismear = 0,
  ibrion = 2,
  algo = 'normal',	
  ediff = 0.0001,
  ediffg = -0.03,
  isif = ISIF,
  lorbit = 11,
  lreal = 'auto',
  nelm = 100,
  nsw = 500,
  prec = 'normal',
  sigma = 0.05,
  kpar = 2,        
  nelmin = 4,
  ncore=ntasks,
  xc = 'rpbe',
  gga = 'RP',
  encut = ENCUT,
  ispin = 1,
  kpts = KPTS,
  gamma = True,
  lasph=True,
  setups='materialsproject',
  dipol = (0.5,0.5,0.5),
  icharg = 2,
  idipol = 3,
  ldipol = True,
  lcharg = True,
  #lvtot = True,
  lwave = True,
  ldau = True,
  ldautype = 2,
  ldau_luj = {'Ni': {'L':2,'U':6.2,'J':0},'Mn': {'L':2,'U':3.9,'J':0}},
  lmaxmix = 4,
  lvhar = True,
                   )
### Start job

# First relaxation to get the correct structure
subprocess.call(["cp", "POSCAR","POSCARBCKP"])
atoms = read('POSCAR')
atoms.set_calculator(calc)
potentialenergy=atoms.get_potential_energy()

# Second relaxation with spin
calc.set(ispin=2,istart=1)
catoms = read('CONTCAR')
catoms = assignMagMom(catoms)
catoms.set_calculator(calc)
potentialenergy=catoms.get_potential_energy()
calc.write_json('OUTPUT.json')

if ID == -1:
  id = db.write(catoms, name=name,path=os.getcwd(),relaxed=False,stable=True)
else:
  db.write(catoms,id=int(ID),name=name,path=os.getcwd(),relaxed=False,stable=True)
  id = int(ID)
if exists('config'):
  f = open("config", "r")
  for line in f:
    flag=line.split(':')
    KEY = flag[0].replace('\n','')
    VALUE = flag[1].replace('\n','')
    db.update(id,data={KEY:VALUE})
  f.close()
f = open("OUTCAR","r")
data = f.readlines()
f.close()
convFlag = 'General timing and accounting informations for this job:'
if len(data)>14:
  if convFlag in data[-14]:
    db.update(id,relaxed=True)
    subprocess.call(["rm", "WAVECAR"])
    subprocess.call(["rm", "CHGCAR"])
    subprocess.call(["rm", "CHG"])
  else:
    db.update(id,relaxed=False)
else:
  db.update(id,relaxed=False)


