from ase.io import read,write
from ase import Atoms
from ase.build import add_adsorbate
import subprocess
import os
from os.path import exists
import operator


def generateCoverage(conf,G,id,config):
  
  occ = str(conf).replace('[','').replace(']','').replace(' ','').replace("'",'')
  
# Retrive paths from config file
  f = open(config)
  lines = f.readlines()
  for line in lines:
    cutLine = line.replace('\n','')
    
    if 'JOBPATH' in cutLine:
      flag = cutLine.split('=')
      jobPath = flag[1] + occ
    elif 'SLABPATH' in cutLine:
      flag = cutLine.split('=')
      slabPath = flag[1]
    elif 'POSPATH' in cutLine:
      flag = cutLine.split('=')
      posPath = flag[1]
    elif 'ISIF' in cutLine:
      isif = cutLine
    elif 'ENCUT' in cutLine:
      encut = cutLine
    elif 'KPTS' in cutLine:
      kpts = cutLine
    elif 'MO' in cutLine:
      mo = cutLine

  f.close()
  if not exists(jobPath):
    os.makedirs(jobPath)

  # Read positions
  f = open(posPath)
  lines = f.readlines()
  sites = {}
  for line in lines:
    flag = line.split(':')
    position = flag[1].replace('\n','').replace('(','').replace(')','').split(',')
    for j in range(0,len(position)):
        position[j] = float(position[j])

    sites[flag[0]]=tuple(position)
  

  atoms = read(slabPath)

  for site in conf:
    siteIdentifier = site.replace('O','').replace('H','').replace('V','')
    pos = sites[siteIdentifier]
    Hpos = tuple(map(operator.add, pos, (0.3,0.3,0.86)))
    if 'V' not in site:
      atoms.append('O')
      atoms.positions[-1]=pos
    if 'OH' in site:
      atoms.append('H')
      atoms.positions[-1] = Hpos
  
  
  write(jobPath+'/POSCAR',atoms)

  fileName = jobPath+'/config'
  File_object = open(fileName,"w")
  File_object.write('NAME:'+occ+'\n')
  File_object.write('OCCUPANCY:'+occ+'\n')
  File_object.write('dG_T:'+str(G)+'\n')
  File_object.write('ID:'+str(id)+'\n')

  File_object.write(isif+'\n')
  File_object.write(encut+'\n')
  File_object.write(kpts+'\n')
  File_object.write(mo+'\n')
  File_object.close()

  subprocess.check_call(['/home/cat/s144109/dataProcessing/scripts/createScripts.sh', occ, config])


