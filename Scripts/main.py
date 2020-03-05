from generateCoverage import generateCoverage
from searchAlgorithm import GPpredictor
from ase.db import connect

def main():
  config='/home/cat/s144109/dataProcessing/NiSb2O6/surfacePourbaix/configuration'
  f = open(config)
  lines = f.readlines()
  
  terminate=True
  for line in lines:
    if 'RUN:TRUE' in line:
      terminate=False
      break 
  if terminate is True:
    print('Search terminated')
    f.close()
    return

  for line in lines:  
    if 'DBPATH' in line:
      flag = line.split('=')
      dbPath = flag[1].replace('\n','')

  f.close()
  [candidates,dG_T,S]=GPpredictor([dbPath],'Ni',2)

  if candidates:
    for rank in range(0,len(candidates)):
      db = connect(dbPath)
      conf = candidates[rank]
      print(conf)
      name = str(conf).replace('[','').replace(']','').replace(' ','').replace("'",'')
      print(name)
      id = db.reserve(name=name)
      
      if id is not None:
        generateCoverage(candidates[rank],dG_T[rank],id,config)
        db.update(id=id,relaxed=False)
        break
