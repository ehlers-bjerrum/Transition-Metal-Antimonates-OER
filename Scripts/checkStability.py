from ase.db import connect
import numpy as np

def checkStability(dbs):
  i=0
  db = connect(dbs[i])
  
  pos3A = np.array([-1.6,0])
  pos3B = np.array([-4.75,0])
  pos5 = np.array([-8.0,0])
  pos7 = np.array([-7.960,3.4])
  pos9A = np.array([-4.6,3.4])
  pos9B = np.array([-1.7,3.4])
  
  positions = [pos3A,pos3B,pos5,pos7,pos9A,pos9B]

  for row in db.select():
    adsorbates = ['O3A','O3B','O5','O7','O9A','O9B']
    conf = row.data.OCCUPANCY.split(',')
    for j in range(0,len(conf)):
      site = conf[j]
      if 'V' in site:
        adsorbates[j] = adsorbates[j].replace('O','V')

    for j in range(0,len(row.symbols)):
      if row.symbols[j] == 'H':
        Hpos = row.positions[j][0:2]
        for i in range(0,len(positions)):
          ads = adsorbates[i]
          pos = positions[i]
          
          dist = []
          dist.append(np.linalg.norm(pos-Hpos))
          dist.append(np.linalg.norm(pos-Hpos+row.cell[0,0:2]))
          dist.append(np.linalg.norm(pos-Hpos-row.cell[0,0:2]))
          dist.append(np.linalg.norm(pos-Hpos+row.cell[1,0:2]))
          dist.append(np.linalg.norm(pos-Hpos-row.cell[1,0:2]))
         
          for d in dist:
            if d < 1.5:
              if 'OH' in ads:
                adsorbates[i] = adsorbates[i].replace('OH','H2O')
              elif 'V' in ads:
                adsorbates[i] = adsorbates[i].replace('V','OH')
              else:
                adsorbates[i] = adsorbates[i].replace('O','OH')
              break


    if conf == adsorbates:
      
      db.update(row.id,stable=True)
    else:
      print(adsorbates)
      print(conf)

      db.update(row.id,stable=False)

