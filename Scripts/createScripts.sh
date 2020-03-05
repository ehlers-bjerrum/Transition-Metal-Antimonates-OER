#!/bin/bash

run=$(grep TARGETPATH= $2)
wd=${run#"TARGETPATH="}

run=$(grep JOBSCRIPTS= $2)
jobScripts=${run#"JOBSCRIPTS="}

run=$(grep WORKFOLDER= $2)
home=${run#"WORKFOLDER="}

run=$(grep DBPATH= $2)
dbName=${run#"DBPATH="}
dbName="'"$dbName"'"
rep1="/"
rep2="\/"
dbName="${dbName//$rep1/$rep2}"

job=$home"jobs/"$1"/"

cd $job
run=$(grep NAME: config)
name=${run#"NAME:"}
run=$(grep ISIF: config)
isif=${run#"ISIF:"}
run=$(grep ENCUT: config)
encut=${run#"ENCUT:"}
run=$(grep KPTS: config)
kpts=${run#"KPTS:"}
run=$(grep ID: config)
id=${run#"ID:"}
run=$(grep MO: config)
mo=${run#"MO:"}
mo="'"$mo"'"

folder=$wd$name

mkdir $folder

cp $jobScripts"run_vasp.py" $folder"/run_vasp.py"
cp $jobScripts"submit_script" $folder"/submit_script"
cp $job"POSCAR" $folder"/POSCAR"
cp $job"config" $folder"/config"


cd $folder
sed -i "s/DBtemp/${dbName}/g" "run_vasp.py"
sed -i "s/temp/${name}/g" "submit_script"
sed -i "s/ISIFtemp/${isif}/g" "run_vasp.py"
sed -i "s/ENCUTtemp/${encut}/g" "run_vasp.py"
sed -i "s/KPTStemp/${kpts}/g" "run_vasp.py"
sed -i "s/MOtemp/${mo}/g" "run_vasp.py"
sed -i "s/IDtemp/${id}/g" "run_vasp.py"

echo "import sys" >> "run_vasp.py"
echo "sys.path.insert(1,'$home')" >> "run_vasp.py"
echo "from main import main" >> "run_vasp.py"
echo "main()" >> "run_vasp.py"


sbatch submit_script

#rm -r $job
cd $home
