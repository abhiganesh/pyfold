# PyFold
> A python package for protein folding and molecular visualization.

Installation tested in a Ubuntu16 system.

## Install CNS Suite
#### 1. Download CNS suite   
Provide your academic profile related information at http://cns-online.org/cns_request/. An email with (a) link to download, (b) login, and (c) password will be sent to you. Follow the link, possibly http://cns-online.org/download/, and download CNS suite "cns_solve_1.3_all_intel-mac_linux.tar.gz".
#### 2. Unzip  
```bash
tar xzvf cns_solve_1.3_all_intel-mac_linux.tar.gz
```
#### 3. Change directory to cns_solve  
```bash
cd cns_solve_1.3
```
#### 4. Unhide the file '.cns_solve_env_sh'  
```bash
mv .cns_solve_env_sh cns_solve_env.sh
```
#### 5. Edit path in 'cns_solve_env.sh'  
Replace '_CNSsolve_location_' with CNS installation directory. For instance, if your CNS installation path is '/home/user/programs/cns_solve_1.3' replace '_CNSsolve_location_' with this path
```bash
vim cns_solve_env.sh
```
#### 6. Test CNS installation  
```bash
source cns_solve_env.sh
cd test 
../bin/run_tests -tidy *.inp
```
#### 7. Change directory  
```bash
cd ../../
```

## Configure and test 'pyfold.py' -- test steps to be added in future release 
#### 1. Update paths for the variables '$cns_suite' and '$program_dssp'  
```bash
vim pyfold.py
```
```python
program_dssp   = "$DIR_BASE/dssp-2.0.4-linux-amd64"
cns_suite      = "/home/badri/DISTFOLD/cns_solve_1.3"
```
#### 2. Test  
##### Example with distances as input (1a3aA)  
a. Build models  
```bash
cd test
rm -r output-1a3aA
python ../pyfold.py -rr ./1a3aA.dist.rr -ss ./1a3aA.ss_sa -o ./output-1a3aA -mcount 20 -selectrr 1.0L
```
b. Evaluate the models
```bash
python ./eval-using-tmscore.py 1guu.pdb output-1guu/ all header
```
##### Example with contacts as input (1guu)  
a. Build models
```bash
rm -r output-1guu
python ../pyfold.py -rr ./1guu.rr -ss ./1guu.ss -o ./output-1guu -mcount 20 -selectrr 1.0L
```
b. Evaluate the models
```bash
python ./eval-using-tmscore.py 1a3aA.pdb output-1a3aA all header
```
