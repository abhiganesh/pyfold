#!/bin/bash                                       
echo "starting cns.."                           
touch iam.running                                 
# CNS-CONFIGURATION                               
source /data/casp13/CONFOLD2/cns_solve_1.3/cns_solve_env.sh                
export KMP_AFFINITY=none                          
export CNS_CUSTOMMODULE=/home/abhi/test7/myout/                 
/data/casp13/CONFOLD2/cns_solve_1.3/intel-x86_64bit-linux/bin/cns_solve < dgsa.inp > dgsa.log 
if [ -f "1a3aA.dist_20.pdb"  ]; then           
   rm iam.running                                 
   echo "trial structures written."             
   rm *embed*                                     
   exit                                           
fi                                                
if [ -f "1a3aA.dista_20.pdb" ]; then 
   rm iam.running                                 
   echo "accepted structures written."          
   rm *embed*                                     
   exit                                           
fi                                                
tail -n 30 dgsa.log                              
echo "ERROR! Final structures not found!"       
echo "CNS FAILED!"                              
mv iam.running iam.failed                         
