import os
import os.path
from os import path
import sys
import getopt
import argparse
import __main__ as main
import datetime

argv = sys.argv[1:]
AA3TO1 = {"ALA":"A", "ASN":"N", "CYS":"C", "GLN":"Q", "HIS":"H", "LEU":"L", "MET":"M", "PRO":"P", "THR":"T", "TYR":"Y", "ARG":"R", "ASP":"D", "GLU":"E", "GLY":"G","ILE":"I", "LYS":"K", "PHE":"F", "SER":"S", "TRP":"W", "VAL":"V"}
CBATOM = {"A":"cb", "N":"cb", "C":"cb", "Q":"cb", "H":"cb", "L":"cb", "M":"cb", "P":"cb", "T":"cb", "Y":"cb", "R":"cb", "D":"cb", "E":"cb", "G":"ca", "I":"cb", "K":"cb", "F":"cb", "S":"cb", "W":"cb", "V":"cb"}
AA1TO3 = dict([[v,k] for k,v in AA3TO1.items()])
cns_suite = "/data/casp13/CONFOLD2/cns_solve_1.3"
cns_executable = cns_suite+"/intel-x86_64bit-linux/bin/cns_solve"

mcount = "20"

__warningDest = sys.stderr

def warn(msg):
    __warningDest.write(msg + '\n')

#line 1321
def seq_fasta(file_fasta):
	#if fasta file does not exist, sys.stderr.write("\n")\\
    seq = ""
    with open(file_fasta) as myfile:
        if '>' in myfile.read():
            with open(file_fasta, 'r') as fin:
                data = fin.read().splitlines(True)
            with open(file_fasta, 'w') as fout:
                fout.writelines(data[1:])
            f = open(file_fasta, "r")
            contents = f.read()
            seq = contents.rstrip()
            f.close()
        else:
            f = open(file_fasta, "r")
            contents = f.read()
            seq = contents.rstrip()
            f.close()        
    return seq

def write_cns_seq(file_rr, file_cns_seq):
    seq = [char for char in seq_rr(file_rr)]
    three_letter_seq = []
    f = open(file_cns_seq, "a+")
    for word in seq:
        letter2three = AA1TO3[word]
        three_letter_seq.append(letter2three)
    s = " ".join(three_letter_seq)
    chunks, chunk_size = int(len(three_letter_seq)/16) + (len(three_letter_seq) % 16 > 0), 64
    lineOfSeq = [ s[i:i+chunk_size] for i in range(0, len(s), chunk_size) ]
    for line in lineOfSeq:
        f.write(line + "\n")
    f.close()            
#    os.system("rm -f file_cns_seq")
#   while three_letter_seq:
#        if len(three_letter_seq) <= 64:
#            print2file(file_cns_seq, three_letter_seq)
#            three_letter_Seq = None
#        else:
#            print2file(file_cns_seq, three_letter_seq[0:64])
#            three_letter_seq = three_letter_seq[64]

#line 616
def load_hbonds(file_hbonds):
    f = open(str(file_hbonds), "r")
    contents = f.read()
    list = contents.split()
    #????
    f.close()
    return hbond_rows    


# subroutine relax_dihedral
def relax_dihedral(ang, dev):
    dev = round(dev +ang *arelaxation, 2)
    return str(ang) + " " + str(dev)
    #or if integers, return (ang, dev)

#subroutine relax_distance
def relax_distance(d, negdevi, posdevi):
    negdevi = round(negdevi + drelaxation, 2)
    posdevi = round(posdevi + drelaxation, 2)
    return str(d) + " " + str(negdevi) + " " + str(posdevi)
    #or if integers, return (d, negdevi, posdevi)

def print2file(file, message):
    f = open(file, "a+")
    f.write(message + "\n")
    f.close()

def count_lines(myfile):
    num_lines = sum(1 for line in open(myfile))
    return num_lines

def seq_rr(file_rr):
    f = open(file_rr, "r")
    myList = []
    with open(file_rr, 'r+') as f:
        lines = f.readlines()
        for text in lines:
            if text.startswith("PFRMAT"): 
                pass
            elif text.startswith("TARGET"):
                pass
            elif text.startswith("AUTHOR"):
                pass
            elif text.startswith("SCORE"):
                pass
            elif text.startswith("REMARK"):
                pass
            elif text.startswith("METHOD"):
                pass
            elif text.startswith("MODEL"):
                pass
            elif text.startswith("PARENT"):
                pass
            elif text.startswith("TER"):
                pass
            elif text.startswith("END"):
                pass
            elif text.startswith(tuple('0123456789')):
                break
            else:
               myList.append(text)
    finalLine = "".join(myList)
    seq = finalLine.strip() 
    return str(seq)

def rr2tbl(file_rr, file_tbl, rrtype):
    r1r2 = {}
    r1a1 = {}
    numb2letter = {}
    r1List = []
    r2List = []
    a1List = []
    r1cb1r2cb2a1 = []
    seq = seq_rr(file_rr)
    with open(file_rr, 'r+') as f:
        lines = f.readlines()
        lines = [line for line in lines if line[0].isdigit()]
        for text in lines:
                splitLine = text.split()
                r1 = splitLine[0].strip()
                r2 = splitLine[1].strip()
                a1 = splitLine[3].strip()
                r1List.append(r1)
                r2List.append(r2)
                a1List.append(a1)
        for r1 in r1List:
                numb2letter.update({r1 : seq[int(r1)-1]})
        for r2 in r2List:
                numb2letter.update({r2 : seq[int(r2)-1]})
        for text in lines:
                splitLine = text.split()
                r1 = splitLine[0].strip()
                r2 = splitLine[1].strip()
                a1 = float(splitLine[3].strip())
                cb1 = CBATOM[numb2letter[r1]]
                cb2 = CBATOM[numb2letter[r2]]
                r1List.append(r1)                
                r2List.append(r2)
                r1cb1r2cb2a1.append((int(r1),cb1,int(r2),cb2,a1))
                sorted_items = sorted(r1cb1r2cb2a1, key=lambda x: (x[0], x[2]))
        if rrtype == "cb":
                for item in sorted_items:
                        print2file(file_tbl,"assign (resid%4d and name %s) (resid%4d and name %s) 3.6 0.01 %s" % (item[0], item[1], item[2], item[3], item[4] - 3.6))
        else:
                for item in sorted_items:
                        print2file(file_tbl,"assign (resid%4d and name %s) (resid%4d and name %s) 3.6 0.01 %s" % (item[0], "ca", item[2], "ca", item[4] - 3.6))
  

    #print(r1r2)
    #print(r1a1)

#line 418 
def build_models():
    tbl_list = {}
    if path.exists("contact.tbl"):
        tbl_list.update({"contact.tbl": count_lines("contact.tbl")-1})
    if path.exists("ssnoe.tbl"):
        tbl_list.update({"ssnoe.tbl": count_lines("ssnoe.tbl")-1}) 
    if path.exists("hbond.tbl"):
        tbl_list.update({"hbond.tbl": count_lines("hbond.tbl")-1})
    if path.exists("dihedral.tbl"): 
        tbl_list.update({"dihedral.tbl": count_lines("dihedral.tbl")-1})  
    # print "\nWARNING! contact/rr restraints not found or empty!" if not defined $tbl_list{"contact.tbl"};
    # print "\nWARNING! secondary structure restraints not found or empty!" if not defined $tbl_list{"dihedral.tbl"};
    # print "\n".seq_fasta("$id.fasta");
    # print "\n".seq_fasta("$id.ss") if -f "$id.ss";
    # for k,v in  tbl_list:
    # 	flag_dihe = 0
    # 	flag_dihe = 1 if k == "dihedral.tbl"
    # 	print "\n".coverage_tbl("${id}.fasta", $tbl, flag_dihe)
    #     print("\n" 0.00 0.00" % (item[0]
    # }
    # prepare CNS task file
    os.system("rm -f iam.*")
    write_cns_customized_modules()
    write_cns_dgsa_file()
    #os.system("sed -i s/contact.tbl//g dgsa.inp")  
    os.system("sed -i s/hbond.tbl//g dgsa.inp")    
    os.system("sed -i s/dihedral.tbl//g dgsa.inp") 
    os.system("sed -i s/ssnoe.tbl//g dgsa.inp")    




    f = open("job.sh", "w+")
    f.write("#!/bin/bash                                       \n")
    f.write("echo \"starting cns..\"                           \n")
    f.write("touch iam.running                                 \n")
    f.write("# CNS-CONFIGURATION                               \n")
    f.write("source " + cns_suite + "/cns_solve_env.sh                \n")
    f.write("export KMP_AFFINITY=none                          \n")
    f.write("export CNS_CUSTOMMODULE="+ dir_out + "/                 \n")
    f.write(cns_executable + " < dgsa.inp > dgsa.log \n")
    f.write("if [ -f " + "\"" + __id + "_" + mcount + ".pdb" + "\"" + "  ]; then           \n")
    f.write("   rm iam.running                                 \n")
    f.write("   echo \"trial structures written.\"             \n")
    f.write("   rm *embed*                                     \n")
    f.write("   exit                                           \n")
    f.write("fi                                                \n")
    f.write("if [ -f " + "\"" + __id + "a" + "_" + mcount + ".pdb" + "\"" + " ]; then \n")
    f.write("   rm iam.running                                 \n")
    f.write("   echo \"accepted structures written.\"          \n")
    f.write("   rm *embed*                                     \n")
    f.write("   exit                                           \n")
    f.write("fi                                                \n")
    f.write("tail -n 30 dgsa.log                              \n")
    f.write("echo \"ERROR! Final structures not found!\"       \n")
    f.write("echo \"CNS FAILED!\"                              \n")
    f.write("mv iam.running iam.failed                         \n")
    f.close()
    print("\n\n" + "Starting job " + dir_out + "/job.sh > job.log")
    os.system("chmod +x " + dir_out + "/job.sh")
    os.system(dir_out + "/job.sh > job.log")

#line 1472
def write_cns_customized_modules():
    # The scalecoolsetup module edited for scaling implementation
    print2file("scalecoolsetupedited", """! module file: scalecoolsetup
! module file: scalecoolsetup
!
! CNS MODULE
! **********
!
! Authors: Gregory L. Warren and Axel T. Brunger
!
! copyright Yale University
!
! version 2/25/98
!
! Function:
!    This module sets up the NMR restraint data weighting 
!    scheme for either the first or second slow-cooling 
!    dynamics sections
!
! Requirements:
!    See the top level script file for the 
!    defined parameter(s). 
! 
!

module { scalecoolsetup }
( 
   &kdani=kdani;             {INPUT: dani force and class numbers}
   &ksani=ksani;             {INPUT: sani force and class numbers}
   &nmr=nmr;                 {INPUT: nmr restraints parameters}
   &input.noe.scale=input.noe.scale;    {INPUT: noe weight}
   &input.cdih.scale=input.cdih.scale;  {INPUT: cdih weight}
   &input.ncycle=input.ncycle;          {INPUT: number of dynamics cycles}
)

checkversion 1.3

set message ? end
evaluate ($message_old=$result)
set echo ? end
evaluate ($echo_old=$result)
if ( $log_level = verbose ) then
  set echo=on message=normal end
else
  set echo=off message=off end
end if

   noe 
      scale * &input.noe.scale
   end

! start: edited by Badri - 4/22/2014     
    evaluate($newScale = &input.noe.scale * 5)
    noe
        scale N2 $newScale
    end
! end: edited by Badri - 4/22/2014     

   restraints dihedral 
      scale = &input.cdih.scale 
   end

   evaluate ($count = 1)
   while (&exist%nmr.dani.file.$count=true) loop nloop
      if (&nmr.dani.file.$count # "" ) then
         if (&kdani.cart.flag=true) then
            evaluate (&kdani.start=1.0)
            evaluate (&kdani.finish=&nmr.dani.force.finl.$count)
         elseif (&kdani.inter.flag=true) then
            evaluate (&kdani.start=&nmr.dani.force.init.$count)
            evaluate (&kdani.finish=1.0)
         else
            evaluate (&kdani.start=&nmr.dani.force.init.$count)
            evaluate (&kdani.finish=&nmr.dani.force.finl.$count)
         end if
         evaluate (&kdani.$count  = &kdani.start)
         evaluate (&kdani.fac.$count = (&kdani.finish/
                                        &kdani.start)^(1/&input.ncycle))
      end if
      evaluate ($count = $count + 1)
   end loop nloop

   evaluate ($count = 1)
   while (&exist%nmr.sani.file.$count=true) loop nloop
      if (&nmr.sani.file.$count # "" ) then
         if (&ksani.cart.flag=true) then
            evaluate (&ksani.start=1.0)
            evaluate (&ksani.finish=&nmr.sani.force.finl.$count)
         elseif (&ksani.inter.flag=true) then
            evaluate (&ksani.start=&nmr.sani.force.init.$count)
            evaluate (&ksani.finish=1.0)
         else
            evaluate (&ksani.start=&nmr.sani.force.init.$count)
            evaluate (&ksani.finish=&nmr.sani.force.finl.$count)
         end if
         evaluate (&ksani.$count  = &ksani.start)
         evaluate (&ksani.fac.$count = (&ksani.finish/
                                        &ksani.start)^(1/&input.ncycle))
      end if
      evaluate ($count = $count + 1)
   end loop nloop

set message=$message_old echo=$echo_old end
    """)

    print2file("scalehotedited", """
! module file: scalehot
!
! CNS MODULE
! **********
!
! Authors: Gregory L. Warren and Axel T. Brunger
!
! copyright Yale University
!
! version 2/25/98
!
! Function:
!    This module sets the NMR restraint data weight 
!    for high temperature dynamics
!
! Requirements:
!    See the top level script file for the 
!    defined parameter(s) nmr and md. 
! 
!

module { scalehot }
( 
   &md=md;                   {INPUT: md parameters}
   &nmr=nmr;                 {INPUT: nmr restraints parameters}
   &input.noe.scale=input.noe.scale;    {INPUT: noe weight}
   &input.cdih.scale=input.cdih.scale;  {INPUT: cdih weight}
)

checkversion 1.3

set message ? end
evaluate ($message_old=$result)
set echo ? end
evaluate ($echo_old=$result)
if ( $log_level = verbose ) then
  set echo=on message=normal end
else
  set echo=off message=off end
end if

   noe 
      scale * &input.noe.scale
   end
! start: edited by Badri - 4/22/2014     
    evaluate($newScale = &input.noe.scale * 5)
    noe
        scale N2 $newScale
    end
! end: edited by Badri - 4/22/2014  
   restraints dihedral
      scale = &input.cdih.scale
   end
   evaluate ($count = 1)
   while (&exist%nmr.dani.file.$count=true) loop nloop
      evaluate ($clsname = "D"+encode($count))
      if (&nmr.dani.file.$count # "" ) then
         dani
            class $$clsname force &nmr.dani.force.init.$count
         end
      end if
      evaluate ($count = $count + 1)
   end loop nloop

   evaluate ($count = 1)
   while (&exist%nmr.sani.file.$count=true) loop nloop
      evaluate ($clsname = "S"+encode($count))
      if (&nmr.sani.file.$count # "" ) then
         sani
            class $$clsname force &nmr.sani.force.init.$count
         end
      end if
      evaluate ($count = $count + 1)
   end loop nloop

set message=$message_old echo=$echo_old end
    """) 

    #line 1653
def write_cns_dgsa_file():
    atomselect = 2
    #dihed_wt1 = sswt * all_noe_wt
    #dihed_wt2 = sswt * all_noe_wt
    print2file("dgsa.inp", """
{+ file: dgsa.inp +}
{+ directory: nmr_calc +}
{+ description: distance geometry, full or substructure, with 
                simulated annealing regularization starting from 
                extended strand or pre-folded structures. +}
{+ authors: Gregory Warren, Michael Nilges, John Kuszewski, 
	    Marius Clore and Axel Brunger +}
{+ copyright: Yale University +}
{+ reference: Clore GM, Gronenborn AM, Tjandra N, Direct structure refinement 
              against residual dipolar couplings in the presence of rhombicity
              of unknown magnitude., J. Magn. Reson., 131, In press, (1998) +}
{+ reference: Clore GM, Gronenborn AM, Bax A, A robust method for determining 
              the magnitude of the fully asymmetric alignment tensor of
              oriented macromolecules in the absence of structural
              information., J. Magn. Reson., In press (1998) +}
{+ reference: Garrett DS, Kuszewski J, Hancock TJ, Lodi PJ, Vuister GW,
              Gronenborn AM, Clore GM, The impact of direct refinement against 
              three-bond HN-C alpha H coupling constants on protein structure
              determination by NMR., J. Magn. Reson. Ser. B, 104(1), 
              99-103, (1994) May +}
{+ reference: Kuszewski J, Nilges M, Brunger AT,   Sampling and efficiency 
              of metric matrix distance geometry:  A novel partial metrization 
              algorithm.  J. Biomol. NMR 2, 33-56, (1992). +} 
{+ reference: Kuszewski J, Qin J, Gronenborn AM, Clore GM, The impact of direct
              refinement against 13C alpha and 13C beta chemical shifts on 
              protein structure determination by NMR., J. Magn. Reson. Ser. B,
              106(1), 92-6, (1995) Jan +}
{+ reference: Kuszewski J, Gronenborn AM, Clore GM, The impact of direct
              refinement against proton chemical shifts on protein structure 
              determination by NMR., J. Magn. Reson. Ser. B, 107(3), 293-7, 
              (1995) Jun +}
{+ reference: Kuszewski J, Gronenborn AM, Clore GM, A potential involving 
              multiple proton chemical-shift restraints for 
              nonstereospecifically assigned methyl and methylene protons.
              J. Magn. Reson. Ser. B, 112(1), 79-81, (1996) Jul. +}
{+ reference: Nilges M, Clore GM, Gronenborn AM, Determination of 
              three-dimensional structures of proteins from interproton 
              distance data by hybrid distance geometry-dynamical simulated 
              annealing calculations. FEBS Lett. 229, 317-324 (1988). +}
{+ reference: Nilges M, Clore GM, Gronenborn AM,  Determination of 
              three-dimensional structures of proteins from interproton 
              distance data by dynamical simulated annealing from a random 
              array of atoms. FEBS LEtt. 239, 129-136 (1988). +}
{+ reference: Nilges M, Kuszewski J, Brunger AT, In: Computational Aspects 
              of the Study of Biological Macromolecules by NMR. 
              (J.C. Hoch, ed.),  New York: Plenum Press, (1991). +}
{+ reference: Tjandra N, Garrett DS, Gronenborn AM, Bax A, Clore GM, Defining
              long range order in NMR structure determination from the 
              dependence of heteronuclear relaxation times on rotational 
              diffusion anisotropy. Nature Struct. Biol., 4(6), 443-9,
              (1997) June +}
{+ reference: Tjandra N, Omichinski JG, Gronenborn AM, Clore GM, Bax A, Use of
              dipolar 1H-15N and 1H-13C couplings in the structure
              determination of magnetically oriented macromolecules in
              solution. Nature Struct. Biol., 4(9), 732-8, (1997) Sept +} 
              
{- Guidelines for using this file:
   - all strings must be quoted by double-quotes
   - logical variables (true/false) are not quoted
   - do not remove any evaluate statements from the file -}
{- begin block parameter definition -} define(
{======================= molecular structure =========================}
{* parameter file(s) *}
{===>} par.1="CNS_TOPPAR:protein.param";
{===>} par.2="";
{===>} par.3="";
{===>} par.4="";
{===>} par.5="";
{* structure file(s) *}
{===>} struct.1="extended.mtf";
{===>} struct.2="";
{===>} struct.3="";
{===>} struct.4="";
{===>} struct.5="";
{* input coordinate file(s) *}
{===>} pdb.in.file.1="extended.pdb";
{===>} pdb.in.file.2="";
{===>} pdb.in.file.3="";
{========================== atom selection ===========================}
{* input "backbone" selection criteria for average structure generation *}
{* for protein      (name n or name ca or name c)
   for nucleic acid (name O5' or name C5' or name C4' or name C3' 
                     or name O3' or name P) *}
{===>} pdb.atom.select=(name n or name ca or name c);
{======================= refinement parameters ========================}
{* distance geometry *}
{+ choice: true false +}
{===>} flg.dg.flag=true;
{* distance geometry/simualted annealing regularization (DGSA) *}
{+ choice: true false +}
{===>} flg.dgsa.flag=true;
{* if only regularizing coordinate files (no DG) then 
   enter the number of coordinate files to be regularized (DGSA) *}
{===>} pdb.dg.count=20;
{* seed for random number generator *}
{* change to get different initial velocities *}
{===>} md.seed=82364;
{* select whether the number of structures will be either trial or 	
   accepted structures and whether to print only the trial, accepted, 	
   both sets of structures. The printing format is as follows:
   trial = pdb.out.name + _#.pdb , accepted = pdb.out.name + a_#.pdb *} 
{* are the number of structures to be trials or accepted? *}
{+ choice: "trial" "accept" +}
{===>} flg.trial.struc="trial";
{* number of trial or accepted structures *}
{===>} pdb.end.count=20;
{* print accepted structures *}
{+ choice: true false +}
{===>} flg.print.accept=true;
{* print trial structures *}
{+ choice: true false +}
{===>} flg.print.trial=true;
{* calculate an average structure for either the trial or 	
   accepted structure.  If calculate accepted average is false then 
   an average for the trial structures will be calculated. *}
{* calculate an average structure? *}
{+ choice: true false +}
{===>} flg.calc.ave.struct=false;
{* calculate an average structure for the accepted structures? *}
{+ choice: true false +}
{===>} flg.calc.ave.accpt=false;
{* minimize average coordinates? *}
{+ choice: true false +}
{===>} flg.min.ave.coor=false;
{============ parameters for the distance geometry stage ==============}
{* shortest path algorithm *}
{+ choice: "auto" "full" "sparse" +}
{===>} md.dg.algo="auto";
{* distance geometry on substructure or complete structure? *}
{* proteins: "sub" or "complete"; dna/rna: "complete" *}
{+ choice: "sub" "complete" +}
{===>} md.dg.type="sub";
{* input atom selection for substructure  *}""")
    if atomselect == 1:
        # as-is
        print2file("dgsa.inp", """
{===>} md.dg.select=(name ca or name ha or name n or name hn
		        or name c or name cb* or name cg*);""")
        
        # print2file("dgsa.inp", "{===>} md.dg.select=(name ca or name ha or name n or name hn")
        # print2file("dgsa.inp", "		        or name c or name cb* or name cg*);")

    if atomselect == 2:
        # add oxygen
        print2file("dgsa.inp", """
{===>} md.dg.select=(name ca or name ha or name n or name hn
		        or name c or name cb* or name cg* or name o);""")       
        # print2file("dgsa.inp", "{===>} md.dg.select=(name ca or name ha or name n or name hn")
        # print2file("dgsa.inp", "		        or name c or name cb* or name cg* or name o);")

    if atomselect == 3:
        # add oxygen and hydrogen
        print2file("dgsa.inp", """
{===>} md.dg.select=(name ca or name ha or name n or name hn or name h
		        or name c or name cb* or name cg* or name o);""")
        # print2file("dgsa.inp", "{===>} md.dg.select=(name ca or name ha or name n or name hn or name h")
        # print2file("dgsa.inp", "		        or name c or name cb* or name cg* or name o);")

    if atomselect == 4:
        # backbone atoms only
        print2file("dgsa.inp", """
{===>} md.dg.select=(name ca or name c or name n or name o);""")
        #print2file("dgsa.inp", "{===>} md.dg.select=(name ca or name c or name n or name o);")

    if atomselect == 5:
        # backbone atoms with cb
        print2file("dgsa.inp", """
{===>} md.dg.select=(name ca or name c or name n or name o
		        or name cb);""")
        # print2file("dgsa.inp", "{===>} md.dg.select=(name ca or name c or name n or name o")
        # print2file("dgsa.inp", "		        or name cb);")

    if atomselect == 6:
        # backbone atoms with cb and hydrogen
        # atom selection according to instructions in the NIH-XPLORE manual
        print2file("dgsa.inp", """
{===>} md.dg.select=(name ca or name c or name n or name o
		        or name cb or name h);""")
        # print2file("dgsa.inp", "{===>} md.dg.select=(name ca or name c or name n or name o")
        # print2file("dgsa.inp", "		        or name cb or name h);")

    if atomselect == 7:
        # backbone atoms with cb, hydrogen and cg
        print2file("dgsa.inp", """
{===>} md.dg.select=(name ca or name c or name n or name o
		        or name cb or name h or name cg*);""")
        # print2file("dgsa.inp", "{===>} md.dg.select=(name ca or name c or name n or name o")
        # print2file("dgsa.inp", "		        or name cb or name h or name cg*);")

    print2file("dgsa.inp", """
{* when using "complete" input the rigid group atom selection *}
{===>} md.dg.group.slct=(known);
{* group interatomic error value in angstroms *}
{===>} md.dg.group.err=0.5;
{* use metrization for complete distance geometry? *}
{+ choice: true false +}
{===>} md.dg.metr.flag=false;
{* ordered or random metrization *}
{+ choice: "ordered" "random" +}
{===>} md.dg.ord="random";
{* input metrization atom selection *}
{===>} md.dg.metr.atom=(all);
{* input number of atoms from selection used during retightening *}
{===>} md.dg.metr.num=4;
{* reference for building the reference data base *}
{+ choice: "parameter" "coordinate" +}
{===>} md.dg.ref="parameter";
{* scale factor for distance geometry restraint term *}
{===>} md.dg.scale=100.;
{* exponent for distance geometry restraint term *}
{===>} md.dg.expo=2;
{* bond length (in angstroms) error value *}
{===>} md.dg.bacc=0.01;
{* angle (in degrees) error value *}
{===>} md.dg.tacc=2.;
{* improper (in degrees) error value *}
{===>} md.dg.iacc=2.;
{* dihedral (in degrees) error value *}
{===>} md.dg.pacc=2.;
{* number of steps of minimization *}
{===>} md.dg.step=200;
{=== parameters for the distance geometry/simulated annealing stage ===}
{* starting temperature *}
{===>} md.hot.temp=2000;
{* number of steps for high temperature dyanmics *}
{===>} md.hot.step=1000;
{* number of steps for slow-cool annealing *}
{===>} md.cool.step=1000;
{* hot molecular dynamics timestep *}
{===>} md.hot.ss=0.003;
{* slow-cool molecular dynamics timestep *}
{===>} md.cool.ss=0.005;
{* initial scale factor for van der Waals (repel) energy term *}
{===>} md.cool.vdw.init=0.003;
{* final scale factor for van der Waals (repel) energy term *}
{===>} md.cool.vdw.finl=4.0;
{* initial van der Waals repel radius *}
{===>} md.cool.init.rad=1;
{* final van der Waals repel radius *}
{===>} md.cool.fina.rad=0.85;
{* scale factor for NOE energy term *}
{===>} md.cool.noe=10;
{* high temperature scale factor for dihedral angle energy term *}
{===>} md.hot.cdih=5;
{* slow-cooling scale factor for dihedral angle energy term *}
{===>} md.cool.cdih=50;
{* slow-cool annealing temperature step *}
{===>} md.cool.tmpstp=25.;
{=============== parameters for final minimization stage ==============}
{* scale factor for NOE energy term *}
{===>} md.pow.noe=10;
{* scale factor for dihedral angle energy term *}
{===>} md.pow.cdih=50;
{* number of minimization steps *}
{===>} md.pow.step=15000;
{* number of cycles of minimization *}
{===>} md.pow.cycl=10;
      
{============================= noe data ===============================}
{- Important - if you do not have a particular data set then
   set the file name to null ("") -}
{* NOE distance restraints files. *}
{* restraint set 1 file *}
{===>} nmr.noe.file.1="contact.tbl";
{* restraint set 2 file *}
{===>} nmr.noe.file.2="ssnoe.tbl";
{* restraint set 3 file *}
{===>} nmr.noe.file.3="";
{* restraint set 4 file *}
{===>} nmr.noe.file.4="";
{* restraint set 5 file *}
{===>} nmr.noe.file.5="";
{* NOE averaging modes *}
{* restraint set 1 *}
{+ choice: "sum" "cent" "R-6" "R-3" "symm" +}
{===>} nmr.noe.ave.mode.1="cent";
{* restraint set 2 *}
{+ choice: "sum" "cent" "R-6" "R-3" "symm" +}
{===>} nmr.noe.ave.mode.2="sum";
{* restraint set 3 *}
{+ choice: "sum" "cent" "R-6" "R-3" "symm" +}
{===>} nmr.noe.ave.mode.3="R-6";
{* restraint set 4 *}
{+ choice: "sum" "cent" "R-6" "R-3" "symm" +}
{===>} nmr.noe.ave.mode.4="";
{* restraint set 5 *}
{+ choice: "sum" "cent" "R-6" "R-3" "symm" +}
{===>} nmr.noe.ave.mode.5="";
{======================== hydrogen bond data ==========================}
{* hydrogen-bond distance restraints file. *}
{===>} nmr.noe.hbnd.file="hbond.tbl";
{* enter hydrogen-bond distance averaging mode *}
{+ choice: "sum" "cent" "R-6" "R-3" "symm" +}
{===>} nmr.noe.ave.mode.hbnd="cent";
{======================= 3-bond J-coupling data =======================}
{* the default setup is for the phi dihedral *}
{* Class 1 *}
{* 3-bond J-coupling non-glycine restraints file *}
{===>} nmr.jcoup.file.1="";
{* 3-bond J-coupling non-glycine potential *}
{+ choice: "harmonic" "square" "multiple" +}
{===>} nmr.jcoup.pot.1="harmonic";
{* 3-bond J-coupling non-glycine force value *}
{===>} nmr.jcoup.force.1.1=1;
{* 3-bond j-coupling multiple class force second value *}
{===>} nmr.jcoup.force.2.1=0;
{* 3-bond j-coupling Karplus coefficients *}
{* the default values are for phi *}
{===>} nmr.jcoup.coef.1.1=6.98;
{===>} nmr.jcoup.coef.2.1=-1.38;
{===>} nmr.jcoup.coef.3.1=1.72;
{===>} nmr.jcoup.coef.4.1=-60.0;
{* Class 2 *}
{* 3-bond j-coupling glycine restraints files *}
{* The potential for the glycine class must be multiple *}
{===>} nmr.jcoup.file.2="";
{* 3-bond J-coupling non-glycine potential *}
{+ choice: "harmonic" "square" "multiple" +}
{===>} nmr.jcoup.pot.2="multiple";
{* 3-bond J-coupling first force value *}
{===>} nmr.jcoup.force.1.2=1;
{* 3-bond j-coupling glycine or multiple force second value *}
{===>} nmr.jcoup.force.2.2=0;
{* 3-bond j-coupling Karplus coefficients *}
{* the default values are for glycine phi *}
{===>} nmr.jcoup.coef.1.2=6.98;
{===>} nmr.jcoup.coef.2.2=-1.38;
{===>} nmr.jcoup.coef.3.2=1.72;
{===>} nmr.jcoup.coef.4.2=0.0;
{================ 1-bond heteronuclear J-coupling data ================}
{* Class 1 *}
{* 1-bond heteronuclear j-coupling file *}
{===>} nmr.oneb.file.1="";
{* 1-bond heteronuclear j-coupling potential *}
{+ choice: "harmonic" "square" +}
{===>} nmr.oneb.pot.1="harmonic";
{* 1-bond heteronuclear j-coupling force value *}
{===>} nmr.oneb.force.1=1.0;
{=============== alpha/beta carbon chemical shift data ================}
{* Class 1 *}
{* carbon, alpha and beta, chemical shift restraints file *}
{===>} nmr.carb.file.1="";
{* carbon, alpha and beta, chemical shift restraint potential *}
{+ choice: "harmonic" "square" +}
{===>} nmr.carb.pot.1="harmonic";
{* carbon, alpha and beta, chemical shift restraint force value *}
{===>} nmr.carb.force.1=0.5;
{===================== proton chemical shift data =====================}
{* Class 1 *}
{* class 1 proton chemical shift restraints file *}
{===>} nmr.prot.file.1="";
{* class 1 proton chemical shift potential *}
{+ choice: "harmonic" "square" "multiple" +}
{===>} nmr.prot.pot.1="harmonic";
{* class 1 proton chemical shift force value *}
{===>} nmr.prot.force.1.1=7.5;
{* 2nd class 1 proton chemical shift force value for multi *}
{===>} nmr.prot.force.2.1=0;
{* class 1 proton chemical shift violation cutoff threshold *}
{===>} nmr.prot.thresh.1=0.3;
{* Class 2 *}
{* class 2 proton chemical shift restraints file *}
{===>} nmr.prot.file.2="";
{* class 2 proton chemical shift potential *}
{+ choice: "harmonic" "square" "multiple" +}
{===>} nmr.prot.pot.2="harmonic";
{* class 2 proton chemical shift force value *}
{===>} nmr.prot.force.1.2=7.5;
{* 2nd class 2 proton chemical shift force value for multi *}
{===>} nmr.prot.force.2.2=0;
{* class 2 proton chemical shift violation cutoff threshold *}
{===>} nmr.prot.thresh.2=0.3;
{* Class 3 *}
{* class 3 proton chemical shift restraints file *}
{===>} nmr.prot.file.3="";
{* class 3 proton chemical shift potential *}
{+ choice: "harmonic" "square" "multiple" +}
{===>} nmr.prot.pot.3="harmonic";
{* class 3 proton chemical shift force value *}
{===>} nmr.prot.force.1.3=7.5;
{* 2nd class 3 proton chemical shift force value for multi *}
{===>} nmr.prot.force.2.3=0;
{* class 3 proton chemical shift violation cutoff threshold *}
{===>} nmr.prot.thresh.3=0.3;
{* Class 4 *}
{* class 4 proton chemical shift restraints file *}
{===>} nmr.prot.file.4="";
{* class 4 proton chemical shift potential *}
{+ choice: "harmonic" "square" "multiple" +}
{===>} nmr.prot.pot.4="multiple";
{* class 4 proton chemical shift force value *}
{===>} nmr.prot.force.1.4=7.5;
{* 2nd class 4 proton chemical shift force value for multi *}
{===>} nmr.prot.force.2.4=0;
{* class 4 proton chemical shift violation cutoff threshold *}
{===>} nmr.prot.thresh.4=0.3;
{================ diffusion anisotropy restraint data =================}
{* fixed or harmonically restrained external axis *}
{+ choice: "fixed" "harm" +}
{===>} nmr.dani.axis="harm";
{* Class 1 *}
{* diffusion anisotropy restraints file *}
{===>} nmr.dani.file.1="";
{* diffusion anisotropy potential *}
{+ choice: "harmonic" "square" +}
{===>} nmr.dani.pot.1="harmonic";
{* diffusion anisotropy initial force value *}
{===>} nmr.dani.force.init.1=0.01;
{* diffusion anisotropy final force value *}
{===>} nmr.dani.force.finl.1=1.0;
{* diffusion anisotropy coefficients *}
{* coef: <Tc> <anis> <rhombicity> <wh> <wn> *}
{* Tc = 1/2(Dx+Dy+Dz) in <ns> *} 
{===>} nmr.dani.coef.1.1=13.1;
{* anis = Dz/0.5*(Dx+Dy) *} 
{===>} nmr.dani.coef.2.1=2.1;
{* rhombicity = 1.5*(Dy-Dx)/(Dz-0.5*(Dy+Dx)) *} 
{===>} nmr.dani.coef.3.1=0.0;
{* wH in <MHz> *} 
{===>} nmr.dani.coef.4.1=600.13;
{* wN in <MHz> *}
{===>} nmr.dani.coef.5.1=60.82;
{============= susceptability anisotropy restraint data ===============}
{* fixed or harmonically restrained external axis *}
{+ choice: "fixed" "harm" +}
{===>} nmr.sani.axis="harm";
{* Class 1 *}
{* susceptability anisotropy restraints file *}
{===>} nmr.sani.file.1="";
{* susceptability anisotropy potential *}
{+ choice: "harmonic" "square" +}
{===>} nmr.sani.pot.1="harmonic";
{* susceptability anisotropy initial force value *}
{===>} nmr.sani.force.init.1=0.01;
{* susceptability anisotropy final force value *}
{===>} nmr.sani.force.finl.1=50.0;
{* susceptability anisotropy coefficients *}
{* coef: <DFS> <axial > <rhombicity>;
   a0+a1*(3*cos(theta)^2-1)+a2*(3/2)*sin(theta)^2*cos(2*phi) *}
{* DFS = a0 *}
{===>} nmr.sani.coef.1.1=-0.0601;
{* axial = a0-a1-3/2*a2 *}
{===>} nmr.sani.coef.2.1=-8.02;
{* rhombicity = a2/a1 *}
{===>} nmr.sani.coef.3.1=0.4;
{======================== other restraint data ========================}
{* dihedral angle restraints file *}
{* Note: the restraint file MUST NOT contain restraints 
         dihedral or end *}
{===>} nmr.cdih.file="dihedral.tbl";
{* DNA-RNA base planarity restraints file *}
{* Note: include weights as $pscale in the restraint file *}
{===>} nmr.plan.file="";
{* input planarity scale factor - this will be written into $pscale *}
{===>} nmr.plan.scale=150;
{* NCS-restraints file *}
{* example is in inputs/xtal_data/eg1_ncs_restrain.dat *}
{===>} nmr.ncs.file="";
{======================== input/output files ==========================}
{* base name for input coordinate files *}
{* used for simulated annealing when distance geometry is not used *}
{===>} pdb.in.name="dg_sub_embed";
{* base name for output coordinate files *}
{===>} pdb.out.name=\"""" + __id + "\";")

    print2file("dgsa.inp", """
{===========================================================================}
{         things below this line do not normally need to be changed         }
{         except for the torsion angle topology setup if you have           }
{         molecules other than protein or nucleic acid                      }
{===========================================================================}
flg.cv.flag=false;
flg.cv.noe=false;
flg.cv.coup=false;
flg.cv.cdih=false;
nmr.cv.numpart=10;
 ) {- end block parameter definition -}
checkversion 1.3
evaluate ($log_level=quiet)
structure 
   if  (&struct.1 # "") then
      @@&struct.1 
   end if
   if  (&struct.2 # "") then
      @@&struct.2 
   end if
   if  (&struct.3 # "") then
      @@&struct.3 
   end if
   if  (&struct.4 # "") then
      @@&struct.4 
   end if
   if  (&struct.5 # "") then
      @@&struct.5 
   end if
end
if ( &BLANK%pdb.in.file.1 = false ) then
   coor @@&pdb.in.file.1
end if
if ( &BLANK%pdb.in.file.2 = false ) then
   coor @@&pdb.in.file.2
end if
if ( &BLANK%pdb.in.file.3 = false ) then
   coor @@&pdb.in.file.3
end if
parameter
   if (&par.1 # "") then
      @@&par.1
   end if
   if (&par.2 # "") then
      @@&par.2
   end if
   if (&par.3 # "") then
      @@&par.3
   end if
   if (&par.4 # "") then
      @@&par.4
   end if
   if (&par.5 # "") then
      @@&par.5
   end if
end
if ( $log_level = verbose ) then
  set message=normal echo=on end
else
  set message=off echo=off end
end if
parameter
   nbonds
      repel=0.5
      rexp=2 irexp=2 rcon=1.
      nbxmod=-2
      wmin=0.01
      cutnb=4.5 ctonnb=2.99 ctofnb=3.
      tolerance=0.5
   end
end
set seed=&md.seed end
{- Read experimental data -}
   @CNS_NMRMODULE:readdata ( nmr=&nmr;
                             flag=&flg;
                             output=$nmr; )
{- Update the way Contact restraints are handled - Badri Adhikari, 3/4/2016  -}  
noe
  ceiling=1000
  asymptote N1 1.0
  sqconstant N1 1.0
  sqexponent N1 2
  soexponent N1 1
  rswitch N1 0.5
end
{- Read and store the number of NMR restraints -}
   @CNS_NMRMODULE:restraintnumber ( num=$num; )
   
{- Set mass and parameter values -}
   
do (fbeta=10) (all)
do (mass=100) (all)
parameter                  
   nbonds  
      repel=0.80  
      rexp=2 irexp=2 rcon=1. 
      nbxmod=3  
      wmin=0.01  
      cutnb=6.0 ctonnb=2.99 ctofnb=3.  
      tolerance=1.5  
   end  
end
evaluate ($nmr.trial.count = 0)    {- Initialize current structure number   -}
evaluate ($nmr.accept.count = 0)    {- Initialize number accepted            -}
evaluate ($nmr.counter 	= 0)
evaluate ($nmr.prev.counter = -1)
@CNS_NMRMODULE:initave  ( ave=$ave;
                          ave2=$ave2;
                          cv=$cv;
                          ener1=$ener1;
                          ener2=$ener2;
                          flag=&flg;
                          nmr.prot=&nmr.prot; )
        
{- Zero the force constant of disulfide bonds. -}
parameter
   bonds ( name SG ) ( name SG ) 0. TOKEN 
end
{- define a distance restraints for each disulfide bond, i.e., 
   treat it as if it were an NOE and break the bond. -}
for $ss_rm_id_1 in id ( name SG ) loop STRM
  for $ss_rm_id_2 in id ( name SG and 
			  bondedto ( id $ss_rm_id_1 )  ) loop STR2
    if ($ss_rm_id_1 > $ss_rm_id_2) then
      pick bond ( id $ss_rm_id_1 ) ( id $ss_rm_id_2 ) equil
      evaluate ($ss_bond=$result) 
      noe 
         assign ( id $ss_rm_id_1 ) ( id $ss_rm_id_2 ) $ss_bond 0.1 0.1
      end 
    end if
  end loop STR2
end loop STRM
{- Count the number of residues and determine molecule type -}
identify (store9) (tag)
evaluate ($nmr.rsn.num = $SELECT)
identify (store9) ( tag and ( resn THY or resn CYT or resn GUA or
                              resn ADE or resn URI ))
evaluate ($nmr.nucl.num = $SELECT)    
if ( &md.dg.ref = "coordinate" ) then
   flag exclude * include bond angl impr vdw end 
   minimize lbfgs nstep=2000 drop=10.  nprint=1000 end
end if
do (refx=x) ( all )
do (refy=y) ( all )
do (refz=z) ( all )
{- generate and store a bounds matrix -}
if (&flg.dg.flag=true) then
   flags exclude * include bond angle dihedral improper vdw noe cdih end
   mmdg
      shortest-path-algorithm=&&md.dg.algo
      scale=&md.dg.scale
      exponent=&md.dg.expo
      baccuracy=&md.dg.bacc
      taccuracy=&md.dg.tacc
      iaccuracy=&md.dg.iacc
      paccuracy=&md.dg.pacc
      if (&md.dg.type="sub") then
   reference=&&md.dg.ref
   storebounds
      else
   reference=&&md.dg.ref
   group &md.dg.group.slct &md.dg.group.err
   storebounds
      end if
   end
      
   {- Begin protocol to generate structures distance geometry structures -}
   while (&pdb.end.count > $nmr.counter) loop dg
      evaluate ($nmr.counter=$nmr.counter + 1)
      evaluate ($embedded=false)
      flags exclude * include dg end
      if (&md.dg.type="sub") then
   igroup interaction=(&md.dg.select) (&md.dg.select) end
      end if
      coor init end
      while ($embedded = false) loop embed
   mmdg
      if (&md.dg.type="sub") then
         recallbounds
         substructure=(&md.dg.select)
         selection=(&md.dg.select)
      else
         recallbounds
         selection=(all)
         if (&md.dg.metr.flag=true) then
   		  &&md.dg.ord
   		  metrization=(&md.dg.metr.atom)=&md.dg.metr.num
         end if
      end if
   end
      end loop embed
      do (x = x * $dgscale) (known)
      do (y = y * $dgscale) (known)
      do (z = z * $dgscale) (known)
      minimize lbfgs
   nstep=&md.dg.step drop=1. nprint=25
      end
      @CNS_NMRMODULE:printdg ( md=&md;
                               output=$nmr;
                               pdb=&pdb; )
   end loop dg
end if
{- initialize and set scaling factors for simulated annealing -}
set seed=&md.seed end
evaluate ($nmr.trial.count = 0)    {- Initialize current structure number   -}
evaluate ($nmr.dg.count = 0)
evaluate ($nmr.accept.count = 0)   {- Initialize number accepted            -}
evaluate ($nmr.counter = 0)
evaluate ($coor_count_init=0.)
evaluate ($coor_input_count=0.)
if (&flg.dg.flag=true) then
   evaluate ($coor_input_count=&pdb.end.count)
else
   evaluate ($coor_input_count=&pdb.dg.count)
end if
@CNS_NMRMODULE:initave  ( flag=&flg;
                          ave=$ave;
                          ave2=$ave2;
                          cv=$cv;
                          ener1=$ener1;
                          ener2=$ener2;
                          nmr.prot=&nmr.prot; )
        
{- scaling of nmr restraint data during regularization -}
@CNS_CUSTOMMODULE:scalehotedited ( md=&md;
                          nmr=&nmr;
                          input.noe.scale=&md.cool.noe;
                          input.cdih.scale=&md.hot.cdih; )
if (&nmr.dani.axis = "harm") then
   do (harmonic=20.0) (resid 500 and name OO)
   do (harmonic=0.0) (resid 500 and name Z )
   do (harmonic=0.0) (resid 500 and name X )
   do (harmonic=0.0) (resid 500 and name Y )
   do (harmonic=0.0) (not (resid 500))
   restraints harmonic exponent=2 end
elseif (&nmr.sani.axis = "harm") then
   do (harmonic=20.0) (resid 500 and name OO)
   do (harmonic=0.0) (resid 500 and name Z )
   do (harmonic=0.0) (resid 500 and name X )
   do (harmonic=0.0) (resid 500 and name Y )
   do (harmonic=0.0) (not (resid 500))
   restraints harmonic exponent=2 end
end if
{- Increase the disulfide bond force constants to their full strength -}
   parameter
      bonds ( name SG ) ( name SG ) 1000. TOKEN 
   end
{- Regularize structures generated by distance geometry - loop until done -}
if (&flg.dgsa.flag=true) then
   while (&pdb.end.count > $nmr.counter) loop dgsa
      {- Set parameter values -}
      parameter
         nbonds
            repel=0.5
            rexp=2 irexp=2 rcon=1.
            nbxmod=-2
            wmin=0.01
            cutnb=4.5 ctonnb=2.99 ctofnb=3.
            tolerance=0.5
         end
      end
      evaluate ($nmr.trial.count = $nmr.trial.count + 1)
      if ($nmr.trial.count <= $coor_input_count) then
         evaluate ($nmr.dg.count=$nmr.dg.count+1)
         evaluate ($coor_count_init=0.)
      else
         evaluate ($coor_count_init=$coor_count_init+1)
         if ($coor_count_init > $coor_input_count ) then
            evaluate ($coor_count_init=1)
         end if
   evaluate ($nmr.dg.count=$coor_count_init)
      end if
   {- $prefix is generated in the macro printdg -}
      if (&flg.dg.flag=true) then
         evaluate ($filename=$nmr.prefix+encode($nmr.dg.count)+".pdb")
      else
         evaluate ($filename=&pdb.in.name+"_"+encode($nmr.dg.count)+".pdb")
      end if
         
      {- Test for correct enantiomer -}
      for $image in ( 1 -1 ) loop imag
         set remarks=reset end 
   coor initialize end
   coor @@$filename
   do (x=x * $image) ( known )
   identity (store1) (not known)
   coor copy end
   do (x=refx) ( all )
   do (y=refy) ( all )
   do (z=refz) ( all )
   for $id in id ( tag ) loop fit
      coordinates
         fit select = ( byresidue (id $id) and not store1 )
      end
     coor copy selection=( byresidue (id $id) ) end
   end loop fit
   coor swap end
         if (&nmr.dani.axis = "fixed" ) then
            fix
               select=(resname ANI)
            end
         elseif (&nmr.sani.axis = "fixed" ) then
            fix
               select=(resname ANI)
            end
         end if
   parameter
      nbonds
         nbxmod=-2
         repel=0.5
      end
   end
   flags exclude * include bond vdw noe cdih coup oneb 
                           carb ncs dani sani harm end
   igroup interaction (all) (all) weights * 1.  vdw 20. end end
   minimize lbfgs nstep=100 nprint=100 end
   flags include angl end
   minimize lbfgs nstep=100 nprint=100 end
   flags include impr dihe end
   evaluate ($nstep1 = int(&md.hot.step/8))
   evaluate ($nstep2 = int(&md.hot.step/2))
   do ( vx = maxwell(0.5) ) ( all )
   do ( vy = maxwell(0.5) ) ( all )
   do ( vz = maxwell(0.5) ) ( all )
   igroup inter (all) (all) weights * 0.1 impr 0.05 vdw 20. end end
   dynamics cartesian
      cmremove=true
      vscaling=false
      tcoupling=true
      timestep=&md.hot.ss
      nstep=$nstep1
      nprint=$nstep1
      temperature=&md.hot.temp
   end
   igroup inter (all) (all) weights * 0.2 impr 0.1  vdw 20. end end
   dynamics cartesian
      cmremove=true
      vscaling=false
      tcoupling=true
      timestep=&md.hot.ss
      nstep=$nstep1
      nprint=$nstep1
      temperature=&md.hot.temp
   end
   parameter  nbonds repel=0.9   end  end
   igroup inter (all) (all) weights * 0.2 impr 0.2 vdw 0.01 end end
   dynamics cartesian
      cmremove=true
      vscaling=false
      tcoupling=true
      timestep=&md.hot.ss
      nstep=$nstep1
      nprint=$nstep1
      temperature=&md.hot.temp
   end
   parameter nbonds nbxmod=-3  end  end
   igroup inter (all) (all) weights * 0.4 impr 0.4 vdw 0.003 end end
   dynamics cartesian
      cmremove=true
      vscaling=false
      tcoupling=true
      timestep=&md.hot.ss
      nstep=$nstep2
      nprint=$nstep2
      temperature=&md.hot.temp
   end
   igroup inter (all) (all) weights * 1.0 impr 1.0 vdw 0.003 end end
   dynamics cartesian
      cmremove=true
      vscaling=false
      tcoupling=true
      timestep=&md.hot.ss
      nstep=$nstep1
      nprint=$nstep1
      temperature=&md.hot.temp
   end
   if ($image = 1) then
      do (store7=x) ( all )
      do (store8=y) ( all )
      do (store9=z) ( all )
      do (store4=vx) ( all )
      do (store5=vy) ( all )
      do (store6=vz) ( all )
   end if
      end loop imag
      {- Establish the correct handedness of the structure -}
      energy end
      evaluate ($e_minus=$ener)
      coor copy end
      do (x=store7) ( all )
      do (y=store8) ( all )
      do (z=store9) ( all )
      energy end
      evaluate ($e_plus=$ener)
      if ( $e_plus > $e_minus ) then
   evaluate ($hand=-1 )
   coor swap end
      else
   evaluate ($hand= 1 )
   do (vx=store4) ( all )
   do (vy=store5) ( all )
   do (vz=store6) ( all )
      end if
   {- Slow-cooling with cartesian dynamics -}
      parameter
   nbonds
      repel=0.80
      rexp=2 irexp=2 rcon=1.
      nbxmod=3
      wmin=0.01
      cutnb=6.0 ctonnb=2.99 ctofnb=3.
      tolerance=0.5
   end
      end
      flags include plan end
      evaluate ($final_t = 0)
      evaluate ($ncycle = int((&md.hot.temp-$final_t)/&md.cool.tmpstp))
      evaluate ($nstep = int(&md.cool.step/$ncycle))
      evaluate ($vdw_step=(&md.cool.vdw.finl/&md.cool.vdw.init)^(1/$ncycle))
      evaluate ($rad_step=(&md.cool.init.rad-&md.cool.fina.rad)/$ncycle)
      evaluate ($radius=&&md.cool.init.rad)
      {- set up nmr restraint scaling -}
      evaluate ($kdani.inter.flag=false)
      evaluate ($ksani.inter.flag=false)
      evaluate ($kdani.cart.flag=false)
      evaluate ($ksani.cart.flag=false)
      @CNS_CUSTOMMODULE:scalecoolsetupedited ( kdani=$kdani;
                                      ksani=$ksani;
                                      nmr=&nmr;
                                      input.noe.scale=&md.cool.noe;
                                      input.cdih.scale=&md.cool.cdih;
                                      input.ncycle=$ncycle; )
      evaluate ($bath=&md.hot.temp)
      evaluate ($k_vdw=&md.cool.vdw.init)
      evaluate ($i_cool = 0)
      while ($i_cool <= $ncycle) loop cool
   evaluate ($i_cool = $i_cool + 1)
   igroup
      interaction (chemical h*) (all) weights * 1 vdw 0. elec 0. end
      interaction (not chemical h*) (not chemical h*) weights * 1 vdw $k_vdw end
   end
   dynamics  cartesian
      cmremove=true
      vscaling = true
      tcoup = false
      timestep = &md.cool.ss
      nstep = $nstep
      nprint = $nstep
      temperature = $bath
   end
   evaluate ($radius=max(&md.cool.fina.rad,$radius-$rad_step))
   parameter  nbonds repel=$radius   end end
   evaluate ($k_vdw=min(&md.cool.vdw.finl,$k_vdw*$vdw_step))
   evaluate ($bath=$bath-&md.cool.tmpstp)
         @CNS_NMRMODULE:scalecool ( kdani=$kdani;
                                    ksani=$ksani;
                                    nmr=&nmr; )
      end loop cool
   {- Final minimization -}
      {- turn on proton chemical shifts -}
      flags include prot end
      
      if ($nmr.nucl.num > 0) then
         flags include elec end
      end if
      noe             
         scale * &md.pow.noe 
      end
        
      restraints dihedral  
         scale = &md.pow.cdih  
      end
 													
      igroup interaction ( all ) ( all ) weights * 1 end end
      evaluate ($count=0 )
      while (&md.pow.cycl > $count) loop pmini
         evaluate ($count=$count + 1)
         minimize lbfgs nstep=&md.pow.step drop=10.0 nprint=25 end
      end loop pmini
      evaluate ($nmr.min.num = $count * &md.pow.step)
      {- translate the geometric center of the structure to the origin -}
      if ($num.dani > 0. ) then
      elseif ($num.sani > 0. ) then
      else
         show ave ( x ) ( all )
         evaluate ($geom_x=-$result)
         show ave ( y ) ( all )
         evaluate ($geom_y=-$result)
         show ave ( z ) ( all )
         evaluate ($geom_z=-$result)
         coor translate vector=( $geom_x $geom_y $geom_z ) selection=( all ) end
      end if
      
      @CNS_NMRMODULE:printaccept ( ave=$ave;
                                   ave2=$ave2;
                                   cv=$cv;
                                   ener1=$ener1;
                                   ener2=$ener2;
                                   flag=&flg;
                                   md=&md;
                                   nmr=&nmr;
                                   num=$num;
                                   output=$nmr;
                                   pdb=&pdb;  )
   end loop dgsa
   @CNS_NMRMODULE:calcave ( ave=$ave;                 
                            ave2=$ave2;               
                            cv=$cv;                   
                            ener1=$ener1;               
                            ener2=$ener2;             
                            flag=&flg;               
                            md=&md;
                            nmr=&nmr;
                            num=$num;                 
                            output=$nmr;           
                            pdb=&pdb;  )
	
      
      
end if
stop""")

def write_cns_generate_seq_file():
    print2file("gseq.inp", "{+ file: generate_seq.inp +}")
    print2file("gseq.inp", "{+ directory: general +}")
    print2file("gseq.inp", "{+ description: Generate structure file for protein, dna/rna, water, ")
    print2file("gseq.inp", "                ligands and/or carbohydrate from sequence information only +}")
    print2file("gseq.inp", "{+ comment: modified by Brian Smith (Edinburgh University) to allow protein")
    print2file("gseq.inp", "            residue renumbering +}")
    print2file("gseq.inp", "{+ authors: Paul Adams, and Axel Brunger +}")
    print2file("gseq.inp", "{+ copyright: Yale University +}")
    print2file("gseq.inp", "{- Guidelines for using this file:")
    print2file("gseq.inp", "   - all strings must be quoted by double-quotes")
    print2file("gseq.inp", "   - logical variables (true/false) are not quoted")
    print2file("gseq.inp", "   - do not remove any evaluate statements from the file -}")
    print2file("gseq.inp", "{- Special patches will have to be entered manually at the relevant points")
    print2file("gseq.inp", "   in the file - see comments throughout the file -}")
    print2file("gseq.inp", "{- begin block parameter definition -} define(")
    print2file("gseq.inp", "{============ protein topology, linkage, and parameter files =============}")
    print2file("gseq.inp", "{* topology files *}")
    print2file("gseq.inp", "{===>} topology_infile_1=\"CNS_TOPPAR:protein.top\";")
    print2file("gseq.inp", "{===>} topology_infile_2=\"CNS_TOPPAR:dna-rna.top\";")
    print2file("gseq.inp", "{===>} topology_infile_3=\"CNS_TOPPAR:water.top\";")
    print2file("gseq.inp", "{===>} topology_infile_4=\"CNS_TOPPAR:ion.top\";")
    print2file("gseq.inp", "{===>} topology_infile_5=\"CNS_TOPPAR:carbohydrate.top\";")
    print2file("gseq.inp", "{===>} topology_infile_6=\"\";")
    print2file("gseq.inp", "{===>} topology_infile_7=\"\";")
    print2file("gseq.inp", "{===>} topology_infile_8=\"\";")
    print2file("gseq.inp", "{* linkage files for linear, continuous polymers (protein, DNA, RNA) *}")
    print2file("gseq.inp", "{===>} link_infile_1=\"CNS_TOPPAR:protein.link\";")
    print2file("gseq.inp", "{===>} link_infile_2=\"CNS_TOPPAR:dna-rna-pho.link\";")
    print2file("gseq.inp", "{===>} link_infile_3=\"\";")
    print2file("gseq.inp", "{* parameter files *}")
    print2file("gseq.inp", "{===>} parameter_infile_1=\"CNS_TOPPAR:protein.param\";")
    print2file("gseq.inp", "{===>} parameter_infile_2=\"CNS_TOPPAR:dna-rna_rep.param\";")
    print2file("gseq.inp", "{===>} parameter_infile_3=\"CNS_TOPPAR:water_rep.param\";")
    print2file("gseq.inp", "{===>} parameter_infile_4=\"CNS_TOPPAR:ion.param\";")
    print2file("gseq.inp", "{===>} parameter_infile_5=\"CNS_TOPPAR:carbohydrate.param\";")
    print2file("gseq.inp", "{===>} parameter_infile_6=\"\";")
    print2file("gseq.inp", "{===>} parameter_infile_7=\"\";")
    print2file("gseq.inp", "{===>} parameter_infile_8=\"\";")
    print2file("gseq.inp", "{====================== other linkages and modifications  ==================}")
    print2file("gseq.inp", "{* extra linkages and modifications by custom patches *}")
    print2file("gseq.inp", "{===>} patch_infile=\"\";")
    print2file("gseq.inp", "{============================= sequence files ==============================}")
    print2file("gseq.inp", "{* multiple sequence files of the same type can be defined by duplicating")
    print2file("gseq.inp", "   the entries below and incrementing the file number *}")
    print2file("gseq.inp", "{* protein sequence file 1 *}")
    print2file("gseq.inp", "{===>} prot_sequence_infile_1=\"input.seq\";")
    print2file("gseq.inp", "{* segid *}")
    print2file("gseq.inp", "{===>} prot_segid_1=\"\";")
    print2file("gseq.inp", "{* start residue numbering at *}")
    print2file("gseq.inp", "{===>} renumber_1=1;")
    print2file("gseq.inp", "{* protein sequence file 2 *}")
    print2file("gseq.inp", "{===>} prot_sequence_infile_2=\"\";")
    print2file("gseq.inp", "{* segid *}")
    print2file("gseq.inp", "{===>} prot_segid_2=\"\";")
    print2file("gseq.inp", "{* start residue numbering at *}")
    print2file("gseq.inp", "{===>} renumber_2=1;")
    print2file("gseq.inp", "{* protein sequence file 3 *}")
    print2file("gseq.inp", "{===>} prot_sequence_infile_3=\"\";")
    print2file("gseq.inp", "{* segid *}")
    print2file("gseq.inp", "{===>} prot_segid_3=\"\";")
    print2file("gseq.inp", "{* start residue numbering at *}")
    print2file("gseq.inp", "{===>} renumber_3=1;")
    print2file("gseq.inp", "{* protein sequence file 4 *}")
    print2file("gseq.inp", "{===>} prot_sequence_infile_4=\"\";")
    print2file("gseq.inp", "{* segid *}")
    print2file("gseq.inp", "{===>} prot_segid_4=\"\";")
    print2file("gseq.inp", "{* start residue numbering at *}")
    print2file("gseq.inp", "{===>} renumber_4=1;")
    print2file("gseq.inp", "{* nucleic acid sequence file 1 *}")
    print2file("gseq.inp", "{===>} nucl_sequence_infile_1=\"\";")
    print2file("gseq.inp", "{* segid *}")
    print2file("gseq.inp", "{===>} nucl_segid_1=\"\";")
    print2file("gseq.inp", "{* nucleic acid sequence file 2 *}")
    print2file("gseq.inp", "{===>} nucl_sequence_infile_2=\"\";")
    print2file("gseq.inp", "{* segid *}")
    print2file("gseq.inp", "{===>} nucl_segid_2=\"\";")
    print2file("gseq.inp", "{* nucleic acid sequence file 3 *}")
    print2file("gseq.inp", "{===>} nucl_sequence_infile_3=\"\";")
    print2file("gseq.inp", "{* segid *}")
    print2file("gseq.inp", "{===>} nucl_segid_3=\"\";")
    print2file("gseq.inp", "{* nucleic acid sequence file 4 *}")
    print2file("gseq.inp", "{===>} nucl_sequence_infile_4=\"\";")
    print2file("gseq.inp", "{* segid *}")
    print2file("gseq.inp", "{===>} nucl_segid_4=\"\";")
    print2file("gseq.inp", "{* water sequence file 1 *}")
    print2file("gseq.inp", "{===>} water_sequence_infile_1=\"\";")
    print2file("gseq.inp", "{* segid *}")
    print2file("gseq.inp", "{===>} water_segid_1=\"\";")
    print2file("gseq.inp", "{* water sequence file 2 *}")
    print2file("gseq.inp", "{===>} water_sequence_infile_2=\"\";")
    print2file("gseq.inp", "{* segid *}")
    print2file("gseq.inp", "{===>} water_segid_2=\"\";")
    print2file("gseq.inp", "{* carbohydrate sequence file 1 *}")
    print2file("gseq.inp", "{===>} carbo_sequence_infile_1=\"\";")
    print2file("gseq.inp", "{* segid *}")
    print2file("gseq.inp", "{===>} carbo_segid_1=\"\";")
    print2file("gseq.inp", "{* carbohydrate sequence file 2 *}")
    print2file("gseq.inp", "{===>} carbo_sequence_infile_2=\"\";")
    print2file("gseq.inp", "{* segid *}")
    print2file("gseq.inp", "{===>} carbo_segid_2=\"\";")
    print2file("gseq.inp", "{* prosthetic group sequence file 1 *}")
    print2file("gseq.inp", "{===>} prost_sequence_infile_1=\"\";")
    print2file("gseq.inp", "{* segid *}")
    print2file("gseq.inp", "{===>} prost_segid_1=\"\";")
    print2file("gseq.inp", "{* prosthetic group sequence file 2 *}")
    print2file("gseq.inp", "{===>} prost_sequence_infile_2=\"\";")
    print2file("gseq.inp", "{* segid *}")
    print2file("gseq.inp", "{===>} prost_segid_2=\"\";")
    print2file("gseq.inp", "{* ligand sequence file 1 *}")
    print2file("gseq.inp", "{===>} lig_sequence_infile_1=\"\";")
    print2file("gseq.inp", "{* segid *}")
    print2file("gseq.inp", "{===>} lig_segid_1=\"\";")
    print2file("gseq.inp", "{* ligand sequence file 2 *}")
    print2file("gseq.inp", "{===>} lig_sequence_infile_2=\"\";")
    print2file("gseq.inp", "{* segid *}")
    print2file("gseq.inp", "{===>} lig_segid_2=\"\";")
    print2file("gseq.inp", "{* ion sequence file 1 *}")
    print2file("gseq.inp", "{===>} ion_sequence_infile_1=\"\";")
    print2file("gseq.inp", "{* segid *}")
    print2file("gseq.inp", "{===>} ion_segid_1=\"\";")
    print2file("gseq.inp", "{* ion sequence file 2 *}")
    print2file("gseq.inp", "{===>} ion_sequence_infile_2=\"\";")
    print2file("gseq.inp", "{* segid *}")
    print2file("gseq.inp", "{===>} ion_segid_2=\"\";")
    print2file("gseq.inp", "{============================= output files ================================}")
    print2file("gseq.inp", "{* output structure file *}")
    print2file("gseq.inp", "{===>} structure_outfile=\"extended.mtf\";")
    print2file("gseq.inp", "{=========================== disulphide bonds ==============================}")
    print2file("gseq.inp", "{* Select pairs of cysteine residues that form disulphide bonds *}")
    print2file("gseq.inp", "{* First 2 entries are the segid and resid of the first cysteine (CYS A). *}")
    print2file("gseq.inp", "{* Second 2 entries are the segid and resid of the second cysteine (CYS B). *}")
    print2file("gseq.inp", "{+ table: rows=8 numbered")
    print2file("gseq.inp", "   cols=5 \"use\" \"segid CYS A\" \"resid CYS A\" \"segid CYS B\" \"resid CYS B\" +}")
    print2file("gseq.inp", "{+ choice: true false +}")
    print2file("gseq.inp", "{===>} ss_use_1=true;")
    print2file("gseq.inp", "{===>} ss_i_segid_1=\"\"; ss_i_resid_1=11;")
    print2file("gseq.inp", "{===>} ss_j_segid_1=\"\"; ss_j_resid_1=27;")
    print2file("gseq.inp", "{+ choice: true false +}")
    print2file("gseq.inp", "{===>} ss_use_2=true;")
    print2file("gseq.inp", "{===>} ss_i_segid_2=\"\"; ss_i_resid_2=45;")
    print2file("gseq.inp", "{===>} ss_j_segid_2=\"\"; ss_j_resid_2=73;")
    print2file("gseq.inp", "{+ choice: true false +}")
    print2file("gseq.inp", "{===>} ss_use_3=false;")
    print2file("gseq.inp", "{===>} ss_i_segid_3=\"\"; ss_i_resid_3=0;")
    print2file("gseq.inp", "{===>} ss_j_segid_3=\"\"; ss_j_resid_3=0;")
    print2file("gseq.inp", "{+ choice: true false +}")
    print2file("gseq.inp", "{===>} ss_use_4=false;")
    print2file("gseq.inp", "{===>} ss_i_segid_4=\"\"; ss_i_resid_4=0;")
    print2file("gseq.inp", "{===>} ss_j_segid_4=\"\"; ss_j_resid_4=0;")
    print2file("gseq.inp", "{+ choice: true false +}")
    print2file("gseq.inp", "{===>} ss_use_5=false;")
    print2file("gseq.inp", "{===>} ss_i_segid_5=\"\"; ss_i_resid_5=0;")
    print2file("gseq.inp", "{===>} ss_j_segid_5=\"\"; ss_j_resid_5=0;")
    print2file("gseq.inp", "{+ choice: true false +}")
    print2file("gseq.inp", "{===>} ss_use_6=false;")
    print2file("gseq.inp", "{===>} ss_i_segid_6=\"\"; ss_i_resid_6=0;")
    print2file("gseq.inp", "{===>} ss_j_segid_6=\"\"; ss_j_resid_6=0;")
    print2file("gseq.inp", "{+ choice: true false +}")
    print2file("gseq.inp", "{===>} ss_use_7=false;")
    print2file("gseq.inp", "{===>} ss_i_segid_7=\"\"; ss_i_resid_7=0;")
    print2file("gseq.inp", "{===>} ss_j_segid_7=\"\"; ss_j_resid_7=0;")
    print2file("gseq.inp", "{+ choice: true false +}")
    print2file("gseq.inp", "{===>} ss_use_8=false;")
    print2file("gseq.inp", "{===>} ss_i_segid_8=\"\"; ss_i_resid_8=0;")
    print2file("gseq.inp", "{===>} ss_j_segid_8=\"\"; ss_j_resid_8=0;")
    print2file("gseq.inp", "{=========================== carbohydrate links  ===========================}")
    print2file("gseq.inp", "{* Select pairs of residues that are linked *}")
    print2file("gseq.inp", "{* First entry is the name of the patch residue. *}")
    print2file("gseq.inp", "{* Second and third entries are the resid and segid for the atoms")
    print2file("gseq.inp", "   referenced by \"-\" in the patch. *}")
    print2file("gseq.inp", "{* Fourth and fifth entries are the resid and segid for the atoms")
    print2file("gseq.inp", "   referenced by \"+\" in the patch *}")
    print2file("gseq.inp", "{+ table: rows=6 numbered")
    print2file("gseq.inp", "          cols=6 \"use\" \"patch name\" \"segid -\" \"resid -\" \"segid +\" \"resid +\" +}")
    print2file("gseq.inp", "{+ choice: true false +}")
    print2file("gseq.inp", "{===>} carbo_use_1=false;")
    print2file("gseq.inp", "{===>} carbo_patch_1=\"B1N\";")
    print2file("gseq.inp", "{===>} carbo_i_segid_1=\"BBBB\"; carbo_i_resid_1=401;")
    print2file("gseq.inp", "{===>} carbo_j_segid_1=\"AAAA\"; carbo_j_resid_1=56;")
    print2file("gseq.inp", "{+ choice: true false +}")
    print2file("gseq.inp", "{===>} carbo_use_2=false;")
    print2file("gseq.inp", "{===>} carbo_patch_2=\"B1N\";")
    print2file("gseq.inp", "{===>} carbo_i_segid_2=\"BBBB\"; carbo_i_resid_2=402;")
    print2file("gseq.inp", "{===>} carbo_j_segid_2=\"AAAA\"; carbo_j_resid_2=182;")
    print2file("gseq.inp", "{+ choice: true false +}")
    print2file("gseq.inp", "{===>} carbo_use_3=false;")
    print2file("gseq.inp", "{===>} carbo_patch_3=\"\";")
    print2file("gseq.inp", "{===>} carbo_i_segid_3=\"\"; carbo_i_resid_3=0;")
    print2file("gseq.inp", "{===>} carbo_j_segid_3=\"\"; carbo_j_resid_3=0;")
    print2file("gseq.inp", "{+ choice: true false +}")
    print2file("gseq.inp", "{===>} carbo_use_4=false;")
    print2file("gseq.inp", "{===>} carbo_patch_4=\"\";")
    print2file("gseq.inp", "{===>} carbo_i_segid_4=\"\"; carbo_i_resid_4=0;")
    print2file("gseq.inp", "{===>} carbo_j_segid_4=\"\"; carbo_j_resid_4=0;")
    print2file("gseq.inp", "{+ choice: true false +}")
    print2file("gseq.inp", "{===>} carbo_use_5=false;")
    print2file("gseq.inp", "{===>} carbo_patch_5=\"\";")
    print2file("gseq.inp", "{===>} carbo_i_segid_5=\"\"; carbo_i_resid_5=0;")
    print2file("gseq.inp", "{===>} carbo_j_segid_5=\"\"; carbo_j_resid_5=0;")
    print2file("gseq.inp", "{+ choice: true false +}")
    print2file("gseq.inp", "{===>} carbo_use_6=false;")
    print2file("gseq.inp", "{===>} carbo_patch_6=\"\";")
    print2file("gseq.inp", "{===>} carbo_i_segid_6=\"\"; carbo_i_resid_6=0;")
    print2file("gseq.inp", "{===>} carbo_j_segid_6=\"\"; carbo_j_resid_6=0;")
    print2file("gseq.inp", "{========================= generate parameters =============================}")
    print2file("gseq.inp", "{* hydrogen flag - determines whether hydrogens will be retained *}")
    print2file("gseq.inp", "{* must be true for NMR, atomic resolution X-ray crystallography ")
    print2file("gseq.inp", "   or modelling.  Set to false for most X-ray crystallographic ")
    print2file("gseq.inp", "   applications at resolution > 1A *}")
    print2file("gseq.inp", "{+ choice: true false +}")
    print2file("gseq.inp", "{===>} hydrogen_flag=false;")
    print2file("gseq.inp", "{* set bfactor flag *}")
    print2file("gseq.inp", "{+ choice: true false +}")
    print2file("gseq.inp", "{===>} set_bfactor=true;")
    print2file("gseq.inp", "{* set bfactor value *}")
    print2file("gseq.inp", "{===>} bfactor=15.0;")
    print2file("gseq.inp", "{* set occupancy flag *}")
    print2file("gseq.inp", "{+ choice: true false +}")
    print2file("gseq.inp", "{===>} set_occupancy=true;")
    print2file("gseq.inp", "{* set occupancy value *}")
    print2file("gseq.inp", "{===>} occupancy=1.0;")
    print2file("gseq.inp", "{===========================================================================}")
    print2file("gseq.inp", "{         things below this line do not need to be changed unless           }")
    print2file("gseq.inp", "{         you need to apply patches - at the appropriate places marked      }")
    print2file("gseq.inp", "{===========================================================================}")
    print2file("gseq.inp", " ) {- end block parameter definition -}")
    print2file("gseq.inp", " checkversion 1.3")
    print2file("gseq.inp", " evaluate ($log_level=quiet)")
    print2file("gseq.inp", " {- read parameter files -}")
    print2file("gseq.inp", " parameter")
    print2file("gseq.inp", "  evaluate ($counter=1)")
    print2file("gseq.inp", "  evaluate ($done=false)")
    print2file("gseq.inp", "  while ( $done = false ) loop read")
    print2file("gseq.inp", "   if ( &exist_parameter_infile_$counter = true ) then")
    print2file("gseq.inp", "      if ( &BLANK%parameter_infile_$counter = false ) then")
    print2file("gseq.inp", "         @@&parameter_infile_$counter")
    print2file("gseq.inp", "      end if")
    print2file("gseq.inp", "   else")
    print2file("gseq.inp", "    evaluate ($done=true)")
    print2file("gseq.inp", "   end if")
    print2file("gseq.inp", "   evaluate ($counter=$counter+1)")
    print2file("gseq.inp", "  end loop read")
    print2file("gseq.inp", " end")
    print2file("gseq.inp", " {- read topology files -}")
    print2file("gseq.inp", " topology")
    print2file("gseq.inp", "  evaluate ($counter=1)")
    print2file("gseq.inp", "  evaluate ($done=false)")
    print2file("gseq.inp", "  while ( $done = false ) loop read")
    print2file("gseq.inp", "   if ( &exist_topology_infile_$counter = true ) then")
    print2file("gseq.inp", "      if ( &BLANK%topology_infile_$counter = false ) then")
    print2file("gseq.inp", "         @@&topology_infile_$counter")
    print2file("gseq.inp", "      end if")
    print2file("gseq.inp", "   else")
    print2file("gseq.inp", "     evaluate ($done=true)")
    print2file("gseq.inp", "   end if")
    print2file("gseq.inp", "   evaluate ($counter=$counter+1)")
    print2file("gseq.inp", "  end loop read")
    print2file("gseq.inp", " end")
    print2file("gseq.inp", " evaluate ($counter=1)")
    print2file("gseq.inp", " evaluate ($done=false)")
    print2file("gseq.inp", " while ( $done = false ) loop prot")
    print2file("gseq.inp", "   if ( &exist_prot_sequence_infile_$counter = true ) then")
    print2file("gseq.inp", "     if ( &BLANK%prot_sequence_infile_$counter = false ) then")
    print2file("gseq.inp", "       do (refx=0) (all)")
    print2file("gseq.inp", "       segment")
    print2file("gseq.inp", "         chain")
    print2file("gseq.inp", "           evaluate ($count=1)")
    print2file("gseq.inp", "           evaluate ($done2=false)")
    print2file("gseq.inp", "           while ( $done2 = false ) loop read")
    print2file("gseq.inp", "             if ( &exist_link_infile_$count = true ) then")
    print2file("gseq.inp", "               if ( &BLANK%link_infile_$count = false ) then")
    print2file("gseq.inp", "                  @@&link_infile_$count")
    print2file("gseq.inp", "               end if")
    print2file("gseq.inp", "             else")
    print2file("gseq.inp", "               evaluate ($done2=true)")
    print2file("gseq.inp", "             end if")
    print2file("gseq.inp", "             evaluate ($count=$count+1)")
    print2file("gseq.inp", "           end loop read")
    print2file("gseq.inp", "           sequence @@&prot_sequence_infile_$counter end")
    print2file("gseq.inp", "         end")
    print2file("gseq.inp", "       end")
    print2file("gseq.inp", "       do (segid=\"T^\" + encode($counter)) (attr refx=9999)")
    print2file("gseq.inp", "     end if")
    print2file("gseq.inp", "     if ( &exist_renumber_$counter = true ) then")
    print2file("gseq.inp", "         if ( &BLANK%renumber_$counter = false ) then")
    print2file("gseq.inp", "           evaluate ($segid=\"T^\" + encode($counter))")
    print2file("gseq.inp", "           do ( resid = adjustl(format(\"I4\",decode(resid) + &renumber_$counter - 1))) ")
    print2file("gseq.inp", "              ( (attr refx=9999) and segid $segid )")
    print2file("gseq.inp", "         end if")
    print2file("gseq.inp", "     end if")
    print2file("gseq.inp", "     evaluate ($counter=$counter+1)")
    print2file("gseq.inp", "   else")
    print2file("gseq.inp", "     evaluate ($done=true)")
    print2file("gseq.inp", "   end if")
    print2file("gseq.inp", " end loop prot")
    print2file("gseq.inp", " evaluate ($counter=1)")
    print2file("gseq.inp", " evaluate ($done=false)")
    print2file("gseq.inp", " while ( $done = false ) loop nseg")
    print2file("gseq.inp", "   if ( &exist_prot_sequence_infile_$counter = true ) then")
    print2file("gseq.inp", "     if ( &BLANK%prot_sequence_infile_$counter = false ) then")
    print2file("gseq.inp", "       evaluate ($segtmp=\"T^\" + encode($counter))")
    print2file("gseq.inp", "       do (segid=capitalize(&prot_segid_$counter)) (segid $segtmp)")
    print2file("gseq.inp", "     end if")
    print2file("gseq.inp", "     evaluate ($counter=$counter+1)")
    print2file("gseq.inp", "   else")
    print2file("gseq.inp", "     evaluate ($done=true)")
    print2file("gseq.inp", "   end if")
    print2file("gseq.inp", " end loop nseg")
    print2file("gseq.inp", " evaluate ($ssc=1)")
    print2file("gseq.inp", " evaluate ($done=false)")
    print2file("gseq.inp", " while ( $done = false ) loop ssbr")
    print2file("gseq.inp", "   if ( &exist_ss_use_$ssc = true ) then")
    print2file("gseq.inp", "     if ( &ss_use_$ssc = true ) then")
    print2file("gseq.inp", "       evaluate ($segidtmp1=capitalize(&ss_i_segid_$ssc))")
    print2file("gseq.inp", "       evaluate ($segidtmp2=capitalize(&ss_j_segid_$ssc))")
    print2file("gseq.inp", "       patch disu")
    print2file("gseq.inp", "         reference=1=(segid $QUOTE%segidtmp1 and resid &ss_i_resid_$ssc)")
    print2file("gseq.inp", "         reference=2=(segid $QUOTE%segidtmp2 and resid &ss_j_resid_$ssc)")
    print2file("gseq.inp", "       end")
    print2file("gseq.inp", "     end if")
    print2file("gseq.inp", "     evaluate ($ssc=$ssc+1)")
    print2file("gseq.inp", "   else")
    print2file("gseq.inp", "     evaluate ($done=true)")
    print2file("gseq.inp", "   end if")
    print2file("gseq.inp", " end loop ssbr")
    print2file("gseq.inp", " {* any special protein patches can be applied here *}")
    print2file("gseq.inp", " {===>}")
    print2file("gseq.inp", " {<===}")
    print2file("gseq.inp", " evaluate ($counter=1)")
    print2file("gseq.inp", " evaluate ($done=false)")
    print2file("gseq.inp", " while ( $done = false ) loop nucl")
    print2file("gseq.inp", "   if ( &exist_nucl_sequence_infile_$counter = true ) then")
    print2file("gseq.inp", "     if ( &BLANK%nucl_sequence_infile_$counter = false ) then")
    print2file("gseq.inp", "       do (refx=0) (all)")
    print2file("gseq.inp", "       segment")
    print2file("gseq.inp", "         chain")
    print2file("gseq.inp", "           evaluate ($count=1)")
    print2file("gseq.inp", "           evaluate ($done2=false)")
    print2file("gseq.inp", "           while ( $done2 = false ) loop read")
    print2file("gseq.inp", "             if ( &exist_link_infile_$count = true ) then")
    print2file("gseq.inp", "               if ( &BLANK%link_infile_$count = false ) then")
    print2file("gseq.inp", "                  @@&link_infile_$count")
    print2file("gseq.inp", "               end if")
    print2file("gseq.inp", "             else")
    print2file("gseq.inp", "               evaluate ($done2=true)")
    print2file("gseq.inp", "             end if")
    print2file("gseq.inp", "             evaluate ($count=$count+1)")
    print2file("gseq.inp", "           end loop read")
    print2file("gseq.inp", "           sequence @@&nucl_sequence_infile_$counter end")
    print2file("gseq.inp", "         end")
    print2file("gseq.inp", "       end")
    print2file("gseq.inp", "       do (segid=capitalize(&nucl_segid_$counter)) (attr refx=9999)")
    print2file("gseq.inp", "     end if")
    print2file("gseq.inp", "     evaluate ($counter=$counter+1)")
    print2file("gseq.inp", "   else")
    print2file("gseq.inp", "     evaluate ($done=true)")
    print2file("gseq.inp", "   end if")
    print2file("gseq.inp", " end loop nucl")
    print2file("gseq.inp", " {* patch rna sugars to dna here if needed - select the residues *}")
    print2file("gseq.inp", " {===>} ")
    print2file("gseq.inp", " for $resid in () loop dna")
    print2file("gseq.inp", "   patch deox reference=nil=(resid $resid) end")
    print2file("gseq.inp", " end loop dna")
    print2file("gseq.inp", " {<===}")
    print2file("gseq.inp", " {* any special nucleic acid patches can be applied here *}")
    print2file("gseq.inp", " {===>}")
    print2file("gseq.inp", " {<===}")
    print2file("gseq.inp", " evaluate ($counter=1)")
    print2file("gseq.inp", " evaluate ($done=false)")
    print2file("gseq.inp", " while ( $done = false ) loop carbo")
    print2file("gseq.inp", "   if ( &exist_carbo_sequence_infile_$counter = true ) then")
    print2file("gseq.inp", "     if ( &BLANK%carbo_sequence_infile_$counter = false ) then")
    print2file("gseq.inp", "       do (refx=0) (all)")
    print2file("gseq.inp", "       segment")
    print2file("gseq.inp", "         chain")
    print2file("gseq.inp", "           sequence @@&carbo_sequence_infile_$counter end")
    print2file("gseq.inp", "         end")
    print2file("gseq.inp", "       end")
    print2file("gseq.inp", "       do (segid=capitalize(&carbo_segid_$counter)) (attr refx=9999)")
    print2file("gseq.inp", "     end if")
    print2file("gseq.inp", "     evaluate ($counter=$counter+1)")
    print2file("gseq.inp", "   else")
    print2file("gseq.inp", "     evaluate ($done=true)")
    print2file("gseq.inp", "   end if")
    print2file("gseq.inp", " end loop carbo")
    print2file("gseq.inp", " evaluate ($carc=1)")
    print2file("gseq.inp", " evaluate ($done=false)")
    print2file("gseq.inp", " while ( $done = false ) loop cabr")
    print2file("gseq.inp", "   if ( &exist_carbo_use_$carc = true ) then")
    print2file("gseq.inp", "     if ( &carbo_use_$carc = true ) then")
    print2file("gseq.inp", "       evaluate ($segidtmp1=capitalize(&carbo_i_segid_$carc))")
    print2file("gseq.inp", "       evaluate ($segidtmp2=capitalize(&carbo_j_segid_$carc))")
    print2file("gseq.inp", "       patch &carbo_patch_$carc")
    print2file("gseq.inp", "         reference=-=(segid $QUOTE%segidtmp1 and")
    print2file("gseq.inp", "                      resid &carbo_i_resid_$carc)")
    print2file("gseq.inp", "         reference=+=(segid $QUOTE%segidtmp2 and")
    print2file("gseq.inp", "                      resid &carbo_j_resid_$carc)")
    print2file("gseq.inp", "       end")
    print2file("gseq.inp", "     end if")
    print2file("gseq.inp", "     evaluate ($carc=$carc+1)")
    print2file("gseq.inp", "   else")
    print2file("gseq.inp", "     evaluate ($done=true)")
    print2file("gseq.inp", "   end if")
    print2file("gseq.inp", " end loop cabr")
    print2file("gseq.inp", " {* any special carbohydrate patches can be applied here *}")
    print2file("gseq.inp", " {===>}")
    print2file("gseq.inp", " {<===}")
    print2file("gseq.inp", " evaluate ($counter=1)")
    print2file("gseq.inp", " evaluate ($done=false)")
    print2file("gseq.inp", " while ( $done = false ) loop prost")
    print2file("gseq.inp", "   if ( &exist_prost_sequence_infile_$counter = true ) then")
    print2file("gseq.inp", "     if ( &BLANK%prost_sequence_infile_$counter = false ) then")
    print2file("gseq.inp", "       do (refx=0) (all)")
    print2file("gseq.inp", "       segment")
    print2file("gseq.inp", "         chain")
    print2file("gseq.inp", "           sequence @@&prost_sequence_infile_$counter end")
    print2file("gseq.inp", "         end")
    print2file("gseq.inp", "       end")
    print2file("gseq.inp", "       do (segid=capitalize(&prost_segid_$counter)) (attr refx=9999)")
    print2file("gseq.inp", "     end if")
    print2file("gseq.inp", "     evaluate ($counter=$counter+1)")
    print2file("gseq.inp", "   else")
    print2file("gseq.inp", "     evaluate ($done=true)")
    print2file("gseq.inp", "   end if")
    print2file("gseq.inp", " end loop prost")
    print2file("gseq.inp", " {* any special prosthetic group patches can be applied here *}")
    print2file("gseq.inp", " {===>}")
    print2file("gseq.inp", " {<===}")
    print2file("gseq.inp", " evaluate ($counter=1)")
    print2file("gseq.inp", " evaluate ($done=false)")
    print2file("gseq.inp", " while ( $done = false ) loop liga")
    print2file("gseq.inp", "   if ( &exist_lig_sequence_infile_$counter = true ) then")
    print2file("gseq.inp", "     if ( &BLANK%lig_sequence_infile_$counter = false ) then")
    print2file("gseq.inp", "       do (refx=0) (all)")
    print2file("gseq.inp", "       segment")
    print2file("gseq.inp", "         chain")
    print2file("gseq.inp", "           sequence @@&lig_sequence_infile_$counter end")
    print2file("gseq.inp", "         end")
    print2file("gseq.inp", "       end")
    print2file("gseq.inp", "       do (segid=capitalize(&lig_segid_$counter)) (attr refx=9999)")
    print2file("gseq.inp", "     end if")
    print2file("gseq.inp", "     evaluate ($counter=$counter+1)")
    print2file("gseq.inp", "   else")
    print2file("gseq.inp", "     evaluate ($done=true)")
    print2file("gseq.inp", "   end if")
    print2file("gseq.inp", " end loop liga")
    print2file("gseq.inp", " {* any special ligand patches can be applied here *}")
    print2file("gseq.inp", " {===>}")
    print2file("gseq.inp", " {<===}")
    print2file("gseq.inp", " evaluate ($counter=1)")
    print2file("gseq.inp", " evaluate ($done=false)")
    print2file("gseq.inp", " while ( $done = false ) loop ion")
    print2file("gseq.inp", "   if ( &exist_ion_sequence_infile_$counter = true ) then")
    print2file("gseq.inp", "     if ( &BLANK%ion_sequence_infile_$counter = false ) then")
    print2file("gseq.inp", "       do (refx=0) (all)")
    print2file("gseq.inp", "       segment")
    print2file("gseq.inp", "         chain")
    print2file("gseq.inp", "           sequence @@&ion_sequence_infile_$counter end")
    print2file("gseq.inp", "         end")
    print2file("gseq.inp", "       end")
    print2file("gseq.inp", "       do (segid=capitalize(&ion_segid_$counter)) (attr refx=9999)")
    print2file("gseq.inp", "     end if")
    print2file("gseq.inp", "     evaluate ($counter=$counter+1)")
    print2file("gseq.inp", "   else")
    print2file("gseq.inp", "     evaluate ($done=true)")
    print2file("gseq.inp", "   end if")
    print2file("gseq.inp", " end loop ion")
    print2file("gseq.inp", " {* any special ion patches can be applied here *}")
    print2file("gseq.inp", " {===>}")
    print2file("gseq.inp", " {<===}")
    print2file("gseq.inp", " evaluate ($counter=1)")
    print2file("gseq.inp", " evaluate ($done=false)")
    print2file("gseq.inp", " while ( $done = false ) loop water")
    print2file("gseq.inp", "   if ( &exist_water_sequence_infile_$counter = true ) then")
    print2file("gseq.inp", "     if ( &BLANK%water_sequence_infile_$counter = false ) then")
    print2file("gseq.inp", "       do (refx=0) (all)")
    print2file("gseq.inp", "       segment")
    print2file("gseq.inp", "         chain")
    print2file("gseq.inp", "           sequence @@&water_sequence_infile_$counter end")
    print2file("gseq.inp", "         end")
    print2file("gseq.inp", "       end")
    print2file("gseq.inp", "       do (segid=capitalize(&water_segid_$counter)) (attr refx=9999)")
    print2file("gseq.inp", "     end if")
    print2file("gseq.inp", "     evaluate ($counter=$counter+1)")
    print2file("gseq.inp", "   else")
    print2file("gseq.inp", "     evaluate ($done=true)")
    print2file("gseq.inp", "   end if")
    print2file("gseq.inp", " end loop water")
    print2file("gseq.inp", " {* any special water patches can be applied here *}")
    print2file("gseq.inp", " {===>}")
    print2file("gseq.inp", " {<===}")
    print2file("gseq.inp", " {* any final patches can be applied here *}")
    print2file("gseq.inp", " {===>}")
    print2file("gseq.inp", " {<===}")
    print2file("gseq.inp", " if (&hydrogen_flag=false) then")
    #	print2file("gseq.inp", "   delete selection=( hydrogen ) end")
    print2file("gseq.inp", " end if")
    print2file("gseq.inp", " if (&set_bfactor=true) then")
    print2file("gseq.inp", "   do (b=&bfactor) ( all )")
    print2file("gseq.inp", " end if")
    print2file("gseq.inp", " if (&set_occupancy=true) then")
    print2file("gseq.inp", "   do (q=&occupancy) ( all )")
    print2file("gseq.inp", " end if")
    print2file("gseq.inp", " write structure output=&structure_outfile end")
    print2file("gseq.inp", " stop")   

def write_cns_generate_extended_file():
    print2file("extn.inp", "{+ file: generate_extended.inp +}")
    print2file("extn.inp", "{+ directory: nmr_calc +}")
    print2file("extn.inp", "{+ description: Generates an extended strand with ideal geometry ")
    print2file("extn.inp", "                for each connected polymer.  ")
    print2file("extn.inp", "                The molecular structure file must not contain any ")
    print2file("extn.inp", "                closed loops except disulfide bonds which are automatically")
    print2file("extn.inp", "                excluded from the generation of the strand conformation.  +}")
    print2file("extn.inp", "{+ authors: Axel T. Brunger +}")
    print2file("extn.inp", "{+ copyright: Yale University +}")
    print2file("extn.inp", "{- begin block parameter definition -} define(")
    print2file("extn.inp", "{======================= molecular structure =========================}")
    print2file("extn.inp", "{* structure file(s) *}")
    print2file("extn.inp", "{===>} structure_file=\"extended.mtf\";")
    print2file("extn.inp", "{* parameter file(s) *}")
    # CNS_TOPPAR:protein-allhdg5-4.para
    print2file("extn.inp", "{===>} par_1=\"CNS_TOPPAR:protein.param\";")
    print2file("extn.inp", "{===>} par_2=\"\";")
    print2file("extn.inp", "{===>} par_3=\"\";")
    print2file("extn.inp", "{===>} par_4=\"\";")
    print2file("extn.inp", "{===>} par_5=\"\";")
    print2file("extn.inp", "{======================= input parameters ============================}")
    print2file("extn.inp", "{* maximum number of trials to generate an acceptable structure *}")
    print2file("extn.inp", "{===>} max_trial=10;")
    print2file("extn.inp", "{=========================== output files ============================}")
    print2file("extn.inp", "{* output coordinates *}")
    print2file("extn.inp", "{===>} output_coor=\"extended.pdb\";")
    print2file("extn.inp", "                                  ")
    print2file("extn.inp", "{===========================================================================}")
    print2file("extn.inp", "{        things below this line do not normally need to be changed          }")
    print2file("extn.inp", "{===========================================================================}")
    print2file("extn.inp", " ) {- end block parameter definition -}")
    print2file("extn.inp", " checkversion 1.3")
    print2file("extn.inp", " evaluate ($log_level=quiet)")
    print2file("extn.inp", " structure @&structure_file end")
    print2file("extn.inp", " parameter")
    print2file("extn.inp", "   if (&par_1 # \" \") then")
    print2file("extn.inp", "      @@&par_1")
    print2file("extn.inp", "   end if")
    print2file("extn.inp", "   if (&par_2 # \" \") then")
    print2file("extn.inp", "      @@&par_2")
    print2file("extn.inp", "   end if")
    print2file("extn.inp", "   if (&par_3 # \" \") then")
    print2file("extn.inp", "      @@&par_3")
    print2file("extn.inp", "   end if")
    print2file("extn.inp", "   if (&par_4 # \" \") then")
    print2file("extn.inp", "      @@&par_4")
    print2file("extn.inp", "   end if")
    print2file("extn.inp", "   if (&par_5 # \" \") then")
    print2file("extn.inp", "      @@&par_5")
    print2file("extn.inp", "   end if")
    print2file("extn.inp", " end")
    print2file("extn.inp", "{ Set force constants for S-S bond lengths and angles to zero  }")
    print2file("extn.inp", "parameter")
    print2file("extn.inp", "   bonds ( name SG ) ( name SG ) 0. 1. ")
    print2file("extn.inp", "end")
    print2file("extn.inp", "igroup interaction=(all) (all) end")
    print2file("extn.inp", "ident (x) ( all )")
    print2file("extn.inp", "do (x=x/5.) ( all )")
    print2file("extn.inp", "do (y=random(0.5) ) ( all )")
    print2file("extn.inp", "do (z=random(0.5) ) ( all )")
    print2file("extn.inp", "flags exclude * include bond angle impr dihe vdw end")
    print2file("extn.inp", "parameter")
    print2file("extn.inp", "   nbonds")
    print2file("extn.inp", "      rcon=50. nbxmod=-3 repel=0.8 cutnb=6. ")
    print2file("extn.inp", "      rexp=2 irexp=2 inhibit=0.0 wmin=0.1 tolerance=0.5")
    print2file("extn.inp", "   end")
    print2file("extn.inp", "end")
    print2file("extn.inp", "evaluate ($count=1) ")
    print2file("extn.inp", "while ($count < 10 ) loop l1")
    print2file("extn.inp", "   do (x=x+gauss(0.1)) ( all ) ")
    print2file("extn.inp", "   do (y=y+gauss(0.1)) ( all ) ")
    print2file("extn.inp", "   do (z=z+gauss(0.1)) ( all ) ")
    print2file("extn.inp", "   minimize lbfgs nstep=200 nprint=200 end")
    print2file("extn.inp", "   evaluate ($count=$count+1)")
    print2file("extn.inp", "end loop l1")
    print2file("extn.inp", "evaluate ($accept=false) ")
    print2file("extn.inp", "evaluate ($trial=1) ")
    print2file("extn.inp", "while ($accept=false) loop accp")
    print2file("extn.inp", "   for $1 in id ( tag ) loop resi")
    print2file("extn.inp", "      igroup ")
    print2file("extn.inp", "         interaction=( byresidue (id $1 ) and not name SG ) ")
    print2file("extn.inp", "                     ( not name SG ) ")
    print2file("extn.inp", "      end")
    print2file("extn.inp", "      evaluate ($accept=true) ")
    print2file("extn.inp", "      print thres=0.1 bonds")
    print2file("extn.inp", "      if ($violations > 0) then")
    print2file("extn.inp", "         evaluate ($accept=false) ")
    print2file("extn.inp", "      end if")
    print2file("extn.inp", "      print thres=10. angles ")
    print2file("extn.inp", "      evaluate ($angles=$result)")
    print2file("extn.inp", "      if ($violations > 0) then")
    print2file("extn.inp", "         evaluate ($accept=false) ")
    print2file("extn.inp", "      end if")
    print2file("extn.inp", "      print thres=10. improper")
    print2file("extn.inp", "      if ($violations > 0) then")
    print2file("extn.inp", "         evaluate ($accept=false) ")
    print2file("extn.inp", "      end if")
    print2file("extn.inp", "      if ($accept=false) then")
    print2file("extn.inp", "         do (x=x+gauss(0.3)) ( byresidue (id $1 ) ) ")
    print2file("extn.inp", "         do (y=y+gauss(0.3)) ( byresidue (id $1 ) ) ")
    print2file("extn.inp", "         do (z=z+gauss(0.3)) ( byresidue (id $1 ) ) ")
    print2file("extn.inp", "      end if")
    print2file("extn.inp", "   end loop resi")
    print2file("extn.inp", "   igroup interaction=( all ) ( all ) end")
    print2file("extn.inp", "   parameter")
    print2file("extn.inp", "      nbonds")
    print2file("extn.inp", "         rcon=50. nbxmod=-3 repel=3. cutnb=10. ")
    print2file("extn.inp", "      end")
    print2file("extn.inp", "   end")
    print2file("extn.inp", "   flags exclude angle improper end")
    print2file("extn.inp", "   ")
    print2file("extn.inp", "   minimize lbfgs nstep=200 nprint=200 end")
    print2file("extn.inp", "   parameter")
    print2file("extn.inp", "      nbonds")
    print2file("extn.inp", "         rcon=50. nbxmod=-3 repel=0.8 cutnb=6. ")
    print2file("extn.inp", "      end")
    print2file("extn.inp", "   end")
    print2file("extn.inp", "   flags include angle improper end")
    print2file("extn.inp", "   ")
    print2file("extn.inp", "   evaluate ($count=1) ")
    print2file("extn.inp", "   while ($count < 5 ) loop l2")
    print2file("extn.inp", "      do (x=x+gauss(0.05)) ( all ) ")
    print2file("extn.inp", "      do (y=y+gauss(0.05)) ( all ) ")
    print2file("extn.inp", "      do (z=z+gauss(0.05)) ( all ) ")
    print2file("extn.inp", "      minimize lbfgs nstep=200 nprint=200 end")
    print2file("extn.inp", "      evaluate ($count=$count+1)")
    print2file("extn.inp", "   end loop l2")
    print2file("extn.inp", "   ")
    print2file("extn.inp", "   parameter")
    print2file("extn.inp", "      nbonds")
    print2file("extn.inp", "         rcon=50. nbxmod=3 repel=0.8 cutnb=6. ")
    print2file("extn.inp", "      end")
    print2file("extn.inp", "   end")
    print2file("extn.inp", "   ")
    print2file("extn.inp", "   minimize lbfgs nstep=300 nprint=300 end   ")
    print2file("extn.inp", "   minimize lbfgs nstep=300 nprint=300 end")
    print2file("extn.inp", "   igroup interaction=( not name SG ) ( not name SG ) end")
    print2file("extn.inp", "   energy end")
    print2file("extn.inp", "   evaluate ($accept=true) ")
    print2file("extn.inp", "   print thres=0.05 bonds")
    print2file("extn.inp", "   evaluate ($bonds=$result)")
    print2file("extn.inp", "   if ($violations > 0) then")
    print2file("extn.inp", "      evaluate ($accept=false) ")
    print2file("extn.inp", "   end if")
    print2file("extn.inp", "   print thres=10. angles ")
    print2file("extn.inp", "   evaluate ($angles=$result)")
    print2file("extn.inp", "   if ($violations > 0) then")
    print2file("extn.inp", "      evaluate ($accept=false) ")
    print2file("extn.inp", "   end if")
    print2file("extn.inp", "   print thres=10. improper")
    print2file("extn.inp", "   evaluate ($impr=$result)")
    print2file("extn.inp", "   if ($violations > 0) then")
    print2file("extn.inp", "      evaluate ($accept=false) ")
    print2file("extn.inp", "   end if")
    print2file("extn.inp", "   print thres=180. dihedral ")
    print2file("extn.inp", "   evaluate ($dihe=$result)")
    print2file("extn.inp", "   evaluate ($trial=$trial + 1) ")
    print2file("extn.inp", "   if ($trial > &max_trial ) then")
    print2file("extn.inp", "      exit loop accp")
    print2file("extn.inp", "   end if")
    print2file("extn.inp", "end loop accp")
    print2file("extn.inp", "remarks extended strand(s) generation")
    print2file("extn.inp", "remarks input molecular structure file=&structure_file ")
    print2file("extn.inp", "remarks final rms deviations (excluding disulfide bonds): ")
    print2file("extn.inp", "remarks    bonds=	 $bonds[F8.4] A  ")
    print2file("extn.inp", "remarks    angles=	 $angles[F8.4] degrees")
    print2file("extn.inp", "remarks    impropers= $impr[F8.4] degrees")
    print2file("extn.inp", "remarks    dihedrals= $dihe[F8.4] degrees (not used in some parameter sets!)")
    print2file("extn.inp", "remarks final van der Waals (repel) energy=$vdw kcal/mole")
    print2file("extn.inp", "write coordinates output=&output_coor format=PDBO end  \n")
    print2file("extn.inp", "stop\n")

#insert lines 398-415

file_rr = "1a3aA.dist.rr"
#lines 1-373 must be written
###############################################################
#116 -- MAKE A COPY OF THE FILES
###############################################################
#FOR TESTING ONLY
print("Start " + str(main.__file__) +": " + str(datetime.datetime.now()) + "\n")
dir_out = os.path.abspath('./myout/')
try:  
    os.mkdir(dir_out)
except OSError:  
    warn("Creation of the file/directory failed")

words = str(os.path.basename(str('1a3aA.dist.rr')))
words = words.split(".")
__id = words[0] + "." + words[1]
os.system("cp " + str(file_rr) + " " + str(dir_out) + "/")
file_rr = __id + ".rr"
os.system("cp " + str(file_rr) + " " + str(dir_out) + "/")
file_ss = __id +".ss"
file_rr = "1a3aA.dist.rr"
try:  
    os.chdir(dir_out)
except OSError:  
    warn("Directory change failed")

##############################################################
#Process Parameters - Line 133
##############################################################
seq = seq_rr(file_rr)
#check errors (if necessary) using chk_errors_seq(seq)
L = len(seq)
#
#lines 372-396
write_cns_seq(file_rr, "input.seq")
write_cns_generate_seq_file()
write_cns_generate_extended_file()
print("\nBuild extended mtf and pdb\n")
#if os.path.isfile("/bash
#377-380 must be written
f = open("job.sh", "w")
f.write("#!/bin/bash \n")
f.write("# CNS-CONFIGURATION\n")
f.write("source " + cns_suite + "/cns_solve_env.sh\n")
f.write("export KMP_AFFINITY=none\n")
f.write("touch extended.running\n")
f.write(cns_executable + " < gseq.inp > gseq.log \n")
f.write(cns_executable + " < extn.inp > extn.log \n")
f.write("rm extended.running\n")
f.close()
os.system("chmod +x job.sh")
os.system("./job.sh")
#if no extended.pdf confess "FAILED! extended.mtf not found!"
#os.system("rm -f gseq.log")
#os.system("rm -f extn.log")
#defined_atoms = xyz_pdb("extended.pdb", "all") #assuming xyz_pdb() alrady exists
#for each keys in residues:

rr2tbl(file_rr, "contact.tbl", "cb")
build_models()

