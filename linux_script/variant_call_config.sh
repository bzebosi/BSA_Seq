#!/bin/bash
# set to "false" to skip downloads
do_download="true"   
sv_call="true" # set to "false" to skip downloads
threads=8
genome_tmp="/nfs5/BPP/Leiboff_Lab/Brian/bsa/others"
ref_dir="${genome_tmp}/genomes"


# Samples
# locations of samples
declare -A sample_loc=(
    [q1]="/nfs5/BPP/Leiboff_Lab/Brian/bsa/bds1_ts"
    [q2]="/nfs5/BPP/Leiboff_Lab/Brian/bsa/S7_6508K"
    [q3]="/nfs5/BPP/Leiboff_Lab/Brian/bsa/tsh3"
    [q4]="/nfs5/BPP/Leiboff_Lab/Brian/bsa/rzl2"
    [q5]="/nfs5/BPP/Leiboff_Lab/Brian/bsa/S4_6508P"
    [q6]="/nfs5/BPP/Leiboff_Lab/Brian/bsa/S2_6609A"
    [q7]="/nfs5/BPP/Leiboff_Lab/Brian/bsa/S3_6508Q"
    [q8]="/nfs5/BPP/Leiboff_Lab/Brian/bsa/S1_6609B"
    [q9]="/nfs5/BPP/Leiboff_Lab/Brian/bsa/S8_6508F"
    [q10]="/nfs5/BPP/Leiboff_Lab/Brian/bsa/S6_6508N"
    [q11]="/nfs5/BPP/Leiboff_Lab/Brian/bsa/S5_6508O"
    [q12]="/nfs5/BPP/Leiboff_Lab/Brian/bsa/S9_5807N"
    [q13]="/nfs5/BPP/Leiboff_Lab/Brian/bsa/rzl1"
    [q14]="/nfs5/BPP/Leiboff_Lab/Brian/bsa/ts3"
    [q17]="/nfs5/BPP/Leiboff_Lab/Brian/bsa/mb6_ts071L"
    [q15]="/nfs5/BPP/Leiboff_Lab/Brian/bsa/mb7_tsN2490"
    [q16]="/nfs5/BPP/Leiboff_Lab/Brian/bsa/mb8_ts1967"
)

# Genomes of interest
declare -A genome_urls=(
  [b73]="https://download.maizegdb.org/Genomes/B73/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz"
  [mo17]="https://download.maizegdb.org/Genomes/Mo17/Zm-Mo17-REFERENCE-CAU-1.0/Zm-Mo17-REFERENCE-CAU-1.0.fa.gz"
  [w22]="https://download.maizegdb.org/Genomes/W22/Zm-W22-REFERENCE-NRGENE-2.0/Zm-W22-REFERENCE-NRGENE-2.0.fa.gz"
  [p39]="https://download.maizegdb.org/Genomes/NAM_Founders/Zm-P39-REFERENCE-NAM-1.0/Zm-P39-REFERENCE-NAM-1.0.fa.gz"
  [ki11]="https://download.maizegdb.org/Genomes/NAM_Founders/Zm-Ki11-REFERENCE-NAM-1.0/Zm-Ki11-REFERENCE-NAM-1.0.fa.gz"
  [a188]="https://download.maizegdb.org/Genomes/A188/Zm-A188-REFERENCE-KSU-1.0/Zm-A188-REFERENCE-KSU-1.0.fa.gz"
  [a632]="https://download.maizegdb.org/Genomes/NAM_Founders/Zm-A632-REFERENCE-CAAS_FIL-1.0/Zm-A632-REFERENCE-CAAS_FIL-1.0.fa.gz"
  [oh43]="https://download.maizegdb.org/Genomes/NAM_Founders/Zm-Oh43-REFERENCE-NAM-1.0/Zm-Oh43-REFERENCE-NAM-1.0.fa.gz"
  [cml247]="https://download.maizegdb.org/Genomes/NAM_Founders/Zm-CML247-REFERENCE-NAM-1.0/Zm-CML247-REFERENCE-NAM-1.0.fa.gz"
  [tx303]="https://download.maizegdb.org/Genomes/NAM_Founders/Zm-Tx303-REFERENCE-NAM-1.0/Zm-Tx303-REFERENCE-NAM-1.0.fa.gz"
  [ky21]="https://download.maizegdb.org/Genomes/NAM_Founders/Zm-Ky21-REFERENCE-NAM-1.0/Zm-Ky21-REFERENCE-NAM-1.0.fa.gz"
  [oh7b]="https://download.maizegdb.org/Genomes/NAM_Founders/Zm-Oh7B-REFERENCE-NAM-1.0/Zm-Oh7B-REFERENCE-NAM-1.0.fa.gz"
  [tzi8]="https://download.maizegdb.org/Genomes/NAM_Founders/Zm-Tzi8-REFERENCE-NAM-1.0/Zm-Tzi8-REFERENCE-NAM-1.0.fa.gz"
  [nc350]="https://download.maizegdb.org/Genomes/NAM_Founders/Zm-NC350-REFERENCE-NAM-1.0/Zm-NC350-REFERENCE-NAM-1.0.fa.gz"
  [cml277]="https://download.maizegdb.org/Genomes/NAM_Founders/Zm-CML277-REFERENCE-NAM-1.0/Zm-CML277-REFERENCE-NAM-1.0.fa.gz"
  [b97]="https://download.maizegdb.org/Genomes/NAM_Founders/Zm-B97-REFERENCE-NAM-1.0/Zm-B97-REFERENCE-NAM-1.0.fa.gz"
)

goi=(b73 w22c a619)

# Collect and sort files
files=(q1)