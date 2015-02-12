#!/bin/bash

DISEASE="${1}";
PATIENT="${2}";
SERVERWORKINGPATH="${3}"; # put: random chars without spaces
FASTQ="${4}";
POOL="${5}"; # put: something useful for me
BARCODELIST="${6}"; # put: random chars without spaces
GENOME="${7}"; # put: assembly complete path
TMPDIR="${8}"; # exlusive temp folder
ASSOCIATIONFILE="${9}";
DBSCHEMA="${10}";
DBTABLE="${11}";
EXPORTPLUGIN="${12}"  # abs path of ExportDataToDB_PipePlugin.py
LTR="${13}";  # put: /opt/applications/scripts/isatk/elements/sequences/LTR.fa
LC="${14}";  # put: /opt/applications/scripts/isatk/elements/sequences/LC.rc.fa
CIGARGENOMEID="${15}"; # put something like 'hg19', 'mm9' ...
VECTORCIGARGENOMEID="${16}";  # put: random chars without spaces
SUBOPTIMALTHRESHOLD="${17}";  # put: 40
TAG="${18}"; # put: content of first 2 cells of AssociationFile
MAXTHREADS="${19}";
FILTERPLUGIN="${20}";

BASENAME="${DISEASE}.${PATIENT}.${POOL}";

VECTORGENOME="/opt/genome/vector/lv/bwa_7/lv.backbone.fa" # gemini set up

# # test echos
# for INPUTVAR in "$@"; do
# 	let INPUTVARNUM++; 
# 	printf -v INPUTNUM '%02d' $INPUTVARNUM;
#     echo "  => Input Variable: Order number = <${INPUTNUM}> ; Var Content = <${INPUTVAR}>";
# done

mkdir ${TMPDIR}
mkdir ${TMPDIR}/reads
mkdir ${TMPDIR}/bam
mkdir ${TMPDIR}/bed
mkdir ${TMPDIR}/sam
# checking vars
RUN_ID=`date +"%Y%m%d%H%M%S"`
RUN_NAME="${DISEASE}|${PATIENT}|${POOL}|${TAG}"


### ------------------ start TRIMMING LTR : GAMMA and LENTI ------------- ###
flexbar2.5 --reads ${FASTQ} --target ${TMPDIR}/reads/${BASENAME}.${TAG}.noLTRLC -f i1.8 -a ${LTR} --threads ${MAXTHREADS} -ae LEFT -at 3.6 -ai -4 -ao 22 -m 18 -q 1 --max-uncalled 800 ## this is for HIV
### ------------------ end TRIMMING LTR ------------- ###

# alignmet
bwa-stable mem -v 0 -k 18 -r 1 -M -T 15 -R "@RG\tID:${DISEASE}.${PATIENT}.${POOL}.${TAG}\tSM:${TAG}\tCN:Andrea.${DISEASE}.${PATIENT}.${POOL}" -t ${MAXTHREADS} ${GENOME} ${TMPDIR}/reads/${BASENAME}.${TAG}.noLTRLC.fastq > ${TMPDIR}/sam/${BASENAME}.${TAG}.noLTRLC.sam

# create BAM and sort them
samtools view -F 260 -uS ${TMPDIR}/sam/${BASENAME}.${TAG}.noLTRLC.sam | samtools sort - ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md ;
samtools index ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.bam

### ------------------ DATA FORMATTING ------------- ###
bamToBed -cigar -i ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.bam > ${TMPDIR}/bed/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.bed

### ------------------ IMPORT DATA INTO DB ------------- ###
NOWIS=`date +'%y-%m-%d %H:%M:%S'`
#echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Import BED data into DB (my script - running only for MY USER! so far)" 
${EXPORTPLUGIN} -b ${TMPDIR}/bed/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.bed -a ${ASSOCIATIONFILE} --patient ${PATIENT} --pool ${POOL} --tag ${TAG} --dbschema ${DBSCHEMA} --dbtable ${DBTABLE}

# Remove TMPDIR
rm -fr ${TMPDIR}
