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
DBHOST="${10}";
DBUSER="${11}";
DBPASSWD="${12}";
DBPORT="${13}";
DBSCHEMA="${14}";
DBTABLE="${15}";
EXPORTPLUGIN="${16}"  # abs path of ExportDataToDB_PipePlugin.py
LTR="${17}";  # put: /opt/applications/scripts/isatk/elements/sequences/LTR.fa
LC="${18}";  # put: /opt/applications/scripts/isatk/elements/sequences/LC.fa
CIGARGENOMEID="${19}"; # put something like 'hg19', 'mm9' ...
VECTORCIGARGENOMEID="${20}";  # put: random chars without spaces
SUBOPTIMALTHRESHOLD="${21}";  # put: 40
TAG="${22}"; # put: content of first 2 cells of AssociationFile
MAXTHREADS="${23}";

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
#echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Trimming LTR"
flexbar2.5 --reads ${FASTQ} --target ${TMPDIR}/reads/${BASENAME}.${TAG}.noLTR -f i1.8 -a ${LTR} --threads ${MAXTHREADS} -ae LEFT -at 2 -ai -4 -ao 15 -m 18 -q 5 ## this is for HIV
### ------------------ end TRIMMING LTR ------------- ###

### ------------------ TRIMMING LC ------------- ###
#echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Trimming LC"
#fastq-mcf ${LC} mld02.454.pool6.${TAG}.noLTR.fastq.gz -o mld02.454.pool6.${TAG}.noLTRLC.fastq.gz -m 12 -p 20 -l 30 -q 10 -P 33
## la precedente istruzione fallisce nel caso di più LC concatenate, esempio la read HHAUOBH02JFYTO
flexbar2.5 --reads ${TMPDIR}/reads/${BASENAME}.${TAG}.noLTR.fastq --target ${TMPDIR}/reads/${BASENAME}.${TAG}.noLTRLC -f i1.8 -a ${LC} --threads ${MAXTHREADS} -ae RIGHT -at 4 -ao 8 -m 2 -q 5 ;


bwa-7.5 mem -r 1 -M -T 15 -R "@RG\tID:${DISEASE}.${PATIENT}.${POOL}.${TAG}\tSM:${TAG}\tCN:Andrea.${DISEASE}.${PATIENT}.${POOL}" -t ${MAXTHREADS} ${GENOME} ${TMPDIR}/reads/${BASENAME}.${TAG}.noLTRLC.fastq > ${TMPDIR}/sam/${BASENAME}.${TAG}.noLTRLC.sam

# create BAM and sort them
#echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Creating BAM and indexes (filter from here the dataset using only valid reads: mapped and primary)"
#echo "samtools view -F 260 -uS ${TMPDIR}/sam/${BASENAME}.${TAG}.noLTRLC.sam > ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.bam;"
# samtools view -F 260 -q 5 -uS ${TMPDIR}/sam/${BASENAME}.${TAG}.noLTRLC.sam > ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.bam;
# samtools sort ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.bam ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted;
samtools view -F 260 -q 5 -uS ${TMPDIR}/sam/${BASENAME}.${TAG}.noLTRLC.sam | samtools sort - ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted;
samtools index ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.bam ;

### ------------------ RECALIBRATION ------------- ###
#echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Recreate/Recalibrate/Fill the MD tag"
samtools fillmd -b ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.bam ${GENOME} > ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.bam 
samtools index ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.bam
rm ${TMPDIR}/bam/${BASENAME}.${TAG}.sorted.bam

### ------------------ FILTERING ------------- ###
#echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Filter by CIGAR and MD tag"
filter_by_cigar_bam --bam ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.bam -o ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.toremovebycigar.list --minStartingMatches_fromNbp 2 --minStartingMatches 5 -p ${MAXTHREADS} -g ${CIGARGENOMEID} --minStartingBasesNoIndels 7 --compareSubOptimal --suboptimalThreshold ${SUBOPTIMALTHRESHOLD} --SAalnAction ignore --XSlikeTag XS --ASlikeTag AS --endClipThreshold 1000 --singleEnd 

CHIMERANROWS=`wc -l ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.toremovebycigar.list | cut -d' ' -f1`
if [ $CHIMERANROWS -gt 0 ] 
	then
	#echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] --- found reads to remove, now removing them all from BAM file"
	FilterSamReads INPUT=${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.bam FILTER=excludeReadList RLF=${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.toremovebycigar.list SO=coordinate O=${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.bam
else
	#echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] --- reads to remove not found, go ahead"
	cp ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.bam ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.bam
fi
samtools index ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.bam
	
#echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Filtering data (Bamtools)"
bamtools filter -in ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.bam -mapQuality ">=12" -out ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.bam

### ------------------ DATA FORMATTING ------------- ###
#echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Convert BAM to BED (returns score as CIGAR)"
bamToBed -cigar -i ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.bam > ${TMPDIR}/bed/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.bed

### ------------------ IMPORT DATA INTO DB ------------- ###
NOWIS=`date +'%y-%m-%d %H:%M:%S'`
#echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Import BED data into DB (my script - running only for MY USER! so far)" 
${EXPORTPLUGIN} -b ${TMPDIR}/bed/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.bed -a ${ASSOCIATIONFILE} --patient ${PATIENT} --pool ${POOL} --tag ${TAG} --host ${DBHOST} --user ${DBUSER} --passwd ${DBPASSWD} --port ${DBPORT} --dbschema ${DBSCHEMA} --dbtable ${DBTABLE}

# Remove TMPDIR
#rm -fr ${TMPDIR}
