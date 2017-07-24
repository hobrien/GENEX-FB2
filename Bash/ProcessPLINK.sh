#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -l h_vmem=30G
#


echo "Starting to process PLINK files for $filename"
BASEDIR=`pwd`
cd ${1%/*}
pwd
filename=${1##*/} 
echo $filename

if [ ! -f $filename.bed  ]
then 
    echo "Converting file to Binary"
    plink --file $filename --out $filename --make-bed
    if [ $? -eq 0 ]
    then
        echo "Finished converting $filename to binary"
    else
        echo "Could not convert $filename to binary"
        exit 1
    fi
fi

if [ ! -f $filename.frq  ]
then 
    echo "Generating frequency file for $filename"
    plink --bfile $filename --out $filename --freq
    if [ $? -eq 0 ]
    then
        echo "Finished generating frequency file for $filename"
    else
        echo "Could not generate frequency file for $filename"
        exit 1
    fi
fi

#if [ ! -d ${filename%/*}/Imputation3 ]
#then
#   mkdir ${filename%/*}/Imputation3
#fi
#cd ${filename%/*}/Imputation3
if [ ! -f ${filename%/*}/Run-plink.sh ]
then
    echo "Running check-bim on $filename"
    perl $BASEDIR/Perl/HRC-1000G-check-bim.pl -b $filename.bim -f $filename.frq -r $BASEDIR/Data/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h
    if [ $? -eq 0 ]
    then
        echo "Finished running check-bim on $filename"
    else
        echo "Could not run check-bim on $filename"
        exit 1
    fi
fi

echo "running plink on $filename"
bash Run-plink.sh

echo "Finished to processing PLINK for $filename"
exit $?

