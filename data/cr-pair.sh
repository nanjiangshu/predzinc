#!/bin/bash 

#create database for zinc-binding prediction
idListFile=$WORKDIR/passe.idlist
modmpath=$CASIODATA3/wk/passe/modm-mtx/
metalProFile=$WORKDIR/passe/passe-mp-3.0/revised-closeMetalPro.dat
ssbondProFile=$WORKDIR/passe/passe.ssbond
min_boundRes=3
max_boundRes=5
encoding_scheme=1
keyMetalList=ZN

rootname=`rootname $idListFile`
train_modmfilelist=${rootname}.modmfilelist
getmodmfilepath -d $modmpath `cat $idListFile` > $train_modmfilelist


commonPara="-l $train_modmfilelist --mode train -t METAL --metal $metalProFile --keymetal $keyMetalList --ssbond $ssbondProFile  --bound-res $min_boundRes $max_boundRes --encode $encoding_scheme --not-use-ic"



#create database for pair-based vector
resList=CHDE
neg_filter=3
cutoff_score2=6.8
trainBaseName=train_pair_${resList}
trainVectorFile=$trainBaseName.vector
trainLabelFile=$trainBaseName.label
echo 
echo "create_svm_vector $commonPara -ns 2 -a $resList -s 0.01 -c $cutoff_score2 --min-hcres 2 --win-pair 100 --win-mini-pair 20 --win3 150 --win4 200 -K 12 -W 5 --neg-filter $neg_filter --vector  $trainVectorFile  --label $trainLabelFile"
create_svm_vector $commonPara -ns 2 -a $resList -s 0.01 -c $cutoff_score2 --min-hcres 2 --win-pair 100 --win-min-pair 20 --win3 150 --win4 200 -K 12 -W 5 --neg-filter $neg_filter --vector  $trainVectorFile  --label $trainLabelFile


    #for pair-based vectors
resList=CHDE
trainBaseName=train_pair_${resList}
./gisttrain.sh $trainBaseName
