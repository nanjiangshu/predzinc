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


#create database for single site vector
resListArray=(C  H  D  E)
neg_filterArray=(2 3 6 8)
for (( i=0; i <4; i++)); do
    resList=${resListArray[$i]}
    neg_filter=${neg_filterArray[$i]};
    trainBaseName=train_single_${resList}
    trainVectorFile=$trainBaseName.vector
    trainLabelFile=$trainBaseName.label
    echo "create_svm_vector $commonPara -ns 1  -s -10 -c -10 -K 12 -P 6 -a $resList --min-hcres 1 --not-use-shcr --neg-filter $neg_filter --vector  $trainVectorFile --label $trainLabelFile" 
    create_svm_vector $commonPara -ns 1  -s -10 -c -10 -K 12 -P 6 -a $resList --min-hcres 1 --not-use-shcr --neg-filter $neg_filter --vector  $trainVectorFile --label $trainLabelFile 
done

#create database for pair-based vector
resList=CHDE
neg_filter=3
cutoff_score2=6.8
trainBaseName=train_pair_${resList}
trainVectorFile=$trainBaseName.vector
trainLabelFile=$trainBaseName.label
echo 
echo "create_svm_vector $commonPara -ns 2 -a $resList -s 0.01 -c $cutoff_score2 --min-hcres 2 --win-pair 100 --win-mini-pair 20 --win3 150 --win4 200 -K 12 -W 5 --neg_filter $neg_filter --vector  $trainVectorFile  --label $trainLabelFile"
create_svm_vector $commonPara -ns 2 -a $resList -s 0.01 -c $cutoff_score2 --min-hcres 2 --win-pair 100 --win-min-pair 20 --win3 150 --win4 200 -K 12 -W 5 --neg_filter $neg_filter --vector  $trainVectorFile  --label $trainLabelFile


#running gist svm to create the weight files
#for single-site vectors
resListArray=(C  H  D  E)
for (( i=0; i <4; i++)); do
    resList=${resListArray[$i]}
    trainBaseName=train_single_${resList}
    ./gisttrain.sh $trainBaseName
done 

    #for pair-based vectors
resList=CHDE
trainBaseName=train_pair_${resList}
./gisttrain.sh $trainBaseName
