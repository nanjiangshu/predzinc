#!/bin/bash -x
function PrintHelp()
{
    echo "usage: gish.sh train_basename test_basename"
}

if [ $#  -lt 2 ]; then
    echo "too few arguments"
    PrintHelp
    exit
fi

TRAIN=$1
TEST=$2

#gist-train-svm -radial -diagfactor 0.1 -power 2.0 -train  $TRAIN.dat -class $TRAIN.label > $TRAIN.weight
#gist-train-svm -radial -diagfactor 0.1  -train  $TRAIN.dat -class $TRAIN.label > $TRAIN.weight
gist-train-svm -radial -power 2.0  -train  $TRAIN.dat -class $TRAIN.label > $TRAIN.weight
#gist-train-svm  -train  $TRAIN.dat -class $TRAIN.label > $TRAIN.weight
gist-classify  -train $TRAIN.dat -learned $TRAIN.weight -test $TEST.dat >$TEST.predict
gist-score-svm -test  $TEST.label $TEST.predict -roc $TEST.roc $TRAIN.weight

