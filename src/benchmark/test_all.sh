#!/bin/bash
# path to binary file for recomended build tree
# (see https://kelidas.gitlab.io/partmod-site/GettingStarted.html)
# NOTE it is a path from one of the subdirectories
path_to_binary_file="../../../../partmod-build/DiscreteModel"
remove_fol="rm -Rf results"
copy_bin="cp $path_to_binary_file DiscreteModel"
run_test="nohup ./DiscreteModel master.inp"
remove_binary="rm DiscreteModel"
remove_out_files="rm *.out"
remove_results="$remove_binary && $remove_fol && $remove_out_files"
save_output=false
while true; do
    read -p $'Do you wish to delete result files after running the test?\nY/N\n' yn
    case $yn in
        [Yy]* ) echo "yes - data will be deleted" ; break;;
        "") echo "yes  - data will be deleted"; break;;
        [SsNn]* ) echo "results will remain in the benchmark folders"; save_output=true; break;;
        * ) echo "Please answer yes or no.";;
    esac
done

for folder in $( ls -d */ ); do
    echo "---------------------------------------------------------"
    cd $folder 
    $copy_bin 
    echo "running test $folder"
    $run_test
    echo "test succesfully finished"
    if [ $save_output = false ]; then
        echo "removing results"
        $remove_results
    else
        echo "results are kept in folder $folder"
    fi
    cd - > /dev/null
done


