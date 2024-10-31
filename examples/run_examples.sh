#!/bin/bash

# Execute from this script's directory
script_dir=$(realpath $(dirname $0))
cd $script_dir
declare -a status_array

# Declase which tests to run
declare -a example_cases=("snapshot" "solidification_mpstats" "solidification_multibeam" "snapshot_T_hist")

# Loop through tests without MPI
echo
echo Running tests without MPI:
echo
executable="$script_dir/../install/bin/3DThesis"
for case_dir in "${example_cases[@]}"; do

    # Change to case directory and clean any old files
    echo "- $case_dir"
    cd $case_dir
    [ -d Data ] && rm -r Data
    [ -f test.log ] && rm test.log
    [ -f 0 ] && rm 0

    # Execute the case and check that output files exist
    status=1
    status_text="Execution Not Successful"
    $executable ./ParamInput.txt > test.log 2>&1
    if [ -d Data ]; then
        file_count="$(find Data -name "$(basename $case_dir)*.csv" -printf "." | wc -m)"
        if [ $(expr $file_count) > 0 ]; then
            status_text="Output File Exists"
            status=0
        fi
    fi
    echo "  - Test Status: $status_text ($status)"
    status_array+=($status)
    cd $script_dir
done

# Loop through tests with MPI
echo
echo Running tests with MPI:
echo
executable="mpirun -n 2 --oversubscribe $script_dir/../install/bin/3DThesis"
for case_dir in "${example_cases[@]}"; do

    # Change to case directory and clean any old files
    echo "- $case_dir"
    cd $case_dir
    [ -d Data ] && rm -r Data
    [ -f test.log ] && rm test.log
    [ -f 0 ] && rm 0 # this file gets created sometimes by the execute command

    # Execute the case and check that output files exist
    status=1
    status_text="Execution Not Successful"
    $executable ./ParamInput.txt > test.log 2>&1
    if [ -d Data ]; then
        file_count="$(find Data -name "$(basename $case_dir)*.csv" -printf "." | wc -m)"
        if [ $(expr $file_count) > 0 ]; then
            status_text="Output File Exists"
            status=0
        fi
    fi
    echo "  - Test Status: $status_text ($status)"
    status_array+=($status)
    cd $script_dir
done

# Output overall summary of tests
echo
echo "Task Status Summary:"
if [[ "${status_array[*]}" =~ 1 ]]; then
    echo Fail
else
    echo Pass
fi