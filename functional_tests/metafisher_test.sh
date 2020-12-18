functional_tests/#!/usr/bin/env bash

# Run a command and check that the output is
# exactly equal the contents of a specified file
# ARG1: command we want to test as a string
# ARG2: a file path containing the expected output
# ARG3: expected exit status
function compare_result_dir {

    exit_status=$?
    output_dir=$2
    expected_output_dir=$3
    expected_exit_status=$4
    rm $output_dir/metafisher_*
    echo $1
    eval $1

    # verbose_message "Testing stdout and exit status: $1"
    #echo "searching for file in $expected_output_dir"
    for expected_output_file in $expected_output_dir/metafisher*; do
        let num_tests+=1
        filename=`basename "$expected_output_file"`
        output=$output_dir/$filename
        if [ ! -f "$output" ]; then
            let num_errors+=1
            echo result file $output not found
            continue
        fi

        difference=$(diff $output $expected_output_file)

        if [ -n "$difference" ]; then
            let num_errors+=1
            echo "Test output failed: $output vs $expected_output_file"
            echo "diff  $expected_output_file $output"
            # echo "Actual output:"
            # head "$output"
            # expected_output=$(cat $expected_output_file)
            # echo "Expected output:"
            # head "$expected_output_file"
            #echo " head Difference:"
            # head "$difference"
        fi

    done
}

function test_exit_status {
    let num_tests+=1
    output=$(eval $1)
    exit_status=$?
    expected_exit_status=$2
    if [ "$exit_status" -ne "$expected_exit_status" ]; then
        let num_errors+=1
        echo "Test exit status failed: $1"
        echo "Actual exit status: $exit_status"
        echo "Expected exit status: $expected_exit_status"
    fi
}

# Number of failed test cases
num_errors=0
# Total number of tests run
num_tests=0

expected_result='functional_tests/output/Desulfovibrio_vulgaris_DP4_expected'
data_dir='data_test/Desulfovibrio_vulgaris_DP4/'

generic_outdir='functional_tests/output/fresh_output_Desulfovibrio_vulgaris_DP4'

echo "regular execution test"
outdir=$generic_outdir
expdir=$expected_result

cmd="./metafisher/metafisher.py \
--gff $data_dir/sequence.gff \
--faa $data_dir/sequence.faa \
-o $outdir --name metafisher -v \
 > ${outdir}/cmd.out"


compare_result_dir "$cmd" $outdir $expdir 0

echo "resize execution test"

outdir=${generic_outdir}_resize
expdir=${expected_result}_resize

cmd="./metafisher/metafisher.py \
--gff $data_dir/sequence.gff \
--faa $data_dir/sequence.faa \
-o $outdir --name metafisher \
--resize --fna $data_dir/sequence.fna \
 > ${outdir}/cmd.out"

compare_result_dir "$cmd" $outdir $expdir 0


echo "exit status test"
test_exit_status "./metafisher/metafisher.py 2>  /dev/null" 2
test_exit_status "./metafisher/metafisher.py -h" 0


echo "rescue execution test"
outdir=${generic_outdir}_rescue
expdir=${expected_result}_rescue

cmd="./metafisher/metafisher.py \
--gff $data_dir/sequence.gff \
--faa $data_dir/sequence.faa \
-o $outdir --name metafisher -v \
--rescue --genomic_seq $data_dir/sequence.fasta \
 > ${outdir}/cmd.out"

compare_result_dir "$cmd" $outdir $expdir 0

# 3. End of testing - check if any errors occurrred
if [ "$num_errors" -gt 0 ]; then
    echo "$test_program failed $num_errors out of $num_tests tests"
    exit 1
else
    echo "$test_program passed all $num_tests successfully"
    exit 0
fi
