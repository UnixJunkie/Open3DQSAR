#!/usr/bin/env bash

# Open3DQSAR self-test
# run after building and installing Open3DQSAR:
# $ ./test.sh
# The output should be "Test completed successfully"
# In case of failure check test/${OUTPUT_FILE}


temp_ref=`mktemp`
temp_test=`mktemp`
TEST_COMPLETED_MSG="Successful completion"
COPY_FILES=\
"binding_data_36_compounds.txt \
ref_e2.sdf \
sample_input_MM.inp"
INPUT_FILE=sample_input_MM.inp
OUTPUT_FILE=sample_input_MM.out
TEST_RESULTS=test_results
REFERENCE_RESULTS=reference_results
OPEN3DTOOL=open3dqsar


# Before leaving, clean up temporary files and cd
# where the user originally was
clean_exit()
{
  rm -f $temp_ref $temp_test
  cd $old_cwd
  exit $1
}

# In case the test is killed before completion
abrupt_exit()
{
  cat << eof
Test aborted
eof
  clean_exit 1
}

# Get the numerical values and round them up to three decimals
get_predicted_val()
{
  sed -n "/BGN COMMAND #00.${1} /,\
/END COMMAND #00.${1} /p" \
    | sed -n '/^PC/,/^$/p' | sed 1,2d | sed '$d' \
    | awk '{printf "%.3f\n", $NF}'
}


# Save the folder where the user originally was
old_cwd=$PWD
# Catch all premature death signals
trap abrupt_exit SIGTSTP SIGINT SIGTERM SIGKILL
# cd into the "test" folder
cwd=`dirname $0`
if [ -z $cwd ]; then
  cwd=.
fi
cd $cwd
# Make a test_results folder, copy the necessary
# input files there, then run the test in there
rm -rf $TEST_RESULTS
mkdir $TEST_RESULTS
cp -R $COPY_FILES $TEST_RESULTS
cd $TEST_RESULTS

# Run the self-test
${OPEN3DTOOL} -i ${INPUT_FILE} -o ${OUTPUT_FILE}
# If the test did not complete successfully,
# print an error message and quit
if (! grep >&/dev/null "$TEST_COMPLETED_MSG" \
  < ${OUTPUT_FILE}); then
  cat << eof
The test case did not complete successfully.
Please check $cwd/${TEST_RESULTS}/${OUTPUT_FILE}
eof
  clean_exit 1
fi
# Find the last PREDICT command ID in the file
predict_cmd=`grep 'BGN COMMAND #00\..... - PREDICT' \
  < ../${REFERENCE_RESULTS}/${OUTPUT_FILE} | tail -n 1 \
  | awk -F\# '{print $2}' | awk '{print $1}' | awk -F. '{print $2}'`
# Get the reference results
get_predicted_val $predict_cmd \
  < ../${REFERENCE_RESULTS}/${OUTPUT_FILE} > $temp_ref
# Get the results obtained in the current test
get_predicted_val $predict_cmd \
  < ${OUTPUT_FILE} > $temp_test
# If they differ, print a warning message and exit
if (! diff >&/dev/null $temp_ref $temp_test); then
  cat << eof
$cwd/${OUTPUT_FILE} presents numerical differences
with respect to $cwd/${REFERENCE_RESULTS}/${OUTPUT_FILE}
Please check $cwd/${TEST_RESULTS}/${OUTPUT_FILE}
eof
  clean_exit 1
fi
# The test was OK, remove the test folder and exit
cat << eof
Test completed successfully
eof
cd ..
rm -rf $TEST_RESULTS
clean_exit 0
