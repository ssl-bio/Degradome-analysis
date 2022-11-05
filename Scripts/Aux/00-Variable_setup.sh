#!/bin/bash

# Script to setup variables for bash and R scripts
# Usage: /bin/bash Var_setting.sh Oliver-2022 Degradome_vars_mint.txt 
#==================================================
# Check the number of arguments
if [[ -z $1 ]]
then
    echo "No basename for input/output was provided. Exiting..."
    exit 1
else
    echo "$1 will be used as basename for input/output"
    ibase=$1
fi

#Set variables
cd Env_variables

VARS="`set -o posix ; set`"

source Env_variables/${2}
echo "ibase=$ibase" > Degradome_${1}.txt
grep -vFe "$VARS" <<<"$(set -o posix ; set)" |
    grep -v ^VARS= | grep -vE "^BASH|^SHLVL" |
    sed 's/\[[0-9]\]=//g' >> Degradome_${1}.txt

unset VARS
