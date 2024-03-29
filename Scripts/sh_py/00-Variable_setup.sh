#!/bin/bash

# Script to setup variables for bash and R scripts

# Usage: /bin/bash Scripts/sh_py/00-Variable_setup.sh <Project_name> <Variable_specification_file> The last located in 'Env_variables'

# Example: /bin/bash ./Scripts/sh_py/00-Variable_setup.sh Oliver-2022 Oliver-2022_vars.txt
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
cd Env_variables || echo "Directory not found";

VARS="$(set -o posix ; set)"

source "$2"
echo "ibase=$ibase" > Degradome_"$1".txt
echo "doc_title='$doc_title'" >> Degradome_"$1".txt
grep -vFe "$VARS" <<<"$(set -o posix ; set)" |
    grep -v ^VARS= | grep -vE "^BASH|^SHLVL" |
    sed 's/\[[0-9]\]=//g' >> Degradome_"$1".txt

unset VARS
