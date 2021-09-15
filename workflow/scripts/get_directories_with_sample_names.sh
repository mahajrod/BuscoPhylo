#!/bin/bash

files_per_dir=20

print_usage() {
  echo "copy names and split files of a directory into multiple sub directories."
  echo "Usage: '-f' files in directory"
  echo "'-i' input directory path"
  echo "'-o' output directory path"
  }

while getopts 'f:i:o:' flag; do
  case "${flag}" in
    f) files_per_dir="${OPTARG}" ;;
    i) inputdir="${OPTARG}" ;;
    o) outdir="${OPTARG}" ;;
    *) print_usage
      exit 1 ;;
  esac
done

N=0 # directory counter
n=0 # file counter

mkdir -p $outdir
dir_name="${outdir}/dir_"

for file in `ls $inputdir/*.fna`; do
  if [ "$(( n % files_per_dir ))" -eq 0 ]; then
    N=$(( N + 1 ))
    dir="${outdir}/dir_$N"
    printf 'Creating directory %s\n' "$dir"
    mkdir "$dir"
  fi
  n=$(( n + 1 ))
  printf 'Moving %s to %s\n' "$file" "$dir"
  filename=echo ${file%.*} | rev | cut -d'_' -f 1 | rev
  touch "${dir}/${filename}"
done