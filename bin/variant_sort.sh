#!/bin/bash

MUTECT1=""
MUTECT2=""
VARSCANSNP=""
VARSCANINDEL=""
TUMORNAME=""
CONTROLNAME=""
usage() {
  echo "USAGE: $0 [-1 MUTECT1] [-2 MUTECT2] [-S VARSCANSNP] [-I VARSCANINDEL] [-t TUMORNAME] [-c CONTROLNAME]" 1>&2
}
exit_abnormal() {
  usage
  exit 1
}

while getopts ":1::2::S::I:t:c:" opt; do
  case "$opt" in
    1) MUTECT1=${OPTARG};;
    2) MUTECT2=${OPTARG};;
    S) VARSCANSNP=${OPTARG};;
    I) VARSCANINDEL=${OPTARG};;
    t) TUMORNAME=${OPTARG};;
    c) CONTROLNAME=${OPTARG};;
    :)
      echo "ERROR: -${OPTARG} requires an argument."
      exit_abnormal
      ;;
    *)
      echo "ERROR: Parameter: -${OPTARG} not known"
      exit_abnormal
      ;;
  esac
done

if [ -z "$VARSCANSNP" ] && [ -z "$VARSCANINDEL"]; then
  for txt in $MUTECT1 $MUTECT2; do
    cat $txt | grep '^\#[A-Z]' > header
    cat $txt | grep -v '#' > $txt.txt
    cat $txt | grep -v '#' | awk '{print $2, $3, $4, $5, $6, $7}'  > $txt.pos
  done
  else
    for txt in $MUTECT1 $MUTECT2 $VARSCANSNP $VARSCANINDEL; do
        cat $txt | grep '^\#[A-Z]' > header
        cat $txt | grep -v '#' > $txt.txt
        cat $txt | grep -v '#' | awk '{print $2, $3, $4, $5, $6, $7}'  > $txt.pos
    done
fi
if [ -z "$CONTROLNAME" ]; then
  sort *.pos | uniq -c | awk '$1 != '1' {print $0}'  > "$TUMORNAME"_least2
  cat *_least2 | awk '{print $2}' > column2
  grep -Ff column2 *.txt > "$TUMORNAME"_final
  cat header *_final > "$TUMORNAME"_final_variants
else
  sort *.pos | uniq -c | awk '$1 != '1' {print $0}'  > "$TUMORNAME"_"$CONTROLNAME"_least2
  cat *_least2 | awk '{print $2}' > column2
  grep -Ff column2 *.txt > "$TUMORNAME"_"$CONTROLNAME"_final
  cat header *_final > "$TUMORNAME"_"$CONTROLNAME"_final_variants
fi

rm *txt.txt
rm *.pos
rm *_least2
rm *_final
rm header
rm column2
