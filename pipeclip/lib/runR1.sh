#!/bin/sh

filename=$1
pvalue=$2
#declare -a epsilon
#declare -a steps
epsilon=(0.01 0.15 0.1)
steps=(0.1 0.08 0.05)
#sleep 500
#status=0

#check if file is generated
#while ["$status" == "0"]
#do
#	if [-p "$filename"]
#		then
#			status = 1
#			filesize = $(stat -c%s "$filename")
#		else
#			sleep 60
#		fi
#done
#
##make sure file does not change
#filestat = 0
#while ["$filestat"=="0"]
#do
#	currentsize = $(stat -c%s "$filename")
#	if ["$filesize" == "$currentsize"]
#		filestat = 1
#	else
#		sleep 60
#	fi
#done

#Call R function
r_status="1"
count=1
SCRIPTPATH="$(
  cd -- "$(dirname "$0")" >/dev/null 2>&1
  pwd -P
)"

for e in "${epsilon[@]}"; do
  for s in "${steps[@]}"; do
    echo "$e,$s"
    if [ -s "$filename.Converge.txt" ]; then
      echo
    else
      #echo "$e,$s"
      Rscript ${SCRIPTPATH}/ZTNB_tryCatch.R $filename $pvalue $e $s
    fi
    #echo "$filename.$count"
    count=$((count + 1))
  done
done
