names=names.tab
reads=reads.tab
aligned=aligned.tab
unique=unique.tab
OUT=stats.tab
ls *.log | xargs -n 1 echo > $names
ls *.log | xargs -n 1 grep -o 'read pairs: [0-9]\+' | cut -f 3 -d ' ' > $reads
for f in *.log
do
	nums=`grep 'aligned' $f | grep -o ' [0-9]\+ ' | awk '{printf "%s",$0} END {print ""}'`
	echo $nums | awk '{print $2+$3}' >> $aligned
	echo $nums | awk '{print $2}' >> $unique
done

printf "File\tPairs\tAligned\tUnique\n" > $OUT
paste $names $reads $aligned $unique >> $OUT
rm $aligned $unique $names $reads
