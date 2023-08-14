# Gene Annotations 


Code is an ATAV query: 

echo "Getting gene annotations..."
resolution=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.resolution)
minsample=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.min_sample)

while read cmmd
	do
	echo $cmmd
    eval $cmmd

    sleep 1

done < $PROJECT/Results/CMH/gene_anno_commands_resolution_"$resolution"_min_sample_"$minsample".txt
exit