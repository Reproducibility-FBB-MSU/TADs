for line in "A549_NA_NA_" "HEK293_siRNA-CTCF_NA_" "HEK293_siRNA-Control_NA_" "HepG2_NA_NA_" "RAD21cv-HEK293_HRV-treated_NA_" "RAD21cv-HEK293_TEV-treated_NA_"
do
	for repl in 1 2
	do
		for chrom in 15 16 17 18 19 20 21 22 "Y"
		do
			for gamma in `seq -5 0.5 5`
			do
				echo $line $repl $chrom $gamma
				./armatus -i matrices/${line}${repl}.20000.chr${chrom}.txt.gz -g $gamma -j -o yielded/armatus_${line}${repl}.20000.chr${chrom}.${gamma} -r 20000
			done
		done
	done
done
	