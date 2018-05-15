for line in "A549_NA_NA_" "HEK293_siRNA-CTCF_NA_" "HEK293_siRNA-Control_NA_" "HepG2_NA_NA_" "RAD21cv-HEK293_HRV-treated_NA_" "RAD21cv-HEK293_TEV-treated_NA_"
do
	for repl in 1 2
	do
		for i in 15 16 17 18 19 20 21 22 "Y"
		do
			wget "https://cb.skoltech.ru/~agalicina/TAD/yielded/hicseg_"${line}${repl}".20000.chr"${i}".txt"
		done
	done
done
	