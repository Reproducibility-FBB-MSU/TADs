import numpy as np
import tadtool.tad as tad

# set initial properties for data location
data_location = "matrices/"
armatus_location = "./armatus"
out_location = "yielded/"
lines = ["A549_NA_NA_", "HEK293_siRNA-CTCF_NA_", "HEK293_siRNA-Control_NA_", "HepG2_NA_NA_", "RAD21cv-HEK293_HRV-treated_NA_", "RAD21cv-HEK293_TEV-treated_NA_"]
chroms = ["15", "16", "17", "18", "19", "20", "21", "22", "Y"]
repls = ["1", "2"]
res = 20000
# set ranges for algorithms
cutoff_values = [i / 10 for i in range(0, 11)]
window_values = [res * i for i in range(1, 5)]
# run tadtool with insulation score
for line in lines:
    for repl in repls:
        for chrom in chroms:
            print(f"\nline={line} repl={repl} chrom={chrom} loading...   ")
            matrix = np.loadtxt(f"{data_location}{line}{repl}.20000.chr{chrom}.txt.gz")
            regions, some_ix = tad.load_regions(f"{data_location}{line}{repl}.20000.chr{chrom}.bed")
            for cutoff in cutoff_values:
                for window in window_values:
                    print(f"\rmethod=ii line={line} repl={repl} chrom={chrom} cutoff={cutoff} window={window} in progress     ", end="")
                    ii = tad.insulation_index(matrix, regions, window_size=window)
                    tads = tad.call_tads_insulation_index(ii, cutoff, regions=regions)
                    with open(f"{out_location}ii_enh_{line}{repl}.20000.chr{chrom}_window{window}_cutoff{cutoff}.txt", "w") as outfile:
                        for some_tad in tads:
                            outfile.write(f"{some_tad.start - 1}\t{some_tad.end}\n")
print("\rFinished." + " " * 100)
