import lavaburst
import numpy as np

# set initial properties for data location
data_location = "matrices/"
armatus_location = "./armatus"
out_location = "yielded/"
lines = ["A549_NA_NA_", "HEK293_siRNA-CTCF_NA_", "HEK293_siRNA-Control_NA_", "HepG2_NA_NA_", "RAD21cv-HEK293_HRV-treated_NA_", "RAD21cv-HEK293_TEV-treated_NA_"]
chroms = ["15", "16", "17", "18", "19", "20", "21", "22", "Y"]
repls = ["1", "2"]
res = 20000
# set ranges for algorithms
corner_gamma = list(range(11))
# run lavaburst with armatus
for line in lines:
    for repl in repls:
        for chrom in chroms:
            print(f"\nline={line} repl={repl} chrom={chrom} loading...   ")
            matrix = np.loadtxt(f"{data_location}{line}{repl}.20000.chr{chrom}.txt.gz")
            good_bins = matrix.astype(bool).sum(axis=0) > 100
            for gamma in corner_gamma:
                print(f"\rmethod=lavaburst_corner line={line} repl={repl} chrom={chrom} gamma={gamma} in progress     ", end="")
                S = lavaburst.scoring.corner_score(matrix, gamma=gamma, binmask=good_bins)
                model = lavaburst.SegModel(S)
                segments = model.optimal_segmentation()
                np.savetxt(f"{out_location}lava_corner_{line}{repl}.20000.chr{chrom}.gamma{gamma}.txt", segments * res, delimiter="\t", fmt="%d")
