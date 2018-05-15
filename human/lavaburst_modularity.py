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
NaN = np.nan
# set ranges for algorithms
armatus_gamma = [i / 2 for i in range(-10, 11)]
modularity_gamma = list(range(41))
potts_gamma = list(range(101))
variance_gamma = list(range(101))
corner_gamma = list(range(11))
cutoff_values = list(range(21))
window_values = [res * i for i in list(np.arange(0.5, 5, 0.5))]
# run lavaburst with armatus
for line in lines:
    for repl in repls:
        for chrom in chroms:
            print(f"\nline={line} repl={repl} chrom={chrom} loading...   ")
            matrix = np.loadtxt(f"{data_location}{line}{repl}.20000.chr{chrom}.txt.gz")
            good_bins = matrix.astype(bool).sum(axis=0) > 100
            for gamma in modularity_gamma:
                print(f"\rmethod=lavaburst_modularity line={line} repl={repl} chrom={chrom} gamma={gamma} in progress     ", end="")
                S = lavaburst.scoring.modularity_score(matrix, gamma=gamma, binmask=good_bins)
                model = lavaburst.SegModel(S)
                segments = model.optimal_segmentation()
                np.savetxt(f"{out_location}lava_modularity_{line}{repl}.20000.chr{chrom}.gamma{gamma}.txt", segments * res, delimiter="\t", fmt="%d")
