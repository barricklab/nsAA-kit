source("R3C-microplate-timecourse.R")

# -------------- Global settings -------------
# Used in the summary graphs
OD600.range = c(0.3, 0.4)

# -------------- Standard plates --------------

analyze_microplate("input-paper/3-aminotyrosine-032615.xlsx",
                   "input-paper/key-3-aminotyrosine-032615.xlsx",
                   "output-paper/3-aminotyrosine-032615")

analyze_microplate("input-paper/3-iodotyrosine-032615.xlsx",
                   "input-paper/key-3-iodotyrosine-032615.xlsx",
                   "output-paper/3-iodotyrosine-032615")

analyze_microplate("input-paper/3-nitrotyrosine-041215.xlsx",
                   "input-paper/key-3-nitrotyrosine-041215.xlsx",
                   "output-paper/3-nitrotyrosine-041215")

analyze_microplate("input-paper/4-cyanophenylalanine-052115.xlsx",
                   "input-paper/key-4-cyanophenylalanine-052115.xlsx",
                   "output-paper/4-cyanophenylalanine-052115")

analyze_microplate("input-paper/5-hydroxytryptophan-081415.xlsx",
                   "input-paper/key-5-hydroxytryptophan-081415.xlsx",
                   "output-paper/5-hydroxytryptophan-081415")

# 4-azidophenylalanine-080115 has clone 1,replicate 3 (=4 total wells) removed

analyze_microplate("input-paper/4-azidophenylalanine-080115.xlsx",
                   "input-paper/key-4-azidophenylalanine-080115.xlsx",
                   "output-paper/4-azidophenylalanine-080115")

analyze_microplate("input-paper/o-nitrobenzyltyrosine-050215.xlsx",
                   "input-paper/key-o-nitrobenzyltyrosine-050215.xlsx",
                   "output-paper/o-nitrobenzyltyrosine-050215")

analyze_microplate("input-paper/tyrosine-081415.xlsx",
                   "input-paper/key-tyrosine-081415.xlsx",
                   "output-paper/tyrosine-081415")

summarize_microplates("output-paper/aaRS_comparison",
                      c(
                        "output-paper/tyrosine-081415",
                        "output-paper/3-aminotyrosine-032615",
                        "output-paper/3-iodotyrosine-032615",
                        "output-paper/3-nitrotyrosine-041215",
                        "output-paper/4-cyanophenylalanine-052115",
                        "output-paper/5-hydroxytryptophan-081415",
                        "output-paper/4-azidophenylalanine-080115",
                        "output-paper/o-nitrobenzyltyrosine-050215"
                      ),
                      OD600.range
                      )

# # -------------- Preconditioned --------------
#
# # 4-azidophenylalanine-preconditioned-080115 has clone 1,replicate 3 (=4 total wells) removed

analyze_microplate("input-paper/4-azidophenylalanine-preconditioned-080115.xlsx",
                   "input-paper/key-4-azidophenylalanine-preconditioned-080115.xlsx",
                   "output-paper/preconditioned-4-azidophenylalanine-080115")

analyze_microplate("input-paper/o-nitrobenzyltyrosine-preconditioned-050215.xlsx",
                   "input-paper/key-o-nitrobenzyltyrosine-preconditioned-050215.xlsx",
                   "output-paper/preconditioned-o-nitrobenzyltyrosine-050215")

# # -------------- Earlier runs --------------
#
# # 5-hydroxytryptophan-041215 had incorrectly transformed amber clones
#
# analyze_microplate("input-paper/5-hydroxytryptophan-041215.xlsx",
#                    "input-paper/key-5-hydroxytryptophan-041215.xlsx",
#                    "output-paper/5-hydroxytryptophan-041215")
#
# # 5-hydroxytryptophan-080315 only ran for 12 hours after starting 2 hours late
#
# analyze_microplate("input-paper/5-hydroxytryptophan-080315.xlsx",
#                    "input-paper/key-5-hydroxytryptophan-080315.xlsx",
#                    "output-paper/5-hydroxytryptophan-080315")
#
# # tyrosine-052115 was rerun later due to lack of confidence in data

analyze_microplate("input-paper/tyrosine-052115.xlsx",
                   "input-paper/key-tyrosine-052115.xlsx",
                   "output-paper/tyrosine-052115")

# # -------------- Including bad wells --------------
#
# analyze_microplate("input-paper/4-azidophenylalanine-080115.xlsx",
#                    "input-paper/key-4-azidophenylalanine-080115-all-replicates-including-bad-well.xlsx",
#                    "output-paper/4-azidophenylalanine-080115-all-replicates-including-bad-well")
#
# analyze_microplate("input-paper/4-azidophenylalanine-preconditioned-080115.xlsx",
#                    "input-paper/key-4-azidophenylalanine-preconditioned-080115-all-replicates-including-bad-well.xlsx",
#                    "output-paper/preconditioned-4-azidophenylalanine-080115-all-replicates-including-bad-well")
#
# # -------------- Highlighted clones --------------
#
# # Clone 1 of 3-aminotyrosine highlighted

analyze_microplate("input-paper/3-aminotyrosine-032615.xlsx",
                   "input-paper/key-3-aminotyrosine-032615-clone-1-highlighted.xlsx",
                   "output-paper/3-aminotyrosine-032615-clone-1-highlighed")

# # Clone 1 of 3-iodotyrosine highlighted

analyze_microplate("input-paper/3-iodotyrosine-032615.xlsx",
                   "input-paper/key-3-iodotyrosine-032615-clone-1-highlighted.xlsx",
                   "output-paper/3-iodotyrosine-032615-clone-1-highlighted")

# # Clone 1 of 3-nitrotyrosine highlighted

analyze_microplate("input-paper/3-nitrotyrosine-041215.xlsx",
                   "input-paper/key-3-nitrotyrosine-041215-clone-1-highlighted.xlsx",
                   "output-paper/3-nitrotyrosine-041215-clone-1-highlighted")

# # Clone 3 of 4-azidophenylalanine highlighted

analyze_microplate("input-paper/4-azidophenylalanine-080115.xlsx",
                   "input-paper/key-4-azidophenylalanine-080115-clone-3-highlight.xlsx",
                   "output-paper/4-azidophenylalanine-080115-clone-3-highlight")

# # -------------- Concentration --------------
#
# # concentration-4-azidophenylalanine-110215 used 5 experimental replicates and 2 blank replicates

analyze_microplate("input-paper/concentration-4-azidophenylalanine-110215.xlsx",
                   "input-paper/key-concentration-1000µM-4-azidophenylalanine-110215.xlsx",
                   "output-paper/concentration-1000µM-4-azidophenylalanine-110215")

analyze_microplate("input-paper/concentration-4-azidophenylalanine-110215.xlsx",
                   "input-paper/key-concentration-320µM-4-azidophenylalanine-110215.xlsx",
                   "output-paper/concentration-320µM-4-azidophenylalanine-110215")

analyze_microplate("input-paper/concentration-4-azidophenylalanine-110215.xlsx",
                   "input-paper/key-concentration-100µM-4-azidophenylalanine-110215.xlsx",
                   "output-paper/concentration-100µM-4-azidophenylalanine-110215")

analyze_microplate("input-paper/concentration-4-azidophenylalanine-110215.xlsx",
                   "input-paper/key-concentration-32µM-4-azidophenylalanine-110215.xlsx",
                   "output-paper/concentration-32µM-4-azidophenylalanine-110215")

analyze_microplate("input-paper/concentration-4-azidophenylalanine-110215.xlsx",
                   "input-paper/key-concentration-10µM-4-azidophenylalanine-110215.xlsx",
                   "output-paper/concentration-10µM-4-azidophenylalanine-110215")

analyze_microplate("input-paper/concentration-4-azidophenylalanine-110215.xlsx",
                   "input-paper/key-concentration-3dot2µM-4-azidophenylalanine-110215.xlsx",
                   "output-paper/concentration-3dot2µM-4-azidophenylalanine-110215")

analyze_microplate("input-paper/concentration-4-azidophenylalanine-110215.xlsx",
                   "input-paper/key-concentration-1µM-4-azidophenylalanine-110215.xlsx",
                   "output-paper/concentration-1µM-4-azidophenylalanine-110215")

summarize_microplates("output-paper/azidophenylalanine_concentration_comparison",
                      c(
                        "output-paper/concentration-1µM-4-azidophenylalanine-110215",
                        "output-paper/concentration-3dot2µM-4-azidophenylalanine-110215",
                        "output-paper/concentration-10µM-4-azidophenylalanine-110215",
                        "output-paper/concentration-32µM-4-azidophenylalanine-110215",
                        "output-paper/concentration-100µM-4-azidophenylalanine-110215",
                        "output-paper/concentration-320µM-4-azidophenylalanine-110215",
                        "output-paper/concentration-1000µM-4-azidophenylalanine-110215"
                      ),
                      OD600.range
)

# # -------------- Polyspecificity --------------

analyze_microplate("input-paper/polyspecificity-iodoYRS-112015.xlsx",
                   "input-paper/key-polyspecificity-aminoy-iodoYRS-112015.xlsx",
                   "output-paper/polyspecificity-aminoy-iodoYRS-112015")

analyze_microplate("input-paper/polyspecificity-iodoYRS-112015.xlsx",
                   "input-paper/key-polyspecificity-azido-iodoYRS-112015.xlsx",
                   "output-paper/polyspecificity-azido-iodoYRS-112015")

analyze_microplate("input-paper/polyspecificity-iodoYRS-112015.xlsx",
                   "input-paper/key-polyspecificity-bromo-iodoYRS-112015.xlsx",
                   "output-paper/polyspecificity-bromo-iodoYRS-112015")

analyze_microplate("input-paper/polyspecificity-iodoYRS-112015.xlsx",
                   "input-paper/key-polyspecificity-chloro-iodoYRS-112015.xlsx",
                   "output-paper/polyspecificity-chloro-iodoYRS-112015")

analyze_microplate("input-paper/polyspecificity-iodoYRS-112015.xlsx",
                   "input-paper/key-polyspecificity-CNF-iodoYRS-112015.xlsx",
                   "output-paper/polyspecificity-CNF-iodoYRS-112015")

analyze_microplate("input-paper/polyspecificity-iodoYRS-112015.xlsx",
                   "input-paper/key-polyspecificity-iodoy-iodoYRS-112015.xlsx",
                   "output-paper/polyspecificity-iodoy-iodoYRS-112015")

analyze_microplate("input-paper/polyspecificity-iodoYRS-112015.xlsx",
                   "input-paper/key-polyspecificity-nitro-iodoYRS-112015.xlsx",
                   "output-paper/polyspecificity-nitro-iodoYRS-112015")

analyze_microplate("input-paper/polyspecificity-iodoYRS-112015.xlsx",
                   "input-paper/key-polyspecificity-tyrosine-iodoYRS-112015.xlsx",
                   "output-paper/polyspecificity-tyrosine-iodoYRS-112015")

summarize_microplates("output-paper/iodoYRS_polyspecificity_comparison",
                      c(
                        "output-paper/polyspecificity-tyrosine-iodoYRS-112015",
                        "output-paper/polyspecificity-iodoy-iodoYRS-112015",
                        "output-paper/polyspecificity-bromo-iodoYRS-112015",
                        "output-paper/polyspecificity-chloro-iodoYRS-112015",
                        "output-paper/polyspecificity-nitro-iodoYRS-112015",
                        "output-paper/polyspecificity-aminoy-iodoYRS-112015",
                        "output-paper/polyspecificity-azido-iodoYRS-112015",
                        "output-paper/polyspecificity-CNF-iodoYRS-112015"
                      ),
                      OD600.range
)

# # -------------- Evolved amberless --------------
#
# # both c321-evolved-without-3-iodotyrosine-081515 and c321-evolved-with-3-iodotyrosine-081515
# # For evolved without, these are really 1 +/-AA replicates of each of 9 transformed clones on the plate but keyed as 3 replicates of 3 clones for continuity

analyze_microplate("input-paper/c321-evolved-without-3-iodotyrosine-081515.xlsx",
"input-paper/key-c321-evolved-without-3-iodotyrosine-081515.xlsx",
"output-paper/c321-evolved-without-3-iodotyrosine-081515")

analyze_microplate("input-paper/c321-evolved-with-3-iodotyrosine-081515.xlsx",
"input-paper/key-c321-evolved-with-3-iodotyrosine-081515.xlsx",
"output-paper/c321-evolved-with-3-iodotyrosine-081515")

summarize_microplates("output-paper/c321_evolved_without_comparison",
                      c(
                        "output-paper/c321-evolved-without-3-iodotyrosine-081515"
                      ),
                      OD600.range
)

summarize_microplates("output-paper/c321_evolved_with_comparison",
                      c(
                        "output-paper/c321-evolved-with-3-iodotyrosine-081515"
                      ),
                      OD600.range
)


