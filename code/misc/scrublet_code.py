#!/usr/bin/env python3

import scrublet as scr
import scipy.io
import sys

patient=sys.argv[1]
doublet_rate=float(sys.argv[2])

print('Patient:', patient, 'Doublet_rate:', doublet_rate)
print(sys.argv)

counts_matrix = scipy.io.mmread(''.join(['data/', patient, '/matrix_', patient, '_raw.mtx'])).T.tocsc()
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=doublet_rate)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=1, min_cells=1, min_gene_variability_pctl=50, n_prin_comps=20)
print(doublet_scores, predicted_doublets)

outFI = open(''.join(['data/', patient, '/doublets_', patient, '_raw.txt']),'w')
outFI.write('predicted_doublets\tdoublet_scores\n')
for i in range(len(predicted_doublets)):
    outFI.write(str(predicted_doublets[i]).upper()+'\t'+str(doublet_scores[i])+'\n')

outFI.close()
