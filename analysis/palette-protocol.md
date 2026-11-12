# DNA-PAINT pair evaluation

The endpoint is NeNA in nanometres on the deposited fixed-cell nuclear pore complex (NPC) measurements. Lower is better. Choose two dyes within each of the 488, 560 and 640 nm channels. All 18 dyes are eligible in their respective channels. This retrospective benchmark was specified before processing the full archive and inspecting the holdout analysis. The paper and its summary workbooks were already public.

## Split and decision rule

1. For each dye and specimen, average its NeNA estimates across the available **training** repeats. Rank dyes separately inside each channel and specimen, best first. Tied measurements receive the mean of their ordinal ranks.
2. A candidate's score is ​`(1 − w) × NPC rank + w × origami rank`. The two smallest scores win. Break remaining ties by the archive dye identifier in ascending order.
3. Using only repeat 1 as input, predict NPC repeat 2 for weights 0, 0.5 and 1. Within a channel compute mean NPC NeNA for the selected pair minus the mean for the two lowest observed dyes in that evaluation repeat. Average regret equally across the three channels. Select the smallest weight among tied minimum regrets, to 12 decimal places.
4. Freeze the weight. Pool repeats 1 and 2 within each dye and specimen and fix the three pairs. The fixed historical baseline uses precisely the same training repeats, dyes and channels with weight 0. Score both policies on NPC repeat 3. The held-out origami measurements and the fourth Cy3B cell acquisition never enter selection or evaluation.
5. Display prospective pair suggestions only if macro-average regret is strictly lower than the baseline and at least two channels strictly improve. Compare regret to a numerical tolerance of 1e−10 nm. Otherwise display the observed results with the recommendation withheld. A missing required measurement prevents scoring the complete panel. The published-sheet replay has no influence on the raw analysis or this decision.

## Interpretation

Each channel contributes one outcome. Replicate labels follow the authors’ directory and sheet order. Neither random field splits nor individual localizations supply independent validation studies. Distinct fluorophores can share a Cy3B reference acquisition in the specificity calculation. The primary NeNA endpoint uses each dye’s own acquisition. The paper's processed localization tables are the starting material. There is no adjustment for instruments, conjugates or sampling conditions outside this experiment.

The test measures how well two choices transfer to one historical NPC repeat. It offers no estimate of a 20 nm event yield for H2B–HaloTag. The workbook replay and per-replicate deviations are available alongside the raw-table evaluation. [Steen et al.](https://doi.org/10.1038/s41592-024-02374-8) · [Deposited localization tables](https://zenodo.org/records/10960858)
