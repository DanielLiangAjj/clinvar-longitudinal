# COSMIC-own-germline sensitivity — cancer vs non-cancer reclassification

Germline-CGC set = 114 CGC v101 genes with GERMLINE=='y' (COSMIC's own annotation).

- Matched to dataset: 114 direct + 0 via synonym = 114/114; unmatched 0 .
- Overlap with paper's 304 CGC: 96/114. Overlap with 106 Hereditary-Cancer: 85/114.

| Gene set | Origin filter | n (exposed) | Crude rate % | Comparator rate % | RR (95% CI) | Crude OR (95% CI) | Adjusted OR (95% CI) |
|---|---|---:|---:|---:|---|---|---|
| COSMIC CGC (304) [baseline] | none (baseline) | 483,822 | 16.619 | 8.191 | 2.029 [2.014, 2.044] | 2.234 [2.215, 2.253] | 1.484 [1.47, 1.499] |
| COSMIC CGC (304) [baseline] | exclude has_somatic | 481,063 | 16.488 | 8.179 | 2.016 [2.001, 2.031] | 2.217 [2.198, 2.236] | 1.479 [1.464, 1.494] |
| Hereditary Cancer (106) | none | 359,743 | 18.182 | 8.336 | 2.181 [2.164, 2.198] | 2.444 [2.421, 2.466] | 1.458 [1.443, 1.474] |
| Hereditary Cancer (106) | exclude has_somatic | 358,054 | 18.046 | 8.319 | 2.169 [2.152, 2.186] | 2.427 [2.404, 2.45] | 1.45 [1.435, 1.466] |
| Germline-CGC (114) | none | 349,547 | 18.329 | 8.351 | 2.195 [2.178, 2.212] | 2.463 [2.44, 2.486] | 1.474 [1.458, 1.49] |
| Germline-CGC (114) | exclude has_somatic | 347,705 | 18.203 | 8.333 | 2.184 [2.167, 2.202] | 2.448 [2.425, 2.471] | 1.468 [1.452, 1.485] |

## Germline-CGC (114) vs TRUE non-cancer (all CGC/Census genes dropped from comparator)

| Gene set | Origin filter | n (exposed) | Crude rate % | Comparator rate % | RR (95% CI) | Crude OR (95% CI) | Adjusted OR (95% CI) |
|---|---|---:|---:|---:|---|---|---|
| Germline-CGC (114) vs true non-cancer | none | 349,547 | 18.329 | 8.19 | 2.238 [2.221, 2.256] | 2.516 [2.492, 2.54] | 1.496 [1.48, 1.513] |
| Germline-CGC (114) vs true non-cancer | exclude has_somatic | 347,705 | 18.203 | 8.178 | 2.226 [2.208, 2.244] | 2.499 [2.475, 2.523] | 1.491 [1.474, 1.507] |

Sanity gate: adjusted 304-CGC baseline reproduced OR=1.484 [1.47, 1.499] (published 1.484 [1.470, 1.499]) — MATCH.
