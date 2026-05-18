# Sequence diversity report

## Within-team pairwise identity (top 5 most similar / 5 most diverse)
```
         team  n_pairs  mean_identity  median_identity  min_identity  max_identity
       DeNovo        1       0.851852         0.851852      0.851852      0.851852
         crow       45       0.703175         0.702381      0.535714      0.845238
       NovoFy       45       0.426285         0.333333      0.290000      0.991379
     BART.bio       45       0.413984         0.212500      0.158333      0.983333
qwen-3.5-plus       45       0.312970         0.300000      0.200000      0.670000
...
         team  n_pairs  mean_identity  median_identity  min_identity  max_identity
     BraiNSEY       28       0.250413         0.197199      0.134454      0.834008
   1000Tokens       45       0.227747         0.213333      0.129412      0.400000
        GLM 5       45       0.220765         0.200000      0.112500      0.812500
grok-4.1-fast       45       0.204236         0.190000      0.060000      0.750000
       Marcel        0            NaN              NaN           NaN           NaN
```

## Cohort-level diversity
- Human cohort: 200 sampled pairs, median pairwise identity 0.214, mean 0.239
- Agent cohort: 200 sampled pairs, median pairwise identity 0.260, mean 0.261
- Mann-Whitney U=15700, p=0.0001999

## K-mer (k=3) Shannon entropy
- Human n=81, median 6.089; agent n=60, median 6.044
- Mann-Whitney U=2832, p=0.09408

## Positive controls (AL002, VHB937)
- Status: PENDING - placeholder sequences only.
- See `data/positive_controls/README.md` for the lookup status; once the patent SEQ IDs are imported the identity table will be regenerated automatically.

