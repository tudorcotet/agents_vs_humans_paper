# Headline statistics - human vs agent (TREM2 hackathon)

**Caveats**
- Selection bias: hit-rate comparisons are conditional on surviving the in-silico ipSAE triage. 41 designs never reached the wet lab. The cohorts had different selection rates: 65 human / 35 agent into the top-100 from 81 / 60 submitted.
- Small N: 65 human, 35 agent in the screened set. Statistical power is limited; report effect sizes alongside p-values.
- Multiple testing: results are exploratory; report unadjusted p-values, note this.
- Replicate noise: for designs with `weird_replicates_flag=True`, the per-design KD has wider effective uncertainty than the SE alone suggests.


## 1. Hit rate x cohort (top-100)
- Human: 25/65 = 38.5% hit
- Agent: 12/35 = 34.3% hit
- Diff (human - agent) = 4.2 pp; bootstrap 95% CI [-14.7, 23.5] pp (n=1000).
- Fisher's exact two-sided: OR=1.20, p=0.8285

## 2. Expression failure x cohort (top-100)
- Human non-expression: 5/65 = 7.7%
- Agent non-expression: 6/35 = 17.1%
- Fisher's exact: OR=0.40, p=0.1857

## 3. pkd by cohort (binders only)
- Human binders n=24, median pkd=7.038 (IQR 6.219 - 7.511)
- Agent binders n=12, median pkd=6.976 (IQR 6.684 - 7.458)
- Mann-Whitney U=134.0, p=0.7499

## 4. Spearman rho: ipSAE (`submitted_ipsae`) vs pkd
- Binders only (n=36): rho=0.397, p=0.01639, 95% CI [0.084, 0.660]
- Top-100 excluding non-binders (n=36): rho=0.397, p=0.01639
- Top-100 left-censored at 10 uM (n=100): rho=0.325, p=0.0009772, 95% CI [0.124, 0.499]

## 5. Sequence length human vs agent (full 141)
- Human n=81, median=85
- Agent n=60, median=80
- Mann-Whitney U=2834.0, p=0.0913

## 6. Method family x cohort (top-100)
```
is_human       False  True 
method_family              
BindCraft          0      2
BoltzGen           0      5
Hallucination      0      8
Mosaic             0      6
Other              7      2
PPIFLOW            0      5
PXDesign          27     11
RFDiffusion        0     26
RFPeptides         1      0
```
chi2=58.807, dof=8, p=7.984e-10.
Cells with chi-square contribution > 1.0:
  - PXDesign x agent: observed 27 vs expected 13.3 (over-represented, contribution 14.11)
  - RFDiffusion x agent: observed 0 vs expected 9.1 (under-represented, contribution 9.10)
  - PXDesign x human: observed 11 vs expected 24.7 (under-represented, contribution 7.60)
  - RFDiffusion x human: observed 26 vs expected 16.9 (over-represented, contribution 4.90)
  - Other x agent: observed 7 vs expected 3.1 (over-represented, contribution 4.71)
  - Hallucination x agent: observed 0 vs expected 2.8 (under-represented, contribution 2.80)
  - Other x human: observed 2 vs expected 5.8 (under-represented, contribution 2.53)
  - Mosaic x agent: observed 0 vs expected 2.1 (under-represented, contribution 2.10)
  - BoltzGen x agent: observed 0 vs expected 1.8 (under-represented, contribution 1.75)
  - PPIFLOW x agent: observed 0 vs expected 1.8 (under-represented, contribution 1.75)
  - Hallucination x human: observed 8 vs expected 5.2 (over-represented, contribution 1.51)
  - RFPeptides x agent: observed 1 vs expected 0.3 (over-represented, contribution 1.21)
  - Mosaic x human: observed 6 vs expected 3.9 (over-represented, contribution 1.13)

## 7. Per-team hit rate (lab-tested designs only)
```
             team  is_human  n_lab  n_binders  n_expressed  hit_rate
            MRAZS      True      8          5            8  0.625000
       1000Tokens      True      7          4            6  0.571429
          GPT 5.2     False      7          4            7  0.571429
claude-sonnet-4.6     False      8          4            8  0.500000
         BraiNSEY      True      5          2            4  0.400000
         EuroBros      True     10          4           10  0.400000
             crow      True     10          4           10  0.400000
    grok-4.1-fast     False      5          2            3  0.400000
         StanFold      True      9          3            9  0.333333
            GLM 5     False      4          1            4  0.250000
         BART.bio      True      5          1            4  0.200000
           NovoFy      True     10          2            9  0.200000
    qwen-3.5-plus     False      5          1            4  0.200000
           DeNovo      True      1          0            0  0.000000
   Gemini 3.1 Pro     False      6          0            3  0.000000
```

