# Method x outcome report

Overall lab hit rate: 37.0% (37/100).

## Methods cross-tab (hit-rate sorted)
```
     design_method_normalized  n_used_total  n_in_top100  n_expressed  n_binders  hit_rate  hit_rate_among_expressed  expression_rate  median_kd_nM_binders  n_human_teams  n_agent_teams
                   Struct_Evo             2            1            1          1  1.000000                  1.000000         1.000000             21.661700              1              0
                PXDesign+MPNN             3            3            3          3  1.000000                  1.000000         1.000000             27.389000              1              0
                       Mosaic             6            6            6          4  0.666667                  0.666667         1.000000             37.062375              1              0
        PXDesign+Protenix+AF2             6            6            5          4  0.666667                  0.800000         0.833333            274.822475              1              0
       AFHall+MPNN+PyR+Boltz2             8            8            8          4  0.500000                  0.500000         1.000000            156.395800              1              0
                    BindCraft             2            2            2          1  0.500000                  0.500000         1.000000              1.908200              1              0
              BoltzGen+Boltz2             2            2            2          1  0.500000                  0.500000         1.000000           1154.571350              2              0
                     PXDesign            34           29           25         13  0.448276                  0.520000         0.862069            137.753367              1              6
      RFDiffusion+MPNN+Boltz2            16           13           13          4  0.307692                  0.307692         1.000000           1073.292600              2              0
 RFDiffusion+LigandMPNN+Boltz            10           10            9          2  0.200000                  0.222222         0.900000            706.005700              1              0
                  MPNN+Boltz2             3            1            0          0  0.000000                       NaN         0.000000                   NaN              1              0
                      Foundry            17            7            6          0  0.000000                  0.000000         0.857143                   NaN              0              5
          PPIFLOW+MPNN+FAMPNN             6            5            5          0  0.000000                  0.000000         1.000000                   NaN              1              0
RFDiffusion3+MPNN+RosettaFold             2            1            0          0  0.000000                       NaN         0.000000                   NaN              1              0
     RFDiffusion2+MPNN+Boltz2             2            2            2          0  0.000000                  0.000000         1.000000                   NaN              1              0
                     BoltzGen            14            3            2          0  0.000000                  0.000000         0.666667                   NaN              4              3
                   RFPeptides             3            1            0          0  0.000000                       NaN         0.000000                   NaN              0              1
                 Protpardelle             1            0            0          0       NaN                       NaN              NaN                   NaN              1              0
        BoltzGen+Boltz2+ipSAE             1            0            0          0       NaN                       NaN              NaN                   NaN              1              0
                      evo+ESM             2            0            0          0       NaN                       NaN              NaN                   NaN              1              0
           AF+PXDesign+Boltz2             1            0            0          0       NaN                       NaN              NaN                   NaN              1              0
```

## Winning methods (hit rate > overall + 1 SE; n_in_top100 >= 3)
```
design_method_normalized  n_in_top100  n_binders  hit_rate  median_kd_nM_binders
           PXDesign+MPNN            3          3  1.000000             27.389000
                  Mosaic            6          4  0.666667             37.062375
   PXDesign+Protenix+AF2            6          4  0.666667            274.822475
```

## Expression-failure-prone methods (expression rate < 80%; n_in_top100 >= 3)
```
design_method_normalized  n_in_top100  n_expressed  expression_rate  n_binders
                BoltzGen            3            2         0.666667          0
```

