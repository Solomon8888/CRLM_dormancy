# GSE114012 ATF3 evidence summary

## Decision rule

- DEG primary cutoff: P.Value < 0.05 and |logFC| >= 0.5.
- DEG strict cutoff: adj.P.Val < 0.05 and |logFC| >= 0.5.
- TF method support: ATF3 appears within rank <= 10 or passes method P/FDR cutoff where available.
- User primary criterion: DEG primary cutoff plus at least two TF methods supporting ATF3, or at least two consensus TF-intersection candidate lists containing ATF3.
- User strict criterion: DEG strict FDR cutoff plus at least two TF methods supporting ATF3.

## Primary matched analysis designs

- DLD1: logFC=1.815, P=1.885e-05, adj.P=0.001767, TF_methods=4, TF_intersections=1, supported_methods=collectri(rank=35,P=0.00577);dorothea(rank=42,P=0.0391);trrust(rank=43,P=0.0206);enrichr(rank=129,adjP=0.0312)
- HCT15: logFC=1.243, P=5.722e-05, adj.P=0.03472, TF_methods=4, TF_intersections=6, supported_methods=chea3(rank=1);collectri(rank=2,P=4.63e-07);enrichr(rank=5,adjP=4.52e-16);dorothea(rank=30,P=0.00958)
- DLD1_HCT15: logFC=1.517, P=0.002052, adj.P=0.1188, TF_methods=5, TF_intersections=11, supported_methods=chea3(rank=1);enrichr(rank=4,adjP=9.99e-05);collectri(rank=8,P=0.00101);dorothea(rank=9,P=0.00137);trrust(rank=13,P=0.0102)
- DLD1_HCT15_SW48: logFC=1.479, P=0.0001821, adj.P=0.1118, TF_methods=4, TF_intersections=10, supported_methods=chea3(rank=1);collectri(rank=5,P=1.4e-05);enrichr(rank=12,adjP=0.000136);dorothea(rank=22,P=0.0269)
- ALL: logFC=1.227, P=0.0006735, adj.P=0.2578, TF_methods=4, TF_intersections=8, supported_methods=chea3(rank=1);collectri(rank=5,P=5e-06);trrust(rank=10,P=0.00121);enrichr(rank=24,adjP=0.00482)
- SW948: logFC=0.9006, P=0.006486, adj.P=0.241, TF_methods=4, TF_intersections=4, supported_methods=chea3(rank=2);collectri(rank=2,P=1.6e-08);enrichr(rank=35,adjP=0.000155);dorothea(rank=41,P=0.0126)
- SW48: logFC=1.371, P=0.0154, adj.P=0.2427, TF_methods=2, TF_intersections=1, supported_methods=enrichr(rank=37,adjP=6.51e-07);collectri(rank=47,P=0.0107)

## Strict FDR matched analysis designs

- DLD1: logFC=1.815, P=1.885e-05, adj.P=0.001767, TF_methods=4, rank1_methods=
- HCT15: logFC=1.243, P=5.722e-05, adj.P=0.03472, TF_methods=4, rank1_methods=chea3

## Output files

- analysis_design_summary.csv
- atf3_deg_summary.csv
- atf3_tf_method_evidence.csv
- atf3_tf_method_summary.csv
- atf3_tf_intersection_evidence.csv
- atf3_tf_intersection_summary.csv
- atf3_integrated_evidence_summary.csv
- project_result_inventory.csv

