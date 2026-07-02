# GSE114012 ATF3 TF排序与交集策略

## 单方法排序检查

`single_method_atf3_rank_by_sorting.csv` 记录了每种TF方法在不同排序口径下ATF3的排名。
推荐主排序口径与 `09_integrate_tf_enrichment_results.R` 中生成method_final结果时使用的排序保持一致。

可以使ATF3排第1的单方法排序口径如下：

- ALL / ChEA3 / current_09_rank: top TF = ATF3
- ALL / ChEA3 / library_count_desc_toprank_asc: top TF = ATF3
- ALL / ChEA3 / score_ascending: top TF = ATF3
- DLD1_HCT15 / ChEA3 / current_09_rank: top TF = ATF3
- DLD1_HCT15 / ChEA3 / library_count_desc_toprank_asc: top TF = ATF3
- DLD1_HCT15 / ChEA3 / score_ascending: top TF = ATF3
- DLD1_HCT15 / ENRICHR / best_adjusted_p_asc: top TF = ATF3
- DLD1_HCT15 / ENRICHR / best_p_value_asc: top TF = ATF3
- DLD1_HCT15 / ENRICHR / library_count_desc_adj_p_asc: top TF = ATF3
- DLD1_HCT15_SW48 / ChEA3 / current_09_rank: top TF = ATF3
- DLD1_HCT15_SW48 / ChEA3 / library_count_desc_toprank_asc: top TF = ATF3
- DLD1_HCT15_SW48 / ChEA3 / score_ascending: top TF = ATF3
- HCT15 / ChEA3 / current_09_rank: top TF = ATF3
- HCT15 / ChEA3 / library_count_desc_toprank_asc: top TF = ATF3
- HCT15 / ChEA3 / score_ascending: top TF = ATF3
- HCT15 / ENRICHR / best_adjusted_p_asc: top TF = ATF3
- HCT15 / ENRICHR / best_p_value_asc: top TF = ATF3
- HCT15 / ENRICHR / library_count_desc_adj_p_asc: top TF = ATF3
- HT55 / ChEA3 / library_count_desc_toprank_asc: top TF = ATF3
- SW948 / ChEA3 / library_count_desc_toprank_asc: top TF = ATF3

## 基于官方method_final排名的最佳交集策略

- 覆盖度最高的推荐方法组合：DoRothEA;ChEA3;CollecTRI。
- ATF3排第1的差异分析方案：5/9（ALL;DLD1_HCT15;DLD1_HCT15_SW48;HCT15;SW948）。
- 在ATF3表达达到主显著阈值的DEG方案中，ATF3排第1的方案数：5（ALL;DLD1_HCT15;DLD1_HCT15_SW48;HCT15;SW948）。
- ATF3进入top10的差异分析方案：7/9。

可复现的交集排序规则：

1. 使用 `09_integrate_tf_enrichment_results.R` 生成的各方法 `method_final` 表。
2. ChEA3使用Integrated--topRank的rank升序。
3. DoRothEA使用ORA P值升序，其次ORA score降序。
4. CollecTRI使用activity P值升序，其次绝对activity score降序。
5. 对入选方法的完整TF列表取交集，不预先限制top-N。
6. 对每个交集TF计算其在所选方法中的rank。
7. 按 Mean_Selected_Rank升序、Best_Selected_Rank升序、Source_Method_Count降序、CheA3_Library_Count降序、TF名称升序排序。

推荐策略逐方案结果：

- ALL: consensus rank= 1, mean selected rank=15, details=DoRothEA(rank=39, P=0.15); ChEA3(rank=1); CollecTRI(rank=5, P=5e-06)
- DLD1_HCT15: consensus rank= 1, mean selected rank=6, details=DoRothEA(rank=9, P=0.00137); ChEA3(rank=1); CollecTRI(rank=8, P=0.00101)
- DLD1_HCT15_SW48: consensus rank= 1, mean selected rank=9.333, details=DoRothEA(rank=22, P=0.0269); ChEA3(rank=1); CollecTRI(rank=5, P=1.4e-05)
- HCT15: consensus rank= 1, mean selected rank=11, details=DoRothEA(rank=30, P=0.00958); ChEA3(rank=1); CollecTRI(rank=2, P=4.63e-07)
- SW948: consensus rank= 1, mean selected rank=15, details=DoRothEA(rank=41, P=0.0126); ChEA3(rank=2); CollecTRI(rank=2, P=1.6e-08)
- HT55: consensus rank= 5, mean selected rank=33.67, details=DoRothEA(rank=28, P=0.0132); ChEA3(rank=4); CollecTRI(rank=69, P=0.00385)
- SW48: consensus rank=10, mean selected rank=81.33, details=DoRothEA(rank=75, P=0.102); ChEA3(rank=122); CollecTRI(rank=47, P=0.0107)
- DLD1: consensus rank=14, mean selected rank=107, details=DoRothEA(rank=42, P=0.0391); ChEA3(rank=244); CollecTRI(rank=35, P=0.00577)
- RKO: consensus rank=28, mean selected rank=164.3, details=DoRothEA(rank=138, P=0.185); ChEA3(rank=171); CollecTRI(rank=184, P=0.0825)

## 同等rank1覆盖度下的最小方法组合

- 最小同覆盖度组合：ChEA3;CollecTRI。
- 该组合同样能使ATF3在 5/9 个方案中排第1，并在 6/9 个方案中进入top10。

## 输出文件

- single_method_atf3_rank_by_sorting.csv
- single_method_atf3_best_sorting.csv
- single_method_atf3_top1_sorting_options.csv
- tf_intersection_combo_search_official_rank_detail.csv
- tf_intersection_combo_search_official_rank_summary.csv
- recommended_tf_intersection_strategy_detail.csv
- minimal_same_top1_intersection_options.csv

