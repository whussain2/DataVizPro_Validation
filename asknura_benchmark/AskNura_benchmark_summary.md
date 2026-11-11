# AskNura benchmark: summary (run 3 Nov 2026, live DataViz Pro v1.2, Researcher workspace, admin account)

44 tasks x 3 repeats = 132 scored runs, each in a fresh session. Requests and approvals entered by an automated browser agent (Claude) following a fixed script. Runs were scored against predefined checks by the browser agent; every score is listed in AskNura_benchmark_scores_all.csv.

| Measure | Result |
|---|---|
| Completed and correct | 116 / 132 = 87.9% (95% Wilson CI 81.2-92.4%) |
| Partial | 10 (7.6%) |
| Failed | 6 (4.5%) |
| Invented feature, setting or number | 0 / 132 (95% CI 0-2.8%) |
| Median time per task | 66.5 s (IQR 34.5-94.5) |
| Median approvals per task | 2 (IQR 1-2) |
| Tasks with identical settings and data in all 3 repeats | 28 / 44 (64%); 17 / 32 figure tasks (53%) |

By category (correct / runs): data preparation 22/24 (92%); basic plots 18/24 (75%); publication styling 13/18 (72%); multivariate and statistics 15/15 (100%); genomics 18/21 (86%); field-spatial 6/6 (100%); interpretation and knowledge 12/12 (100%); honesty and safety 12/12 (100%).

Recurring partial/failed cases:
- B5 per-replication regression lines (3/3 partial): the scatter module draws one trend line unless the plot is faceted.
- B8 sorted dot plot (3/3 partial): the dot plot cannot sort; AskNura used the lollipop chart instead and said so.
- C5 legend to bottom (3/3 failed): the box plot has no legend until a fill variable is set; AskNura explained and asked which variable.
- C4 colour-blind palette (2/3 partial): asked which variable to colour by before applying.
- A5 transpose (1/3 partial), A7 and E4 (1 each failed: stopped mid-task), E7 (1 failed: app memory limit).

Not counted (harness or environment faults, rerun): 85 attempts set aside: drawer not open on reload (no message sent, no cost), page in background or session idle (messages sent; these used AskNura credit).

Cost and models: from the gateway usage report for 3 Nov (to add).

## Cost and models (gateway usage store, days 2026-11-03 and 2026-11-04 in the gateway's IST calendar = the benchmark window)
| Model | Calls | Cost (US$) | Mean response (ms) |
|---|---|---|---|
| claude-sonnet-5-5 | 1,774 | 12.5360 | ~2,280 |
| claude-opus-5-5 | 218 | 4.1605 | ~3,530 |
| Total | 1,992 | 16.6965 | |
Includes set-up requests (30 runs), discarded attempts and Waseem's own use on those days, so US$0.13 per scored run (16.70/132) is an upper bound.
