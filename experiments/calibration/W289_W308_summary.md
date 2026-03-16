# W289-W308 Sweep Summary

Baseline: W285 (RMSLE=0.666, AK-PWS=19.5%)
Target: RMSLE < 0.60, AK-PWS >= 18%, OR <= 1.0%, CA-N <= 0.5%

| Config | Changes vs W285 | Block |
|--------|-----------------|-------|
| W289 | v_max_warm=0.85 | A:Isolate |
| W290 | v_max_warm=1.0 | A:Isolate |
| W291 | pathogen_adapt_rate=0.0003 | A:Isolate |
| W292 | k_vbnc=4.0 | A:Isolate |
| W293 | T_vbnc_min=11.5 | A:Isolate |
| W294 | wavefront_D_P=150.0 | A:Isolate |
| W295 | P_env_max=4000.0 | A:Isolate |
| W296 | pathogen_adapt_rate=0.0005, v_max_warm=0.85 | B:Pair |
| W297 | k_vbnc=4.0, v_max_warm=0.85 | B:Pair |
| W298 | T_vbnc_min=11.5, pathogen_adapt_rate=0.0003 | B:Pair |
| W299 | T_vbnc_min=12.0, k_vbnc=4.0 | B:Pair |
| W300 | v_max_warm=0.85, wavefront_D_P=200.0 | B:Pair |
| W301 | cumulative_dose_threshold=2000.0, wavefront_D_P=150.0 | B:Pair |
| W302 | k_vbnc=4.0, pathogen_adapt_rate=0.0003, v_max_warm=0.85 | C:Triple |
| W303 | T_vbnc_min=11.5, k_vbnc=4.0, v_max_warm=0.85 | C:Triple |
| W304 | pathogen_adapt_rate=0.0003, v_max_warm=0.85, wavefront_D_P=200.0 | C:Triple |
| W305 | T_vbnc_min=11.5, v_max_warm=0.85, fw_lat_min=46.0, fw_lat_max=63.0 | C:Triple |
| W306 | T_vbnc_min=11.5, k_vbnc=4.0, pathogen_adapt_rate=0.0003, v_max_warm=0.85 | D:Best |
| W307 | T_vbnc_min=12.0, k_vbnc=4.5, pathogen_adapt_rate=0.0002, v_max_warm=0.9 | D:Best |
| W308 | K_half=1800000, T_vbnc_min=11.5, k_vbnc=4.0, pathogen_adapt_rate=0.0003, v_max_warm=0.85, wavefront_D_P=200.0, fw_lat_min=46.0, fw_lat_max=63.0 | D:Best |
