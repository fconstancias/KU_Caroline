---
title: " Study II - Importing metaphlan profiles into phyloseq"
author: "Florentin Constancias"
date: "November 20, 2024"
output: 
  html_document: 
    toc: yes
    keep_md: yes
---









# Import and get ready:


``` r
here::here("../../data/processed_data/metaphlan/merged_abundance_table.txt") %>% 
  metaphlan_2phyloseq(merged_metaphlan = .) -> ps

ps %>% 
  sample_names() %>% 
  head()
```

```
## [1] "X11563A_trimmed_HF" "X11564A_trimmed_HF" "X11566A_trimmed_HF"
## [4] "X11567A_trimmed_HF" "X11568A_trimmed_HF" "X11576A_trimmed_HF"
```


``` r
ps %>% 
  sample_names() %>% 
  str_extract(.,"[^_]+") %>% 
  str_remove_all(., "[A-Z]")  -> sample_names(ps)
```


``` r
here::here("../../data/processed_data/metaphlan/Metadata_Deerland_Study.xlsx") %>% 
  read_excel(sheet = "Trial2") -> meta

meta %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item" id="htmlwidget-39d925b953787dfc4931" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-39d925b953787dfc4931">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108","109","110","111","112","113","114","115","116","117","118","119","120","121","122","123","124","125","126","127","128","129","130","131","132","133","134","135","136","137","138","139","140","141","142","143","144","145","146","147","148","149","150","151","152","153","154","155","156","157","158","159","160","161","162","163","164","165","166","167","168","169","170","171","172","173","174","175","176","177","178","179","180","181","182","183","184","185","186","187","188","189","190","191","192","193","194","195","196","197","198","199","200","201","202","203","204","205","206","207","208","209","210","211","212","213","214","215","216","217","218","219","220","221","222","223","224","225","226","227","228","229","230","231","232","233","234","235","236","237","238","239"],[11563,11564,11566,11567,11568,11576,11578,11580,11581,11582,11583,11585,11587,11591,11593,11597,11598,11600,11601,11602,11603,11604,11606,11610,11612,11613,11614,11616,11617,11619,11623,11625,11627,11629,11636,11637,11638,11639,11640,11642,11643,11644,11646,11647,11648,11656,11658,11660,11661,11662,11663,11665,11667,11671,11673,11677,11678,11680,11681,11682,11683,11684,11686,11690,11692,11693,11694,11696,11697,11699,11703,11705,11707,11709,11716,11717,11718,11719,11720,11722,11723,11724,11726,11727,11728,11736,11738,11740,11741,11742,11743,11745,11747,11752,11754,11758,11759,11761,11763,11764,11765,11767,11771,11773,11774,11775,11777,11778,11780,11784,11786,11788,11790,11797,11798,11799,11800,11803,11805,11806,11807,11809,11810,11811,11819,11821,11823,11824,11825,11826,11828,11830,11834,11836,11840,11841,11843,11844,11845,11846,11847,11849,11853,11855,11856,11857,11859,11860,11862,11866,11868,11870,11872,11879,11880,11881,11882,11883,11885,11886,11887,11889,11890,11891,11899,11901,11903,11904,11905,11906,11908,11910,11914,11916,11920,11921,11923,11924,11925,11926,11927,11929,11933,11935,11936,11937,11939,11940,11942,11946,11948,11950,11952,11959,11960,11961,11962,11963,11965,11966,11967,11969,11970,11971,11979,11981,11983,11984,11985,11986,11988,11990,11994,11996,12000,12001,12003,12004,12005,12006,12007,12009,12013,12015,12016,12017,12019,12020,12022,12026,12028,12030,12032,12039,12040,12041,12042,12043,12045],["21B1","21B2","21B4","21B5","21B6","21B14","21B16","21B18","21B19","21B20","21B21","21B23","21B25","21B29","21B31","21B35","21B36","21B38","21B39","21B40","21B41","21B42","21B44","21B48","21B50","21B51","21B52","21B54","21B55","21B57","21B61","21B63","21B65","21B67","21B74","21B75","21B76","21B77","21B78","21B80","22B1","22B2","22B4","22B5","22B6","22B14","22B16","22B18","22B19","22B20","22B21","22B23","22B25","22B29","22B31","22B35","22B36","22B38","22B39","22B40","22B41","22B42","22B44","22B48","22B50","22B51","22B52","22B54","22B55","22B57","22B61","22B63","22B65","22B67","22B74","22B75","22B76","22B77","22B78","22B80","23B1","23B2","23B4","23B5","23B6","23B14","23B16","23B18","23B19","23B20","23B21","23B23","23B25","23B29","23B31","23B35","23B36","23B38","23B40","23B41","23B42","23B44","23B48","23B50","23B51","23B52","23B54","23B55","23B57","23B61","23B63","23B65","23B67","23B74","23B75","23B76","23B77","23B78","23B80","21F1","21F2","21F4","21F5","21F6","21F14","21F16","21F18","21F19","21F20","21F21","21F23","21F25","21F29","21F31","21F35","21F36","21F38","21F39","21F40","21F41","21F42","21F44","21F48","21F50","21F51","21F52","21F54","21F55","21F57","21F61","21F63","21F65","21F67","21F74","21F75","21F76","21F77","21F78","21F80","22F1","22F2","22F4","22F5","22F6","22F14","22F16","22F18","22F19","22F20","22F21","22F23","22F25","22F29","22F31","22F35","22F36","22F38","22F39","22F40","22F41","22F42","22F44","22F48","22F50","22F51","22F52","22F54","22F55","22F57","22F61","22F63","22F65","22F67","22F74","22F75","22F76","22F77","22F78","22F80","23F1","23F2","23F4","23F5","23F6","23F14","23F16","23F18","23F19","23F20","23F21","23F23","23F25","23F29","23F31","23F35","23F36","23F38","23F39","23F40","23F41","23F42","23F44","23F48","23F50","23F51","23F52","23F54","23F55","23F57","23F61","23F63","23F65","23F67","23F74","23F75","23F76","23F77","23F78","23F80"],["Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva"],["TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3"],["D01","D02","D04","D05","D06","D14","D16","D18","D19","D20","D21","D23","D25","D29","D31","D35","D36","D38","D39","D40","D41","D42","D44","D48","D50","D51","D52","D54","D55","D57","D61","D63","D65","D67","D74","D75","D76","D77","D78","D80","D01","D02","D04","D05","D06","D14","D16","D18","D19","D20","D21","D23","D25","D29","D31","D35","D36","D38","D39","D40","D41","D42","D44","D48","D50","D51","D52","D54","D55","D57","D61","D63","D65","D67","D74","D75","D76","D77","D78","D80","D01","D02","D04","D05","D06","D14","D16","D18","D19","D20","D21","D23","D25","D29","D31","D35","D36","D38","D40","D41","D42","D44","D48","D50","D51","D52","D54","D55","D57","D61","D63","D65","D67","D74","D75","D76","D77","D78","D80","D01","D02","D04","D05","D06","D14","D16","D18","D19","D20","D21","D23","D25","D29","D31","D35","D36","D38","D39","D40","D41","D42","D44","D48","D50","D51","D52","D54","D55","D57","D61","D63","D65","D67","D74","D75","D76","D77","D78","D80","D01","D02","D04","D05","D06","D14","D16","D18","D19","D20","D21","D23","D25","D29","D31","D35","D36","D38","D39","D40","D41","D42","D44","D48","D50","D51","D52","D54","D55","D57","D61","D63","D65","D67","D74","D75","D76","D77","D78","D80","D01","D02","D04","D05","D06","D14","D16","D18","D19","D20","D21","D23","D25","D29","D31","D35","D36","D38","D39","D40","D41","D42","D44","D48","D50","D51","D52","D54","D55","D57","D61","D63","D65","D67","D74","D75","D76","D77","D78","D80"],["Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo"],[1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500],["21B1_NatriumChlorid0.9Study2_Plaque_TP1_P1Placebo_1000","21B2_NatriumChlorid0.9Study2_Plaque_TP1_P2Placebo_1000","21B4_NatriumChlorid0.9Study2_Plaque_TP1_P4Placebo_1000","21B5_NatriumChlorid0.9Study2_Plaque_TP1_P5Placebo_1000","21B6_NatriumChlorid0.9Study2_Plaque_TP1_P6Placebo_1000","21B14_NatriumChlorid0.9Study2_Plaque_TP1_P14Placebo_1000","21B16_NatriumChlorid0.9Study2_Plaque_TP1_P16Placebo_1000","21B18_NatriumChlorid0.9Study2_Plaque_TP1_P18Placebo_1000","21B19_NatriumChlorid0.9Study2_Plaque_TP1_P19Placebo_1000","21B20_NatriumChlorid0.9Study2_Plaque_TP1_P20Placebo_1000","21B21_NatriumChlorid0.9Study2_Plaque_TP1_P21Placebo_1000","21B23_NatriumChlorid0.9Study2_Plaque_TP1_P23Placebo_1000","21B25_NatriumChlorid0.9Study2_Plaque_TP1_P25Placebo_1000","21B29_NatriumChlorid0.9Study2_Plaque_TP1_P29Placebo_1000","21B31_NatriumChlorid0.9Study2_Plaque_TP1_P31Placebo_1000","21B35_NatriumChlorid0.9Study2_Plaque_TP1_P35Placebo_1000","21B36_NatriumChlorid0.9Study2_Plaque_TP1_P36Placebo_1000","21B38_NatriumChlorid0.9Study2_Plaque_TP1_P38Placebo_1000","21B39_NatriumChlorid0.9Study2_Plaque_TP1_P39Placebo_1000","21B40_NatriumChlorid0.9Study2_Plaque_TP1_P40Placebo_1000","21B41_NatriumChlorid0.9Study2_Plaque_TP1_P41Placebo_1000","21B42_NatriumChlorid0.9Study2_Plaque_TP1_P42Placebo_1000","21B44_NatriumChlorid0.9Study2_Plaque_TP1_P44Placebo_1000","21B48_NatriumChlorid0.9Study2_Plaque_TP1_P48Placebo_1000","21B50_NatriumChlorid0.9Study2_Plaque_TP1_P50Placebo_1000","21B51_NatriumChlorid0.9Study2_Plaque_TP1_P51Placebo_1000","21B52_NatriumChlorid0.9Study2_Plaque_TP1_P52Placebo_1000","21B54_NatriumChlorid0.9Study2_Plaque_TP1_P54Placebo_1000","21B55_NatriumChlorid0.9Study2_Plaque_TP1_P55Placebo_1000","21B57_NatriumChlorid0.9Study2_Plaque_TP1_P57Placebo_1000","21B61_NatriumChlorid0.9Study2_Plaque_TP1_P61Placebo_1000","21B63_NatriumChlorid0.9Study2_Plaque_TP1_P63Placebo_1000","21B65_NatriumChlorid0.9Study2_Plaque_TP1_P65Placebo_1000","21B67_NatriumChlorid0.9Study2_Plaque_TP1_P67Placebo_1000","21B74_NatriumChlorid0.9Study2_Plaque_TP1_P74Placebo_1000","21B75_NatriumChlorid0.9Study2_Plaque_TP1_P75Placebo_1000","21B76_NatriumChlorid0.9Study2_Plaque_TP1_P76Placebo_1000","21B77_NatriumChlorid0.9Study2_Plaque_TP1_P77Placebo_1000","21B78_NatriumChlorid0.9Study2_Plaque_TP1_P78Placebo_1000","21B80_NatriumChlorid0.9Study2_Plaque_TP1_P80Placebo_1000","22B1_NatriumChlorid0.9Study2_Plaque_TP2_P1Placebo_1000","22B2_NatriumChlorid0.9Study2_Plaque_TP2_P2Placebo_1000","22B4_NatriumChlorid0.9Study2_Plaque_TP2_P4Placebo_1000","22B5_NatriumChlorid0.9Study2_Plaque_TP2_P5Placebo_1000","22B6_NatriumChlorid0.9Study2_Plaque_TP2_P6Placebo_1000","22B14_NatriumChlorid0.9Study2_Plaque_TP2_P14Placebo_1000","22B16_NatriumChlorid0.9Study2_Plaque_TP2_P16Placebo_1000","22B18_NatriumChlorid0.9Study2_Plaque_TP2_P18Placebo_1000","22B19_NatriumChlorid0.9Study2_Plaque_TP2_P19Placebo_1000","22B20_NatriumChlorid0.9Study2_Plaque_TP2_P20Placebo_1000","22B21_NatriumChlorid0.9Study2_Plaque_TP2_P21Placebo_1000","22B23_NatriumChlorid0.9Study2_Plaque_TP2_P23Placebo_1000","22B25_NatriumChlorid0.9Study2_Plaque_TP2_P25Placebo_1000","22B29_NatriumChlorid0.9Study2_Plaque_TP2_P29Placebo_1000","22B31_NatriumChlorid0.9Study2_Plaque_TP2_P31Placebo_1000","22B35_NatriumChlorid0.9Study2_Plaque_TP2_P35Placebo_1000","22B36_NatriumChlorid0.9Study2_Plaque_TP2_P36Placebo_1000","22B38_NatriumChlorid0.9Study2_Plaque_TP2_P38Placebo_1000","22B39_NatriumChlorid0.9Study2_Plaque_TP2_P39Placebo_1000","22B40_NatriumChlorid0.9Study2_Plaque_TP2_P40Placebo_1000","22B41_NatriumChlorid0.9Study2_Plaque_TP2_P41Placebo_1000","22B42_NatriumChlorid0.9Study2_Plaque_TP2_P42Placebo_1000","22B44_NatriumChlorid0.9Study2_Plaque_TP2_P44Placebo_1000","22B48_NatriumChlorid0.9Study2_Plaque_TP2_P48Placebo_1000","22B50_NatriumChlorid0.9Study2_Plaque_TP2_P50Placebo_1000","22B51_NatriumChlorid0.9Study2_Plaque_TP2_P51Placebo_1000","22B52_NatriumChlorid0.9Study2_Plaque_TP2_P52Placebo_1000","22B54_NatriumChlorid0.9Study2_Plaque_TP2_P54Placebo_1000","22B55_NatriumChlorid0.9Study2_Plaque_TP2_P55Placebo_1000","22B57_NatriumChlorid0.9Study2_Plaque_TP2_P57Placebo_1000","22B61_NatriumChlorid0.9Study2_Plaque_TP2_P61Placebo_1000","22B63_NatriumChlorid0.9Study2_Plaque_TP2_P63Placebo_1000","22B65_NatriumChlorid0.9Study2_Plaque_TP2_P65Placebo_1000","22B67_NatriumChlorid0.9Study2_Plaque_TP2_P67Placebo_1000","22B74_NatriumChlorid0.9Study2_Plaque_TP2_P74Placebo_1000","22B75_NatriumChlorid0.9Study2_Plaque_TP2_P75Placebo_1000","22B76_NatriumChlorid0.9Study2_Plaque_TP2_P76Placebo_1000","22B77_NatriumChlorid0.9Study2_Plaque_TP2_P77Placebo_1000","22B78_NatriumChlorid0.9Study2_Plaque_TP2_P78Placebo_1000","22B80_NatriumChlorid0.9Study2_Plaque_TP2_P80Placebo_1000","23B1_NatriumChlorid0.9Study2_Plaque_TP3_P1Placebo_1000","23B2_NatriumChlorid0.9Study2_Plaque_TP3_P2Placebo_1000","23B4_NatriumChlorid0.9Study2_Plaque_TP3_P4Placebo_1000","23B5_NatriumChlorid0.9Study2_Plaque_TP3_P5Placebo_1000","23B6_NatriumChlorid0.9Study2_Plaque_TP3_P6Placebo_1000","23B14_NatriumChlorid0.9Study2_Plaque_TP3_P14Placebo_1000","23B16_NatriumChlorid0.9Study2_Plaque_TP3_P16Placebo_1000","23B18_NatriumChlorid0.9Study2_Plaque_TP3_P18Placebo_1000","23B19_NatriumChlorid0.9Study2_Plaque_TP3_P19Placebo_1000","23B20_NatriumChlorid0.9Study2_Plaque_TP3_P20Placebo_1000","23B21_NatriumChlorid0.9Study2_Plaque_TP3_P21Placebo_1000","23B23_NatriumChlorid0.9Study2_Plaque_TP3_P23Placebo_1000","23B25_NatriumChlorid0.9Study2_Plaque_TP3_P25Placebo_1000","23B29_NatriumChlorid0.9Study2_Plaque_TP3_P29Placebo_1000","23B31_NatriumChlorid0.9Study2_Plaque_TP3_P31Placebo_1000","23B35_NatriumChlorid0.9Study2_Plaque_TP3_P35Placebo_1000","23B36_NatriumChlorid0.9Study2_Plaque_TP3_P36Placebo_1000","23B38_NatriumChlorid0.9Study2_Plaque_TP3_P38Placebo_1000","23B40_NatriumChlorid0.9Study2_Plaque_TP3_P40Placebo_1000","23B41_NatriumChlorid0.9Study2_Plaque_TP3_P41Placebo_1000","23B42_NatriumChlorid0.9Study2_Plaque_TP3_P42Placebo_1000","23B44_NatriumChlorid0.9Study2_Plaque_TP3_P44Placebo_1000","23B48_NatriumChlorid0.9Study2_Plaque_TP3_P48Placebo_1000","23B50_NatriumChlorid0.9Study2_Plaque_TP3_P50Placebo_1000","23B51_NatriumChlorid0.9Study2_Plaque_TP3_P51Placebo_1000","23B52_NatriumChlorid0.9Study2_Plaque_TP3_P52Placebo_1000","23B54_NatriumChlorid0.9Study2_Plaque_TP3_P54Placebo_1000","23B55_NatriumChlorid0.9Study2_Plaque_TP3_P55Placebo_1000","23B57_NatriumChlorid0.9Study2_Plaque_TP3_P57Placebo_1000","23B61_NatriumChlorid0.9Study2_Plaque_TP3_P61Placebo_1000","23B63_NatriumChlorid0.9Study2_Plaque_TP3_P63Placebo_1000","23B65_NatriumChlorid0.9Study2_Plaque_TP3_P65Placebo_1000","23B67_NatriumChlorid0.9Study2_Plaque_TP3_P67Placebo_1000","23B74_NatriumChlorid0.9Study2_Plaque_TP3_P74Placebo_1000","23B75_NatriumChlorid0.9Study2_Plaque_TP3_P75Placebo_1000","23B76_NatriumChlorid0.9Study2_Plaque_TP3_P76Placebo_1000","23B77_NatriumChlorid0.9Study2_Plaque_TP3_P77Placebo_1000","23B78_NatriumChlorid0.9Study2_Plaque_TP3_P78Placebo_1000","23B80_NatriumChlorid0.9Study2_Plaque_TP3_P80Placebo_1000","21F1_Study2_Saliva_TP1_P1Placebo_500","21F2_Study2_Saliva_TP1_P2Placebo_500","21F4_Study2_Saliva_TP1_P4Placebo_500","21F5_Study2_Saliva_TP1_P5Placebo_500","21F6_Study2_Saliva_TP1_P6Placebo_500","21F14_Study2_Saliva_TP1_P14Placebo_500","21F16_Study2_Saliva_TP1_P16Placebo_500","21F18_Study2_Saliva_TP1_P18Placebo_500","21F19_Study2_Saliva_TP1_P19Placebo_500","21F20_Study2_Saliva_TP1_P20Placebo_500","21F21_Study2_Saliva_TP1_P21Placebo_500","21F23_Study2_Saliva_TP1_P23Placebo_500","21F25_Study2_Saliva_TP1_P25Placebo_500","21F29_Study2_Saliva_TP1_P29Placebo_500","21F31_Study2_Saliva_TP1_P31Placebo_500","21F35_Study2_Saliva_TP1_P35Placebo_500","21F36_Study2_Saliva_TP1_P36Placebo_500","21F38_Study2_Saliva_TP1_P38Placebo_500","21F39_Study2_Saliva_TP1_P39Placebo_500","21F40_Study2_Saliva_TP1_P40Placebo_500","21F41_Study2_Saliva_TP1_P41Placebo_500","21F42_Study2_Saliva_TP1_P42Placebo_500","21F44_Study2_Saliva_TP1_P44Placebo_500","21F48_Study2_Saliva_TP1_P48Placebo_500","21F50_Study2_Saliva_TP1_P50Placebo_500","21F51_Study2_Saliva_TP1_P51Placebo_500","21F52_Study2_Saliva_TP1_P52Placebo_500","21F54_Study2_Saliva_TP1_P54Placebo_500","21F55_Study2_Saliva_TP1_P55Placebo_500","21F57_Study2_Saliva_TP1_P57Placebo_500","21F61_Study2_Saliva_TP1_P61Placebo_500","21F63_Study2_Saliva_TP1_P63Placebo_500","21F65_Study2_Saliva_TP1_P65Placebo_500","21F67_Study2_Saliva_TP1_P67Placebo_500","21F74_Study2_Saliva_TP1_P74Placebo_500","21F75_Study2_Saliva_TP1_P75Placebo_500","21F76_Study2_Saliva_TP1_P76Placebo_500","21F77_Study2_Saliva_TP1_P77Placebo_500","21F78_Study2_Saliva_TP1_P78Placebo_500","21F80_Study2_Saliva_TP1_P80Placebo_500","22F1_Study2_Saliva_TP2_P1Placebo_500","22F2_Study2_Saliva_TP2_P2Placebo_500","22F4_Study2_Saliva_TP2_P4Placebo_500","22F5_Study2_Saliva_TP2_P5Placebo_500","22F6_Study2_Saliva_TP2_P6Placebo_500","22F14_Study2_Saliva_TP2_P14Placebo_500","22F16_Study2_Saliva_TP2_P16Placebo_500","22F18_Study2_Saliva_TP2_P18Placebo_500","22F19_Study2_Saliva_TP2_P19Placebo_500","22F20_Study2_Saliva_TP2_P20Placebo_500","22F21_Study2_Saliva_TP2_P21Placebo_500","22F23_Study2_Saliva_TP2_P23Placebo_500","22F25_Study2_Saliva_TP2_P25Placebo_500","22F29_Study2_Saliva_TP2_P29Placebo_500","22F31_Study2_Saliva_TP2_P31Placebo_500","22F35_Study2_Saliva_TP2_P35Placebo_500","22F36_Study2_Saliva_TP2_P36Placebo_500","22F38_Study2_Saliva_TP2_P38Placebo_500","22F39_Study2_Saliva_TP2_P39Placebo_500","22F40_Study2_Saliva_TP2_P40Placebo_500","22F41_Study2_Saliva_TP2_P41Placebo_500","22F42_Study2_Saliva_TP2_P42Placebo_500","22F44_Study2_Saliva_TP2_P44Placebo_500","22F48_Study2_Saliva_TP2_P48Placebo_500","22F50_Study2_Saliva_TP2_P50Placebo_500","22F51_Study2_Saliva_TP2_P51Placebo_500","22F52_Study2_Saliva_TP2_P52Placebo_500","22F54_Study2_Saliva_TP2_P54Placebo_500","22F55_Study2_Saliva_TP2_P55Placebo_500","22F57_Study2_Saliva_TP2_P57Placebo_500","22F61_Study2_Saliva_TP2_P61Placebo_500","22F63_Study2_Saliva_TP2_P63Placebo_500","22F65_Study2_Saliva_TP2_P65Placebo_500","22F67_Study2_Saliva_TP2_P67Placebo_500","22F74_Study2_Saliva_TP2_P74Placebo_500","22F75_Study2_Saliva_TP2_P75Placebo_500","22F76_Study2_Saliva_TP2_P76Placebo_500","22F77_Study2_Saliva_TP2_P77Placebo_500","22F78_Study2_Saliva_TP2_P78Placebo_500","22F80_Study2_Saliva_TP2_P80Placebo_500","23F1_Study2_Saliva_TP3_P1Placebo_500","23F2_Study2_Saliva_TP3_P2Placebo_500","23F4_Study2_Saliva_TP3_P4Placebo_500","23F5_Study2_Saliva_TP3_P5Placebo_500","23F6_Study2_Saliva_TP3_P6Placebo_500","23F14_Study2_Saliva_TP3_P14Placebo_500","23F16_Study2_Saliva_TP3_P16Placebo_500","23F18_Study2_Saliva_TP3_P18Placebo_500","23F19_Study2_Saliva_TP3_P19Placebo_500","23F20_Study2_Saliva_TP3_P20Placebo_500","23F21_Study2_Saliva_TP3_P21Placebo_500","23F23_Study2_Saliva_TP3_P23Placebo_500","23F25_Study2_Saliva_TP3_P25Placebo_500","23F29_Study2_Saliva_TP3_P29Placebo_500","23F31_Study2_Saliva_TP3_P31Placebo_500","23F35_Study2_Saliva_TP3_P35Placebo_500","23F36_Study2_Saliva_TP3_P36Placebo_500","23F38_Study2_Saliva_TP3_P38Placebo_500","23F39_Study2_Saliva_TP3_P39Placebo_500","23F40_Study2_Saliva_TP3_P40Placebo_500","23F41_Study2_Saliva_TP3_P41Placebo_500","23F42_Study2_Saliva_TP3_P42Placebo_500","23F44_Study2_Saliva_TP3_P44Placebo_500","23F48_Study2_Saliva_TP3_P48Placebo_500","23F50_Study2_Saliva_TP3_P50Placebo_500","23F51_Study2_Saliva_TP3_P51Placebo_500","23F52_Study2_Saliva_TP3_P52Placebo_500","23F54_Study2_Saliva_TP3_P54Placebo_500","23F55_Study2_Saliva_TP3_P55Placebo_500","23F57_Study2_Saliva_TP3_P57Placebo_500","23F61_Study2_Saliva_TP3_P61Placebo_500","23F63_Study2_Saliva_TP3_P63Placebo_500","23F65_Study2_Saliva_TP3_P65Placebo_500","23F67_Study2_Saliva_TP3_P67Placebo_500","23F74_Study2_Saliva_TP3_P74Placebo_500","23F75_Study2_Saliva_TP3_P75Placebo_500","23F76_Study2_Saliva_TP3_P76Placebo_500","23F77_Study2_Saliva_TP3_P77Placebo_500","23F78_Study2_Saliva_TP3_P78Placebo_500","23F80_Study2_Saliva_TP3_P80Placebo_500"],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,2,1,2,2,1,2,2,2,1,2,1,1,2,2,2,2,1,2,1,2,1,1,2,2,1,1,2,2,2,2,2,1,1,2,1,1,1,2,2,1,2,1,2,2,1,2,2,2,1,2,1,1,2,2,2,2,1,2,1,2,1,1,2,2,1,1,2,2,2,2,2,1,1,2,1,1,1,2,2,1,2,1,2,2,1,2,2,2,1,2,1,1,2,2,2,2,1,1,2,1,2,1,2,1,1,2,2,1,1,2,2,2,2,1,1,1,2,1,1,2,1,2,2,1,2,2,2,1,2,1,1,2,2,2,2,1,2,1,2,1,1,2,2,1,1,2,2,2,2,2,1,1,2,1,1,1,2,2,1,2,1,2,2,1,2,2,2,1,2,1,1,2,2,2,2,1,2,1,2,1,1,2,2,1,1,2,2,2,2,2,1,1,2,1,1,1,2,2,1,2,1,2,2,1,2,2,2,1,2,1,1,2,2,2,2,1,2,1,2,1,1,2,2,1,1,2,2,2,2,2,1,1,2,1,1,1,2,2],[1.01785714285714,0.857142857142857,0.4,0.541666666666667,1.03571428571429,0.821428571428571,1.16071428571429,0.854166666666667,0.869047619047619,1.25,1.13095238095238,0.68452380952381,0.511904761904762,1.13690476190476,1.13095238095238,0.518518518518518,0.732142857142857,1.29761904761905,0.357142857142857,1.22222222222222,1.41071428571429,0.767857142857143,0.273809523809524,0.791666666666667,0.845238095238095,1.25,1.01785714285714,1.17857142857143,0.595238095238095,1.54761904761905,0.375,0.869047619047619,0.839285714285714,0.404761904761905,1.27380952380952,1.41666666666667,0.755952380952381,1.20238095238095,1.5448717948718,1.04166666666667,1.8452380952381,2.47619047619048,2.32666666666667,2.76785714285714,1.97619047619048,2.53571428571429,2.33928571428571,2.38888888888889,2.54761904761905,2.92857142857143,2.89880952380952,2.70238095238095,2.5952380952381,2.56547619047619,3.2202380952381,2.5679012345679,1.86309523809524,2.89880952380952,2.11904761904762,2.91666666666667,2.875,2.625,2.45833333333333,1.95833333333333,2.375,2.80357142857143,3.19047619047619,2.42261904761905,3.0952380952381,2.69642857142857,2.57738095238095,2.5,2.19047619047619,1.95833333333333,3.00595238095238,3.01190476190476,2.01190476190476,2.5,2.53205128205128,2.64880952380952,1.61904761904762,0.744047619047619,0.566666666666667,0.7738095238095239,0.5238095238095239,0.666666666666667,1.00595238095238,1.28472222222222,0.8988095238095239,1.75,0.869047619047619,0.755952380952381,1.14880952380952,0.636904761904762,1.0952380952381,0.944444444444444,0.732142857142857,0.886904761904762,0.861111111111111,1.61309523809524,0.839285714285714,0.261904761904762,0.482142857142857,0.529761904761905,1,1.125,0.410714285714286,0.732142857142857,1.19642857142857,0.916666666666667,0.386904761904762,0.369047619047619,1.38690476190476,1.44047619047619,1.27380952380952,1.33333333333333,1.48214285714286,1.71153846153846,1,1.01785714285714,0.857142857142857,0.4,0.541666666666667,1.03571428571429,0.821428571428571,1.16071428571429,0.854166666666667,0.869047619047619,1.25,1.13095238095238,0.68452380952381,0.511904761904762,1.13690476190476,1.13095238095238,0.518518518518518,0.732142857142857,1.29761904761905,0.357142857142857,1.22222222222222,1.41071428571429,0.767857142857143,0.273809523809524,0.791666666666667,0.845238095238095,1.25,1.01785714285714,1.17857142857143,0.595238095238095,1.54761904761905,0.375,0.869047619047619,0.839285714285714,0.404761904761905,1.27380952380952,1.41666666666667,0.755952380952381,1.20238095238095,1.5448717948718,1.04166666666667,1.8452380952381,2.47619047619048,2.32666666666667,2.76785714285714,1.97619047619048,2.53571428571429,2.33928571428571,2.38888888888889,2.54761904761905,2.92857142857143,2.89880952380952,2.70238095238095,2.5952380952381,2.56547619047619,3.2202380952381,2.5679012345679,1.86309523809524,2.89880952380952,2.11904761904762,2.91666666666667,2.875,2.625,2.45833333333333,1.95833333333333,2.375,2.80357142857143,3.19047619047619,2.42261904761905,3.0952380952381,2.69642857142857,2.57738095238095,2.5,2.19047619047619,1.95833333333333,3.00595238095238,3.01190476190476,2.01190476190476,2.5,2.53205128205128,2.64880952380952,1.61904761904762,0.744047619047619,0.566666666666667,0.7738095238095239,0.5238095238095239,0.666666666666667,1.00595238095238,1.28472222222222,0.8988095238095239,1.75,0.869047619047619,0.755952380952381,1.14880952380952,0.636904761904762,1.0952380952381,0.944444444444444,0.732142857142857,0.886904761904762,0.303571428571429,0.861111111111111,1.61309523809524,0.839285714285714,0.261904761904762,0.482142857142857,0.529761904761905,1,1.125,0.410714285714286,0.732142857142857,1.19642857142857,0.916666666666667,0.386904761904762,0.369047619047619,1.38690476190476,1.44047619047619,1.27380952380952,1.33333333333333,1.48214285714286,1.71153846153846,1],[7.73809523809524,1.78571428571429,0.666666666666667,3.57142857142857,5.95238095238095,0,1.78571428571429,2.77777777777778,5.35714285714286,10.1190476190476,10.7142857142857,8.928571428571431,4.16666666666667,1.19047619047619,9.52380952380952,1.85185185185185,4.76190476190476,8.928571428571431,2.38095238095238,0.694444444444444,2.38095238095238,2.38095238095238,1.19047619047619,1.19047619047619,4.16666666666667,2.97619047619048,16.6666666666667,10.1190476190476,6.54761904761905,3.57142857142857,7.73809523809524,4.76190476190476,1.19047619047619,4.76190476190476,7.14285714285714,10.7142857142857,3.57142857142857,12.5,7.69230769230769,10.1190476190476,25.5952380952381,10.1190476190476,26,11.9047619047619,10.1190476190476,14.8809523809524,28.5714285714286,22.9166666666667,27.3809523809524,26.7857142857143,40.4761904761905,22.6190476190476,20.8333333333333,16.0714285714286,39.8809523809524,33.9506172839506,13.0952380952381,17.2619047619048,15.4761904761905,19.4444444444444,18.452380952381,13.6904761904762,17.2619047619048,11.9047619047619,24.4047619047619,30.952380952381,27.3809523809524,35.1190476190476,35.7142857142857,27.9761904761905,33.3333333333333,20.8333333333333,30.3571428571429,23.2142857142857,26.1904761904762,39.8809523809524,32.7380952380952,36.3095238095238,27.5641025641026,25.5952380952381,15.4761904761905,5.35714285714286,5.33333333333333,4.76190476190476,8.928571428571431,4.76190476190476,7.73809523809524,16.6666666666667,2.38095238095238,13.0952380952381,16.6666666666667,10.1190476190476,5.35714285714286,3.57142857142857,11.9047619047619,1.23456790123457,1.78571428571429,8.33333333333333,8.33333333333333,11.9047619047619,2.97619047619048,1.78571428571429,10.1190476190476,10.7142857142857,8.928571428571431,20.8333333333333,5.95238095238095,8.928571428571431,11.9047619047619,1.19047619047619,6.54761904761905,6.54761904761905,5.35714285714286,7.73809523809524,13.6904761904762,8.928571428571431,10.1190476190476,12.8205128205128,8.928571428571431,7.73809523809524,1.78571428571429,0.666666666666667,3.57142857142857,5.95238095238095,0,1.78571428571429,2.77777777777778,5.35714285714286,10.1190476190476,10.7142857142857,8.928571428571431,4.16666666666667,1.19047619047619,9.52380952380952,1.85185185185185,4.76190476190476,8.928571428571431,2.38095238095238,0.694444444444444,2.38095238095238,2.38095238095238,1.19047619047619,1.19047619047619,4.16666666666667,2.97619047619048,16.6666666666667,10.1190476190476,6.54761904761905,3.57142857142857,7.73809523809524,4.76190476190476,1.19047619047619,4.76190476190476,7.14285714285714,10.7142857142857,3.57142857142857,12.5,7.69230769230769,10.1190476190476,25.5952380952381,10.1190476190476,26,11.9047619047619,10.1190476190476,14.8809523809524,28.5714285714286,22.9166666666667,27.3809523809524,26.7857142857143,40.4761904761905,22.6190476190476,20.8333333333333,16.0714285714286,39.8809523809524,33.9506172839506,13.0952380952381,17.2619047619048,15.4761904761905,19.4444444444444,18.452380952381,13.6904761904762,17.2619047619048,11.9047619047619,24.4047619047619,30.952380952381,27.3809523809524,35.1190476190476,35.7142857142857,27.9761904761905,33.3333333333333,20.8333333333333,30.3571428571429,23.2142857142857,26.1904761904762,39.8809523809524,32.7380952380952,36.3095238095238,27.5641025641026,25.5952380952381,15.4761904761905,5.35714285714286,5.33333333333333,4.76190476190476,8.928571428571431,4.76190476190476,7.73809523809524,16.6666666666667,2.38095238095238,13.0952380952381,16.6666666666667,10.1190476190476,5.35714285714286,3.57142857142857,11.9047619047619,1.23456790123457,1.78571428571429,8.33333333333333,7.14285714285714,8.33333333333333,11.9047619047619,2.97619047619048,6.54761904761905,4.16666666666667,3.57142857142857,8.928571428571431,20.8333333333333,8.33333333333333,8.928571428571431,5.95238095238095,4.16666666666667,2.38095238095238,17.8571428571429,10.1190476190476,20.8333333333333,13.6904761904762,8.928571428571431,10.1190476190476,12.8205128205128,8.928571428571431],[24,23,25,21,27,23,28,23,27,25,21,26,21,23,22,26,26,24,25,33,29,22,25,23,20,21,25,24,24,29,24,24,19,21,25,23,27,24,27,22,24,23,25,21,27,23,28,23,27,25,21,26,21,23,22,26,26,24,25,33,29,22,25,23,20,21,25,24,24,29,24,24,19,21,25,23,27,24,27,22,24,23,25,21,27,23,28,23,27,25,21,26,21,23,22,26,26,24,33,29,22,22,23,22,21,25,30,24,20,25,24,24,26,26,23,27,24,27,22,24,23,25,21,27,23,28,23,27,25,21,26,21,23,22,26,26,24,25,33,29,22,25,23,20,21,25,24,24,29,24,24,19,21,25,23,27,24,27,22,24,23,25,21,27,23,28,23,27,25,21,26,21,23,22,26,26,24,25,33,29,22,25,23,20,21,25,24,24,29,24,24,19,21,25,23,27,24,27,22,24,23,25,21,27,23,28,23,27,25,21,26,21,23,22,26,26,24,25,33,29,22,25,23,20,21,25,24,24,29,24,24,19,21,25,23,27,24,27,22]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ID<\/th>\n      <th>ID_external<\/th>\n      <th>Sample<\/th>\n      <th>Time<\/th>\n      <th>Subject<\/th>\n      <th>Group<\/th>\n      <th>Quant<\/th>\n      <th>Full_name<\/th>\n      <th>group<\/th>\n      <th>sex<\/th>\n      <th>mean_plaque<\/th>\n      <th>mean_bleeding<\/th>\n      <th>age<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,7,9,10,11,12,13]},{"orderable":false,"targets":0},{"name":" ","targets":0},{"name":"ID","targets":1},{"name":"ID_external","targets":2},{"name":"Sample","targets":3},{"name":"Time","targets":4},{"name":"Subject","targets":5},{"name":"Group","targets":6},{"name":"Quant","targets":7},{"name":"Full_name","targets":8},{"name":"group","targets":9},{"name":"sex","targets":10},{"name":"mean_plaque","targets":11},{"name":"mean_bleeding","targets":12},{"name":"age","targets":13}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


``` r
ps %>% 
  physeq_add_metadata(metadata = meta,
                      sample_column = "ID") -> ps
```

phyloseq doenst like numbers in sample names, just add S_ to those:


``` r
ps %>% 
  sample_names() %>% 
  paste0("S_", .) -> sample_names(ps)
```


``` r
ps
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 603 taxa and 239 samples ]
## sample_data() Sample Data:       [ 239 samples by 13 sample variables ]
## tax_table()   Taxonomy Table:    [ 603 taxa by 8 taxonomic ranks ]
```


``` r
ps %>% 
  microViz::ps_mutate(Time = factor(Time, levels = c("TP1", "TP2", "TP3")),
                      Period = Time,
                      time = parse_number(as.character(Time)),
                      Sample_Type =  Sample,
                      Period = recode(Period, TP1  = "Baseline"),
                      Sex = recode(sex, Male = 1, Female = 2),
                      Sample = factor(Sample, levels = c("Saliva", "Plaque")),
                      Subject = fct_relevel(Subject, "D01"),
                      Sample_Time = paste0(Sample,"_", Time) 
  ) -> ps_up

ps_up %>% 
  sample_data(.) %>% 
  data.frame() -> dfdf

dfdf %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item" id="htmlwidget-74e53db2a95ab5e3d88c" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-74e53db2a95ab5e3d88c">{"x":{"filter":"none","vertical":false,"data":[["S_11563","S_11564","S_11566","S_11567","S_11568","S_11576","S_11578","S_11580","S_11581","S_11582","S_11583","S_11585","S_11587","S_11591","S_11593","S_11597","S_11598","S_11600","S_11601","S_11602","S_11603","S_11604","S_11606","S_11610","S_11612","S_11613","S_11614","S_11616","S_11617","S_11619","S_11623","S_11625","S_11627","S_11629","S_11636","S_11637","S_11638","S_11639","S_11640","S_11642","S_11643","S_11644","S_11646","S_11647","S_11648","S_11656","S_11658","S_11660","S_11661","S_11662","S_11663","S_11665","S_11667","S_11671","S_11673","S_11677","S_11678","S_11680","S_11681","S_11682","S_11683","S_11684","S_11686","S_11690","S_11692","S_11693","S_11694","S_11696","S_11697","S_11699","S_11703","S_11705","S_11707","S_11709","S_11716","S_11717","S_11718","S_11719","S_11720","S_11722","S_11723","S_11724","S_11726","S_11727","S_11728","S_11736","S_11738","S_11740","S_11741","S_11742","S_11743","S_11745","S_11747","S_11752","S_11754","S_11758","S_11759","S_11761","S_11763","S_11764","S_11765","S_11767","S_11771","S_11773","S_11774","S_11775","S_11777","S_11778","S_11780","S_11784","S_11786","S_11788","S_11790","S_11797","S_11798","S_11799","S_11800","S_11803","S_11805","S_11806","S_11807","S_11809","S_11810","S_11811","S_11819","S_11821","S_11823","S_11824","S_11825","S_11826","S_11828","S_11830","S_11834","S_11836","S_11840","S_11841","S_11843","S_11844","S_11845","S_11846","S_11847","S_11849","S_11853","S_11855","S_11856","S_11857","S_11859","S_11860","S_11862","S_11866","S_11868","S_11870","S_11872","S_11879","S_11880","S_11881","S_11882","S_11883","S_11885","S_11886","S_11887","S_11889","S_11890","S_11891","S_11899","S_11901","S_11903","S_11904","S_11905","S_11906","S_11908","S_11910","S_11914","S_11916","S_11920","S_11921","S_11923","S_11924","S_11925","S_11926","S_11927","S_11929","S_11933","S_11935","S_11936","S_11937","S_11939","S_11940","S_11942","S_11946","S_11948","S_11950","S_11952","S_11959","S_11960","S_11961","S_11962","S_11963","S_11965","S_11966","S_11967","S_11969","S_11970","S_11971","S_11979","S_11981","S_11983","S_11984","S_11985","S_11986","S_11988","S_11990","S_11994","S_11996","S_12000","S_12001","S_12003","S_12004","S_12005","S_12006","S_12007","S_12009","S_12013","S_12015","S_12016","S_12017","S_12019","S_12020","S_12022","S_12026","S_12028","S_12030","S_12032","S_12039","S_12040","S_12041","S_12042","S_12043","S_12045"],[11563,11564,11566,11567,11568,11576,11578,11580,11581,11582,11583,11585,11587,11591,11593,11597,11598,11600,11601,11602,11603,11604,11606,11610,11612,11613,11614,11616,11617,11619,11623,11625,11627,11629,11636,11637,11638,11639,11640,11642,11643,11644,11646,11647,11648,11656,11658,11660,11661,11662,11663,11665,11667,11671,11673,11677,11678,11680,11681,11682,11683,11684,11686,11690,11692,11693,11694,11696,11697,11699,11703,11705,11707,11709,11716,11717,11718,11719,11720,11722,11723,11724,11726,11727,11728,11736,11738,11740,11741,11742,11743,11745,11747,11752,11754,11758,11759,11761,11763,11764,11765,11767,11771,11773,11774,11775,11777,11778,11780,11784,11786,11788,11790,11797,11798,11799,11800,11803,11805,11806,11807,11809,11810,11811,11819,11821,11823,11824,11825,11826,11828,11830,11834,11836,11840,11841,11843,11844,11845,11846,11847,11849,11853,11855,11856,11857,11859,11860,11862,11866,11868,11870,11872,11879,11880,11881,11882,11883,11885,11886,11887,11889,11890,11891,11899,11901,11903,11904,11905,11906,11908,11910,11914,11916,11920,11921,11923,11924,11925,11926,11927,11929,11933,11935,11936,11937,11939,11940,11942,11946,11948,11950,11952,11959,11960,11961,11962,11963,11965,11966,11967,11969,11970,11971,11979,11981,11983,11984,11985,11986,11988,11990,11994,11996,12000,12001,12003,12004,12005,12006,12007,12009,12013,12015,12016,12017,12019,12020,12022,12026,12028,12030,12032,12039,12040,12041,12042,12043,12045],["21B1","21B2","21B4","21B5","21B6","21B14","21B16","21B18","21B19","21B20","21B21","21B23","21B25","21B29","21B31","21B35","21B36","21B38","21B39","21B40","21B41","21B42","21B44","21B48","21B50","21B51","21B52","21B54","21B55","21B57","21B61","21B63","21B65","21B67","21B74","21B75","21B76","21B77","21B78","21B80","22B1","22B2","22B4","22B5","22B6","22B14","22B16","22B18","22B19","22B20","22B21","22B23","22B25","22B29","22B31","22B35","22B36","22B38","22B39","22B40","22B41","22B42","22B44","22B48","22B50","22B51","22B52","22B54","22B55","22B57","22B61","22B63","22B65","22B67","22B74","22B75","22B76","22B77","22B78","22B80","23B1","23B2","23B4","23B5","23B6","23B14","23B16","23B18","23B19","23B20","23B21","23B23","23B25","23B29","23B31","23B35","23B36","23B38","23B40","23B41","23B42","23B44","23B48","23B50","23B51","23B52","23B54","23B55","23B57","23B61","23B63","23B65","23B67","23B74","23B75","23B76","23B77","23B78","23B80","21F1","21F2","21F4","21F5","21F6","21F14","21F16","21F18","21F19","21F20","21F21","21F23","21F25","21F29","21F31","21F35","21F36","21F38","21F39","21F40","21F41","21F42","21F44","21F48","21F50","21F51","21F52","21F54","21F55","21F57","21F61","21F63","21F65","21F67","21F74","21F75","21F76","21F77","21F78","21F80","22F1","22F2","22F4","22F5","22F6","22F14","22F16","22F18","22F19","22F20","22F21","22F23","22F25","22F29","22F31","22F35","22F36","22F38","22F39","22F40","22F41","22F42","22F44","22F48","22F50","22F51","22F52","22F54","22F55","22F57","22F61","22F63","22F65","22F67","22F74","22F75","22F76","22F77","22F78","22F80","23F1","23F2","23F4","23F5","23F6","23F14","23F16","23F18","23F19","23F20","23F21","23F23","23F25","23F29","23F31","23F35","23F36","23F38","23F39","23F40","23F41","23F42","23F44","23F48","23F50","23F51","23F52","23F54","23F55","23F57","23F61","23F63","23F65","23F67","23F74","23F75","23F76","23F77","23F78","23F80"],["Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva"],["TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3"],["D01","D02","D04","D05","D06","D14","D16","D18","D19","D20","D21","D23","D25","D29","D31","D35","D36","D38","D39","D40","D41","D42","D44","D48","D50","D51","D52","D54","D55","D57","D61","D63","D65","D67","D74","D75","D76","D77","D78","D80","D01","D02","D04","D05","D06","D14","D16","D18","D19","D20","D21","D23","D25","D29","D31","D35","D36","D38","D39","D40","D41","D42","D44","D48","D50","D51","D52","D54","D55","D57","D61","D63","D65","D67","D74","D75","D76","D77","D78","D80","D01","D02","D04","D05","D06","D14","D16","D18","D19","D20","D21","D23","D25","D29","D31","D35","D36","D38","D40","D41","D42","D44","D48","D50","D51","D52","D54","D55","D57","D61","D63","D65","D67","D74","D75","D76","D77","D78","D80","D01","D02","D04","D05","D06","D14","D16","D18","D19","D20","D21","D23","D25","D29","D31","D35","D36","D38","D39","D40","D41","D42","D44","D48","D50","D51","D52","D54","D55","D57","D61","D63","D65","D67","D74","D75","D76","D77","D78","D80","D01","D02","D04","D05","D06","D14","D16","D18","D19","D20","D21","D23","D25","D29","D31","D35","D36","D38","D39","D40","D41","D42","D44","D48","D50","D51","D52","D54","D55","D57","D61","D63","D65","D67","D74","D75","D76","D77","D78","D80","D01","D02","D04","D05","D06","D14","D16","D18","D19","D20","D21","D23","D25","D29","D31","D35","D36","D38","D39","D40","D41","D42","D44","D48","D50","D51","D52","D54","D55","D57","D61","D63","D65","D67","D74","D75","D76","D77","D78","D80"],["Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo","Placebo"],[1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500],["21B1_NatriumChlorid0.9Study2_Plaque_TP1_P1Placebo_1000","21B2_NatriumChlorid0.9Study2_Plaque_TP1_P2Placebo_1000","21B4_NatriumChlorid0.9Study2_Plaque_TP1_P4Placebo_1000","21B5_NatriumChlorid0.9Study2_Plaque_TP1_P5Placebo_1000","21B6_NatriumChlorid0.9Study2_Plaque_TP1_P6Placebo_1000","21B14_NatriumChlorid0.9Study2_Plaque_TP1_P14Placebo_1000","21B16_NatriumChlorid0.9Study2_Plaque_TP1_P16Placebo_1000","21B18_NatriumChlorid0.9Study2_Plaque_TP1_P18Placebo_1000","21B19_NatriumChlorid0.9Study2_Plaque_TP1_P19Placebo_1000","21B20_NatriumChlorid0.9Study2_Plaque_TP1_P20Placebo_1000","21B21_NatriumChlorid0.9Study2_Plaque_TP1_P21Placebo_1000","21B23_NatriumChlorid0.9Study2_Plaque_TP1_P23Placebo_1000","21B25_NatriumChlorid0.9Study2_Plaque_TP1_P25Placebo_1000","21B29_NatriumChlorid0.9Study2_Plaque_TP1_P29Placebo_1000","21B31_NatriumChlorid0.9Study2_Plaque_TP1_P31Placebo_1000","21B35_NatriumChlorid0.9Study2_Plaque_TP1_P35Placebo_1000","21B36_NatriumChlorid0.9Study2_Plaque_TP1_P36Placebo_1000","21B38_NatriumChlorid0.9Study2_Plaque_TP1_P38Placebo_1000","21B39_NatriumChlorid0.9Study2_Plaque_TP1_P39Placebo_1000","21B40_NatriumChlorid0.9Study2_Plaque_TP1_P40Placebo_1000","21B41_NatriumChlorid0.9Study2_Plaque_TP1_P41Placebo_1000","21B42_NatriumChlorid0.9Study2_Plaque_TP1_P42Placebo_1000","21B44_NatriumChlorid0.9Study2_Plaque_TP1_P44Placebo_1000","21B48_NatriumChlorid0.9Study2_Plaque_TP1_P48Placebo_1000","21B50_NatriumChlorid0.9Study2_Plaque_TP1_P50Placebo_1000","21B51_NatriumChlorid0.9Study2_Plaque_TP1_P51Placebo_1000","21B52_NatriumChlorid0.9Study2_Plaque_TP1_P52Placebo_1000","21B54_NatriumChlorid0.9Study2_Plaque_TP1_P54Placebo_1000","21B55_NatriumChlorid0.9Study2_Plaque_TP1_P55Placebo_1000","21B57_NatriumChlorid0.9Study2_Plaque_TP1_P57Placebo_1000","21B61_NatriumChlorid0.9Study2_Plaque_TP1_P61Placebo_1000","21B63_NatriumChlorid0.9Study2_Plaque_TP1_P63Placebo_1000","21B65_NatriumChlorid0.9Study2_Plaque_TP1_P65Placebo_1000","21B67_NatriumChlorid0.9Study2_Plaque_TP1_P67Placebo_1000","21B74_NatriumChlorid0.9Study2_Plaque_TP1_P74Placebo_1000","21B75_NatriumChlorid0.9Study2_Plaque_TP1_P75Placebo_1000","21B76_NatriumChlorid0.9Study2_Plaque_TP1_P76Placebo_1000","21B77_NatriumChlorid0.9Study2_Plaque_TP1_P77Placebo_1000","21B78_NatriumChlorid0.9Study2_Plaque_TP1_P78Placebo_1000","21B80_NatriumChlorid0.9Study2_Plaque_TP1_P80Placebo_1000","22B1_NatriumChlorid0.9Study2_Plaque_TP2_P1Placebo_1000","22B2_NatriumChlorid0.9Study2_Plaque_TP2_P2Placebo_1000","22B4_NatriumChlorid0.9Study2_Plaque_TP2_P4Placebo_1000","22B5_NatriumChlorid0.9Study2_Plaque_TP2_P5Placebo_1000","22B6_NatriumChlorid0.9Study2_Plaque_TP2_P6Placebo_1000","22B14_NatriumChlorid0.9Study2_Plaque_TP2_P14Placebo_1000","22B16_NatriumChlorid0.9Study2_Plaque_TP2_P16Placebo_1000","22B18_NatriumChlorid0.9Study2_Plaque_TP2_P18Placebo_1000","22B19_NatriumChlorid0.9Study2_Plaque_TP2_P19Placebo_1000","22B20_NatriumChlorid0.9Study2_Plaque_TP2_P20Placebo_1000","22B21_NatriumChlorid0.9Study2_Plaque_TP2_P21Placebo_1000","22B23_NatriumChlorid0.9Study2_Plaque_TP2_P23Placebo_1000","22B25_NatriumChlorid0.9Study2_Plaque_TP2_P25Placebo_1000","22B29_NatriumChlorid0.9Study2_Plaque_TP2_P29Placebo_1000","22B31_NatriumChlorid0.9Study2_Plaque_TP2_P31Placebo_1000","22B35_NatriumChlorid0.9Study2_Plaque_TP2_P35Placebo_1000","22B36_NatriumChlorid0.9Study2_Plaque_TP2_P36Placebo_1000","22B38_NatriumChlorid0.9Study2_Plaque_TP2_P38Placebo_1000","22B39_NatriumChlorid0.9Study2_Plaque_TP2_P39Placebo_1000","22B40_NatriumChlorid0.9Study2_Plaque_TP2_P40Placebo_1000","22B41_NatriumChlorid0.9Study2_Plaque_TP2_P41Placebo_1000","22B42_NatriumChlorid0.9Study2_Plaque_TP2_P42Placebo_1000","22B44_NatriumChlorid0.9Study2_Plaque_TP2_P44Placebo_1000","22B48_NatriumChlorid0.9Study2_Plaque_TP2_P48Placebo_1000","22B50_NatriumChlorid0.9Study2_Plaque_TP2_P50Placebo_1000","22B51_NatriumChlorid0.9Study2_Plaque_TP2_P51Placebo_1000","22B52_NatriumChlorid0.9Study2_Plaque_TP2_P52Placebo_1000","22B54_NatriumChlorid0.9Study2_Plaque_TP2_P54Placebo_1000","22B55_NatriumChlorid0.9Study2_Plaque_TP2_P55Placebo_1000","22B57_NatriumChlorid0.9Study2_Plaque_TP2_P57Placebo_1000","22B61_NatriumChlorid0.9Study2_Plaque_TP2_P61Placebo_1000","22B63_NatriumChlorid0.9Study2_Plaque_TP2_P63Placebo_1000","22B65_NatriumChlorid0.9Study2_Plaque_TP2_P65Placebo_1000","22B67_NatriumChlorid0.9Study2_Plaque_TP2_P67Placebo_1000","22B74_NatriumChlorid0.9Study2_Plaque_TP2_P74Placebo_1000","22B75_NatriumChlorid0.9Study2_Plaque_TP2_P75Placebo_1000","22B76_NatriumChlorid0.9Study2_Plaque_TP2_P76Placebo_1000","22B77_NatriumChlorid0.9Study2_Plaque_TP2_P77Placebo_1000","22B78_NatriumChlorid0.9Study2_Plaque_TP2_P78Placebo_1000","22B80_NatriumChlorid0.9Study2_Plaque_TP2_P80Placebo_1000","23B1_NatriumChlorid0.9Study2_Plaque_TP3_P1Placebo_1000","23B2_NatriumChlorid0.9Study2_Plaque_TP3_P2Placebo_1000","23B4_NatriumChlorid0.9Study2_Plaque_TP3_P4Placebo_1000","23B5_NatriumChlorid0.9Study2_Plaque_TP3_P5Placebo_1000","23B6_NatriumChlorid0.9Study2_Plaque_TP3_P6Placebo_1000","23B14_NatriumChlorid0.9Study2_Plaque_TP3_P14Placebo_1000","23B16_NatriumChlorid0.9Study2_Plaque_TP3_P16Placebo_1000","23B18_NatriumChlorid0.9Study2_Plaque_TP3_P18Placebo_1000","23B19_NatriumChlorid0.9Study2_Plaque_TP3_P19Placebo_1000","23B20_NatriumChlorid0.9Study2_Plaque_TP3_P20Placebo_1000","23B21_NatriumChlorid0.9Study2_Plaque_TP3_P21Placebo_1000","23B23_NatriumChlorid0.9Study2_Plaque_TP3_P23Placebo_1000","23B25_NatriumChlorid0.9Study2_Plaque_TP3_P25Placebo_1000","23B29_NatriumChlorid0.9Study2_Plaque_TP3_P29Placebo_1000","23B31_NatriumChlorid0.9Study2_Plaque_TP3_P31Placebo_1000","23B35_NatriumChlorid0.9Study2_Plaque_TP3_P35Placebo_1000","23B36_NatriumChlorid0.9Study2_Plaque_TP3_P36Placebo_1000","23B38_NatriumChlorid0.9Study2_Plaque_TP3_P38Placebo_1000","23B40_NatriumChlorid0.9Study2_Plaque_TP3_P40Placebo_1000","23B41_NatriumChlorid0.9Study2_Plaque_TP3_P41Placebo_1000","23B42_NatriumChlorid0.9Study2_Plaque_TP3_P42Placebo_1000","23B44_NatriumChlorid0.9Study2_Plaque_TP3_P44Placebo_1000","23B48_NatriumChlorid0.9Study2_Plaque_TP3_P48Placebo_1000","23B50_NatriumChlorid0.9Study2_Plaque_TP3_P50Placebo_1000","23B51_NatriumChlorid0.9Study2_Plaque_TP3_P51Placebo_1000","23B52_NatriumChlorid0.9Study2_Plaque_TP3_P52Placebo_1000","23B54_NatriumChlorid0.9Study2_Plaque_TP3_P54Placebo_1000","23B55_NatriumChlorid0.9Study2_Plaque_TP3_P55Placebo_1000","23B57_NatriumChlorid0.9Study2_Plaque_TP3_P57Placebo_1000","23B61_NatriumChlorid0.9Study2_Plaque_TP3_P61Placebo_1000","23B63_NatriumChlorid0.9Study2_Plaque_TP3_P63Placebo_1000","23B65_NatriumChlorid0.9Study2_Plaque_TP3_P65Placebo_1000","23B67_NatriumChlorid0.9Study2_Plaque_TP3_P67Placebo_1000","23B74_NatriumChlorid0.9Study2_Plaque_TP3_P74Placebo_1000","23B75_NatriumChlorid0.9Study2_Plaque_TP3_P75Placebo_1000","23B76_NatriumChlorid0.9Study2_Plaque_TP3_P76Placebo_1000","23B77_NatriumChlorid0.9Study2_Plaque_TP3_P77Placebo_1000","23B78_NatriumChlorid0.9Study2_Plaque_TP3_P78Placebo_1000","23B80_NatriumChlorid0.9Study2_Plaque_TP3_P80Placebo_1000","21F1_Study2_Saliva_TP1_P1Placebo_500","21F2_Study2_Saliva_TP1_P2Placebo_500","21F4_Study2_Saliva_TP1_P4Placebo_500","21F5_Study2_Saliva_TP1_P5Placebo_500","21F6_Study2_Saliva_TP1_P6Placebo_500","21F14_Study2_Saliva_TP1_P14Placebo_500","21F16_Study2_Saliva_TP1_P16Placebo_500","21F18_Study2_Saliva_TP1_P18Placebo_500","21F19_Study2_Saliva_TP1_P19Placebo_500","21F20_Study2_Saliva_TP1_P20Placebo_500","21F21_Study2_Saliva_TP1_P21Placebo_500","21F23_Study2_Saliva_TP1_P23Placebo_500","21F25_Study2_Saliva_TP1_P25Placebo_500","21F29_Study2_Saliva_TP1_P29Placebo_500","21F31_Study2_Saliva_TP1_P31Placebo_500","21F35_Study2_Saliva_TP1_P35Placebo_500","21F36_Study2_Saliva_TP1_P36Placebo_500","21F38_Study2_Saliva_TP1_P38Placebo_500","21F39_Study2_Saliva_TP1_P39Placebo_500","21F40_Study2_Saliva_TP1_P40Placebo_500","21F41_Study2_Saliva_TP1_P41Placebo_500","21F42_Study2_Saliva_TP1_P42Placebo_500","21F44_Study2_Saliva_TP1_P44Placebo_500","21F48_Study2_Saliva_TP1_P48Placebo_500","21F50_Study2_Saliva_TP1_P50Placebo_500","21F51_Study2_Saliva_TP1_P51Placebo_500","21F52_Study2_Saliva_TP1_P52Placebo_500","21F54_Study2_Saliva_TP1_P54Placebo_500","21F55_Study2_Saliva_TP1_P55Placebo_500","21F57_Study2_Saliva_TP1_P57Placebo_500","21F61_Study2_Saliva_TP1_P61Placebo_500","21F63_Study2_Saliva_TP1_P63Placebo_500","21F65_Study2_Saliva_TP1_P65Placebo_500","21F67_Study2_Saliva_TP1_P67Placebo_500","21F74_Study2_Saliva_TP1_P74Placebo_500","21F75_Study2_Saliva_TP1_P75Placebo_500","21F76_Study2_Saliva_TP1_P76Placebo_500","21F77_Study2_Saliva_TP1_P77Placebo_500","21F78_Study2_Saliva_TP1_P78Placebo_500","21F80_Study2_Saliva_TP1_P80Placebo_500","22F1_Study2_Saliva_TP2_P1Placebo_500","22F2_Study2_Saliva_TP2_P2Placebo_500","22F4_Study2_Saliva_TP2_P4Placebo_500","22F5_Study2_Saliva_TP2_P5Placebo_500","22F6_Study2_Saliva_TP2_P6Placebo_500","22F14_Study2_Saliva_TP2_P14Placebo_500","22F16_Study2_Saliva_TP2_P16Placebo_500","22F18_Study2_Saliva_TP2_P18Placebo_500","22F19_Study2_Saliva_TP2_P19Placebo_500","22F20_Study2_Saliva_TP2_P20Placebo_500","22F21_Study2_Saliva_TP2_P21Placebo_500","22F23_Study2_Saliva_TP2_P23Placebo_500","22F25_Study2_Saliva_TP2_P25Placebo_500","22F29_Study2_Saliva_TP2_P29Placebo_500","22F31_Study2_Saliva_TP2_P31Placebo_500","22F35_Study2_Saliva_TP2_P35Placebo_500","22F36_Study2_Saliva_TP2_P36Placebo_500","22F38_Study2_Saliva_TP2_P38Placebo_500","22F39_Study2_Saliva_TP2_P39Placebo_500","22F40_Study2_Saliva_TP2_P40Placebo_500","22F41_Study2_Saliva_TP2_P41Placebo_500","22F42_Study2_Saliva_TP2_P42Placebo_500","22F44_Study2_Saliva_TP2_P44Placebo_500","22F48_Study2_Saliva_TP2_P48Placebo_500","22F50_Study2_Saliva_TP2_P50Placebo_500","22F51_Study2_Saliva_TP2_P51Placebo_500","22F52_Study2_Saliva_TP2_P52Placebo_500","22F54_Study2_Saliva_TP2_P54Placebo_500","22F55_Study2_Saliva_TP2_P55Placebo_500","22F57_Study2_Saliva_TP2_P57Placebo_500","22F61_Study2_Saliva_TP2_P61Placebo_500","22F63_Study2_Saliva_TP2_P63Placebo_500","22F65_Study2_Saliva_TP2_P65Placebo_500","22F67_Study2_Saliva_TP2_P67Placebo_500","22F74_Study2_Saliva_TP2_P74Placebo_500","22F75_Study2_Saliva_TP2_P75Placebo_500","22F76_Study2_Saliva_TP2_P76Placebo_500","22F77_Study2_Saliva_TP2_P77Placebo_500","22F78_Study2_Saliva_TP2_P78Placebo_500","22F80_Study2_Saliva_TP2_P80Placebo_500","23F1_Study2_Saliva_TP3_P1Placebo_500","23F2_Study2_Saliva_TP3_P2Placebo_500","23F4_Study2_Saliva_TP3_P4Placebo_500","23F5_Study2_Saliva_TP3_P5Placebo_500","23F6_Study2_Saliva_TP3_P6Placebo_500","23F14_Study2_Saliva_TP3_P14Placebo_500","23F16_Study2_Saliva_TP3_P16Placebo_500","23F18_Study2_Saliva_TP3_P18Placebo_500","23F19_Study2_Saliva_TP3_P19Placebo_500","23F20_Study2_Saliva_TP3_P20Placebo_500","23F21_Study2_Saliva_TP3_P21Placebo_500","23F23_Study2_Saliva_TP3_P23Placebo_500","23F25_Study2_Saliva_TP3_P25Placebo_500","23F29_Study2_Saliva_TP3_P29Placebo_500","23F31_Study2_Saliva_TP3_P31Placebo_500","23F35_Study2_Saliva_TP3_P35Placebo_500","23F36_Study2_Saliva_TP3_P36Placebo_500","23F38_Study2_Saliva_TP3_P38Placebo_500","23F39_Study2_Saliva_TP3_P39Placebo_500","23F40_Study2_Saliva_TP3_P40Placebo_500","23F41_Study2_Saliva_TP3_P41Placebo_500","23F42_Study2_Saliva_TP3_P42Placebo_500","23F44_Study2_Saliva_TP3_P44Placebo_500","23F48_Study2_Saliva_TP3_P48Placebo_500","23F50_Study2_Saliva_TP3_P50Placebo_500","23F51_Study2_Saliva_TP3_P51Placebo_500","23F52_Study2_Saliva_TP3_P52Placebo_500","23F54_Study2_Saliva_TP3_P54Placebo_500","23F55_Study2_Saliva_TP3_P55Placebo_500","23F57_Study2_Saliva_TP3_P57Placebo_500","23F61_Study2_Saliva_TP3_P61Placebo_500","23F63_Study2_Saliva_TP3_P63Placebo_500","23F65_Study2_Saliva_TP3_P65Placebo_500","23F67_Study2_Saliva_TP3_P67Placebo_500","23F74_Study2_Saliva_TP3_P74Placebo_500","23F75_Study2_Saliva_TP3_P75Placebo_500","23F76_Study2_Saliva_TP3_P76Placebo_500","23F77_Study2_Saliva_TP3_P77Placebo_500","23F78_Study2_Saliva_TP3_P78Placebo_500","23F80_Study2_Saliva_TP3_P80Placebo_500"],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,2,1,2,2,1,2,2,2,1,2,1,1,2,2,2,2,1,2,1,2,1,1,2,2,1,1,2,2,2,2,2,1,1,2,1,1,1,2,2,1,2,1,2,2,1,2,2,2,1,2,1,1,2,2,2,2,1,2,1,2,1,1,2,2,1,1,2,2,2,2,2,1,1,2,1,1,1,2,2,1,2,1,2,2,1,2,2,2,1,2,1,1,2,2,2,2,1,1,2,1,2,1,2,1,1,2,2,1,1,2,2,2,2,1,1,1,2,1,1,2,1,2,2,1,2,2,2,1,2,1,1,2,2,2,2,1,2,1,2,1,1,2,2,1,1,2,2,2,2,2,1,1,2,1,1,1,2,2,1,2,1,2,2,1,2,2,2,1,2,1,1,2,2,2,2,1,2,1,2,1,1,2,2,1,1,2,2,2,2,2,1,1,2,1,1,1,2,2,1,2,1,2,2,1,2,2,2,1,2,1,1,2,2,2,2,1,2,1,2,1,1,2,2,1,1,2,2,2,2,2,1,1,2,1,1,1,2,2],[1.01785714285714,0.857142857142857,0.4,0.541666666666667,1.03571428571429,0.821428571428571,1.16071428571429,0.854166666666667,0.869047619047619,1.25,1.13095238095238,0.68452380952381,0.511904761904762,1.13690476190476,1.13095238095238,0.518518518518518,0.732142857142857,1.29761904761905,0.357142857142857,1.22222222222222,1.41071428571429,0.767857142857143,0.273809523809524,0.791666666666667,0.845238095238095,1.25,1.01785714285714,1.17857142857143,0.595238095238095,1.54761904761905,0.375,0.869047619047619,0.839285714285714,0.404761904761905,1.27380952380952,1.41666666666667,0.755952380952381,1.20238095238095,1.5448717948718,1.04166666666667,1.8452380952381,2.47619047619048,2.32666666666667,2.76785714285714,1.97619047619048,2.53571428571429,2.33928571428571,2.38888888888889,2.54761904761905,2.92857142857143,2.89880952380952,2.70238095238095,2.5952380952381,2.56547619047619,3.2202380952381,2.5679012345679,1.86309523809524,2.89880952380952,2.11904761904762,2.91666666666667,2.875,2.625,2.45833333333333,1.95833333333333,2.375,2.80357142857143,3.19047619047619,2.42261904761905,3.0952380952381,2.69642857142857,2.57738095238095,2.5,2.19047619047619,1.95833333333333,3.00595238095238,3.01190476190476,2.01190476190476,2.5,2.53205128205128,2.64880952380952,1.61904761904762,0.744047619047619,0.566666666666667,0.7738095238095239,0.5238095238095239,0.666666666666667,1.00595238095238,1.28472222222222,0.8988095238095239,1.75,0.869047619047619,0.755952380952381,1.14880952380952,0.636904761904762,1.0952380952381,0.944444444444444,0.732142857142857,0.886904761904762,0.861111111111111,1.61309523809524,0.839285714285714,0.261904761904762,0.482142857142857,0.529761904761905,1,1.125,0.410714285714286,0.732142857142857,1.19642857142857,0.916666666666667,0.386904761904762,0.369047619047619,1.38690476190476,1.44047619047619,1.27380952380952,1.33333333333333,1.48214285714286,1.71153846153846,1,1.01785714285714,0.857142857142857,0.4,0.541666666666667,1.03571428571429,0.821428571428571,1.16071428571429,0.854166666666667,0.869047619047619,1.25,1.13095238095238,0.68452380952381,0.511904761904762,1.13690476190476,1.13095238095238,0.518518518518518,0.732142857142857,1.29761904761905,0.357142857142857,1.22222222222222,1.41071428571429,0.767857142857143,0.273809523809524,0.791666666666667,0.845238095238095,1.25,1.01785714285714,1.17857142857143,0.595238095238095,1.54761904761905,0.375,0.869047619047619,0.839285714285714,0.404761904761905,1.27380952380952,1.41666666666667,0.755952380952381,1.20238095238095,1.5448717948718,1.04166666666667,1.8452380952381,2.47619047619048,2.32666666666667,2.76785714285714,1.97619047619048,2.53571428571429,2.33928571428571,2.38888888888889,2.54761904761905,2.92857142857143,2.89880952380952,2.70238095238095,2.5952380952381,2.56547619047619,3.2202380952381,2.5679012345679,1.86309523809524,2.89880952380952,2.11904761904762,2.91666666666667,2.875,2.625,2.45833333333333,1.95833333333333,2.375,2.80357142857143,3.19047619047619,2.42261904761905,3.0952380952381,2.69642857142857,2.57738095238095,2.5,2.19047619047619,1.95833333333333,3.00595238095238,3.01190476190476,2.01190476190476,2.5,2.53205128205128,2.64880952380952,1.61904761904762,0.744047619047619,0.566666666666667,0.7738095238095239,0.5238095238095239,0.666666666666667,1.00595238095238,1.28472222222222,0.8988095238095239,1.75,0.869047619047619,0.755952380952381,1.14880952380952,0.636904761904762,1.0952380952381,0.944444444444444,0.732142857142857,0.886904761904762,0.303571428571429,0.861111111111111,1.61309523809524,0.839285714285714,0.261904761904762,0.482142857142857,0.529761904761905,1,1.125,0.410714285714286,0.732142857142857,1.19642857142857,0.916666666666667,0.386904761904762,0.369047619047619,1.38690476190476,1.44047619047619,1.27380952380952,1.33333333333333,1.48214285714286,1.71153846153846,1],[7.73809523809524,1.78571428571429,0.666666666666667,3.57142857142857,5.95238095238095,0,1.78571428571429,2.77777777777778,5.35714285714286,10.1190476190476,10.7142857142857,8.928571428571431,4.16666666666667,1.19047619047619,9.52380952380952,1.85185185185185,4.76190476190476,8.928571428571431,2.38095238095238,0.694444444444444,2.38095238095238,2.38095238095238,1.19047619047619,1.19047619047619,4.16666666666667,2.97619047619048,16.6666666666667,10.1190476190476,6.54761904761905,3.57142857142857,7.73809523809524,4.76190476190476,1.19047619047619,4.76190476190476,7.14285714285714,10.7142857142857,3.57142857142857,12.5,7.69230769230769,10.1190476190476,25.5952380952381,10.1190476190476,26,11.9047619047619,10.1190476190476,14.8809523809524,28.5714285714286,22.9166666666667,27.3809523809524,26.7857142857143,40.4761904761905,22.6190476190476,20.8333333333333,16.0714285714286,39.8809523809524,33.9506172839506,13.0952380952381,17.2619047619048,15.4761904761905,19.4444444444444,18.452380952381,13.6904761904762,17.2619047619048,11.9047619047619,24.4047619047619,30.952380952381,27.3809523809524,35.1190476190476,35.7142857142857,27.9761904761905,33.3333333333333,20.8333333333333,30.3571428571429,23.2142857142857,26.1904761904762,39.8809523809524,32.7380952380952,36.3095238095238,27.5641025641026,25.5952380952381,15.4761904761905,5.35714285714286,5.33333333333333,4.76190476190476,8.928571428571431,4.76190476190476,7.73809523809524,16.6666666666667,2.38095238095238,13.0952380952381,16.6666666666667,10.1190476190476,5.35714285714286,3.57142857142857,11.9047619047619,1.23456790123457,1.78571428571429,8.33333333333333,8.33333333333333,11.9047619047619,2.97619047619048,1.78571428571429,10.1190476190476,10.7142857142857,8.928571428571431,20.8333333333333,5.95238095238095,8.928571428571431,11.9047619047619,1.19047619047619,6.54761904761905,6.54761904761905,5.35714285714286,7.73809523809524,13.6904761904762,8.928571428571431,10.1190476190476,12.8205128205128,8.928571428571431,7.73809523809524,1.78571428571429,0.666666666666667,3.57142857142857,5.95238095238095,0,1.78571428571429,2.77777777777778,5.35714285714286,10.1190476190476,10.7142857142857,8.928571428571431,4.16666666666667,1.19047619047619,9.52380952380952,1.85185185185185,4.76190476190476,8.928571428571431,2.38095238095238,0.694444444444444,2.38095238095238,2.38095238095238,1.19047619047619,1.19047619047619,4.16666666666667,2.97619047619048,16.6666666666667,10.1190476190476,6.54761904761905,3.57142857142857,7.73809523809524,4.76190476190476,1.19047619047619,4.76190476190476,7.14285714285714,10.7142857142857,3.57142857142857,12.5,7.69230769230769,10.1190476190476,25.5952380952381,10.1190476190476,26,11.9047619047619,10.1190476190476,14.8809523809524,28.5714285714286,22.9166666666667,27.3809523809524,26.7857142857143,40.4761904761905,22.6190476190476,20.8333333333333,16.0714285714286,39.8809523809524,33.9506172839506,13.0952380952381,17.2619047619048,15.4761904761905,19.4444444444444,18.452380952381,13.6904761904762,17.2619047619048,11.9047619047619,24.4047619047619,30.952380952381,27.3809523809524,35.1190476190476,35.7142857142857,27.9761904761905,33.3333333333333,20.8333333333333,30.3571428571429,23.2142857142857,26.1904761904762,39.8809523809524,32.7380952380952,36.3095238095238,27.5641025641026,25.5952380952381,15.4761904761905,5.35714285714286,5.33333333333333,4.76190476190476,8.928571428571431,4.76190476190476,7.73809523809524,16.6666666666667,2.38095238095238,13.0952380952381,16.6666666666667,10.1190476190476,5.35714285714286,3.57142857142857,11.9047619047619,1.23456790123457,1.78571428571429,8.33333333333333,7.14285714285714,8.33333333333333,11.9047619047619,2.97619047619048,6.54761904761905,4.16666666666667,3.57142857142857,8.928571428571431,20.8333333333333,8.33333333333333,8.928571428571431,5.95238095238095,4.16666666666667,2.38095238095238,17.8571428571429,10.1190476190476,20.8333333333333,13.6904761904762,8.928571428571431,10.1190476190476,12.8205128205128,8.928571428571431],[24,23,25,21,27,23,28,23,27,25,21,26,21,23,22,26,26,24,25,33,29,22,25,23,20,21,25,24,24,29,24,24,19,21,25,23,27,24,27,22,24,23,25,21,27,23,28,23,27,25,21,26,21,23,22,26,26,24,25,33,29,22,25,23,20,21,25,24,24,29,24,24,19,21,25,23,27,24,27,22,24,23,25,21,27,23,28,23,27,25,21,26,21,23,22,26,26,24,33,29,22,22,23,22,21,25,30,24,20,25,24,24,26,26,23,27,24,27,22,24,23,25,21,27,23,28,23,27,25,21,26,21,23,22,26,26,24,25,33,29,22,25,23,20,21,25,24,24,29,24,24,19,21,25,23,27,24,27,22,24,23,25,21,27,23,28,23,27,25,21,26,21,23,22,26,26,24,25,33,29,22,25,23,20,21,25,24,24,29,24,24,19,21,25,23,27,24,27,22,24,23,25,21,27,23,28,23,27,25,21,26,21,23,22,26,26,24,25,33,29,22,25,23,20,21,25,24,24,29,24,24,19,21,25,23,27,24,27,22],["Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","Baseline","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3"],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3],["Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva","Saliva"],[1,2,1,2,2,1,2,2,2,1,2,1,1,2,2,2,2,1,2,1,2,1,1,2,2,1,1,2,2,2,2,2,1,1,2,1,1,1,2,2,1,2,1,2,2,1,2,2,2,1,2,1,1,2,2,2,2,1,2,1,2,1,1,2,2,1,1,2,2,2,2,2,1,1,2,1,1,1,2,2,1,2,1,2,2,1,2,2,2,1,2,1,1,2,2,2,2,1,1,2,1,2,1,2,1,1,2,2,1,1,2,2,2,2,1,1,1,2,1,1,2,1,2,2,1,2,2,2,1,2,1,1,2,2,2,2,1,2,1,2,1,1,2,2,1,1,2,2,2,2,2,1,1,2,1,1,1,2,2,1,2,1,2,2,1,2,2,2,1,2,1,1,2,2,2,2,1,2,1,2,1,1,2,2,1,1,2,2,2,2,2,1,1,2,1,1,1,2,2,1,2,1,2,2,1,2,2,2,1,2,1,1,2,2,2,2,1,2,1,2,1,1,2,2,1,1,2,2,2,2,2,1,1,2,1,1,1,2,2],["Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP1","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP2","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Plaque_TP3","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP1","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP2","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3","Saliva_TP3"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ID<\/th>\n      <th>ID_external<\/th>\n      <th>Sample<\/th>\n      <th>Time<\/th>\n      <th>Subject<\/th>\n      <th>Group<\/th>\n      <th>Quant<\/th>\n      <th>Full_name<\/th>\n      <th>group<\/th>\n      <th>sex<\/th>\n      <th>mean_plaque<\/th>\n      <th>mean_bleeding<\/th>\n      <th>age<\/th>\n      <th>Period<\/th>\n      <th>time<\/th>\n      <th>Sample_Type<\/th>\n      <th>Sex<\/th>\n      <th>Sample_Time<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,7,9,10,11,12,13,15,17]},{"orderable":false,"targets":0},{"name":" ","targets":0},{"name":"ID","targets":1},{"name":"ID_external","targets":2},{"name":"Sample","targets":3},{"name":"Time","targets":4},{"name":"Subject","targets":5},{"name":"Group","targets":6},{"name":"Quant","targets":7},{"name":"Full_name","targets":8},{"name":"group","targets":9},{"name":"sex","targets":10},{"name":"mean_plaque","targets":11},{"name":"mean_bleeding","targets":12},{"name":"age","targets":13},{"name":"Period","targets":14},{"name":"time","targets":15},{"name":"Sample_Type","targets":16},{"name":"Sex","targets":17},{"name":"Sample_Time","targets":18}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```



``` r
ps_up %>% 
  generate_color_palette(var = "Time",
                         pal = "jama",
                         print = TRUE) -> time_pal
```

<img src="01_metaphlan_import_files/figure-html/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

``` r
ps_up %>% 
  generate_color_palette(var = "Period",
                         pal = "jama",
                         print = TRUE) -> period_pal

ps_up %>% 
  generate_color_palette(var = "Sample_Type",
                         pal = "jco",
                         print = TRUE) -> sample_pal
```

<img src="01_metaphlan_import_files/figure-html/unnamed-chunk-8-2.png" style="display: block; margin: auto;" />

``` r
ps_up %>% 
  generate_color_palette(var = "sex",
                         pal = "uchicago",
                         print = TRUE) -> sex_pal
```

<img src="01_metaphlan_import_files/figure-html/unnamed-chunk-8-3.png" style="display: block; margin: auto;" />

``` r
ps_up %>% 
  generate_color_palette(var = "Subject",
                         pal = "randomcoloR",
                         print = TRUE, runTsne = TRUE) -> sub_pal
```

<img src="01_metaphlan_import_files/figure-html/unnamed-chunk-8-4.png" style="display: block; margin: auto;" />

# Explore metadata -> to add cluster:


``` r
pd <- position_dodge(0)


ps_up %>% 
  sample_data() %>% 
  data.frame() %>% 
  rownames_to_column("sample_id") -> meta

meta %>% 
  select(sample_id, Sample, Subject, Time, time, starts_with("mean")) %>% 
  pivot_longer(cols = starts_with("mean")) -> meta_long


meta_long %>% 
  filter(Sample == "Plaque") %>% 
  ggplot(data = ., aes(x=Time, y=value)) +
  geom_point(aes(color = Subject),  position = pd) +
  geom_line(aes(group=Subject), position = pd, linetype = "dashed", color = "black", linewidth = 0.08) +
  geom_boxplot(aes(group=Time, fill = NULL, color = NULL), outlier.shape = NA, alpha = 0.4) +
  # ggtitle("Plaque development during sugar rinsing  \n & oral hygiene interuption") +
  facet_grid(name ~ Sample, scales = "free_y", switch = "y") +
  theme(strip.placement = "outside") +
  # ylab("Mean plaque (Quigley Hein plaque index)") + xlab("Time of sampling") + 
  scale_x_discrete(breaks = c("TP1", "TP2", "TP3"), labels = c("Baseline", 
                                                               "2 weeks treatment", "4 weeks")) + theme_linedraw() + theme(strip.placement = "outside") + xlab(NULL) + ylab(NULL) + theme(legend.position = "none") + scale_color_manual(values = sub_pal) -> p_meta

p_meta
```

<img src="01_metaphlan_import_files/figure-html/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />


``` r
meta_long %>% 
  group_by(name, Subject) %>% 
  arrange(time) %>%
  mutate(delta = value - first(value)) -> meta_long_delta

meta_long_delta %>% 
  filter(Sample == "Plaque") %>% 
  ggplot(data = ., aes(x=Time, y=delta)) +
  geom_point(aes(color = Subject),  position = pd) +
  geom_line(aes(group=Subject), position = pd, linetype = "dashed", color = "black", linewidth = 0.08) +
  geom_boxplot(aes(group=Time, fill = NULL, color = NULL), outlier.shape = NA, alpha = 0.4) +
  # ggtitle("Plaque development during sugar rinsing  \n & oral hygiene interuption") +
  facet_grid(name ~ Sample, scales = "free_y", switch = "y") +
  theme(strip.placement = "outside") +
  # ylab("Mean plaque (Quigley Hein plaque index)") + xlab("Time of sampling") + 
  scale_x_discrete(breaks = c("TP1", "TP2", "TP3"), labels = c("Baseline", 
                                                               "2 weeks treatment", "4 weeks")) + theme_linedraw() + theme(strip.placement = "outside") + xlab(NULL) + ylab(NULL) + theme(legend.position = "none") + scale_color_manual(values = sub_pal) -> p_meta_delta

p_meta_delta
```

<img src="01_metaphlan_import_files/figure-html/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />


``` r
meta_long %>% 
  group_by(name, Subject) %>% 
  arrange(time) %>%
  mutate(ratio = value / first(value)) -> meta_long_ratio

meta_long_ratio %>% 
  filter(Sample == "Plaque") %>% 
  ggplot(data = ., aes(x=Time, y=ratio)) +
  geom_point(aes(color = Subject),  position = pd) +
  geom_line(aes(group=Subject), position = pd, linetype = "dashed", color = "black", linewidth = 0.08) +
  geom_boxplot(aes(group=Time, fill = NULL, color = NULL), outlier.shape = NA, alpha = 0.4) +
  # ggtitle("Plaque development during sugar rinsing  \n & oral hygiene interuption") +
  facet_grid(name ~ Sample, scales = "free_y", switch = "y") +
  theme(strip.placement = "outside") +
  # ylab("Mean plaque (Quigley Hein plaque index)") + xlab("Time of sampling") + 
  scale_x_discrete(breaks = c("TP1", "TP2", "TP3"), labels = c("Baseline", 
                                                               "2 weeks treatment", "4 weeks")) + theme_linedraw() + theme(strip.placement = "outside") + xlab(NULL) + ylab(NULL) + theme(legend.position = "none") + scale_color_manual(values = sub_pal) -> p_meta_ratio

p_meta_ratio
```

<img src="01_metaphlan_import_files/figure-html/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />


``` r
meta %>% 
  filter(Sample == "Plaque") %>% 
  select(Subject, Time, time, mean_plaque, mean_bleeding) %>% 
  pivot_longer(cols = mean_plaque:mean_bleeding) %>% 
  group_by(name, Subject) %>% 
  arrange(time, Subject) %>%
  mutate(delta = value- first(value)) %>% 
  select(-value) %>% 
  group_by(Subject) %>% 
  arrange(time) %>%
  pivot_wider(names_from = "name", values_from = "delta") %>% 
  ggplot(data = ., aes(x=mean_plaque, y=mean_bleeding)) +
  geom_point(aes(color = Time)) +
  geom_path(aes(group=Subject),  linetype = "dashed", color = "black", linewidth = 0.08, 
            arrow = arrow(
              angle = 30, length = unit(0.15, "inches"),
              ends = "last", type = "open"
            )) -> clin_xy_d

clin_xy_d
```

<img src="01_metaphlan_import_files/figure-html/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />


``` r
# meta %>% 
#   filter(Sample == "Plaque") %>% 
#   select(Subject, Time, time, mean_plaque, mean_bleeding) %>% 
#   pivot_longer(cols = mean_plaque:mean_bleeding) %>% 
#   group_by(name, Subject) %>% 
#   arrange(time, Subject) %>%
#   mutate(delta = value- first(value)) %>% 
#   

# meta_long_delta %>% 
#   dplyr::filter(name == "mean_bleeding") -> bleeding_only
# 
# bleeding_only %>% 
#   pull(delta) %>% 
#   cut_number(., n = 3, labels = c("Low", "Medium", "High")) %>% 
#   cbind(.,
#         bleeding_only)
```


``` r
meta %>% 
  filter(Sample == "Plaque") %>% 
  select(Subject, Time, time, mean_plaque, mean_bleeding) %>% 
  pivot_longer(cols = mean_plaque:mean_bleeding) %>% 
  group_by(name, Subject) %>% 
  arrange(time, Subject) %>%
  pivot_wider(names_from = "name", values_from = "value") %>% 
  ggplot(data = ., aes(x=mean_plaque, y=mean_bleeding)) +
  geom_point(aes(color = Time)) +
  geom_path(aes(group=Subject),  linetype = "dashed", color = "black", linewidth = 0.08, 
            arrow = arrow(
              angle = 30, length = unit(0.15, "inches"),
              ends = "last", type = "open"
            )) -> clin_xy
clin_xy
```

<img src="01_metaphlan_import_files/figure-html/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />


``` r
meta %>% 
  filter(Sample == "Plaque") %>% 
  select(mean_plaque, mean_bleeding) -> tmp

tmp %>% 
  kmeans(.,3,  nstart = 25) -> km


factoextra::fviz_cluster(km, data = tmp,
                         # palette = c("#2E9FDF"), 
                         geom = "point",
                         ellipse.type = "convex", 
                         ggtheme = theme_bw()
) -> kall

kall$data %>% 
  rename(cluster_all = cluster) -> kall$data

kall
```

<img src="01_metaphlan_import_files/figure-html/unnamed-chunk-15-1.png" style="display: block; margin: auto;" />


``` r
# clin_xy_d$data %>% 
#   left_join(., 
#             result_clust,
#             by = c("Subject" = "name")) %>% 
#   # pivot_wider(names_from = "name", values_from = "delta") %>% 
#   ggplot(data = ., aes(x=mean_plaque, y=mean_bleeding)) +
#   geom_point(aes(color = cluster_Dtp3, shape = cluster_Dtp3)) +
#   geom_path(aes(group=Subject),  linetype = "dashed", color = "black", linewidth = 0.08, 
#             arrow = arrow(
#               angle = 30, length = unit(0.15, "inches"),
#               ends = "last", type = "open"
#             )) -> pclust_3
# 
# pclust_3
```


``` r
meta %>% 
  filter(Sample == "Plaque", 
         Period == "Baseline") %>% 
  rownames_to_column("toto") %>% 
  column_to_rownames("Subject") %>% 
  select(mean_plaque, mean_bleeding) -> tmp

tmp %>% 
  kmeans(.,3,  nstart = 25) -> km


factoextra::fviz_cluster(km, data = tmp ,
                         # palette = c("#2E9FDF"), 
                         geom = "text", repel = TRUE, show.clust.cent = FALSE,
                         ellipse.type = "convex", 
                         ggtheme = theme_bw()
) + ggtitle("Baseline") -> kbase

kbase$data %>% 
  rename(cluster_base = cluster) -> kbase$data

kbase
```

<img src="01_metaphlan_import_files/figure-html/unnamed-chunk-17-1.png" style="display: block; margin: auto;" />


``` r
meta_long_delta %>% 
  filter(Sample == "Plaque", 
         Time == "TP2") %>% 
  select(-value) %>% 
  pivot_wider(names_from = "name", values_from = "delta") %>% 
  column_to_rownames("Subject") %>% 
  select(mean_plaque, mean_bleeding) -> tmp

tmp %>% 
  kmeans(.,3,  nstart = 25) -> km


factoextra::fviz_cluster(km, data = tmp ,
                         # palette = c("#2E9FDF"), 
                         geom = "text", repel = TRUE, show.clust.cent = FALSE,
                         ellipse.type = "convex", 
                         ggtheme = theme_bw()
) + ggtitle("Delta TP2") -> kdtp2

kdtp2$data %>% 
  rename(cluster_Dtp2 = cluster) -> kdtp2$data
```


``` r
meta_long_ratio %>% 
  filter(Sample == "Plaque", 
         Time == "TP2") %>% 
  select(-value) %>% 
  pivot_wider(names_from = "name", values_from = "ratio") %>% 
  mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  column_to_rownames("Subject") %>% 
  select(mean_plaque, mean_bleeding) %>% 
  drop_na(.) -> tmp

tmp %>% 
  kmeans(.,3,  nstart = 25) -> km


factoextra::fviz_cluster(km, data = tmp ,
                         # palette = c("#2E9FDF"), 
                         geom = "text", repel = TRUE, show.clust.cent = FALSE,
                         ellipse.type = "convex", 
                         ggtheme = theme_bw()
) + ggtitle("Ratio TP2") -> krtp2

krtp2$data %>% 
  rename(cluster_Rtp2 = cluster) -> krtp2$data

kdtp2
```

<img src="01_metaphlan_import_files/figure-html/unnamed-chunk-19-1.png" style="display: block; margin: auto;" />


``` r
meta %>% 
  filter(Sample == "Plaque", 
         Period == "TP2") %>% 
  rownames_to_column("toto") %>% 
  column_to_rownames("Subject") %>% 
  select(mean_plaque, mean_bleeding) -> tmp

tmp %>% 
  kmeans(.,3,  nstart = 25) -> km


factoextra::fviz_cluster(km, data = tmp ,
                         # palette = c("#2E9FDF"), 
                         geom = "text", repel = TRUE, show.clust.cent = FALSE,
                         ellipse.type = "convex", 
                         ggtheme = theme_bw()
) + ggtitle("TP2") -> ktp2

ktp2$data %>% 
  rename(cluster_tp2 = cluster) -> ktp2$data

ktp2
```

<img src="01_metaphlan_import_files/figure-html/unnamed-chunk-20-1.png" style="display: block; margin: auto;" />


``` r
meta_long_delta %>% 
  filter(Sample == "Plaque", 
         Time == "TP3") %>% 
  select(-value) %>% 
  pivot_wider(names_from = "name", values_from = "delta") %>% 
  column_to_rownames("Subject") %>% 
  select(mean_plaque, mean_bleeding) -> tmp

tmp %>% 
  kmeans(.,3,  nstart = 25) -> km


factoextra::fviz_cluster(km, data = tmp ,
                         # palette = c("#2E9FDF"), 
                         geom = "text", repel = TRUE, show.clust.cent = FALSE,
                         ellipse.type = "convex", 
                         ggtheme = theme_bw()
) + ggtitle("TP3")  -> kdtp3

kdtp3$data %>% 
  rename(cluster_Dtp3 = cluster) -> kdtp3$data

kdtp3
```

<img src="01_metaphlan_import_files/figure-html/unnamed-chunk-21-1.png" style="display: block; margin: auto;" />


``` r
kdtp3$data %>% 
  select(-x,-y) %>% 
  full_join(., kdtp2$data %>% 
              select(-x,-y) ) -> result_clust

clin_xy_d$data %>% 
  left_join(., 
            result_clust,
            by = c("Subject" = "name")) %>% 
  # pivot_wider(names_from = "name", values_from = "delta") %>% 
  ggplot(data = ., aes(x=mean_plaque, y=mean_bleeding)) +
  geom_point(aes(color = cluster_Dtp3, shape = cluster_Dtp3)) +
  geom_path(aes(group=Subject),  linetype = "dashed", color = "black", linewidth = 0.08, 
            arrow = arrow(
              angle = 30, length = unit(0.15, "inches"),
              ends = "last", type = "open"
            )) -> pclust_3

pclust_3
```

<img src="01_metaphlan_import_files/figure-html/unnamed-chunk-22-1.png" style="display: block; margin: auto;" />


``` r
clin_xy_d$data %>% 
  left_join(., 
            result_clust,
            by = c("Subject" = "name")) %>% 
  # pivot_wider(names_from = "name", values_from = "delta") %>% 
  ggplot(data = ., aes(x=mean_plaque, y=mean_bleeding)) +
  geom_point(aes(color = cluster_Dtp2, shape = cluster_Dtp2)) +
  geom_path(aes(group=Subject),  linetype = "dashed", color = "black", linewidth = 0.08, 
            arrow = arrow(
              angle = 30, length = unit(0.15, "inches"),
              ends = "last", type = "open"
            )) -> pclust_2

pclust_2
```

<img src="01_metaphlan_import_files/figure-html/unnamed-chunk-23-1.png" style="display: block; margin: auto;" />


``` r
meta_long_delta %>% 
  select(-value) %>% 
  pivot_wider(names_from = "name", values_from = "delta") %>% 
  rename_at(vars(starts_with('mean')), funs(paste0('delta_', .))) %>% 
  select(sample_id, delta_mean_plaque, delta_mean_bleeding) %>% 
  full_join(., 
            meta_long_ratio %>% 
              select(-value) %>% 
              pivot_wider(names_from = "name", values_from = "ratio") %>% 
              rename_at(vars(starts_with('mean')), funs(paste0('ratio_', .))) %>% 
              select(sample_id, ratio_mean_plaque, ratio_mean_bleeding)) -> tmp


ps_up %>% 
  sample_data() %>% 
  data.frame() %>% 
  rownames_to_column("sample_id") %>% 
  left_join(.,
            tmp) %>% 
  column_to_rownames("sample_id") -> sample_data(ps_up)


ps_up %>% 
  sample_data() %>% 
  data.frame() %>% 
  rownames_to_column("tmp") %>% 
  left_join(.,
            result_clust,
            by = c("Subject" = "name")) %>% 
  column_to_rownames("tmp") -> sample_data(ps_up)
```

Maybe also add the delta and ratio to all samples.


``` r
# ggscatter(
#   ind.coord, x = "Dim.1", y = "Dim.2", 
#   color = "cluster", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
#   shape = "Species", size = 1.5,  legend = "right", ggtheme = theme_bw(),
#   xlab = paste0("Dim 1 (", variance.percent[1], "% )" ),
#   ylab = paste0("Dim 2 (", variance.percent[2], "% )" )
# ) +
#   stat_mean(aes(color = cluster), size = 4)
```


``` r
ps_up %>% 
   microViz::ps_mutate(cluster_Dtp2 = as.factor(cluster_Dtp2)) -> ps_up
```



``` r
save(period_pal,
     sample_pal,
     sex_pal,
     sub_pal,
     time_pal,
     p_meta_delta, p_meta,p_meta_ratio, 
     ps_up, clin_xy, clin_xy_d, pclust_3, pclust_2, file = here::here("../../data/processed_data/metaphlan/01_data.Rdata"))
```



``` r
saveRDS(object = ps_up, file = here::here("../../data/processed_data/metaphlan/metaphlan4.1.1_phyloseq.RDS"))
```




``` r
sessionInfo()
```

```
## R version 4.4.0 (2024-04-24)
## Platform: aarch64-apple-darwin20
## Running under: macOS Sonoma 14.6.1
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## time zone: Europe/Zurich
## tzcode source: internal
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] phyloseq_1.50.0 readxl_1.4.3    lubridate_1.9.3 forcats_1.0.0  
##  [5] stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2     readr_2.1.5    
##  [9] tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0
## 
## loaded via a namespace (and not attached):
##   [1] permute_0.9-7           rlang_1.1.4             magrittr_2.0.3         
##   [4] ade4_1.7-22             compiler_4.4.0          mgcv_1.9-1             
##   [7] vctrs_0.6.5             reshape2_1.4.4          pkgconfig_2.0.3        
##  [10] crayon_1.5.3            fastmap_1.2.0           backports_1.5.0        
##  [13] XVector_0.45.0          labeling_0.4.3          utf8_1.2.4             
##  [16] rmarkdown_2.29          tzdb_0.4.0              UCSC.utils_1.2.0       
##  [19] microViz_0.12.4         xfun_0.49               zlibbioc_1.51.2        
##  [22] cachem_1.1.0            GenomeInfoDb_1.42.0     jsonlite_1.8.9         
##  [25] biomformat_1.34.0       highr_0.11              rhdf5filters_1.17.0    
##  [28] Rhdf5lib_1.27.0         broom_1.0.7             parallel_4.4.0         
##  [31] cluster_2.1.6           R6_2.5.1                bslib_0.8.0            
##  [34] stringi_1.8.4           car_3.1-3               jquerylib_0.1.4        
##  [37] cellranger_1.1.0        Rcpp_1.0.13-1           iterators_1.0.14       
##  [40] knitr_1.48              IRanges_2.39.2          Matrix_1.7-1           
##  [43] splines_4.4.0           igraph_2.1.1            timechange_0.3.0       
##  [46] tidyselect_1.2.1        rstudioapi_0.17.1       abind_1.4-8            
##  [49] yaml_2.3.10             vegan_2.6-8             codetools_0.2-20       
##  [52] curl_6.0.0              lattice_0.22-6          plyr_1.8.9             
##  [55] Biobase_2.66.0          withr_3.0.2             Rtsne_0.17             
##  [58] evaluate_1.0.1          survival_3.7-0          Biostrings_2.73.2      
##  [61] ggpubr_0.6.0            pillar_1.9.0            carData_3.0-5          
##  [64] DT_0.33                 foreach_1.5.2           stats4_4.4.0           
##  [67] speedyseq_0.5.3.9021    generics_0.1.3          rprojroot_2.0.4        
##  [70] S4Vectors_0.43.2        hms_1.1.3               munsell_0.5.1          
##  [73] scales_1.3.0            glue_1.8.0              tools_4.4.0            
##  [76] data.table_1.16.2       ggsignif_0.6.4          rhdf5_2.49.0           
##  [79] grid_4.4.0              ape_5.8                 crosstalk_1.2.1        
##  [82] colorspace_2.1-1        nlme_3.1-166            GenomeInfoDbData_1.2.13
##  [85] Formula_1.2-5           cli_3.6.3               fansi_1.0.6            
##  [88] V8_6.0.0                gtable_0.3.6            ggsci_3.2.0            
##  [91] rstatix_0.7.2           sass_0.4.9              digest_0.6.37          
##  [94] BiocGenerics_0.52.0     ggrepel_0.9.6           farver_2.1.2           
##  [97] htmlwidgets_1.6.4       factoextra_1.0.7        htmltools_0.5.8.1      
## [100] multtest_2.61.0         lifecycle_1.0.4         httr_1.4.7             
## [103] here_1.0.1              MASS_7.3-61
```
