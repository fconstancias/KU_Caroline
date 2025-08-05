---
title: " Study II - Beta-diversity analyses of metaphlan profiles"
author: "Florentin Constancias"
date: "December 30, 2024"
output: 
  html_document: 
    toc: yes
    keep_md: yes
---










``` r
source("https://raw.githubusercontent.com/fconstancias/KU_Caroline/refs/heads/main/code/functions/phyloseq_functions.R")

source_dir = "https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/"
# source_dir = "~/Documents/GitHub/DivComAnalyses/R/"

source(paste0(source_dir,"phyloseq_normalisation.R"))
source(paste0(source_dir,"phyloseq_varia.R"))
source(paste0(source_dir,"phyloseq_beta.R"))
```



``` r
load(here::here("../../data/processed_data/metaphlan/01_data.Rdata"))
load(here::here("../../data/processed_data/metaphlan/02_data.Rdata"))
```


``` r
pd <- position_dodge(0.3)
```


# Compute the distance:


``` r
ps_up %>%
  subset_taxa(Class != "UNCLASSIFIED") %>% 
  transform_sample_counts(function(x) x/sum(x) * 1) %>%
  phyloseq_compute_bdiv() -> beta


ps_up %>% 
  subset_taxa(Class != "UNCLASSIFIED") %>% 
  transform_sample_counts(function(x) x/sum(x) * 1) %>%
  microViz::dist_calc(., dist = "robust.aitchison") %>% 
  microViz::dist_get() %>% 
  magrittr::divide_by(100) -> beta$rAitchison

beta$bray <- NULL
beta$sorensen <- NULL
```


# Example1 - all samples:


``` r
ps_up %>%
  subset_taxa(Class != "UNCLASSIFIED") %>%
  # subset_samples(Sample == "Saliva" &
  # Time == "TP1") %>% 
  transform_sample_counts(function(x) x/sum(x) * 1) %>% 
  phyloseq_explore_beta(ps_up = .,
                        beta = beta,
                        color_group = "Time",
                        shape_group = "Sample",
                        distance_for_more = "rAitchison", #distance for detailed plot (path, and envfit vectors)
                        tax_rank_fit  = "Genus", # taxonomic level for envfit vector fitting
                        metadata_sel = c("delta_mean_plaque", "delta_mean_bleeding"), # metadata selection
                        pval_cutoff = 0.05, # pvalue cutoff 
                        top_r = 12, # max number of taxa or metadata to be displayed as vector of the ordination
                        alpha = NULL,
                        col_pal = time_pal,
                        fill_pal = time_pal,
                        path_group = "interaction(Sample,Subject)",
                        facet_formula = "Sample ~ .",
                        axis1 = 1,
                        axis2 = 2,
                        seed = 123,
                        perm = 999,
                        permanova_terms = c("Subject","Time" ,"cluster_Dtp2"),
                        metadata_dist_boxplot = c("Subject", "Time", "Sample"),
                        strata = "none") -> out
```

Focus on Robust Aitchinson distance, clear discrimination between site, display trajectories.


``` r
out$PCOA$p + facet_null() + stat_ellipse(data=out$PCOA$p$data,
                                         segments = 21, linetype = 2, color = "grey80",
                                         aes(group = interaction(time, Sample)),
)
```

<img src="emvit_demo_files/figure-html/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />
Corresponding legend:


``` r
out$PCOA$p_leg
```

<img src="emvit_demo_files/figure-html/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

Genus with diff in proportions fitting ordination:


``` r
out$env_fit_tax$envfit %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item" id="htmlwidget-9eceb20dc2a931752ae4" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-9eceb20dc2a931752ae4">{"x":{"filter":"none","vertical":false,"data":[["Corynebacterium","Veillonella","Actinomyces","Selenomonas","Candidatus_Nanosynsacchari","Candidatus_Saccharibacteria_unclassified","Porphyromonas","Capnocytophaga","GGB1472","Pauljensenia","Prevotella","GGB3883","Tannerella","Schaalia","Lachnoanaerobaculum","Aggregatibacter","Streptococcus","Cardiobacterium","Actinobaculum","Kingella","GGB4308","Rothia","Campylobacter","Granulicatella","Leptotrichia","Olsenella","Gemella","Ottowia","Megasphaera","Eikenella","Candidatus_Nanoperiomorbus","Arachnia","Fusobacterium","GGB12798","GGB6813","GGB1144","Neisseria","Abiotrophia","GGB4303","Haemophilus","Eubacteriales_Family_XIII_Incertae_Sedis_unclassified","Mogibacterium","Lancefieldella","Oribacterium","Centipeda","GGB1202","Lautropia","Dialister","Peptostreptococcus","Treponema","Parvimonas","GGB12785","Alloprevotella","Catonella","Delftia","GGB6679","GGB3887","Ralstonia","GGB1022","Candidatus_Nanosynbacter","GGB1838","Anaeroglobus","Candidatus_Absconditabacteria_unclassified","Bacteroidota_unclassified","GGB4299","GGB12794","GGB3385","GGB12443","GGB3059","Johnsonella","GGB74353","GGB4400","GGB1833","GGB2676","Fretibacterium","GGB71270","GGB1025","GGB10852","GGB12783","GGB12788","Hallella","GGB2663","Shuttleworthia","Bifidobacterium","GGB3388","GGB2671","GGB96534","GGB4937","GGB4538","GGB2666","Solobacterium","GGB49434","Candidatus_Gracilibacteria_unclassified","GGB12789","Filifactor","GGB1843","GGB4936","Lachnospiraceae_unclassified","Peptidiphaga","Metamycoplasma","GGB2964","GGB49219","Peptostreptococcaceae_unclassified","GGB18703","GGB4733","GGB1460","GGB2675","GGB45697","GGB1088","GGB97311","GGB56136","GGB72530","GGB71076","Slackia","GGB3386","GGB73508","Colibacter","GGB71456","GGB6673","Desulfobulbus","Bulleidia","GGB49397","GGB4300","Cryptobacterium","GGB3387","Scardovia","GGB12761","GGB4533","GGB12787","GGB1473","GGB1024","Lactococcus","Cutibacterium","Propionibacterium","GGB4964","Candida","GGB97303","GGB12786","GGB1186","Peptoanaerobacter","GGB2674","Caldibacillus","Ureibacillus","GGB16822","Prevotellaceae_unclassified","GGB1832","GGB1611","GGB96297","GGB49499","GGB38873","GGB49401","GGB1844","GGB3389","GGB12441","GGB3886","GGB1203","Pseudoleptotrichia","GGB1840","Candidatus_Nanogingivalis","Simonsiella","GGB72444","Paracoccus","Lacticaseibacillus","Lactiplantibacillus","Escherichia","Bilophila","Desulfovibrio","Phocaeicola","GGB49229","Geobacillus","Dolosigranulum","GGB12763","GGB9835","GGB1026","Alloscardovia","GGB6688","Parascardovia","Sphingomonas","GGB3008","GGB49400","GGB49528","GGB4333","GGB71303","GGB2672","GGB70946","GGB72533","Staphylococcus","GGB4786","GGB4783","GGB4721","Limosilactobacillus","Lactobacillus","GGB1845","GGB1188","Pseudomonas","Sneathia","Acinetobacter","Klebsiella","GGB3390","Arcanobacterium","GGB4393","Latilactobacillus","GGB1021","Leuconostoc","Weissella","GGB49504","GGB12440","Methylobacterium","Ligilactobacillus","GGB96582","Moraxella"],[-0.5540780293977819,0.2279045799327089,-0.1753061452773274,-0.2952457621763102,-0.1523548787183606,0.01835683465677965,0.03476371246150621,-0.4467987895890816,-0.3107250372603324,-0.2298766491795372,0.3591110221717918,0.3588909871592711,-0.3405372694362842,0.708288283137362,-0.1691273602258908,-0.1772648420862203,0.5943250478396928,-0.3231606918478073,-0.3012549784404946,-0.2332533926014322,-0.05584031002948301,0.5335029507014287,-0.1444491544614101,0.7499999999999999,-0.3334240650753781,-0.1935958942945634,0.2708097714061311,-0.1111069658468919,0.4033976849511314,-0.2596401953417009,-0.1839634212155043,-0.1986888742949907,0.2510171657100136,-0.2015613790801842,-0.1942659491611694,0.3919555383613507,0.03865248077732793,0.07148861431418227,-0.2149337377028147,-0.01113536571327803,0.5961261509388655,0.5885196735605661,0.4283107086480816,0.6455098338629327,-0.1967633273098689,-0.1637837388676677,-0.1445916466563306,-0.1633259667282654,0.5266180442359145,-0.145398907614585,0.09341496744893135,0.06349651210009606,0.04430379788199677,0.3134096717209407,-0.1330693040760057,-0.1641286474969529,0.2196253951427577,-0.09185153138706449,-0.09583166620895606,0.06699334552624121,-0.04250609195877373,-0.1241075338951873,-0.1423824378091206,-0.0800402337756668,-0.1581300382942582,-0.2153030686471691,-0.3349045639115487,-0.0658838233490025,-0.04145670654397284,-0.1145688206240866,-0.07505382445652994,-0.106441781853049,-0.08641432737275306,-0.02288559246750531,-0.08388755247762279,-0.08780295468216649,-0.09409694740168691,-0.05390334148836298,-0.04951916624747599,-0.2133272521851012,-0.07537617280016021,-0.01184315738397235,-0.03047954746634298,0.0289753122232106,-0.09758031719159758,-0.05690808414733417,-0.04872660115892658,-0.07853879609678494,0.01319697552607622,0.01641790661856786,0.6502981593236344,-0.0722428436351103,-0.0674796241793265,-0.008889611505176106,0.0627593443081197,0.04726068499148139,0.208324398369602,0.08698316018599841,-0.3121399119302097,0.01480409945215325,-0.01941954191901596,-0.01703893160711058,-0.1011214424112836,-0.09749650145782604,0.0251742487801677,-0.04096073949402398,-0.004526307989118575,0.002815764403123371,0.1888111000792236,-0.1597242170775257,0.09273046216872827,-0.0525291832350651,-0.007059972419767796,0.00313328505136239,-0.0317929967901438,0.07276363933507346,-0.05591942424670988,-0.001159779847812128,-0.0888213107864592,0.02943789410501283,0.00505969427407308,-0.005566088484075875,-0.1011583045477178,-0.02003215473774271,-0.07265568776651671,-0.01610032254618015,-0.1120419761638661,-0.1245322233089615,-0.1129348342068335,-0.1346125025227591,0.1082080385423013,-0.0002546408260763283,0.0438934190965634,0.07412878106435788,0.2318320581987141,-0.04531275678606486,-0.02674970942832647,-0.06038482328235933,0.2111884096104969,-0.07415201578577031,-0.02003294707015352,-0.05414964804360995,-0.05414964804361009,-0.05414964804361042,0.2029247685705373,-0.03336247739983037,0.07258850294579566,-0.08063221302299937,-0.04779981756734467,-0.04779981756734482,0.1519686769175542,0.2615206581636716,-0.03625343624739909,0.2737812753871666,0.424011153507092,0.4276111735737686,-0.05303770437086722,0.02800221883405998,0.1722860799463662,-0.01282689433856336,-0.00519121046937164,-0.04464794942351025,0.0300936009964267,-0.04464794942350973,-0.07104756818614313,-0.03547134141344766,-0.03134660933927803,-0.06515881936951129,-0.06515881936951173,-0.01576584073742421,-0.04186317068108089,0.3236072481786635,0.2585111411035133,0.3182413530750004,0.1473913313048979,0.2330689782395448,0.07307557587815969,0.2822820140110451,0.1099988995372886,0.09187020046053504,0.06938516504593596,0.04895359846194557,0.1168594475393705,0.1635512440643989,0.131277415950104,0.1514352612583295,0.08305129647581441,0.1935583514946369,0.1742375619252017,0.09200875152891314,0.03774612283909622,0.04963781491162002,0.1842492831670363,0.08662446389489797,0.1006612643526089,0.06041003288224615,0.06041003288224581,0.03494767023185641,0.0914918699022768,0.07533953244654773,0.07214751732582196,0.08667881733013307,0.06116586250747289,0.04968035581755229,0.04812997829502174,0.03996026215918255,0.01650411841273085,0.04147454639498212,0.01753035009339722,0.06394395983836712,0.04467984199592583],[0.4711863320528276,0.4906191839555807,0.5642316952634435,0.4553436322587415,0.2686441740195496,0.1343650306536324,-0.1461293669991595,-0.4928884066462598,-0.1904092224342959,0.4991376688387579,0.4081294777891542,0.3783066090021979,0.3481290290924874,0.2612964353009416,0.4821267283772689,-0.3557451341450321,0.4100707201960312,-0.6575724394380765,0.2975833364841413,-0.3307165466134115,0.2042375323170274,-0.001350776602479835,-0.03800970649249735,0.05721059486902509,0.1829629709219408,0.1750388238396452,-0.2011365550542536,-0.2969100747592124,0.1508524973937061,-0.5518323924343039,-0.05832788611145563,-0.4158108107205086,-0.182749027217783,0.0446351532744713,-0.09223188421793724,0.03085743434834639,-0.593006480469307,-0.4664182737460034,0.1878537115355549,-0.2425726301764939,0.05209599810366202,0.2435890164103104,0.3719184730125327,0.2068738004824989,0.144120401431004,0.1984074138281985,-0.6088563662784848,0.111807934272951,-0.001745306308014921,-0.2353869973861842,-0.09829141925121367,0.0991028108613296,-0.1277793060072671,-0.2246150394845588,0.1956077577085638,-0.189867455331623,0.1817880120942792,0.1789876889826796,-0.2014041600592489,0.009412969988073233,-0.2330966275769888,-0.08009623999015282,-0.1301374246190254,-0.2546159987005876,-0.1524907399535991,-0.1434964706564296,0.01981619814930185,-0.2529718061124538,-0.1493130383664934,-0.2790258213259428,-0.200327224645371,-0.4759187088334421,-0.3549978579805098,-0.1064649591628752,-0.2329932519128339,-0.1637078017902702,0.01640159536088194,-0.1622600186932777,-0.1987101614331573,0.0007904716508395281,-0.08751233014933453,0.02267629116902297,-0.07297183977203561,0.05105210807533572,0.0884839389005153,-0.1612742573387365,-0.1595416446495342,-0.2584734156772157,0.08934187207402064,-0.1625178382393545,0.137536005418119,0.04595129018468012,-0.255342349012613,-0.1551215556830175,-0.1163167282313417,-0.3762116553638982,-0.1311738812458536,0.08362665595008763,-0.06129012054307697,0.0399772591784356,-0.1450699618576538,-0.09789002817765782,-0.3031922419750324,0.2133653337453016,-0.1031646654812753,-0.1169468677217729,-0.07356943969748876,0.03137547458931722,-0.07280236199517401,0.06010951959081456,0.02576587040839278,-0.154380214591762,-0.116830330847095,-0.07502521667249397,-0.1338192099317631,-0.06633598552845645,0.08736803392546777,-0.08887813997403449,-0.1652390426286107,-0.06809930396029612,-0.1458313948038896,-0.0861487407509048,-0.1117723677417937,0.1810104314072985,-0.1254640176594165,0.1573934269780023,0.07016152453686866,0.2109774137742828,-0.02009307821092239,-0.03422891261635936,-0.1624222477313268,-0.0171507462068981,-0.0136287591107806,0.1272724681591282,0.1858470269314777,0.1302040428630728,-0.1737837255362128,0.1188819963781385,0.02220798555300208,-0.1564510570850972,-0.1719556571679074,0.01590108009906631,0.01590108009906635,0.01590108009906644,-0.1392669452153945,-0.2037483843919676,0.05744823313532246,-0.1317951520618028,-0.1057434933681645,-0.1057434933681648,0.001174739353334076,0.1157609053547614,-0.2163172908497207,-0.04691812088248297,0.06642076132889631,-0.05896969094681616,-0.1207384577476907,-0.1019410935006941,-0.1078473266952964,-0.08625066407523437,-0.09992989649849746,0.02396572130666675,0.06550554102427869,0.02396572130666651,0.03920424558602968,-0.107803586266912,-0.01129795579804804,0.01714663059147487,0.01714663059147497,0.04779290314282051,0.09749371489980765,-0.06854633453198303,-0.1235197245587874,-0.1065068236572709,0.2201744994945991,-0.1677415429643651,-0.1128702552676697,0.1215033622121946,-0.02388725450789073,0.003351008951838237,-0.1201895690404989,-0.09100190281171772,0.02970103640004664,-0.08833470057011251,-0.1860484583581348,0.000193582933729295,0.03643889951028069,0.007366237639833279,0.03576852917380724,-0.1414580743821295,0.06251122290732432,0.0974479127879386,-0.01562269879680228,0.01359550595566045,0.0897663386604956,0.09219499301912971,0.09219499301912905,0.03855528337679811,-0.0009482427282687247,0.0006548836790358151,0.06190579855389794,0.001703854999083525,0.05337528584410019,0.08975110971080216,0.08959680418240276,-0.08802653709720618,-0.09335069357717275,-0.1086341060587932,0.0006141119650412863,-0.04150790252488342,0.04253564534938736],[0.001,0.001,0.001,0.001,0.001,0.151,0.096,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.003,0.001,0.09,0.001,0.001,0.002,0.001,0.001,0.001,0.001,0.01,0.001,0.001,0.005,0.005,0.001,0.001,0.001,0.001,0.002,0.001,0.001,0.001,0.001,0.002,0.001,0.001,0.011,0.001,0.001,0.101,0.223,0.135,0.001,0.005,0.002,0.001,0.011,0.006,0.628,0.002,0.075,0.019,0.002,0.007,0.002,0.001,0.001,0.053,0.001,0.006,0.001,0.001,0.228,0.002,0.014,0.36,0.021,0.012,0.007,0.203,0.9350000000000001,0.535,0.891,0.118,0.026,0.03,0.003,0.433,0.043,0.001,0.532,0.001,0.063,0.131,0.001,0.004,0.196,0.001,0.912,0.08799999999999999,0.363,0.001,0.001,0.308,0.126,0.638,0.907,0.008,0.029,0.388,0.033,0.242,0.586,0.074,0.345,0.324,0.435,0.019,0.5590000000000001,0.105,0.526,0.051,0.021,0.068,0.054,0.105,0.005,0.23,0.102,0.011,0.968,0.8129999999999999,0.079,0.001,0.096,0.016,0.12,0.003,0.031,0.025,0.897,0.897,0.897,0.002,0.007,0.42,0.049,0.17,0.17,0.083,0.001,0.005,0.001,0.001,0.001,0.109,0.309,0.007,0.478,0.35,0.927,0.618,0.927,0.535,0.181,0.913,0.772,0.772,0.93,0.279,0.001,0.001,0.001,0.002,0.001,0.116,0.001,0.242,0.441,0.109,0.355,0.183,0.025,0.006,0.06900000000000001,0.45,0.014,0.027,0.033,0.631,0.278,0.019,0.49,0.104,0.203,0.203,0.902,0.436,0.674,0.392,0.517,0.556,0.284,0.309,0.384,0.345,0.175,1,0.65,0.843],[0.4818750030024196,0.2665681125918806,0.3179802036172416,0.2682625412463161,0.08688167282390799,0.01675200968270299,0.02055164670557825,0.4031282709093615,0.1209706162000868,0.2750703091740158,0.2691939323867558,0.247686350915904,0.2160247689142898,0.5191566656493409,0.2377865121953303,0.1438991204747338,0.474916974405329,0.4889938568969323,0.1633309917691838,0.1491851216395175,0.04083594403784235,0.2592624049246249,0.02032209336008122,0.5153536636358991,0.1317566998055302,0.06204756515465942,0.1036530022877473,0.09154416873162757,0.1689563968567425,0.3387869545768293,0.033925588321764,0.1934498245194148,0.08781544115512176,0.03882122933689369,0.04212472626086106,0.1408056676110996,0.3216793490618965,0.2028143768332248,0.07422387034162195,0.05371071426634642,0.3261697167119702,0.3695374852643495,0.2930982785089344,0.4185327554120531,0.05418530290741479,0.06029197714071562,0.3567138806752811,0.03568514576005472,0.2526151260033857,0.06972624571107328,0.01674893490308667,0.01261863574787147,0.01666041364791596,0.1354280020324981,0.0509820310189227,0.05737464184109058,0.07403865754278267,0.03686648138428106,0.04531406339743065,0.0041688543972399,0.05113776585156479,0.01987374858288467,0.0338926373083908,0.06488750648022346,0.04395793084647925,0.06098063599853656,0.1025234172863013,0.06224563938179744,0.02187309751387206,0.08287352995857286,0.04168576578664523,0.2166341663648842,0.1215947488473217,0.01080175330794645,0.05585812312603372,0.03143424802739055,0.008310221176084891,0.02662867930407082,0.03820054348926789,0.04145355216961652,0.01215117111935851,0.0005961504527883462,0.005696570827932708,0.003138803100433526,0.0158050463135509,0.0266414610403144,0.0253479163102251,0.06647344973838903,0.00742928918051764,0.02430384036556389,0.4024320904279217,0.006677280828936869,0.06353709265672144,0.02199031033771513,0.01591161138653956,0.1309566940755195,0.05520471426119775,0.01326200553286233,0.09217035885210739,0.001655388445236377,0.01951333435466914,0.008992960640954903,0.09304779363691024,0.0501262544076003,0.01027175826267805,0.0139860347737531,0.004948787567092556,0.0009039148422459286,0.03730053477065857,0.02652948036330499,0.008437354041892012,0.02422274467057728,0.01247835527071598,0.005136111506258827,0.01723244706316212,0.008831030390840523,0.009801254168767906,0.007196593446397634,0.03205689489804928,0.005013601421292521,0.01939490464562467,0.006788442705730119,0.02070081050026107,0.03021043909490701,0.0191468458401782,0.0228011693328052,0.01591865013721892,0.0546710288188606,0.01198542047274523,0.01757290702239978,0.03469554361927065,0.0002679938795243351,0.001924129036076227,0.01976013374187118,0.08041760158313886,0.0173125712627492,0.02816119287033069,0.01619484111741859,0.04107516874381439,0.02730416435330348,0.02729925677608709,0.002901191951994045,0.00290119195199406,0.002901191951994095,0.05517565856385273,0.03882776274041642,0.007805721116810325,0.02174418367978063,0.01226642855608374,0.0122664285560838,0.02103764997254013,0.074504525853132,0.04382033923612891,0.07028153516832192,0.1677822742424822,0.1697238822938303,0.0158409784292099,0.01018014135857031,0.03763184498199739,0.006926094958715713,0.009120621111015838,0.002338963743917687,0.00473350084853413,0.002338963743917634,0.005997925532730193,0.01173203138217932,0.001011312435533419,0.004135121020139665,0.00413512102013972,0.002307017344689014,0.01025432257416895,0.09966916423389088,0.07477003997591725,0.1025849143992467,0.06394498644819954,0.07511001964537757,0.01646854207828265,0.08602954182693529,0.01154122808829599,0.007698211968258047,0.01754347575696458,0.009726239000843276,0.01324268959942762,0.03147289277790986,0.04722732517069587,0.02088900922207204,0.007492305029079311,0.03417554788924808,0.02881865574198059,0.02593833235591886,0.004857218993674003,0.01089418117398956,0.031144827225683,0.007003455665084785,0.01656960592664784,0.01106659370779045,0.0110665937077903,0.002466537016889093,0.007625613076279317,0.005170609792253488,0.008232198772238903,0.006846314469939133,0.00600289354739018,0.009585596979075435,0.009422259689534144,0.008512661467764669,0.008185876459440293,0.01231652540282432,0.0002802702016931409,0.005293817748481487,0.003466432570680042],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Eukaryota","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Actinobacteria","Firmicutes","Actinobacteria","Firmicutes","Candidatus_Saccharibacteria","Candidatus_Saccharibacteria","Bacteroidota","Bacteroidota","Bacteroidota","Actinobacteria","Bacteroidota","Firmicutes","Bacteroidota","Actinobacteria","Firmicutes","Proteobacteria","Firmicutes","Proteobacteria","Actinobacteria","Proteobacteria","Firmicutes","Actinobacteria","Proteobacteria","Firmicutes","Fusobacteria","Actinobacteria","Firmicutes","Proteobacteria","Firmicutes","Proteobacteria","Candidatus_Saccharibacteria","Actinobacteria","Fusobacteria","Candidatus_Saccharibacteria","Proteobacteria","Bacteroidota","Proteobacteria","Firmicutes","Firmicutes","Proteobacteria","Firmicutes","Firmicutes","Actinobacteria","Firmicutes","Firmicutes","Bacteroidota","Proteobacteria","Firmicutes","Firmicutes","Spirochaetes","Firmicutes","Candidatus_Saccharibacteria","Bacteroidota","Firmicutes","Proteobacteria","Proteobacteria","Firmicutes","Proteobacteria","Bacteroidota","Candidatus_Saccharibacteria","Bacteroidota","Firmicutes","Candidatus_Absconditabacteria","Bacteroidota","Firmicutes","Candidatus_Saccharibacteria","Firmicutes","Proteobacteria","Firmicutes","Firmicutes","Bacteroidota","Fusobacteria","Bacteroidota","Spirochaetes","Synergistetes","Firmicutes","Bacteroidota","Synergistetes","Candidatus_Saccharibacteria","Candidatus_Saccharibacteria","Bacteroidota","Spirochaetes","Firmicutes","Actinobacteria","Firmicutes","Spirochaetes","Firmicutes","Candidatus_Absconditabacteria","Firmicutes","Spirochaetes","Firmicutes","Candidatus_Saccharibacteria","Candidatus_Gracilibacteria","Candidatus_Saccharibacteria","Firmicutes","Bacteroidota","Candidatus_Absconditabacteria","Firmicutes","Actinobacteria","Tenericutes","Firmicutes","Bacteroidota","Firmicutes","Candidatus_Saccharibacteria","Firmicutes","Bacteroidota","Spirochaetes","Bacteroidota","Bacteroidota","Firmicutes","Fusobacteria","Proteobacteria","Spirochaetes","Actinobacteria","Firmicutes","Spirochaetes","Firmicutes","Bacteroidota","Proteobacteria","Proteobacteria","Firmicutes","Firmicutes","Firmicutes","Actinobacteria","Firmicutes","Actinobacteria","Candidatus_Saccharibacteria","Firmicutes","Candidatus_Saccharibacteria","Bacteroidota","Bacteroidota","Firmicutes","Actinobacteria","Actinobacteria","Firmicutes","Ascomycota","Proteobacteria","Candidatus_Saccharibacteria","Bacteroidota","Firmicutes","Spirochaetes","Firmicutes","Firmicutes","Firmicutes","Bacteroidota","Bacteroidota","Bacteroidota","Firmicutes","Spirochaetes","Firmicutes","Firmicutes","Bacteroidota","Firmicutes","Proteobacteria","Firmicutes","Bacteroidota","Fusobacteria","Bacteroidota","Candidatus_Saccharibacteria","Proteobacteria","Proteobacteria","Proteobacteria","Firmicutes","Firmicutes","Proteobacteria","Proteobacteria","Proteobacteria","Bacteroidota","Bacteroidota","Firmicutes","Firmicutes","Candidatus_Saccharibacteria","Actinobacteria","Bacteroidota","Actinobacteria","Proteobacteria","Actinobacteria","Proteobacteria","Firmicutes","Firmicutes","Fusobacteria","Tenericutes","Firmicutes","Spirochaetes","Firmicutes","Proteobacteria","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Bacteroidota","Bacteroidota","Proteobacteria","Fusobacteria","Proteobacteria","Proteobacteria","Firmicutes","Actinobacteria","Fusobacteria","Firmicutes","Bacteroidota","Firmicutes","Firmicutes","Spirochaetes","Proteobacteria","Proteobacteria","Firmicutes","Firmicutes","Proteobacteria"],["Actinomycetia","Negativicutes","Actinomycetia","Negativicutes","Candidatus_Saccharimonadia","Candidatus_Saccharibacteria_unclassified","Bacteroidia","Flavobacteriia","Bacteroidia","Actinomycetia","Bacteroidia","Clostridia","Bacteroidia","Actinomycetia","Clostridia","Gammaproteobacteria","Bacilli","Gammaproteobacteria","Actinomycetia","Betaproteobacteria","CFGB82399","Actinomycetia","Epsilonproteobacteria","Bacilli","Fusobacteriia","Coriobacteriia","Bacilli","Betaproteobacteria","Negativicutes","Betaproteobacteria","Candidatus_Nanoperiomorbia","Actinomycetia","Fusobacteriia","Candidatus_Saccharibacteria_unclassified","Gammaproteobacteria","Bacteroidia","Betaproteobacteria","Bacilli","Negativicutes","Gammaproteobacteria","Clostridia","Clostridia","Coriobacteriia","Clostridia","Negativicutes","CFGB570","Betaproteobacteria","Negativicutes","Clostridia","Spirochaetia","Tissierellia","Candidatus_Nanosyncoccalia","Bacteroidia","Clostridia","Betaproteobacteria","Betaproteobacteria","Clostridia","Betaproteobacteria","CFGB494","Candidatus_Saccharimonadia","Flavobacteriia","Negativicutes","Candidatus_Absconditabacteria_unclassified","Bacteroidota_unclassified","Negativicutes","Candidatus_Saccharibacteria_unclassified","Clostridia","Epsilonproteobacteria","CFGB1261","Clostridia","Bacteroidia","Fusobacteriia","Flavobacteriia","Spirochaetia","Synergistia","CFGB67123","CFGB496","Synergistia","CFGB12783","CFGB4355","Bacteroidia","Spirochaetia","Clostridia","Actinomycetia","Clostridia","CFGB1021","CFGB67123","Candidatus_Absconditabacteria_unclassified","CFGB1764","CFGB1017","Erysipelotrichia","Candidatus_Saccharibacteria_unclassified","Candidatus_Gracilibacteria_unclassified","CFGB4355","Clostridia","Flavobacteriia","Candidatus_Absconditabacteria_unclassified","Clostridia","Actinomycetia","Tenericutes_unclassified","CFGB1207","CFGB77157","Clostridia","Candidatus_Saccharibacteria_unclassified","CFGB1781","Bacteroidia","Spirochaetia","Bacteroidia","Bacteroidia","Negativicutes","Fusobacteriia","Gammaproteobacteria","Spirochaetia","Coriobacteriia","Clostridia","CFGB1017","Negativicutes","CFGB72483","CFGB73473","Deltaproteobacteria","Erysipelotrichia","Clostridia","Negativicutes","Coriobacteriia","Clostridia","Actinomycetia","Candidatus_Saccharibacteria_unclassified","CFGB1762","CFGB4355","Bacteroidia","CFGB81629","Bacilli","Actinomycetia","Actinomycetia","CFGB1872","Saccharomycetes","Betaproteobacteria","Candidatus_Nanosyncoccalia","Bacteroidia","Clostridia","Spirochaetia","Bacilli","Bacilli","Bacilli","Bacteroidia","Flavobacteriia","CFGB651","CFGB1230","CFGB1021","CFGB14224","CFGB47647","Flavobacteriia","Clostridia","CFGB4196","CFGB1534","CFGB570","Fusobacteriia","Flavobacteriia","Candidatus_Nanosyncoccalia","Betaproteobacteria","Gammaproteobacteria","Alphaproteobacteria","Bacilli","Bacilli","Gammaproteobacteria","Deltaproteobacteria","Deltaproteobacteria","Bacteroidia","CFGB47554","Bacilli","Bacilli","CFGB9314","Actinomycetia","CFGB497","Actinomycetia","Betaproteobacteria","Actinomycetia","Alphaproteobacteria","CFGB79381","CFGB4967","Fusobacteriia","Tenericutes_unclassified","CFGB1762","CFGB1023","Clostridia","Betaproteobacteria","Bacilli","Clostridia","Clostridia","CFGB1781","Bacilli","Bacilli","Flavobacteriia","CFGB72586","Gammaproteobacteria","Fusobacteriia","Gammaproteobacteria","Gammaproteobacteria","Clostridia","Actinomycetia","Fusobacteriia","Bacilli","CFGB494","Bacilli","Bacilli","Spirochaetia","Epsilonproteobacteria","Alphaproteobacteria","Bacilli","Clostridia","Gammaproteobacteria"],["Corynebacteriales","Veillonellales","Actinomycetales","Selenomonadales","Candidatus_Saccharimonadales","Candidatus_Saccharibacteria_unclassified","Bacteroidales","Flavobacteriales","Bacteroidales","Actinomycetales","Bacteroidales","Eubacteriales","Bacteroidales","Actinomycetales","Eubacteriales","Pasteurellales","Lactobacillales","Cardiobacteriales","Actinomycetales","Neisseriales","OFGB82399","Micrococcales","Campylobacterales","Lactobacillales","Fusobacteriales","Coriobacteriales","Bacillales","Burkholderiales","Veillonellales","Neisseriales","Candidatus_Nanoperiomorbales","Propionibacteriales","Fusobacteriales","Candidatus_Saccharibacteria_unclassified","Pasteurellales","Bacteroidales","Neisseriales","Lactobacillales","Selenomonadales","Pasteurellales","Eubacteriales","Eubacteriales","Coriobacteriales","Eubacteriales","Selenomonadales","OFGB570","Burkholderiales","Veillonellales","Eubacteriales","Spirochaetales","Tissierellales","Candidatus_Nanogingivales","Bacteroidales","Eubacteriales","Burkholderiales","Neisseriales","Eubacteriales","Burkholderiales","OFGB494","Candidatus_Nanosynbacterales","Flavobacteriales","Veillonellales","Candidatus_Absconditabacteria_unclassified","Bacteroidota_unclassified","Selenomonadales","Candidatus_Saccharibacteria_unclassified","Eubacteriales","Campylobacterales","OFGB1261","Eubacteriales","Bacteroidales","Fusobacteriales","Flavobacteriales","Spirochaetales","Synergistales","OFGB67123","OFGB496","Synergistales","OFGB12783","OFGB4355","Bacteroidales","Spirochaetales","Eubacteriales","Bifidobacteriales","Eubacteriales","OFGB1021","OFGB67123","Candidatus_Absconditabacteria_unclassified","OFGB1764","OFGB1017","Erysipelotrichales","Candidatus_Saccharibacteria_unclassified","Candidatus_Gracilibacteria_unclassified","OFGB4355","Eubacteriales","Flavobacteriales","Candidatus_Absconditabacteria_unclassified","Eubacteriales","Actinomycetales","Mycoplasmoidales","OFGB1207","OFGB77157","Eubacteriales","Candidatus_Saccharibacteria_unclassified","OFGB1781","Bacteroidales","Spirochaetales","Bacteroidales","Bacteroidales","Selenomonadales","Fusobacteriales","Cardiobacteriales","Spirochaetales","Eggerthellales","Eubacteriales","OFGB1017","Veillonellales","OFGB72483","OFGB73473","Desulfobacterales","Erysipelotrichales","Eubacteriales","Selenomonadales","Eggerthellales","Eubacteriales","Bifidobacteriales","Candidatus_Saccharibacteria_unclassified","OFGB1762","OFGB4355","Bacteroidales","OFGB81629","Lactobacillales","Propionibacteriales","Propionibacteriales","OFGB1872","Saccharomycetales","Burkholderiales","Candidatus_Nanogingivales","Bacteroidales","Eubacteriales","Spirochaetales","Bacillales","Bacillales","Bacillales","Bacteroidales","Flavobacteriales","OFGB651","OFGB1230","OFGB1021","OFGB14224","OFGB47647","Flavobacteriales","Eubacteriales","OFGB4196","OFGB1534","OFGB570","Fusobacteriales","Flavobacteriales","Candidatus_Nanogingivales","Neisseriales","Moraxellales","Rhodobacterales","Lactobacillales","Lactobacillales","Enterobacterales","Desulfovibrionales","Desulfovibrionales","Bacteroidales","OFGB47554","Bacillales","Lactobacillales","OFGB9314","Actinomycetales","OFGB497","Bifidobacteriales","Neisseriales","Bifidobacteriales","Sphingomonadales","OFGB79381","OFGB4967","Fusobacteriales","Mycoplasmoidales","OFGB1762","OFGB1023","Eubacteriales","Neisseriales","Bacillales","Eubacteriales","Eubacteriales","OFGB1781","Lactobacillales","Lactobacillales","Flavobacteriales","OFGB72586","Pseudomonadales","Fusobacteriales","Moraxellales","Enterobacterales","Eubacteriales","Actinomycetales","Fusobacteriales","Lactobacillales","OFGB494","Lactobacillales","Lactobacillales","Spirochaetales","Campylobacterales","Hyphomicrobiales","Lactobacillales","Eubacteriales","Moraxellales"],["Corynebacteriaceae","Veillonellaceae","Actinomycetaceae","Selenomonadaceae","Candidatus_Saccharimonadaceae","Candidatus_Saccharibacteria_unclassified","Porphyromonadaceae","Flavobacteriaceae","Porphyromonadaceae","Actinomycetaceae","Prevotellaceae","Lachnospiraceae","Tannerellaceae","Actinomycetaceae","Lachnospiraceae","Pasteurellaceae","Streptococcaceae","Cardiobacteriaceae","Actinomycetaceae","Neisseriaceae","FGB82399","Micrococcaceae","Campylobacteraceae","Carnobacteriaceae","Leptotrichiaceae","Atopobiaceae","Bacillales_unclassified","Comamonadaceae","Veillonellaceae","Neisseriaceae","Candidatus_Nanoperiomorbaceae","Propionibacteriaceae","Fusobacteriaceae","Candidatus_Saccharibacteria_unclassified","Pasteurellaceae","Prevotellaceae","Neisseriaceae","Aerococcaceae","Selenomonadaceae","Pasteurellaceae","Eubacteriales_Family_XIII_Incertae_Sedis","Eubacteriales_Family_XIII_Incertae_Sedis","Atopobiaceae","Lachnospiraceae","Selenomonadaceae","FGB570","Burkholderiaceae","Veillonellaceae","Peptostreptococcaceae","Treponemataceae","Peptoniphilaceae","Candidatus_Nanogingivalaceae","Prevotellaceae","Lachnospiraceae","Comamonadaceae","Neisseriaceae","Lachnospiraceae","Burkholderiaceae","FGB494","Candidatus_Nanosynbacteraceae","Flavobacteriaceae","Veillonellaceae","Candidatus_Absconditabacteria_unclassified","Bacteroidota_unclassified","Selenomonadaceae","Candidatus_Saccharibacteria_unclassified","Lachnospiraceae","Campylobacteraceae","FGB1261","Lachnospiraceae","Prevotellaceae","Leptotrichiaceae","Flavobacteriaceae","Treponemataceae","Synergistaceae","FGB67123","FGB496","Synergistaceae","FGB12783","FGB4355","Prevotellaceae","Treponemataceae","Lachnospiraceae","Bifidobacteriaceae","Lachnospiraceae","FGB1021","FGB67123","Candidatus_Absconditabacteria_unclassified","FGB1764","FGB1017","Erysipelotrichaceae","Candidatus_Saccharibacteria_unclassified","Candidatus_Gracilibacteria_unclassified","FGB4355","Peptostreptococcaceae","Weeksellaceae","Candidatus_Absconditabacteria_unclassified","Lachnospiraceae","Actinomycetaceae","Metamycoplasmataceae","FGB1207","FGB77157","Peptostreptococcaceae","Candidatus_Saccharibacteria_unclassified","FGB1781","Muribaculaceae","Treponemataceae","Porphyromonadaceae","Porphyromonadaceae","Selenomonadaceae","Fusobacteriaceae","Cardiobacteriaceae","Treponemataceae","Eggerthellaceae","Lachnospiraceae","FGB1017","Veillonellaceae","FGB72483","FGB73473","Desulfobulbaceae","Erysipelotrichaceae","Lachnospiraceae","Selenomonadaceae","Eggerthellaceae","Lachnospiraceae","Bifidobacteriaceae","Candidatus_Saccharibacteria_unclassified","FGB1762","FGB4355","Porphyromonadaceae","FGB81629","Streptococcaceae","Propionibacteriaceae","Propionibacteriaceae","FGB1872","Debaryomycetaceae","Comamonadaceae","Candidatus_Nanogingivalaceae","Prevotellaceae","Peptostreptococcaceae","Treponemataceae","Bacillaceae","Planococcaceae","Paenibacillaceae","Prevotellaceae","Flavobacteriaceae","FGB651","FGB1230","FGB1021","FGB14224","FGB47647","Weeksellaceae","Lachnospiraceae","FGB4196","FGB1534","FGB570","Leptotrichiaceae","Flavobacteriaceae","Candidatus_Nanogingivalaceae","Neisseriaceae","Moraxellaceae","Rhodobacteraceae","Lactobacillaceae","Lactobacillaceae","Enterobacteriaceae","Desulfovibrionaceae","Desulfovibrionaceae","Bacteroidaceae","FGB47554","Bacillaceae","Carnobacteriaceae","FGB9314","Actinomycetaceae","FGB497","Bifidobacteriaceae","Neisseriaceae","Bifidobacteriaceae","Sphingomonadaceae","FGB79381","FGB4967","Leptotrichiaceae","Metamycoplasmataceae","FGB1762","FGB1023","Eubacteriales_unclassified","Neisseriaceae","Staphylococcaceae","Peptostreptococcaceae","Peptostreptococcaceae","FGB1781","Lactobacillaceae","Lactobacillaceae","Weeksellaceae","FGB72586","Pseudomonadaceae","Leptotrichiaceae","Moraxellaceae","Enterobacteriaceae","Lachnospiraceae","Actinomycetaceae","Leptotrichiaceae","Lactobacillaceae","FGB494","Lactobacillaceae","Lactobacillaceae","Treponemataceae","Campylobacteraceae","Methylobacteriaceae","Lactobacillaceae","Lachnospiraceae","Moraxellaceae"],["Corynebacterium","Veillonella","Actinomyces","Selenomonas","Candidatus_Nanosynsacchari","Candidatus_Saccharibacteria_unclassified","Porphyromonas","Capnocytophaga","GGB1472","Pauljensenia","Prevotella","GGB3883","Tannerella","Schaalia","Lachnoanaerobaculum","Aggregatibacter","Streptococcus","Cardiobacterium","Actinobaculum","Kingella","GGB4308","Rothia","Campylobacter","Granulicatella","Leptotrichia","Olsenella","Gemella","Ottowia","Megasphaera","Eikenella","Candidatus_Nanoperiomorbus","Arachnia","Fusobacterium","GGB12798","GGB6813","GGB1144","Neisseria","Abiotrophia","GGB4303","Haemophilus","Eubacteriales_Family_XIII_Incertae_Sedis_unclassified","Mogibacterium","Lancefieldella","Oribacterium","Centipeda","GGB1202","Lautropia","Dialister","Peptostreptococcus","Treponema","Parvimonas","GGB12785","Alloprevotella","Catonella","Delftia","GGB6679","GGB3887","Ralstonia","GGB1022","Candidatus_Nanosynbacter","GGB1838","Anaeroglobus","Candidatus_Absconditabacteria_unclassified","Bacteroidota_unclassified","GGB4299","GGB12794","GGB3385","GGB12443","GGB3059","Johnsonella","GGB74353","GGB4400","GGB1833","GGB2676","Fretibacterium","GGB71270","GGB1025","GGB10852","GGB12783","GGB12788","Hallella","GGB2663","Shuttleworthia","Bifidobacterium","GGB3388","GGB2671","GGB96534","GGB4937","GGB4538","GGB2666","Solobacterium","GGB49434","Candidatus_Gracilibacteria_unclassified","GGB12789","Filifactor","GGB1843","GGB4936","Lachnospiraceae_unclassified","Peptidiphaga","Metamycoplasma","GGB2964","GGB49219","Peptostreptococcaceae_unclassified","GGB18703","GGB4733","GGB1460","GGB2675","GGB45697","GGB1088","GGB97311","GGB56136","GGB72530","GGB71076","Slackia","GGB3386","GGB73508","Colibacter","GGB71456","GGB6673","Desulfobulbus","Bulleidia","GGB49397","GGB4300","Cryptobacterium","GGB3387","Scardovia","GGB12761","GGB4533","GGB12787","GGB1473","GGB1024","Lactococcus","Cutibacterium","Propionibacterium","GGB4964","Candida","GGB97303","GGB12786","GGB1186","Peptoanaerobacter","GGB2674","Caldibacillus","Ureibacillus","GGB16822","Prevotellaceae_unclassified","GGB1832","GGB1611","GGB96297","GGB49499","GGB38873","GGB49401","GGB1844","GGB3389","GGB12441","GGB3886","GGB1203","Pseudoleptotrichia","GGB1840","Candidatus_Nanogingivalis","Simonsiella","GGB72444","Paracoccus","Lacticaseibacillus","Lactiplantibacillus","Escherichia","Bilophila","Desulfovibrio","Phocaeicola","GGB49229","Geobacillus","Dolosigranulum","GGB12763","GGB9835","GGB1026","Alloscardovia","GGB6688","Parascardovia","Sphingomonas","GGB3008","GGB49400","GGB49528","GGB4333","GGB71303","GGB2672","GGB70946","GGB72533","Staphylococcus","GGB4786","GGB4783","GGB4721","Limosilactobacillus","Lactobacillus","GGB1845","GGB1188","Pseudomonas","Sneathia","Acinetobacter","Klebsiella","GGB3390","Arcanobacterium","GGB4393","Latilactobacillus","GGB1021","Leuconostoc","Weissella","GGB49504","GGB12440","Methylobacterium","Ligilactobacillus","GGB96582","Moraxella"],[null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null],[0.003403225806451613,0.003403225806451613,0.003403225806451613,0.003403225806451613,0.003403225806451613,0.2259645390070922,0.159496062992126,0.003403225806451613,0.003403225806451613,0.003403225806451613,0.003403225806451613,0.003403225806451613,0.003403225806451613,0.003403225806451613,0.003403225806451613,0.003403225806451613,0.003403225806451613,0.003403225806451613,0.003403225806451613,0.003403225806451613,0.008440000000000001,0.003403225806451613,0.15192,0.003403225806451613,0.003403225806451613,0.005861111111111111,0.003403225806451613,0.003403225806451613,0.003403225806451613,0.003403225806451613,0.02344444444444444,0.003403225806451613,0.003403225806451613,0.01302469135802469,0.01302469135802469,0.003403225806451613,0.003403225806451613,0.003403225806451613,0.003403225806451613,0.005861111111111111,0.003403225806451613,0.003403225806451613,0.003403225806451613,0.003403225806451613,0.005861111111111111,0.003403225806451613,0.003403225806451613,0.02495698924731183,0.003403225806451613,0.003403225806451613,0.1664921875,0.3116092715231788,0.2034642857142857,0.003403225806451613,0.01302469135802469,0.005861111111111111,0.003403225806451613,0.02495698924731183,0.01507142857142857,0.7007421052631578,0.005861111111111111,0.1307851239669421,0.04008999999999999,0.005861111111111111,0.01678409090909091,0.005861111111111111,0.003403225806451613,0.003403225806451613,0.09724347826086957,0.003403225806451613,0.01507142857142857,0.003403225806451613,0.003403225806451613,0.3165,0.005861111111111111,0.03077083333333333,0.4548502994011976,0.04344117647058824,0.02693617021276596,0.01678409090909091,0.2855533333333334,0.9439473684210528,0.6135054347826088,0.9397219512195123,0.1830735294117647,0.05224761904761904,0.05861111111111111,0.008440000000000001,0.5256914285714286,0.08100892857142857,0.003403225806451613,0.6135054347826088,0.003403225806451613,0.1136153846153846,0.1988561151079137,0.003403225806451613,0.01110526315789474,0.2813333333333333,0.003403225806451613,0.9397219512195123,0.149741935483871,0.4559107142857143,0.003403225806451613,0.003403225806451613,0.4049627329192546,0.1926521739130435,0.7048062827225132,0.9397219512195123,0.01896629213483146,0.05718691588785047,0.4815764705882353,0.06272972972972973,0.3294322580645161,0.6612085561497326,0.1301166666666667,0.4438719512195122,0.422,0.5256914285714286,0.04008999999999999,0.6341344086021506,0.1678409090909091,0.6131823204419889,0.09439473684210525,0.04344117647058824,0.1215932203389831,0.09822413793103447,0.1678409090909091,0.01302469135802469,0.3171895424836602,0.1668372093023256,0.02495698924731183,0.9726095238095238,0.8752193877551019,0.1366311475409836,0.003403225806451613,0.159496062992126,0.0348041237113402,0.1848175182481752,0.008440000000000001,0.0600091743119266,0.05072115384615385,0.9397219512195123,0.9397219512195123,0.9397219512195123,0.005861111111111111,0.01678409090909091,0.5152325581395348,0.09149557522123894,0.2508391608391609,0.2508391608391609,0.1423821138211382,0.003403225806451613,0.01302469135802469,0.003403225806451613,0.003403225806451613,0.003403225806451613,0.171634328358209,0.4049627329192546,0.01678409090909091,0.5666179775280898,0.4475757575757576,0.9434134615384615,0.6936063829787235,0.9434134615384615,0.6135054347826088,0.2633862068965517,0.9397219512195123,0.8353435897435897,0.8353435897435897,0.9434134615384615,0.3749617834394904,0.003403225806451613,0.003403225806451613,0.003403225806451613,0.005861111111111111,0.003403225806451613,0.1813037037037037,0.003403225806451613,0.3294322580645161,0.5286988636363636,0.171634328358209,0.4512349397590361,0.264472602739726,0.05072115384615385,0.01507142857142857,0.1223445378151261,0.5364406779661017,0.03077083333333333,0.05374528301886792,0.06272972972972973,0.7007421052631578,0.3749617834394904,0.04008999999999999,0.5775977653631285,0.1678409090909091,0.2855533333333334,0.2855533333333334,0.9397219512195123,0.5256914285714286,0.7368601036269431,0.4836959064327486,0.6060388888888889,0.6341344086021506,0.3792658227848101,0.4049627329192546,0.4794319526627219,0.4438719512195122,0.2564236111111111,1,0.7143229166666667,0.9029086294416243]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>Family<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0},{"name":" ","targets":0},{"name":"Axis.1","targets":1},{"name":"Axis.2","targets":2},{"name":"pval","targets":3},{"name":"r","targets":4},{"name":"Kingdom","targets":5},{"name":"Phylum","targets":6},{"name":"Class","targets":7},{"name":"Order","targets":8},{"name":"Family","targets":9},{"name":"tax_rank_plot","targets":10},{"name":"Species","targets":11},{"name":"Strain","targets":12},{"name":"pval.adj","targets":13}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

Plot:


``` r
out$env_fit_tax$plot 
```

<img src="emvit_demo_files/figure-html/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

To be able to fit the vectors we have to compute the ordination again. And sometimes,  ordiantions of the same data can be mirrored (sample are still similarly distant from each other but it is a problem for vector fitting !) 
see: https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/visualizing-and-interpreting-ordinations/
It shouldn't be the case because we are setting the same seed before but you never know...


``` r
out$env_fit_tax$ploted_ordination_used_for_envfit 
```

<img src="emvit_demo_files/figure-html/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />
Same plot!  

Same with metadata:


``` r
out$env_fit_meta$plot
```

<img src="emvit_demo_files/figure-html/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

# Run fonction on TP2 Plaque:


``` r
ps_up %>%
  subset_taxa(Class != "UNCLASSIFIED") %>%
  subset_samples(Sample == "Plaque" &
                   Time == "TP2") %>%
  transform_sample_counts(function(x) x/sum(x) * 1) %>% 
  phyloseq_explore_beta(ps_up = .,
                        beta = beta,
                        color_group = "Time",
                        shape_group = "Sample",
                        distance_for_more = "rAitchison", #distance for detailed plot (path, and envfit vectors)
                        metadata_sel = c("delta_mean_plaque", "delta_mean_bleeding"), # metadata selection
                        col_pal = time_pal,
                        fill_pal = time_pal,
                        path_group = "interaction(Sample,Subject)",
                        facet_formula = "Sample ~ .",
                        permanova_terms = c("cluster_Dtp2"),
                        strata = "none") -> plaque_tp2
```

Error in phyloseq_add_metadata_vector_fix(dist = beta[[distance_for_more]],  : 
No significant metadata vectors after adjustment.

No assocaition between   c("delta_mean_plaque", "delta_mean_bleeding") and microbiota at TP2.

We will run it without:

We can specify the taxonomic level for testing association between ordination and proportion of taxa: i.e., Phylum


``` r
ps_up %>%
  subset_taxa(Class != "UNCLASSIFIED") %>%
  subset_samples(Sample == "Plaque" &
                   Time == "TP2") %>%
  transform_sample_counts(function(x) x/sum(x) * 1) %>% 
  phyloseq_explore_beta(ps_up = .,
                        beta = beta,
                        color_group = "Time",
                        shape_group = "Sample",
                        distance_for_more = "rAitchison", #distance for detailed plot (path, and envfit vectors)
                        tax_rank_fit  = "Phylum", # taxonomic level for envfit vector fitting
                        metadata_sel = NULL,
                        col_pal = time_pal,
                        fill_pal = time_pal,
                        path_group = "interaction(Sample,Subject)",
                        facet_formula = "Sample ~ .",
                        permanova_terms = c("cluster_Dtp2"),
                        strata = "none") -> plaque_tp2


plaque_tp2$env_fit_tax$plot
```

<img src="emvit_demo_files/figure-html/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />


Example with `Species`:



``` r
ps_up %>%
  subset_taxa(Class != "UNCLASSIFIED") %>%
  subset_samples(Sample == "Plaque" &
                   Time == "TP2") %>%
  transform_sample_counts(function(x) x/sum(x) * 1) %>% 
  phyloseq_explore_beta(ps_up = .,
                        beta = beta,
                        color_group = "Time",
                        shape_group = "Sample",
                        distance_for_more = "rAitchison", #distance for detailed plot (path, and envfit vectors)
                        tax_rank_fit  = "Species", # taxonomic level for envfit vector fitting
                        metadata_sel = NULL,
                        col_pal = time_pal,
                        fill_pal = time_pal,
                        path_group = "interaction(Sample,Subject)",
                        facet_formula = "Sample ~ .",
                        permanova_terms = c("cluster_Dtp2"),
                        strata = "none") -> plaque_tp2




plaque_tp2$PCOA$p
```

<img src="emvit_demo_files/figure-html/unnamed-chunk-13-1.png" style="display: block; margin: auto;" />



``` r
plaque_tp2$env_fit_tax$plot
```

<img src="emvit_demo_files/figure-html/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />



``` r
sessionInfo()
```

```
## R version 4.4.0 (2024-04-24)
## Platform: aarch64-apple-darwin20
## Running under: macOS 15.2
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## time zone: Europe/Paris
## tzcode source: internal
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] microbiome_1.28.0    gginnards_0.2.0-1    vegan_2.6-8         
##  [4] lattice_0.22-6       permute_0.9-7        GUniFrac_1.8        
##  [7] ape_5.8-1            reshape2_1.4.4       scales_1.3.0        
## [10] microViz_0.12.4      speedyseq_0.5.3.9021 phyloseq_1.50.0     
## [13] readxl_1.4.3         lubridate_1.9.4      forcats_1.0.0       
## [16] stringr_1.5.1        dplyr_1.1.4          purrr_1.0.2         
## [19] readr_2.1.5          tidyr_1.3.1          tibble_3.2.1        
## [22] ggplot2_3.5.1        tidyverse_2.0.0     
## 
## loaded via a namespace (and not attached):
##   [1] inline_0.3.20           rlang_1.1.4             magrittr_2.0.3         
##   [4] clue_0.3-66             ade4_1.7-22             matrixStats_1.4.1      
##   [7] compiler_4.4.0          mgcv_1.9-1              vctrs_0.6.5            
##  [10] rmutil_1.1.10           pkgconfig_2.0.3         crayon_1.5.3           
##  [13] fastmap_1.2.0           backports_1.5.0         XVector_0.46.0         
##  [16] labeling_0.4.3          modeest_2.4.0           rmarkdown_2.29         
##  [19] tzdb_0.4.0              UCSC.utils_1.2.0        xfun_0.49              
##  [22] zlibbioc_1.52.0         cachem_1.1.0            GenomeInfoDb_1.42.1    
##  [25] jsonlite_1.8.9          biomformat_1.34.0       rhdf5filters_1.18.0    
##  [28] Rhdf5lib_1.28.0         broom_1.0.7             parallel_4.4.0         
##  [31] cluster_2.1.8           R6_2.5.1                bslib_0.8.0            
##  [34] stringi_1.8.4           RColorBrewer_1.1-3      GGally_2.2.1           
##  [37] rpart_4.1.23            car_3.1-3               jquerylib_0.1.4        
##  [40] cellranger_1.1.0        Rcpp_1.0.13-1           iterators_1.0.14       
##  [43] knitr_1.49              IRanges_2.40.1          Matrix_1.7-1           
##  [46] splines_4.4.0           igraph_2.1.2            timechange_0.3.0       
##  [49] tidyselect_1.2.1        rstudioapi_0.17.1       abind_1.4-8            
##  [52] yaml_2.3.10             timeDate_4041.110       codetools_0.2-20       
##  [55] plyr_1.8.9              Biobase_2.66.0          withr_3.0.2            
##  [58] stable_1.1.6            evaluate_1.0.1          Rtsne_0.17             
##  [61] survival_3.7-0          ggstats_0.7.0           Biostrings_2.74.0      
##  [64] pillar_1.10.0           ggpubr_0.6.0            carData_3.0-5          
##  [67] DT_0.33                 foreach_1.5.2           stats4_4.4.0           
##  [70] generics_0.1.3          rprojroot_2.0.4         S4Vectors_0.44.0       
##  [73] hms_1.1.3               munsell_0.5.1           timeSeries_4041.111    
##  [76] glue_1.8.0              statip_0.2.3            tools_4.4.0            
##  [79] data.table_1.16.4       spatial_7.3-17          fBasics_4041.97        
##  [82] ggsignif_0.6.4          cowplot_1.1.3           rhdf5_2.50.1           
##  [85] grid_4.4.0              crosstalk_1.2.1         colorspace_2.1-1       
##  [88] nlme_3.1-166            GenomeInfoDbData_1.2.13 Formula_1.2-5          
##  [91] cli_3.6.3               gtable_0.3.6            stabledist_0.7-2       
##  [94] rstatix_0.7.2           sass_0.4.9              digest_0.6.37          
##  [97] BiocGenerics_0.52.0     ggrepel_0.9.6           htmlwidgets_1.6.4      
## [100] farver_2.1.2            htmltools_0.5.8.1       multtest_2.62.0        
## [103] lifecycle_1.0.4         httr_1.4.7              here_1.0.1             
## [106] statmod_1.5.0           MASS_7.3-61
```
