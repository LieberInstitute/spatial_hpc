library(scater)
library(scran)
library(SingleCellExperiment)

plotColData(sce,
    x = "k_50_label", y = "detected",
    colour_by = "k_50_label", point_size = 2.5, point_alpha = 0.5
)

###first who are neurons?
sce$neuron<-ifelse(sce$k_50_label %in% c(11,1,5,22,28,14,17,15,12,4),F,T)

table(sce$neuron)
#FALSE  TRUE
#23498 53753

##how do these fall across sorts?
#        HPC-PI+ HPC-PI+NeuN+
#  FALSE   17471         6027
#  TRUE    20343        33410

##what about sample #?


##ok now let's make a "broad" category with excit, inhib, and glial/non-neuron subtypes
broadTab<-data.frame('cluster'=seq(1,38,1),'broadType'=rep(NA,38))
broadTab$broadType[broadTab$cluster %in%c(37,17,31,27,1,22,10,34)]<-paste0('Inhib_',c(1:8))
broadTab$broadType[broadTab$cluster %in%38]<-'CR'
broadTab$broadType[broadTab$cluster %in%c(9,13,21,25,28,33,6,16,32,35,19,24,8,26,36,7,23)]<-paste0('Excit_',c(1:17))
broadTab$broadType[broadTab$cluster %in% c(15,5,30)]<-paste0('CP_',c(1:3))
broadTab$broadType[broadTab$cluster %in% c(2,14,20)]<-paste0('Astro_',c(1:3))
broadTab$broadType[broadTab$cluster %in% c(11,18)]<-paste0('Oligo_',c(1:2))
broadTab$broadType[broadTab$cluster %in% c(3)]<-'Endo'
broadTab$broadType[broadTab$cluster %in% c(29)]<-'PC'
broadTab$broadType[broadTab$cluster %in% c(12)]<-'OPC'
broadTab$broadType[broadTab$cluster %in% c(4)]<-'Micro'

c(paste0('Inhib_',c(1:8)),'CR', paste0('Excit_',c(1:17)),paste0('CP_',c(1:3)),paste0('Astro_',c(1:3)),paste0('Oligo_',c(1:2)),'Endo','PC','OPC','Micro')

levels(broadTab$broadType)[c(37,17,31,27,1,22,10,34,38,9,13,21,25,28,33,6,16,32,35,19,24,8,26,36,7,23,15,5,30,2,14,20,11,18,3,29,12,4)]



table(sce$broadType)
#          Astro              CP           Excit           Inhib           Micro
#           4040            4059           45914            7839            3101
#          Oligo             OPC Vascular/Immune
#           8456            1280            2562


table(sce$broadType,sce$BrNum)

#                  Br2720 Br2743 Br3942 Br6423 Br6432 Br6471 Br6522 Br8325
#  Astro              165    241    400    286   1035    666    256    260
#  CP                  55    144    377     74    191   1808    712     50
#  Excit              517   4549   4305    838   5426  12036   3161   3657
#  Inhib              115    500    789    245   1735   1277    461    410
#  Micro              139    198    221    155    249    531    396    330
#  Oligo              257   1053   1788    660    661    719    471    632
#  OPC                 87    116    123     65    324    172     90     60
#  Vascular/Immune     48    151    140     53    242    701    164    139
#
#                  Br8492 Br8667
#  Astro              296    435
#  CP                 533    115
#  Excit             5942   5483
#  Inhib             1378    929
#  Micro              696    186
#  Oligo              362   1853
#  OPC                 88    155
#  Vascular/Immune    765    159

table(sce$broadType,sce$Sample)
#                  10c-scp 11c-scp 12c-scp 13c-scp 14c-scp 15c-scp 16c-scp
#  Astro               348      87     286     114     192      49     448
#  CP                   69      46     335      42     109      35    1415
#  Excit              2408    3075    1702    2603    1099    3450    3445
#  Inhib               357     572     192     597     217     283     303
#  Micro               171      15     138      83     192       6     327
#  Oligo              1714     139    1626     162    1039      14     463
#  OPC                 107      48      92      31     102      14     127
#  Vascular/Immune     134      25      93      47      58      93     353
#
#                  17c-scp 18c-scp 19c-scp 1c-scp 20c-scp 21c-scp 22c-scp
#  Astro               218     717     318    177     192     104     203
#  CP                  393     117      74     34     291     242      48
#  Excit              8591    2848    2578   1180    2669    3273     330
#  Inhib               974     539    1196    128     506     872      89
#  Micro               204     169      80    223     442     254     108
#  Oligo               256     509     152    521     245     117     470
#  OPC                  45     237      87     49      57      31      45
#  Vascular/Immune     348     142     100     93     400     365      33
#
#                  23c-scp 24c-scp 25c-scp 26c-scp 27c-scp 2c-scp
#  Astro                83     129      36     174      82     83
#  CP                   26      52       3     456     256     16
#  Excit               508     383     134    1691    1470   2477
#  Inhib               156      97      18     160     301    282
#  Micro                47     107      32     245     151    107
#  Oligo               190     207      50     369     102    111
#  OPC                  20      70      17      69      21     11
#  Vascular/Immune      20      39       9      94      70     46


