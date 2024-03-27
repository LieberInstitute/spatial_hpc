sce$cell.class<-factor(ifelse(sce$fine.type %in% c('L6.1','L6.2','L6b'),'L6/6b',
    ifelse(sce$fine.type %in% 'Sub.1','Sub.1',
                       ifelse(sce$fine.type %in% 'Sub.2','Sub.2',
                              ifelse(sce$fine.type %in% c('CA1','ProS'),'CA1/ProS',
                              ifelse(sce$cell.type %in% c('Choroid'),'Choroid',
                            ifelse(sce$cell.type %in% c('Ependy'),'Ependy',
                                     ifelse(sce$fine.type %in% c('L2/3.1','L2/3.5'),'L2/3.PrS.PaS',
                                ifelse(sce$fine.type %in% c('L2/3.2','L2/3.4','L2/3.6','L2/3.3'),'L2/3.Prs.Ent',
                                                     ifelse(sce$fine.type %in% c('HATA'),'HATA/Amy',
                                                            ifelse(sce$fine.type %in% c('AHi.1','AHi.2','AHi.3','AHi.4'),'HATA/Amy',
                                                                          ifelse(sce$fine.type %in% c('PENK'),'GABA.Amy',
                                                                                 ifelse(sce$fine.type %in% c('LAMP5.CGE','LAMP5.MGE','CXCL14'),'GABA.LAMP5',
                                                                                        ifelse(sce$fine.type %in% c('HTR3A','VIP'),'GABA.CGE',
                                                                                               ifelse(sce$fine.type %in% c('PV.FS','SST','CORT','CRABP1','C1QL1'),'GABA.MGE',
                                                                                                      as.character(sce$cell.type))))))))))))))))



sce$broad.class<-factor(ifelse(sce$fine.type %in% c('OPC','COP'),'OPC',
                              ifelse(sce$cell.type %in% c('Oligo'),'Oligo',
                                     ifelse(sce$cell.type %in% c('Astro'),'Astro',
                                            ifelse(sce$cell.type %in% c('Micro/Macro/T'),'Micro/Macro/T',as.character(sce$broad.type))))))

sce$cart.group<-factor(
    ifelse(sce$broad.class %in% c('ExcN','InhN'),'Neuron',
           ifelse(sce$broad.class %in% c('CSF','Micro/Macro/T','OPC', 'Vascular'),'Other',
                  as.character(sce$broad.class))))
