spe_pseudo$broad<-factor(ifelse(spe_pseudo$cluster %in% c(1,4,8,12,13,14,17,20,21,24,3),'Neuron',
                  ifelse(spe_pseudo$cluster %in% c(5,10,18),'Neuropil',
                         ifelse(spe_pseudo$cluster %in% c(15,16),'unsure',
                         ifelse(spe_pseudo$cluster %in% c(2,6,9,19,22),'WM',
                          'CSF/Vascular')))))
spe$broad.class<-factor(ifelse(spe$cluster %in% c(1,4,8,12,14,17,19),'Neuron',
                                      ifelse(spe$cluster %in% c(3,13,15,16),'Neuron_low',
                                             ifelse(spe$cluster %in% c(5,10),'Neuropil',
                                                    ifelse(spe$cluster %in% c(2,6,9,18,20),'WM',
                                                           ifelse(spe$cluster %in% c(7),'Vascular',
                                                                  'Choroid'))))))
spe$broad.class2<-factor(ifelse(spe$cluster %in% 12,'GCL',
                                ifelse(spe$cluster %in% 1,'GABA',
                                ifelse(spe$cluster %in% c(4,8,14,17,19),'Neuron',
                                ifelse(spe$cluster %in% c(3,13,15,16),'Neuron_low',
                                       ifelse(spe$cluster %in% c(5,10),'Neuropil',
                                              ifelse(spe$cluster %in% c(2,6,9,18,20),'WM',
                                                     ifelse(spe$cluster %in% c(7),'Vascular',
                                                            'Choroid'))))))
spe$cluster<-factor(spe$cluster,levels=c(1,4,8,12,14,17,19,
                                         3,13,15,5,10,16,2,6,9,18,20,
                                         7,11))

spe_pseudo$broad.class<-factor(ifelse(spe_pseudo$cluster %in% c(1,3,4,7,8,9,10,15,18),'Neuron',
                                ifelse(spe_pseudo$cluster %in% c(2,5,13,14,16),'Neuropil',
                                              ifelse(spe_pseudo$cluster %in% c(6,19,17,11),'WM',
                                                     ifelse(spe_pseudo$cluster %in% c(20),'Vascular',
                                                     'Choroid')))))

