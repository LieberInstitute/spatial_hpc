
refine <- function(sample_id, pred, dis, shape = "hexagon"){
  refined.pred = NULL
  pred = as.data.frame(pred, row.names = sample_id)
  dis_df = as.data.frame(dis, row.names = sample_id)
  colnames(dis_df) = sample_id
  if (shape == "hexagon"){
    num_nbs = 6
  } else if (shape == "square"){
    num_obs = 4
  } else {
    cat("Shape not recognized, shape = 'hexagon' for Visium data, 'square' for ST data")
  }
  for (i in 1:length(sample_id)){
    index = sample_id[i]
    dis_tmp = as.data.frame(cbind(dis_df[index],sample_id))
    dis_tmp = dis_tmp[order(dis_tmp[,1]),]
    nbs = dis_tmp[1:num_nbs+1,]
    nbs_pred = pred[nbs$sample_id,]
    self_pred = pred[index,]
    v_c = as.data.frame(table(nbs_pred))
    if (((v_c[v_c$nbs_pred == self_pred, "Freq"] < num_nbs/2) && v_c$nbs_pred == self_pred) && (max(v_c$Freq) > num_nbs/2)){
      refined.pred = c(refined.pred, v_c[v_c$Freq == max(v_c$Freq),"nbs_pred"])
    } else {
      refined.pred = c(refined.pred, self_pred)
    }
  }
  return(refined.pred)
}




