library(CytobankAPI)
library(CytobankBridgeR)


#authenticate cytobank premium access

cyto_session <- authenticate(site="premium", username="emilsinclair", password="Laila1tov!123")

#auth_token = cyto_session@auth_token

#draw automatic gates around each cluster in the cluster channel

CytobankBridgeR::gates.apply_cluster_gates(UserSession=cyto_session, experiment_id=172038, clusters=c(1:15), name="CLUSTER_", channel_name="cluster", grouped = FALSE)

# list all fcs files from a given experiment for API based download

exp<-read.csv("exp.csv", header=T)

############################    retrieve list of experiment IDs and associated FCS file IDs

for (i in exp$ID){
  files<-data.frame(fcs_files.list(UserSession=cyto_session, experiment_id=i))
  results<-as.data.frame(lapply(rbind(files, results), unlist))
  }

fcs<-results$id


res<-data.frame(results$id, results$experimentId)
colnames(res)<-c("FCS", "ExpID")
  
################## Download FCS files for each experiment

for (i in unique(res$ExpID)){
  fcs_files.download_fcs_files_stable(cyto_session, experiment_id=i, fcs_files=res$FCS[res$ExpID==i])
}



################ experimental area ##################



for (i in unique(res$ExpID)){
     print(res$FCS[res$ExpID==i])
   }









