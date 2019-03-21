library(CytobankAPI)
library(CytobankBridgeR)


#authenticate cytobank premium access

cyto_session <- authenticate(site="premium", username="emilsinclair", password="Laila1tov!123")

#auth_token = cyto_session@auth_token

#draw automatic gates around each cluster in the cluster channel

CytobankBridgeR::gates.apply_cluster_gates(UserSession=cyto_session, experiment_id=173406, clusters=c(1:30), name="CLUSTER_", channel_name="cluster", grouped = FALSE)


#make a new experiment

new_experiment <- experiments.new(cyto_session, experiment_name="API Test3", purpose="Let's try this", output="raw")


#upload fcs files

fcs_files.upload(cyto_session, 166074, file_path="C:/Users/Florian/PhD/flowdata/re_vac063_T_cell_sorts/by_person/818/818c-1.fcs", timeout = 300)
fcs_files.upload(cyto_session, 166074, file_path="C:/Users/Florian/PhD/flowdata/re_vac063_T_cell_sorts/by_person/818/818c+8.fcs", timeout = 300)
fcs_files.upload(cyto_session, 166074, file_path="C:/Users/Florian/PhD/flowdata/re_vac063_T_cell_sorts/by_person/818/818c+9.fcs", timeout = 300)

#download sample tags from clean baseline v infected experiment

sample_tags.download(cyto_session, 160746, directory = getwd())

#upload downloaded file to new experiment

sample_tags.upload(cyto_session, experiment_id=166074, "C:/Users/Florian/PhD/flowdata/re_vac063_T_cell_sorts/by_person/818/experiment_160746_annotations.tsv")

liste<-fcs_files.list(cyto_session, experiment_id=166074, output = "default")
View(liste)


spade.copy_settings(cyto_session, spade,
                    output = "default", timeout = cyto_session@short_timeout)

cyto_spade<-spade.new(cyto_session, 172509, "bigspade15",
          timeout = cyto_session@long_timeout)

cyto_spade$Value[cyto_spade@Name==down_sampled_events_target]<-15



CytobankAPI::experiments.show(UserSession=cyto_session, experiment_id=172509)
