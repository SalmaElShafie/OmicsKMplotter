#' Draws Kaplan Meier survival curves from all TCGA dataset types (mRNA, miRNA, rnaseq, mutations, RPPA, methylation) classifying patients into groups using any number of user-identified targets from the dataset; where groups are divided as patients with high/low expression of this genes with the mean being the cutoff (higher or lower than mean in TCGA cohort for that cancertype= low or high groups in the KM plot). For example if 2 genes were selected as targets, there would be 4 groups: 1.high expression of both genes 2. high expression of first but low expression of second 3. low expression of first gene but high expression of second 4. low expression of both genes.
#' @param cancertype cancer type using TCGA format (ex. BRCA)
#' @param datasettype omic type using RTCGA format (ex. mRNA, miRNASeq, methylation, rnaseq, mutations, RPPA)
#' @param UserList vector containing any number of target genes/proteins/probes of interest using TCGA/RTCGA format for the corresponding datasettypes (ex. when using miRNASeq datasettype: UserList=c("hsa-let-7d","hsa-let-7e", "hsa-let-7a3"), rnaseq: UserList=c("AADAC|13", "AADAT|51166", "AAGAB|79719"), RPPA:UserList=c("ACC1", "AR", "ACVRL1"), mRNA: UserList=c("ELMO2", "CREB3L1", "PNMA1"), methylation: UserList= c("cg00000292","cg00002426","cg00003994"), mutations: UserList= c("TP53", "PIK3CA", "SOX15", "MSH3")
#' @return graph of Kaplan Meier survival curve with n categories of patients based on selected number of target features and two possible categories for each (high and low expression)
#' @export export this function
#' @examples OmicsKMplotter(cancertype= "OV", datasettype = "miRNASeq", UserList=c("hsa-let-7d","hsa-let-7e", "hsa-let-7a3")) or OmicsKMplotter(cancertype= "OV", datasettype = "rnaseq", UserList=c("AADAC|13", "AADAT|51166", "AAGAB|79719")) or OmicsKMplotter(cancertype= "OV", datasettype = "RPPA", UserList=c("ACC1", "AR", "ACVRL1")) or OmicsKMplotter(cancertype= "OV", datasettype = "mRNA", UserList=c("ELMO2", "CREB3L1", "PNMA1")) or OmicsKMplotter(cancertype= "BRCA", datasettype = "methylation", UserList= c("cg00000292","cg00002426","cg00003994")) or OmicsKMplotter(cancertype= "BRCA", datasettype = "mutations", UserList= c("TP53", "PIK3CA", "SOX15", "MSH3"))
#' @importFrom magrittr %>%

OmicsKMplotter<- function(cancertype, datasettype, UserList) {

  j=paste("RTCGA.",datasettype, "::", sep ="", rlang::parse_expr(paste(cancertype, datasettype, sep=".")))
  c=rlang::parse_expr(j)
  SelectedDataset=data.frame()

  d=paste0("RTCGA.clinical::",rlang::parse_expr(paste0(cancertype, ".clinical")))
  e=rlang::parse_expr(d)
  as.data.frame(eval(e)) ->>clinical_Set


  if (datasettype == "miRNASeq") {
    as.data.frame(eval(c)) -> SelectedDataset
    sd= dplyr::filter(SelectedDataset, miRNA_ID =="reads_per_million_miRNA_mapped")
    SelectedDataset= sd %>% dplyr::mutate(bcr_patient_barcode= tolower(substr(rownames(sd), 1, 12)))
    colnames(SelectedDataset) = gsub("-","",colnames(SelectedDataset))
    UserList= gsub("-","",UserList)
    jointdataset <-merge (clinical_Set,SelectedDataset , by.x = 'patient.bcr_patient_barcode', by.y ='bcr_patient_barcode')
    Surv_features<- RTCGA::survivalTCGA(jointdataset, extract.cols=UserList)
    for(i in 1: length(UserList))
    {
      Surv_features_Per_miRNA<- cut(as.numeric(Surv_features[,UserList[i]]), breaks=c(-Inf, mean(as.numeric(Surv_features[,UserList[i]])), Inf), labels=c("Low","High"))

      Surv_features <-cbind(Surv_features,Surv_features_Per_miRNA )
      names(Surv_features)[names(Surv_features) == "Surv_features_Per_miRNA"] <-paste(UserList[i],"Cat",sep="")
    }
    UserList= paste(UserList,"Cat", sep = "")
    d=rlang::parse_expr(paste0("survival::Surv(times, patient.vital_status)~",(paste(UserList, collapse=" + "))))
    ffit <- survminer::surv_fit(as.formula(d) , data=Surv_features) #####works toooooo!!!
    graph=  survminer::ggsurvplot(ffit, conf.int=TRUE, pval=TRUE)
    return(graph)
  }

  else if (datasettype == "mutations") {
    as.data.frame(eval(c)) -> SelectedDataset
    ee=rlang::parse_expr(paste0(cancertype, ".clinical"))
    #RTCGA::survivalTCGA(eval(ee)) -> data.surv
    RTCGA::survivalTCGA(clinical_Set) -> data.surv
    SelectedDataset %>%
      dplyr::mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) %>%
      dplyr::select(bcr_patient_barcode) %>%
      unique -> patients_with_mutations_information
    data.surv %>%
      dplyr::filter(bcr_patient_barcode %in%
                      patients_with_mutations_information$bcr_patient_barcode) -> patients_with_survival_and_mutations_info
    SelectedDataset %>%
      dplyr::filter(Hugo_Symbol %in% UserList) %>%
      dplyr::filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% # cancer tissue
      dplyr::mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12))  -> data.mutations

    patients_with_survival_and_mutations_info %>%
      dplyr::left_join(data.mutations,
                       by = "bcr_patient_barcode") %>%
      dplyr::select(times, bcr_patient_barcode, patient.vital_status, Hugo_Symbol) -> data.clinical_mutations
    slimZ=data.clinical_mutations %>%
      count(bcr_patient_barcode, Hugo_Symbol) %>%
      tidyr::spread(Hugo_Symbol, n, fill = 0)
    jointdataset1 <-merge (data.surv,slimZ , by.x = 'bcr_patient_barcode', by.y ='bcr_patient_barcode')
    d=rlang::parse_expr(paste0("survival::Surv(times, patient.vital_status)~",(paste(UserList, collapse=" + "))))
    ffit <- survminer::surv_fit(as.formula(d) , data=jointdataset1) #####works toooooo!!!
    graph=  survminer::ggsurvplot(ffit, conf.int=TRUE, pval=TRUE)
    return(graph)
  }
  else if (datasettype== "methylation") {
    as.data.frame(eval(c)) ->SelectedDataset
    SelectedDataset=SelectedDataset %>% dplyr::filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% dplyr::mutate(patient.bcr_patient_barcode = tolower(substr(bcr_patient_barcode, 1, 12)))
    jointdataset <-merge (clinical_Set,SelectedDataset , by.x = 'patient.bcr_patient_barcode', by.y ='patient.bcr_patient_barcode')
    Surv_features<- RTCGA::survivalTCGA(jointdataset, extract.cols=UserList)
    for(i in 1: length(UserList))
    {
      Surv_features_Per_mRNA<- cut(as.numeric(Surv_features[,UserList[i]]), breaks=c(-Inf, mean(as.numeric(Surv_features[,UserList[i]])), Inf), labels=c("Low","High"))

      Surv_features <-cbind(Surv_features,Surv_features_Per_mRNA )
      names(Surv_features)[names(Surv_features) == "Surv_features_Per_mRNA"] <-paste(UserList[i],"Cat",sep="")
    }
    UserList= paste(UserList,"Cat", sep = "")
    d=rlang::parse_expr(paste0("survival::Surv(times, patient.vital_status)~",(paste(UserList, collapse=" + "))))
    ffit <- survminer::surv_fit(as.formula(d) , data=Surv_features) #####works toooooo!!!
    graph=  survminer::ggsurvplot(ffit, conf.int=TRUE, pval=TRUE)
    return(graph)
  }
  else if (datasettype == "mRNA") {
    as.data.frame(eval(c)) ->SelectedDataset
    SelectedDataset=SelectedDataset %>% dplyr::filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% dplyr::mutate(patient.bcr_patient_barcode = tolower(substr(bcr_patient_barcode, 1, 12)))
    jointdataset <-merge (clinical_Set,SelectedDataset , by.x = 'patient.bcr_patient_barcode', by.y ='patient.bcr_patient_barcode')
    Surv_features<- RTCGA::survivalTCGA(jointdataset, extract.cols=UserList)
    for(i in 1: length(UserList))
    {
      Surv_features_Per_mRNA<- cut(as.numeric(Surv_features[,UserList[i]]), breaks=c(-Inf, mean(as.numeric(Surv_features[,UserList[i]])), Inf), labels=c("Low","High"))
      Surv_features <-cbind(Surv_features,Surv_features_Per_mRNA )
      names(Surv_features)[names(Surv_features) == "Surv_features_Per_mRNA"] <-paste(UserList[i],"Cat",sep="")
    }
    UserList= paste(UserList,"Cat", sep = "")
    d=rlang::parse_expr(paste0("survival::Surv(times, patient.vital_status)~",(paste(UserList, collapse=" + "))))
    ffit <- survminer::surv_fit(as.formula(d) , data=Surv_features) #####works toooooo!!!
    graph=  survminer::ggsurvplot(ffit, conf.int=TRUE, pval=TRUE)
    return(graph)
  }
  else if (datasettype == "RPPA") {
    as.data.frame(eval(c)) ->SelectedDataset
    #return(RPPA_Set)
    SelectedDataset=SelectedDataset %>% dplyr::filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% dplyr::mutate(patient.bcr_patient_barcode = tolower(substr(bcr_patient_barcode, 1, 12)))
    colnames(SelectedDataset)=gsub("-", "", colnames(SelectedDataset))
    colnames(SelectedDataset)=gsub("_", "", colnames(SelectedDataset))
    lista=colnames(SelectedDataset)
    colnamesmod=c(lista[1], paste("p", lista[2:(length(lista)-1)],sep=""), lista[length(lista)])
    colnames(SelectedDataset)= colnamesmod
    UserList = gsub("-", "",UserList)
    UserList = gsub("_", "",UserList)
    UserList= paste("p", UserList, sep="")
    jointdataset <-merge (clinical_Set,SelectedDataset , by.x = 'patient.bcr_patient_barcode', by.y ='patient.bcrpatientbarcode')
    Surv_features<- RTCGA::survivalTCGA(jointdataset, extract.cols=UserList)
    for(i in 1: length(UserList))
    {
      Surv_features_Per_RPPA<- cut(as.numeric(Surv_features[,UserList[i]]), breaks=c(-Inf, mean(as.numeric(Surv_features[,UserList[i]])), Inf), labels=c("Low","High"))
      Surv_features <-cbind(Surv_features,Surv_features_Per_RPPA )
      names(Surv_features)[names(Surv_features) == "Surv_features_Per_RPPA"] <-paste(UserList[i],"Cat",sep="")
    }
    UserList= paste(UserList,"Cat", sep = "")
    d=rlang::parse_expr(paste0("survival::Surv(times, patient.vital_status)~",(paste(UserList, collapse=" + "))))
    ffit <- survminer::surv_fit(as.formula(d) , data=Surv_features) #####works toooooo!!!
    graph=  survminer::ggsurvplot(ffit, conf.int=TRUE, pval=TRUE)
    return(graph)
  }
  else if (datasettype =="rnaseq") {
    as.data.frame(eval(c)) ->SelectedDataset
    SelectedDataset=SelectedDataset %>% dplyr::filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% dplyr::mutate(patient.bcr_patient_barcode = tolower(substr(bcr_patient_barcode, 1, 12)))
    colnames(SelectedDataset)=gsub("[?][|]", "g", colnames(SelectedDataset))
    colnames(SelectedDataset)=gsub("[|]", "", colnames(SelectedDataset))
    UserList = gsub("[?][|]","g",UserList)
    UserList = gsub("[|]", "",UserList)
    jointdataset <-merge (clinical_Set,SelectedDataset , by.x = 'patient.bcr_patient_barcode', by.y ='patient.bcr_patient_barcode')
    Surv_features<- RTCGA::survivalTCGA(jointdataset, extract.cols=UserList)
    for(i in 1: length(UserList))
    {
      Surv_features_Per_rnaseq<- cut(as.numeric(Surv_features[,UserList[i]]), breaks=c(-Inf, mean(as.numeric(Surv_features[,UserList[i]])), Inf), labels=c("Low","High"))
      Surv_features <-cbind(Surv_features,Surv_features_Per_rnaseq )
      names(Surv_features)[names(Surv_features) == "Surv_features_Per_rnaseq"] <-paste(UserList[i],"Cat",sep="")
    }
    UserList= paste(UserList,"Cat", sep = "")
    d=rlang::parse_expr(paste0("survival::Surv(times, patient.vital_status)~",(paste(UserList, collapse=" + "))))
    ffit <- survminer::surv_fit(as.formula(d) , data=Surv_features) #####works toooooo!!!
    shortest= c(paste0(gene1, "=0, ", gene2, "=0"), paste0(gene1, "=0, ", gene2, "=1"), paste0(gene1, "=1, ", gene2, "=0"), paste0(gene1, "=1, ", gene2, "=1"))
    names(sfit$strata) = shortest
    graph=  survminer::ggsurvplot(ffit, conf.int=TRUE, pval=TRUE)
    return(graph)
  }

}

####examples of using function and it works
#OmicsKMplotter(cancertype= "OV", datasettype = "miRNASeq", UserList=c("hsa-let-7d","hsa-let-7e", "hsa-let-7a3"))
#OmicsKMplotter(cancertype= "OV", datasettype = "rnaseq", UserList=c("AADAC|13", "AADAT|51166", "AAGAB|79719"))
#OmicsKMplotter(cancertype= "OV", datasettype = "RPPA", UserList=c("ACC1", "AR", "ACVRL1"))
#OmicsKMplotter(cancertype= "OV", datasettype = "mRNA", UserList=c("ELMO2", "CREB3L1", "PNMA1"))
#OmicsKMplotter(cancertype= "BRCA", datasettype = "methylation", UserList= c("cg00000292","cg00002426","cg00003994"))
#OmicsKMplotter(cancertype= "BRCA", datasettype = "mutations", UserList= c("TP53", "PIK3CA", "SOX15", "MSH3"))

