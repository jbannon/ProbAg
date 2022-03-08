setwd("~/Desktop/ProbAgony/")
library(TRONCO)

run.model.experiments = function(data,name){
  for(k in seq(20,500,by=10)){
    temp = data[sample(nrow(data), k), ]
    model= import.genotypes(temp)
    model.fit = tronco.capri(model,silent = T)
    print("***should plot here*****")
    print(paste0("./figs/patient_exp/brca/",name,"/random_",k,".pdf"))
    pdf(paste0("./figs/patient_exp/brca/",name,"/random_",k,".pdf"))
    tronco.plot(model.fit,disconnected = T,legend = F)
    dev.off()
  }
}




mutations = read.csv(paste0("./data/_pancandata/brca/point_mutations.csv"),row.names = 1)
options = names(sort(colSums(mutations),decreasing = T))[1:100]
options_short = options[1:10]
options_medium = options[1:30]
options_long = options[1:60]


mutations_s = mutations[,options_short]
mutations_m = mutations[,options_medium]
mutations_l = mutations[,options_long]
mutations_all = mutations[,options]

model = import.genotypes(mutations_s)
model.fit = tronco.capri(model)
pdf("./figs/patient_exp/brca/short/all_patients.pdf")
tronco.plot(model.fit,disconnected = T,legend = F)
dev.off()

model = import.genotypes(mutations_m)
model.fit = tronco.capri(model)
pdf("./figs/patient_exp/brca/medium/all_patients.pdf")
tronco.plot(model.fit,disconnected = T,legend = F)
dev.off()

model = import.genotypes(mutations_l)
model.fit = tronco.capri(model)
pdf("./figs/patient_exp/brca/long/all_patients.pdf")
tronco.plot(model.fit,disconnected = T,legend = F)
dev.off()


run.model.experiments(mutations_s,"short")
run.model.experiments(mutations_m,"medium")
run.model.experiments(mutations_l,"long")


sum(model.fit$adj.matrix.prima.facie)
mutations = read.csv(paste0("./data/_pancandata/ov/point_mutations.csv"),row.names = 1)
options = names(sort(colSums(mutations),decreasing = T))[1:100]
options_short = options[1:10]

mutations_s = mutations[,options_short]
model = import.genotypes(mutations_s)
model.fit = tronco.capri(model)
tronco.plot(model.fit,disconnected = T,legend = F)

library(bnlearn)
model.fit
net = as.bnlearn.network(model.fit)
