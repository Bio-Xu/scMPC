source("./Function/1.FunctionalGeneSet.R")

G0051988<-integrationGeneSet(GOTermsID="GO:0051988",AmiGO="GOterms\\0051988.txt",
                             N_pdf="G0051988.pdf")#18

G0090224<-integrationGeneSet(GOTermsID="GO:0090224",AmiGO="GOterms\\0090224.txt",
                             N_pdf="G0090224.pdf")#49

G0030997<-integrationGeneSet(GOTermsID="GO:0030997",AmiGO="GOterms\\0030997.txt",
                             N_pdf="G0030997.pdf")#3

G0046605<-integrationGeneSet(GOTermsID="GO:0046605",AmiGO="GOterms\\0046605.txt",
                             N_pdf="G0046605.pdf")#55

G0051276<-integrationGeneSet(GOTermsID="GO:0051276",AmiGO="GOterms\\0051276.txt",
                             N_pdf="G0051276.pdf")#1101

G0051383<-integrationGeneSet(GOTermsID="GO:0051383",AmiGO="GOterms\\0051383.txt",
                             N_pdf="G0051383.pdf")#24

G0000819<-integrationGeneSet(GOTermsID="GO:0000819",AmiGO="GOterms\\0000819.txt",
                             N_pdf="G0000819.pdf")#239

G0016477<-integrationGeneSet(GOTermsID="GO:0016477",AmiGO="GOterms\\0016477.txt",
                             N_pdf="G0016477.pdf")#1599

G0048870<-integrationGeneSet(GOTermsID="GO:0048870",AmiGO="GOterms\\0048870.txt",
                             N_pdf="G0048870.pdf")#1810

G0034330<-integrationGeneSet(GOTermsID="GO:0034330",AmiGO="GOterms\\0034330.txt",
                             N_pdf="G0034330.pdf")#743

G0031577<-integrationGeneSet(GOTermsID="GO:0031577",AmiGO="GOterms\\0031577.txt",
                             N_pdf="G0031577.pdf")#45

G0045321<-integrationGeneSet(GOTermsID="GO:0045321",AmiGO="GOterms\\0045321.txt",
                             N_pdf="G0045321.pdf")#1017

G0001837<-integrationGeneSet(GOTermsID="GO:0001837",AmiGO="GOterms\\0001837.txt",
                             N_pdf="G0001837.pdf")#172

G0010718<-integrationGeneSet(GOTermsID="GO:0010718",AmiGO="GOterms\\0010718.txt",
                             N_pdf="G0010718.pdf")#64

G0045765<-integrationGeneSet(GOTermsID="GO:0045765",AmiGO="GOterms\\0045765.txt",
                             N_pdf="G0045765.pdf")#430

G0045766<-integrationGeneSet(GOTermsID="GO:0045766",AmiGO="GOterms\\0045766.txt",
                             N_pdf="G0045766.pdf")##220

G0001570<-integrationGeneSet(GOTermsID="GO:0001570",AmiGO="GOterms\\0001570.txt",
                             N_pdf="G0001570.pdf")#84

G0001525<-integrationGeneSet(GOTermsID="GO:0001525",AmiGO="GOterms\\0001525.txt",
                             N_pdf="G0001525.pdf")#624

G0002857<-integrationGeneSet(GOTermsID="GO:0002857",AmiGO="GOterms\\0002857.txt",
                             N_pdf="G0002857.pdf")#7

G0002842<-integrationGeneSet(GOTermsID="GO:0002842",AmiGO="GOterms\\0002842.txt",
                             N_pdf="G0002842.pdf")#6

G0002419<-integrationGeneSet(GOTermsID="GO:0002419",AmiGO="GOterms\\0002419.txt",
                             N_pdf="G0002419.pdf")#5

G0002420<-integrationGeneSet(GOTermsID="GO:0002420",AmiGO="GOterms\\0002420.txt",
                             N_pdf="G0002420.pdf")#9

G0002367<-integrationGeneSet(GOTermsID="GO:0002367",AmiGO="GOterms\\0002367.txt",
                             N_pdf="G0002367.pdf")#101

G0050776<-integrationGeneSet(GOTermsID="GO:0050776",AmiGO="GOterms\\0050776.txt",
                             N_pdf="G0050776.pdf")#1025

G0007155<-integrationGeneSet(GOTermsID="GO:0007155",AmiGO="GOterms\\0007155.txt",
                             N_pdf="G0007155.pdf")#1540

G0033631<-integrationGeneSet(GOTermsID="GO:0033631",AmiGO="GOterms\\0033631.txt",
                             N_pdf="G0033631.pdf")#17

G0044331<-integrationGeneSet(GOTermsID="GO:0044331",AmiGO="GOterms\\0044331.txt",
                             N_pdf="G0044331.pdf")#48

G0006096<-integrationGeneSet(GOTermsID="GO:0006096",AmiGO="GOterms\\0006096.txt",
                             N_pdf="G0006096.pdf")#84

G0006091<-integrationGeneSet(GOTermsID="GO:0006091",AmiGO="GOterms\\0006091.txt",
                             N_pdf="G0006091.pdf")#527

G0045005<-integrationGeneSet(GOTermsID="GO:0045005",AmiGO="GOterms\\0045005.txt",
                             N_pdf="G0045005.pdf")#61

G0006281<-integrationGeneSet(GOTermsID="GO:0006281",AmiGO="GOterms\\0006281.txt",
                             N_pdf="G0006281.pdf")#582

G0071456<-integrationGeneSet(GOTermsID="GO:0071456",AmiGO="GOterms\\0071456.txt",
                             N_pdf="G0071456.pdf")#170

G0042060<-integrationGeneSet(GOTermsID="GO:0042060",AmiGO="GOterms\\0042060.txt",
                             N_pdf="G0042060.pdf")#454

G0006954<-integrationGeneSet(GOTermsID="GO:0006954",AmiGO="GOterms\\0006954.txt",
                             N_pdf="G0006954.pdf")#874

G0007162<-integrationGeneSet(GOTermsID="GO:0007162",AmiGO="GOterms\\0007162.txt",
                             N_pdf="G0007162.pdf")#343


G0007062<-integrationGeneSet(GOTermsID="GO:0007062",AmiGO="GOterms\\0007062.txt",
                             N_pdf="G0007062.pdf")#64

####HallMarks
m_df = msigdbr(species = "Homo sapiens", category = "H")

H<-m_df[grep(pattern="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",m_df$gs_name),]
H1<-H[,"gene_symbol"]%>% unlist()#204
H1<-unique(H1)#200

H<-m_df[grep(pattern="HALLMARK_ANGIOGENESIS",m_df$gs_name),]
H2<-H[,"gene_symbol"]%>% unlist()#36
H2<-unique(H2)#36

H<-m_df[grep(pattern="HALLMARK_DNA_REPAIR",m_df$gs_name),]
H3<-H[,"gene_symbol"]%>% unlist()#170
H3<-unique(H3)#150

H<-m_df[grep(pattern="HALLMARK_GLYCOLYSIS",m_df$gs_name),]
H4<-H[,"gene_symbol"]%>% unlist()#215
H4<-unique(H4)#200

H<-m_df[grep(pattern="HALLMARK_INFLAMMATORY_RESPONSE",m_df$gs_name),]
H5<-H[,"gene_symbol"]%>% unlist()#222
H5<-unique(H5)#200

Cluster1=G0042060
Cluster2=c(G0007162,G0033631,G0044331,G0007155)
Cluster3=G0001837
Cluster4=c(G0016477,G0048870)
Cluster5=c(G0034330,G0051276,G0051383,G0007062,G0000819)
Cluster6=G0010718
Cluster7=c(G0045765,G0045766,G0001570,G0001525)
Cluster8=G0030949
Cluster9=c(G0045005,G0006281)
Cluster10=c(G0051988,G0030997,G0046605,G0090224)
Cluster11=G0031577
Cluster12=G0006096
Cluster13=G0071456
Cluster14=G0006091
Cluster15=c(G0002419,G0002420)
Cluster16=c(G0002857,G0002842)
Cluster17=G0002367
Cluster18=G0050776
Cluster19=G0006954
Cluster20=G0045321
GeneSets<-list(Cluster1=Cluster1,Cluster2=Cluster2,Cluster3=Cluster3,Cluster4=Cluster4,
              Cluster5=Cluster5,Cluster6=Cluster6,Cluster7=Cluster7,Cluster8=Cluster8,
              Cluster9=Cluster9,Cluster10=Cluster10,Cluster11=Cluster11,Cluster12=Cluster12,
              Cluster13=Cluster13,Cluster14=Cluster14,Cluster15=Cluster15,Cluster16=Cluster16,
              Cluster17=Cluster17,Cluster18=Cluster18,Cluster19=Cluster19,Cluster20=Cluster20,
              H1=H1,H2=H2,H3=H3,H4=H4,H5=H5)
saveRDS(GeneSets, file = "GeneSets.rds")
