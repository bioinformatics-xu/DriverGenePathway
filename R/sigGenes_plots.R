#creating plots

sigGenes_plot <- function(q_result){
gene = lgq.btBinom = lgq.fisher = lgq.lrt = lgq.ct = lgq.projection = gene_projection = gene_beta = gene_fisher = gene_lrt = gene_ct = NULL
dir.create("plots")
setwd("plots")
q_result1 <- q_result[1:100,]

q_result2 <- q_result1[,c(1,2)]
q_result2$lgq.btBinom <- -log10(q_result2$q.btBinom+0.00001)
b <- ggplot2::ggplot(q_result2, ggplot2::aes(x=gene,y=lgq.btBinom))+
  ggplot2::geom_point(ggplot2::aes(size=2*(lgq.btBinom^2)), color = "#00AFBB") +
  ggrepel::geom_text_repel(ggplot2::aes(label = gene), size = 2.5)+
  ggplot2::theme_minimal()+
  ggplot2::theme(legend.background = ggplot2::element_rect(fill="lightblue"))+
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45,size = 4))+
  ggplot2::labs(title="Plot of result genes tested by beta-binomial", x ="Genes", y = "-log(q.beta-binomial)")
ggplot2::ggsave(b, file="beta-binomial.pdf", width = 20)

q_result2 <- q_result1[,c(1,3)]
q_result2$lgq.fisher <- -log10(q_result2$q.fisher+0.00001)
b <- ggplot2::ggplot(q_result2, ggplot2::aes(x=gene,y=lgq.fisher))+
  ggplot2::geom_point(ggplot2::aes(size=2*(lgq.fisher^2)), color = "#00AFBB") +
  ggrepel::geom_text_repel(ggplot2::aes(label = gene), size = 2.5)+
  ggplot2::theme_minimal()+
  ggplot2::theme(legend.background = ggplot2::element_rect(fill="lightblue"))+
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45,size = 4))+
  ggplot2::labs(title="Plot of result genes tested by fisher", x ="Genes", y = "-log(q.fisher)")
ggplot2::ggsave(b, file="fisher.pdf", width = 20)

q_result2 <- q_result1[,c(1,4)]
q_result2$lgq.lrt <- -log10(q_result2$q.lrt+0.00001)
b <- ggplot2::ggplot(q_result2, ggplot2::aes(x=gene,y=lgq.lrt))+
  ggplot2::geom_point(ggplot2::aes(size=2*(lgq.lrt^2)), color = "#00AFBB") +
  ggrepel::geom_text_repel(ggplot2::aes(label = gene), size = 2.5)+
  ggplot2::theme_minimal()+
  ggplot2::theme(legend.background = ggplot2::element_rect(fill="lightblue"))+
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45,size = 4))+
  ggplot2::labs(title="Plot of result genes tested by lrt", x ="Genes", y = "-log(q.lrt)")
ggplot2::ggsave(b, file="lrt.pdf", width = 20)

q_result2 <- q_result1[,c(1,5)]
q_result2$lgq.ct <- -log10(q_result2$q.ct+0.00001)
b <- ggplot2::ggplot(q_result2,ggplot2::aes(x=gene,y=lgq.ct))+
  ggplot2::geom_point(ggplot2::aes(size=2*(lgq.ct^2)), color = "#00AFBB") +
  ggrepel::geom_text_repel(ggplot2::aes(label = gene), size = 2.5)+
  ggplot2::theme_minimal()+
  ggplot2::theme(legend.background = ggplot2::element_rect(fill="lightblue"))+
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45,size = 4))+
  ggplot2::labs(title="Plot of result genes tested by ct", x ="Genes", y = "-log(q.ct)")
ggplot2::ggsave(b, file="ct.pdf", width=20)

q_result2 <- q_result1[,c(1,6)]
q_result2$lgq.projection <- -log10(q_result2$q.projection+0.00001)
b <- ggplot2::ggplot(q_result2,ggplot2::aes(x=gene,y=lgq.projection))+
  ggplot2::geom_point(ggplot2::aes(size=2*(lgq.projection^2)), color = "#00AFBB") +
  ggrepel::geom_text_repel(ggplot2::aes(label = gene), size = 2.5)+
  ggplot2::theme_minimal()+
  ggplot2::theme(legend.background = ggplot2::element_rect(fill="lightblue"))+
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45,size = 4))+
  ggplot2::labs(title="Plot of result genes tested by 2D projection", x ="Genes", y = "-log(q.projection)")
ggplot2::ggsave(b, file="projection.pdf", width=20)

VennDiagram::venn.diagram(x=list(gene_projection=gene_projection,gene_beta=gene_beta,gene_fisher=gene_fisher,gene_lrt=gene_lrt,gene_ct=gene_ct),
             filename = "vennplot",
             fill =c("cornflowerblue","green","yellow","darkorchid1","red"),
             cat.col =c("darkblue", "darkgreen", "orange","darkorchid4","black"),cat.cex = 1,cat.pos = 0, cat.dist = 0.07,cat.fontfamily = "serif")

setwd("..")
}
