# library(ggplot2)

plotEffect <- function(preOutM){
  effect = geneNumber = NULL
  gene_effect <- as.data.frame(tapply(preOutM$gene, preOutM$effect, function(t) length(unique(t))))
  gene_effect$effect <- rownames(gene_effect)
  colnames(gene_effect) <- c("geneNumber","effect")
  gene_effect <- subset(gene_effect,select = c("effect","geneNumber"))
  options(scipen=3)


  pic2<-ggplot2::ggplot(gene_effect,ggplot2::aes(x=effect,y=geneNumber,fill=effect)) +
    ggplot2::geom_bar(stat='identity')+
    ggplot2::theme_bw()+
    ggplot2::geom_text(ggplot2::aes(label = geneNumber, vjust = -0.2, hjust = 0.5), size = 3)+
    ggplot2::theme(legend.position='none',axis.text.x = ggplot2::element_text(angle = 30, vjust = 0.85, hjust = 0.75))
  ggplot2::ggsave(pic2, file="EffectPlot.pdf", width=4, height=6)
  ggplot2::ggsave(pic2, file="EffectPlot.png", width=4, height=6)
}

plotCategory <- function(preOutM){
  category = geneNumber = NULL
  gene_category <- as.data.frame(tapply(preOutM$gene, preOutM$categ, function(t) length(unique(t))))
  gene_category$category <- rownames(gene_category)
  colnames(gene_category) <- c("geneNumber","category")
  gene_category <- subset(gene_category,select = c("category","geneNumber"))
  gene_category <- gene_category[!gene_category$categ == "---",]

  pic1<-ggplot2::ggplot(gene_category,ggplot2::aes(x=category,y=geneNumber,fill=category)) +
    ggplot2::geom_bar(stat='identity')+
    ggplot2::theme_bw()+
    ggplot2::geom_text(ggplot2::aes(label = geneNumber, vjust = -0.2, hjust = 0.5), size = 3)+
    ggplot2::theme(legend.position='none',axis.text.x = ggplot2::element_text(angle = 30, vjust = 0.85, hjust = 0.75))
  ggplot2::ggsave(pic1, file="CategoryPlot.pdf", width=4, height=6)
  ggplot2::ggsave(pic1, file="CategoryPlot.png", width=4, height=6)
}

