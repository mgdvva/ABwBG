---
title: "Zadanie 13"
author: "Magdalena Jackiewicz"
date: "2025-01-14"
output: html_document
---

# Proste wykresy eksploracyjne

## Instalacja i załadowanie potrzebnych pakietów do zobrazowania danych

```{r}
install.packages("ggplot2") 
library(ggplot2)
```

## Załadowanie przykładowego zbioru danych: iris

```{r}
data(iris)
```

## Boxplot

```{r}
ggplot(iris, aes(x = Species, y = Sepal.Length, fill=Species)) +
  geom_boxplot() +
  labs(title="Boxplot - Długość działki kielicha (Sepal.Length) wg gatunków irysa")
```

### Na wykresie typu Boxplot widać zależność długości działki kielicha od gatunku Irysa. Najdłuższe są u gatunku virginica.

## 2. Histogram

```{r}
ggplot(iris, aes(x = Sepal.Width)) +
  geom_histogram(binwidth = 0.2, fill="yellow", color="red") +
  labs(title="Histogram - Szerokość działki kielicha (Sepal.Width)")
```

### Na podanym histogramie widać, że najczęstsza szerokość działki kielicha oscyluje około 3.

## Wykres punktowy

```{r}
ggplot(iris, aes(x=Sepal.Length, y=Petal.Length, color=Species)) +
  geom_point() +
  labs(title="Scatter plot - Iris")
```

### Ten wykres punktowy przedstawia zależność długości płatków od długości działki kielicha u różnych gatunków Irysa

## Wykres skrzypcowy + boxplot

```{r}
ggplot(iris, aes(x=Species, y=Sepal.Width, fill=Species)) +
  geom_violin(alpha=0.5) +
  geom_boxplot(width=0.1, color="black", outlier.shape=NA) +
  labs(title="Violin + Boxplot - Iris")
```

### Połączenie tych 2 wykresów przedstawia dane dotyczące szerokości działki kielicha dla 3 gatunków Irysa.

# Załadowanie kolejnych danych przykładowych - StackedBar Plot

## Skumulowany wykres słupkowy

```{r}
df_bar <- data.frame(
  Sample = rep(c("S1","S2","S3"), each=3),
  Category = rep(c("A","B","C"), times=3),
  Count = c(10, 5, 15, 8, 12, 6, 20, 10, 5)
)
ggplot(df_bar, aes(x=Sample, y=Count, fill=Category)) +
  geom_bar(stat="identity") +
  labs(title="Skumulowany wykres słupkowy")
```

# Waffle Plot

## Instalacja i załadowanie potrzebnych pakietów

```{r}
install.packages("waffle")
library(waffle)
```

## Załadowanie przykładowych danych

```{r}
parts <- c(`Cat A (50)`=50, `Cat B (30)`=30, `Cat C (20)`=20)
```

```{r}
waffle(parts, rows=5, 
       title="Przykładowy Waffle Plot",
       legend_pos = "bottom")
```

# Time Series Plot

## Załadowanie przykładowych danych

```{r}
df_time <- data.frame(
  Time = rep(1:5, times=2),
  Expression = c(1,2,3,2.5,4, 2,3.5,3,4,5),
  Gene = rep(c("GeneA","GeneB"), each=5)
)
```

```{r}
ggplot(df_time, aes(x=Time, y=Expression, color=Gene)) +
  geom_line() +
  geom_point() +
  labs(title="Analiza czasowa ekspresji genów")
```

### Wykres przedstawia analizę ekspresji 2 genów w czasie

# Waterfall Plot

# Klasyczny Waterfall

## Załadowanie sztucznych danych: zmiana wielkości guza w %

```{r}
set.seed(123)
df_wf <- data.frame(
  Pacjent = paste0("P", 1:20),
  Zmiana = sample(seq(-100, 100, by=10), 20)
)
```

## Sortowanie według wartości

```{r}
df_wf <- df_wf[order(df_wf$Zmiana), ]
```

# Prosty "waterfall" z ggplot2 - barplot z uporządkowanymi słupkami

```{r}
df_wf$Pacjent <- factor(df_wf$Pacjent, levels=df_wf$Pacjent)
```

```{r}
ggplot(df_wf, aes(x=Pacjent, y=Zmiana, fill=Zmiana>0)) +
  geom_bar(stat="identity") +
  coord_flip() +
  scale_fill_manual(values=c("pink","black"), name="Zmiana dodatnia?") +
  labs(title="Klasyczny Waterfall Plot - Zmiana wielkości",
       x="Pacjent", y="Zmiana (%)")
```

### Na podanym wykresie widać jak zmienia się wielkość guza w % u poszczególnych pacjentów

# Volcano Plot

## Załadowanie przykładowych danych

```{r}
set.seed(123)
df_volcano <- data.frame(
  gene = paste0("Gene", 1:100),
  log2FC = rnorm(100, 0, 1),
  pval = runif(100, 0, 0.05)
)
df_volcano$negLogP <- -log10(df_volcano$pval)
```

```{r}
plot(df_volcano$log2FC, df_volcano$negLogP,
     pch=20, col="grey50",
     xlab="log2 Fold Change", ylab="-log10(p-value)",
     main="Volcano Plot (base R)")

abline(h=-log10(0.05), col="red", lty=2)
abline(v=c(-1, 1), col="blue", lty=2)
```

### Na tym wykresie czerwone linie przerwane określa próg istotności statystycznej a niebieskie mogą być ustalane według danego badania

# HeatMap

## Instalacja i załadowanie potrzebnych pakietów

```{r}
install.packages("pheatmap")
library(pheatmap)
```

## Losowa macierz 10 genów x 5 próbek

```{r}
set.seed(123)
mat <- matrix(rnorm(50), nrow=10, ncol=5)
rownames(mat) <- paste0("Gene", 1:10)
colnames(mat) <- paste0("Sample", 1:5)

pheatmap(mat, 
         scale="row", 
         cluster_rows=TRUE, 
         cluster_cols=TRUE,
         main="Heatmap - 10 genów x 5 próbek")
```

### Dzięki Heatmapie można zaobserwować ugrupowanie podanych elementów

# Wykresy redukcji wymiarów

## PCA

```{r}
data(iris)
pca_result <- prcomp(iris[,1:4], center = TRUE, scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  Species = iris$Species
)

ggplot(pca_df, aes(x=PC1, y=PC2, color=Species)) +
  geom_point() +
  labs(title="PCA - Iris")
```

### Na podanym wykresie widać jak zredukowano wielowymiarowe dane do 2 wymiarów i można porównać według tych wymiarów 3 gatunki Irysa

## t-SNE

## Instalacja i załadowanie potrzebnych pakietów

```{r}
install.packages("Rtsne")
library(Rtsne)
```

## Usunięcie duplikatów względem kolumn 1:4 (Sepal.Length, Sepal.Width, Petal.Length, Petal.Width)

```{r}
iris_nodup <- iris[!duplicated(iris[,1:4]), ]
```

## Wywołanie wykresu Rtsne

```{r}
tsne_out <- Rtsne(iris_nodup[,1:4], pca=FALSE, perplexity=20, max_iter=500)

library(ggplot2)

# Tworzymy data.frame z wynikami t-SNE
tsne_df <- data.frame(
  X = tsne_out$Y[,1],
  Y = tsne_out$Y[,2],
  Species = iris_nodup$Species  # bo usunęliśmy te same wiersze
)

ggplot(tsne_df, aes(x=X, y=Y, color=Species)) +
  geom_point() +
  labs(title="t-SNE - Iris (bez duplikatów)")
```

### Na tym wykresie widać 2 czynniki wpływające na pewne cechy u 3 gatunków Irysa już bez duplikacji próbek

# Manhattan Plot

## Instalacja i załadowanie potrzebnych pakietów

```{r}
install.packages("qqman")
library(qqman)
```

# Wygenerowanie przykładu: 500 SNP w 5 chromosomach

```{r}
set.seed(123)
chrom <- rep(1:5, each=100)
bp <- rep(1:100, times=5)
pval <- runif(500, min=0, max=0.1)
df_gwas <- data.frame(CHR=chrom, BP=bp, P=pval, SNP=paste0("rs",1:500))

manhattan(df_gwas,
          genomewideline = -log10(5e-8),
          suggestiveline = -log10(1e-3),
          main="Przykładowy Manhattan Plot")
```

### Na tym wykresie widać istotność danych SNP i na którym chromosomie się znajdują, jeśli punkty (SNP) są ponad tą niebieską linią, znaczy że są istotnie statystycznie

# Diagram Venna

## Instalacja i załadowanie potrzebnych pakietów

```{r}
install.packages("VennDiagram")
library(VennDiagram)
library(grid)
```

```{r}
setA <- paste0("Gene", 1:10)
setB <- paste0("Gene", 6:15)

venn <- venn.diagram(
  x = list(A=setA, B=setB),
  filename = NULL,
  fill = c("skyblue", "pink"),
  alpha = c(0.5, 0.5),
  cex = 2,
  cat.cex = 2
)
grid.newpage()
grid.draw(venn)
```

### Diagram Venna pozwala na porównanie zbiorów i znalezienie ich części wspólnej

# UpSet Plot

## Instalacja i załadowanie potrzebnych pakietów

```{r}
install.packages("UpSetR")
library(UpSetR)
```

```{r}
listInput <- list(
  SetA = setA,
  SetB = setB,
  SetC = paste0("Gene", 8:12)
)

upset(fromList(listInput), 
      order.by = "freq", 
      main.bar.color = "steelblue",
      sets.bar.color = "tomato")
```

### Wykres UpSet Plot pozwala na analizę dużych danych z wieloma zbiorami i czy dane zbiory mają wspólne cechy

# Pathway and Annotation Plots

## Instalacja i załadowanie potrzebnych pakietów

```{r}
BiocManager::install("pathview")
library(pathview)
```

## Załadowanie przykładowych sztucznych danych

```{r}
genelist <- c("1840"=2, "4609"=-1, "7124"=1)  # Entrez ID
```

## Wywołanie wykresu

```{r}
pv.out <- pathview(gene.data = genelist,
                   pathway.id = "hsa04110",
                   species = "hsa",
                   out.suffix="example")
```

### Ten wykres zapisuje sie w pliku .png, wykres pokazuje różne ścieżki biologiczne i ich powiązania z innymi danymi

# Drzewo filogenetyczne

## Instalacja i załadowanie potrzebnych pakietów

```{r}
install.packages("ape")
library(ape)
```

## Wywołanie wykresu

```{r}
tree <- rtree(10)  # losowe drzewo z 10 taksonami
plot(tree, main="Random Phylogenetic Tree")
```

### Drzewo filogenetyczne pokazuje powiązania między danymi, pokazując ich wspólne pochodzenie

# Synteny Plots

## Instalacja i załadowanie potrzebnych pakietów

```{r}
install.packages("genoPlotR")
library(genoPlotR)
```

## Przykładowe dane

```{r}
data("barto", package="genoPlotR")

plot_gene_map(dna_segs = barto$dna_segs,
              comparisons = barto$comparisons,
              main = "Synteny plot - Bartonella genomes (genoPlotR)")
```

### Ten wykres przedstawia porównanie i powiązania między genomami bakterii *Bartonella*

# Circos Plot

## Instalacja i załadowanie potrzebnych pakietów

```{r}
install.packages("circlize")
library(circlize)

library(dplyr)
library(circlize)
```

## Przygotowanie zakresów sektorów

```{r}
bed <- data.frame(
  chr   = c("chr1","chr1","chr2","chr2"),
  start = c(1, 50, 1, 50),
  end   = c(25, 75, 25, 75),
  value = c(10, 20, 5, 15)
)
```

## Grupowanie, żeby wyliczyć minimalny start i maksymalny end dla każdego chromosomu

```{r}
chr_ranges <- bed %>%
  group_by(chr) %>%
  summarise(
    min_start = min(start),
    max_end   = max(end)
  )
```

```{r}
circos.clear()
circos.initialize(
  factors = chr_ranges$chr, 
  xlim    = cbind(chr_ranges$min_start, chr_ranges$max_end)
)

circos.trackPlotRegion(
  ylim = c(0, 1),
  panel.fun = function(x, y) {
    # Odczytujemy informację o sektorze
    sector.name = CELL_META$sector.index
    # Rysujemy napis na środku sektora
    circos.text(
      CELL_META$xcenter,
      0.5,
      sector.name,
      facing = "bending.inside"
    )
  }
)
```

```{r}
for(i in seq_len(nrow(bed))) {
  # Wyciągamy chrom, start, end
  chr   = bed$chr[i]
  start = bed$start[i]
  end   = bed$end[i]
  val   = bed$value[i]
   circos.rect(
    xleft       = start, 
    ybottom     = 0, 
    xright      = end, 
    ytop        = val/20, 
    sector.index= chr,
    col         = "skyblue", 
    border      = "black"
  )
}
circos.clear()
```

# Ideograms

## Instalacja i załadowanie potrzebnych pakietów

```{r}
BiocManager::install("karyoploteR")
library(karyoploteR)
```

## Wywołanie wykresu

```{r}
kp <- plotKaryotype(genome="hg19")
```

# Przykładowo można zaznaczyć region chr1

```{r}
plot.new()  
region <- toGRanges(data.frame(chr="chr1", start=1e6, end=2e6))
kpRect(kp, data=region, y0=0, y1=1, col="green", border=NA)
```

### Ten typ wykresu pozwala na przedstawienie struktury chromosomów, pokazuje ich rozmieszczenie i inne cechy, można dodatkowo zaznczyć konkrenty region, który nas interesuje
