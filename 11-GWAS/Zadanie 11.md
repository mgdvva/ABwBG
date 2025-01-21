---
title: "GWAS"
author: "Magdalena Jackiewicz"
date: "2025-01-07"
output: html_document
---

# Analiza GWAS

## Trzeba zacząć od wczytywania i załadowania potrzebnych pakietów

```{r}
packages <- c("rrBLUP"
   , "BGLR"
   , "DT"
   , "SNPRelate"
   , "dplyr"
   , "qqman"
   , "poolr")

{for (pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    }
  }
}

library(pkg, character.only = TRUE)
```

## Wczytwanie pliku z danymi genotypowymi, informacjami o osobnikach i imformacje o mapowaniu

```{r}
Geno <- read_ped("C:/Users/Magda/OneDrive - Szkoła Główna Gospodarstwa Wiejskiego/Pulpit/ABWG-GWAS/sativas413.ped")
```

## Kolumny są wczytywane jako osobne wartości

```{r}
p = Geno$p
n = Geno$n
Geno = Geno$x
```

## Podpatrzenie pliku

```{r}
head(Geno)
Geno
```

### Po tej komendzie pokazują się linijki danych zapisane w postaci "0", "2" i "3". Danych jest 985

```{r}
FAM <- read.table("C:/Users/Magda/OneDrive - Szkoła Główna Gospodarstwa Wiejskiego/Pulpit/ABWG-GWAS/sativas413.fam")
head(FAM)
```

```{r}
MAP <- read.table("C:/Users/Magda/OneDrive - Szkoła Główna Gospodarstwa Wiejskiego/Pulpit/ABWG-GWAS/sativas413.map")
head(MAP)
```

### Po wcześniejszych komendach, wczytano pliki .fam i .map, i widać co jest w nich zawarte

## W następnym kroku, trzeba przekodować wartości markerów: 2- NA (brakujące dane); 0-0 (homozygota dla allelu głównego); 1-1 (heterozygota); 3-2 (homozygota dla allelu mniejszościowego). Dane będą bardziej przejrzyste i my będziemy wiedzieć co jest czym

```{r}
Geno[Geno == 2] <- NA
Geno[Geno == 0] <- 0
Geno[Geno == 1] <- 1
Geno[Geno == 3] <- 2
```

## Dalej trzeba przekonwertować dane na macierz i wykonać jej transpozycję

```{r}
Geno <- matrix(Geno, nrow = p, ncol = n, byrow = TRUE)
Geno <- t(Geno)
```

## Dalej trzeba podać wymiary macierzy - liczbę osobników i markery SNP.

```{r}
dim(Geno)
```

### Wynik - 413 wierszy i 36901 kolumn

## Dalej będą wczytywane dane fenotypowe i trzeba sprawdzić ich zgodność z danymi genotypowymi

```{r}
rice.pheno <- read.table("C:/Users/Magda/OneDrive - Szkoła Główna Gospodarstwa Wiejskiego/Pulpit/ABWG-GWAS/RiceDiversity_44K_Phenotypes_34traits_PLINK.txt",
                         header = TRUE, stringsAsFactors = FALSE, sep = "\t")
head(rice.pheno)
```

```{r}
dim(rice.pheno)
```

### Wynik - 413 wierszy i 38 kolumn

## Dalej trzeba przypisać nazwy wierszy dla macierzy Geno na podstawie drugiej kolumny z ramki FAM - V2. Ta kolumna ma identyfikatory próbek

```{r}
rownames(Geno) <- FAM$V2
```

## I na koniec trzeba sprawdzić, czy się wszystko ze sobą zgadza

```{r}
table(rownames(Geno) == rice.pheno$NSFTVID)
```

### Wynik - TRUE i nadal jest 413 wierszy

## Dalej trzeba wyodrębnić pierwszą cechę

```{r}
y <- matrix(rice.pheno$Flowering.time.at.Arkansas)
rownames(y) <- rice.pheno$NSFTVID
index <- !is.na(y)
y <- y[index, 1, drop = FALSE]
Geno <- Geno[index, ]
table(rownames(Geno) == rownames(y))
```

### Tworzona jest tabela częstości z wektora logicznego,a wynik pokazuje ile razy nazwy wierszy są zgodne lub czy się różnią

## Dalej jest przprowadzana kontrola jakości (QC) danych markerowych. Brakujące dane markerowe są zastąpione średnią wartością dla każdego markera

```{r}
for (j in 1:ncol(Geno)){
  Geno[, j] <- ifelse(is.na(Geno[, j]), mean(Geno[, j], nar.rm = TRUE), Geno[, j])
}
```

### Pętla idzie przez wszystkie kolumny macierzy Geno, obliczana jest średnia wartość kolumny, ignorując przy tym NA. Wartość NA jest zastępowana właśnie tą obliczoną średnią i kolumna jest aktualizowana

## Dalej trzeba odflitrować markery z MAF \< 5%

```{r}
p <- colSums(Geno)/(2 * nrow(Geno))
maf <- ifelse(p > 0.5, 1-p, p)
maf.index <- which(maf < 0.05)
Geno1 <- Geno[, -maf.index]
dim(Geno1)
```

### Na początku są obliczane frekewencje allelu mniejszościowego dla każdego SNP, a po tym jest definiowany MAP i trzeba sprawdzić wymiary nowej macierzy. Wynik - 347 wierszy i 36762 kolumn

## Po tym trzeba zaaktualizować plik .map i podać nowe wymiary danych genotypowych i informacje o markerach

```{r}
MAP <- read.table("C:/Users/Magda/OneDrive - Szkoła Główna Gospodarstwa Wiejskiego/Pulpit/ABWG-GWAS/sativas413.map")
dim(MAP)
MAP1 <- MAP[-maf.index, ]
dim(MAP1)
```

### Wynik MAP -36901 wierszy i 4 kolumny, MAP1- 36762 wierszy i 4 kolumny

## Kolejny krok to analiza PCA - tworzenie macierz markerów

```{r}
Geno1 <- as.matrix(Geno1)
sample <- row.names(Geno1)
length(sample)
```

## Obiekt Geno został przekonwertowany na macierz - Geno1 i następnie przypisano do wektora. Wynik - 374 elemntów wektora

```{r}
colnames(Geno1) <- MAP1$V2
snp.id <- colnames(Geno1)
length(snp.id)
```

### Analiza SNP przypisanych do wektora. Wynik 36763 elementów wektora

## Dalej trzeba stworzyć plik GDS

```{r}
snpgdsCreateGeno("44k.gds", genmat = Geno1, sample.id = sample, snp.id = snp.id, 
                 snp.chromosome = MAP1$V1, snp.position = MAP1$V4, snpfirstdim = FALSE)

geno_44k <- snpgdsOpen("44k.gds")
snpgdsSummary("44k.gds")
```

### Wynik - całkowita liczba próbek: 374, całkowita liczba SNP: 36762

## Następnie trzeba przeprowadzić PCA (Analiza Składowych Głównych)

```{r}
pca <- snpgdsPCA(geno_44k, snp.id = colnames(Geno1))
pca <- data.frame(sample.id = row.names(Geno1), 
                  EV1 = pca$eigenvect[, 1], 
                  EV2 = pca$eigenvect[, 2], 
                  EV3 = pca$eigenvect[, 3], 
                  EV4 = pca$eigenvect[, 4], 
                  stringsAsFactors = FALSE)

plot(pca$EV2, pca$EV1, xlab = "PC2", ylab = "PC1")
```

### Ta komenda wykonuje PCA na danych genotypowych , oblicza pierwsze cztery składowe główne, a dalej tworzy ramkę danych zawierającą te składowe dla każdej próbki.Później tworzy wykres 2D z pierwszą składową główną na osi Y i drugą składową główną na osi X, co umożliwia wizualizację struktury danych w przestrzeni tych dwóch składowych.

## W kolejnm etapie trzeba wczytać dodatkowe informacje o próbkach

```{r}
pca_1 <- read.csv("C:/Users/Magda/OneDrive - Szkoła Główna Gospodarstwa Wiejskiego/Pulpit/ABWG-GWAS/RiceDiversity.44K.germplasm.csv", 
                  header = TRUE, skip = 1, stringsAsFactors = FALSE)
pca_2 <- pca_1[match(pca$sample.id, pca_1$NSFTV.ID), ]

pca_population <- cbind(pca_2$Sub.population, pca)
colnames(pca_population)[1] <- "population"

plot(pca_population$EV1, pca_population$EV2, xlab = "PC1", ylab = "PC2", 
     col = c(1:6)[factor(pca_population$population)])
legend(x = "topright", legend = levels(factor(pca_population$population)), 
       col = c(1:6), pch = 1, cex = 0.6)
```

## Dalej dane trzeba przygotować do analizy GWAS - dane genotypowe i fenotypowe

```{r}
geno_final <- data.frame(marker = MAP1[, 2], chrom = MAP1[, 1], pos = MAP1[, 4], 
                         t(Geno1 - 1), check.names = FALSE)

pheno_final <- data.frame(NSFTV_ID = rownames(y), y = y)
```

### W wyniku tej komendy powstaje ramka, która posiada dane kolumny: marker - identyfikatory markerów; chrom - chromosom; pos - pozycja markera; przekształcone dane genotypowe z Geno1 po transpozycji i przekształceniu wartości

## Przeprowadzana jest analiza GWAS

```{r}
GWAS <- GWAS(pheno_final, geno_final, min.MAF = 0.05, P3D = TRUE, plot = FALSE)
```

\##\@ Po tej analizie oszacowano komponenty wariancji

## Dalej wyodrębniono istotne statystyczne markery SNP

```{r}
GWAS_1 <- GWAS %>% filter(y != "0")
GWAS_1 %>% filter(y < 1e-04)
```

## Można zobaczyć, które markery SNP spełniły istotność statystyczną

```{r}
head(GWAS_1)
```

## Żeby zwizualizować dane będzie tworzony wykres Manhattan

```{r}
manhattan(x = GWAS_1, chr = "chrom", bp = "pos", p = "y", snp = "marker", 
          col = c("violet", "green"), suggestiveline = -log10(1e-04), logp = TRUE)
```

### Na wykresie widać wyniki asocjacji genotypów z fenotypami w formie punktów rozłożonych na osiach, reprezentujących różne markery genetyczne na chromosomach.Na osi X są markery genetyczne posortowane według chromosomów, a na osi Y jest logarytmiczna transformacja wartości p z testu asocjacji, im wyższa wartość, tym silniejsza asocjacja. Każdy punkt na osi, to jedno SNP. Widać, które SNP wykazują najsilniejszą statystyczną asocjację z fenotypem.

## *Wnioski*: Dzięki analizie PCA można zobaczyć zróżnicowanie próbek, czyli jak one są rozmieszczone w przestrzeni genotypowej, a dzięki analizie GWAS można zaanalizować SNP, np. przesrotować tylko te, które są istotnie statystycznie i skupić się na nich w kontekście badanej cechy. Połączenie tych dwóch metod analiz pozwala na eliminację błędów oraz zwiększyść dokładność wyników i wyciągnąć tylko te próbki, które są potrzebne do badania.
