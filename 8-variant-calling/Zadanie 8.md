---
title: "Zadanie 8"
author: "Magdalena Jackiewicz"
date: "2024-12-03"
output: html_document
---

# Trzeba zacząć od instalacji potrzebnych pakietów i ich załadowania

```{r}
BiocManager::install("VariantTools")
BiocManager::install("Rsamtools")
BiocManager::install("GenomicRanges")
BiocManager::install("VariantAnnotation")
BiocManager::install("GenomicFeatures")
BiocManager::install("BiocParallel")
```

```{r}
library(VariantTools)
library(Rsamtools)
library(GenomicRanges)
library(GenomicFeatures)
library(VariantAnnotation)
library(BiocParallel)
```

# Następnie trzeba ustawić katalog roboczy i sprawdzić dostępność tych danych

```{r}
setwd("C:/Users/Magda/OneDrive - Szkoła Główna Gospodarstwa Wiejskiego/Pulpit/ABWG-GWAS")
```

# Można sprawdzić czy pliki są dostępne

```{r}
list.files()
```

## Po tej komendzie pokazuje się jakie mamy pliki w tym katalogu roboczym

# Dalej trzeba wczytać plik BAM i genom referencyjny, który jest w formacie FASTA

```{r}
bamfile <- "C:/Users/Magda/OneDrive - Szkoła Główna Gospodarstwa Wiejskiego/Pulpit/ABWG-GWAS/aligned_sample.BAM"
bam <- BamFile(bamfile)

ref_genome <- "C:/Users/Magda/OneDrive - Szkoła Główna Gospodarstwa Wiejskiego/Pulpit/ABWG-GWAS/ecoli_reference.fasta"
fa <- FaFile(ref_genome)
```

## Plik BAM to binarny format pliku używany w bioinformatyce do przechowywania informacji o sekwencjach DNA lub RNA wyrównanych do odniesienia genomowego. Jest plikiem skompresowanym, dzięki czemu można w nim przechowywać duże biory danych genomowych

# Następnie trzeba plik BAM przesortować według współrzędnych

```{r}
input_bam <- "C:/Users/Magda/OneDrive - Szkoła Główna Gospodarstwa Wiejskiego/Pulpit/ABWG-GWAS/aligned_sample.BAM"
output_bam <- "C:/Users/Magda/OneDrive - Szkoła Główna Gospodarstwa Wiejskiego/Pulpit/ABWG-GWAS/sorted_aligned_sample.BAM"
```

## Zdefiniowano ściężkę wejściową i wyjściową pliku BAM

```{r}
sortBam(file = input_bam, destination = output_bam, overwrite = TRUE)
sorted_bam <- "C:/Users/Magda/OneDrive - Szkoła Główna Gospodarstwa Wiejskiego/Pulpit/ABWG-GWAS/sorted_aligned_sample.BAM.bam"
```

## Przesortowano plik BAM i zdefiniowano już przesortowany plik

# W kolejnym kroku trzeba zaindeksować plik FASTA i przesortowany plik BAM

```{r}
indexFa(ref_genome)
indexBam(sorted_bam)
```

# Następnie trzeba przeporwadzić kontrolę jakości danych sekwencyjnych

```{r}
scanBamHeader(bam)
```

## Ta komenda pokazuje nagłówek pliku BAM

```{r}
idxstats <- idxstatsBam(sorted_bam)
print(idxstats)
```

## Tą komendą można sprawdzić podstawowe statystyki pliku BAM

# Następnie trzeba obliczyć pokrycie genomu

```{r}
coverage_data <- coverage(sorted_bam)
summary(coverage_data[[1]]) # dla genomów prokariota
```

## Wyniki: wartość minimalna = 0; pierwszy kwantyl = poniżej 25%; mediana = 36; średnia arytmetyczna = 32,02; trzeci kwartyl = poniżej 75%; wartość maksymalna = 393

```{r}
plot(coverage_data[[1]], main="Pokrycie genomu dla sekwencji U00096.3", ylab="Pokrycie", xlab="Pozycja w genomie")
```

## Ta komenda generuje wykres typu Manhatan. Na osi X jest pozycja w genomie, która tutaj jest pokazana w numerycznej skali pozycji genomowej, a na osi Y jest wartość pokrycia genomu, jako liczna odczytów na danej pozycji w genomie. Na wykresie widać, to pokrycie w genomie jest dosyć niskie, głównie wartości oscylują poniżej 100 na osi Y, natomiast są regiony o bardzo wysokim pokryciu, widać w okolicach 4e+06 a 4.5e+06 (oś X). Te niskie wartości mogą być spowodowane np. niedostateczną ilością danych, a wysokie wartości mogą być efektem zduplikowanymi regionami w genomie lub fragmentami wysoce konserwatywnymi. Z wykresu można wyczytać, że pokrycie nie jest równomierne, przez nierównomierne fragmentowanie DNA lub przez różne poziomy amplifikacji PCR.

# Dalej trzeba wykryć warianty

```{r}
pileup_param <- PileupParam(
    distinguish_strands = FALSE,
    distinguish_nucleotides = TRUE,
    min_base_quality = 20
)

pile <- pileup(sorted_bam, scanBamParam = ScanBamParam(), pileupParam = pileup_param)
```

## Trzeba zacząć od zdefiniowana parametrów skanowania i wykonać pileup

## Następnie, trzeba przekonwertować dane pileup do ramki danych z uzgodnieniem nazw sekwencji

```{r}
pile_df<-as.data.frame(pile)
class(pile_df)
pile_df <- pile_df %>%
    mutate(seqnames = as.character(seqnames)) %>%
    mutate(seqnames = ifelse(seqnames == "U00096.3", "NC_000913.3", seqnames))
```

# Dalej trzeba pogrupować dane według pozycji

```{r}
variant_candidates <- pile_df %>%
    group_by(seqnames, pos) %>%
    summarise(
        total = sum(count),
        A = sum(count[nucleotide == "A"]),
        C = sum(count[nucleotide == "C"]),
        G = sum(count[nucleotide == "G"]),
        T = sum(count[nucleotide == "T"]),
        .groups = 'drop'
    ) %>%
    mutate(
        ref = as.character(getSeq(fa, GRanges(seqnames, IRanges(pos, pos))))
    ) %>%
    rowwise() %>%
    mutate(
        # Obliczanie alternatywnych alleli
        alt_alleles = list(setdiff(c("A", "C", "G", "T"), ref)),
        # Liczenie odczytów dla referencyjnego i alternatywnych alleli
        ref_count = sum(c_across(c("A", "C", "G", "T"))[ref]),
        alt_count = sum(c_across(c("A", "C", "G", "T"))[alt_alleles])
    ) %>%
    ungroup() %>%
    # Filtracja na podstawie minimalnej liczby odczytów dla wariantu
    filter(alt_count >= 5) %>%
    # Opcjonalne filtrowanie na podstawie proporcji
    filter((alt_count / total) >= 0.2)

```

## Ta komenda na początku grupuje dane i odczytuje liczbę odczytów dla każdego nukleotydu. Dalej dodaje sekwencję referencyjną. Dalej oblicza allele alternatywne i liczy odczyty dla skewencji referencyjnej i alleli alternatywnych. Na koniec dane są filtrowane - usuwane są ugrupowania i uwzględniane są pozycje z co najmniej 5 odczytami alternatywnych alleli i te, gdzie alternatywne allele stanowią co najmniej 20% wszystkich odczytów

```{r}
head(variant_candidates)
```

# Możemy podejrzeć jak wyświetlają się warianty

# I trzeba przefiltrować i eksportować wyniki do pliku

```{r}
# Filtracja wariantów na podstawie jakości i głębokości pokrycia
filtered_variants <- variant_candidates %>%
    filter(total >= 10, alt_count / total >= 0.2, alt_count >= 5)

# Wyświetlenie liczby wariantów przed i po filtrowaniu
cat("Liczba wariantów przed filtrowaniem:", nrow(variant_candidates), "\n")
cat("Liczba wariantów po filtrowaniu:", nrow(filtered_variants), "\n")

# Konwersja do data.frame dla eksportu
df_variants <- as.data.frame(filtered_variants)

# Eksport do pliku CSV
write.csv(df_variants, "C:/Users/Magda/OneDrive - Szkoła Główna Gospodarstwa Wiejskiego/Pulpit/ABWG-GWAS/wyniki_wariantow.csv", row.names = FALSE)
```

## Na początku uwzględnio tylko warianty, które mają całkowitą liczbę odczytów większą lub równą 10; warianty, w których proporcja odczytów wariantu alternatywnego do całkowitej liczby odczytów wynosi co najmniej 20; i te, które mają co najmniej 5 odczytów alternatywnego wariantu. Później zwracano liczbę wierszy i przekonwertowano do ramki danych, a na koniec wyeksportowano dane do pliku w Excelu.

# Dzięki *Variant Calling* jest możliwe indetyfikacja różnic między genomem analizowanym a genomem referencyjnym, Dzięki temu można np. zidentyfikować mutacje, wykrywać polimorfizmy czy analizować różnorodność genetyczną między populacjami
