---
title: "Zadanie 12"
author: "Magdalena Jackiewicz"
date: "2024-12-17"
output: html_document
---

# Zaczynamy od zainstalowania i załadowania potrzebnych pakietów

```{r}
BiocManager::install("VariantAnnotation")
BiocManager::install("GenomicRanges")
BiocManager::install("AnnotationHub")
```

```{r}
library(VariantAnnotation)
library(GenomicRanges)
library(AnnotationHub)
```

# Następnie trzeba znaleźć przykładowy plik VCF z pakietu VariantAnnotation

```{r}
vcf_file <- system.file("extdata", "chr22.vcf.gz", package = "VariantAnnotation")
vcf_file
```

## Ten przykładowy plik chr22.vcf.gz zawiera dane z fragmentami wariantów genetycznych ma chromosomie 22. Zawiera informacje: *CHROM* - chromosom (w tym przypadku zawsze 22); *POS* - pozycja w chromosomie; *ID* - identyfikator wariantu; *REF* - nukleotyd w sekwencji referencyjnej; *ALT* - nukleotyd lub sekwencja alternatywna; *QUAL* - jakość wariantu (liczba związana z pewnością detekcji); *FILTER* - status filtru; *INFO* - dodatkowe dane o wariancie.

# Można zobaczyć ogólnie jak wyglądają dane w tym pliku

```{r}
vcf_data <- readVcf(vcf_file, genome = "hg19")
show(vcf_data)
```

## Najpierw plik jest wczytywany do obiektu VCF, a dalej można podejrzeć dane

# Dalej trzeba wyświetlić podstawowe informacje o zawartości pliku

```{r}
vcf_header <- header(vcf_data)
vcf_header
```

## Dzięki tej komendzie są wyświetlanie nagłówki i metadane. Wyświetlane naglowki: *class* - dane są przechowywane w obiekcie klasy VCFHeader; *samples* - jest 5 próbek, a każda próbka ma dane genotypowe; *meta* - mówi o specyfikacji formatu; *fixed* - ma 2 kolumny: FILTER - zawiera informacje o filtrach zastosowanych do każdego wariantu i ALT - inaczej alternatywna sekwencja nukleotydów lub allel w pozycji wariantu; *info* - ma 22 różne podkategorie opisujące szczegóły każdego wariantu; *geno* - opisuje format danych genotypowych dla próbek: GT - genotyp dla każdej próbki, DS - liczba kopii allelu alternatywnego i GL -logarytmiczne prawdopodobieństwa dla genotypów

# W kolejnym kroku trzeba sprawdzić kolumny INFO i FORMAT

```{r}
info(vcf_data)
geno(vcf_data)
```

## W pliku jest 10376 wierszy o informajcach wariantów genetycznych i 22 kolumny o różnych podkategoriach opisujących warianty w kolumnie. A dane genotypowe są zorganizowane w 3 pola: GT, DS, GL.

# Dalej trzeba policzyć liczbę wariantów

```{r}
num_variants <- length(rowRanges(vcf_data))
num_variants
```

## Liczba wariantów wynosi 10376

# Następnie trzeba przeprowadzić filtrację i analizę jakości.Trzeba zacząć od sprawdzenia kolumny QUALw pliku VCF. Ale przed tym musiałam zainstalowaći załadować kolejny potrzebny pakiet

```{r}
install.packages("vcfR")
library(vcfR)
```

```{r}
vcf_file <- read.vcfR("C:/Users/Magda/AppData/Local/R/win-library/4.4/VariantAnnotation/extdata/chr22.vcf.gz")
str(vcf_file) # Sprawdzenie struktury badanego pliku

qual_values <- vcf_file@fix[, "QUAL"] # Indeksowanie po nazwie kolumny
head(qual_values) # Wyświetlenie pierwszych kilku wartości
names(vcf_file@fix) # Sprawdzenie jak kolumny są nazwane
summary(qual_values) # Podsumowanie
```

## Wynik komendy o sprawdzenie nazw kolumn daje wynik NULL, to oznacza, że dany plik nie jest tym typem co obsługuje funckję "names()". Jeśli chodzi o podsumowanie, to wychodzi, że w wektorze jest 10376 wartości, klasa obiektu to "character", co oznacza, że wartości w qual_values są zapisane jako ciągi znaków, a mode obiektu to także "character", co oznacza, że dane są przechowywane w formie tekstu. Dane liczbowe można przekonwertować z formy tekstu na liczbowe.

# Dalej należy ustalić kryterium jakości i odfiltrować warianty

```{r}
qual_values_numeric <- as.numeric(qual_values) # Konwersja danych tekstowych na dane numeryczne
quality_threshold <- 70 # Ustalenie progu jakości, tutaj wynosi 70
filtered_vcf <- vcf_file[qual_values_numeric > quality_threshold, ] # Filtracja wariantów, które mają QUAL > 70
head(filtered_vcf) #Sprawdzenie wyników
```

## Pokazują się dane liczbowe pierwszych sześciu kolumn.

# Można jeszcze zobaczyć ile wariantów zostało odfiltrowanych

```{r}
n_variants_original <- nrow(vcf_file)  # Liczba wariantów przed filtracją
n_variants_filtered <- nrow(filtered_vcf)  # Liczba wariantów po filtracji
cat("Liczba wariantów przed filtracją:", n_variants_original, "\n")
cat("Liczba wariantów po filtracji:", n_variants_filtered, "\n")
```

## Wynik:

### Liczba wariantów przed filtracją: 10376

### Liczba wariantów po filtracji: 10340

# Po filtracji można zapisać nowe dane

```{r}
write.vcf(filtered_vcf, file = "filtered_vcf_file.vcf")
```

# Dalej trzeba wczytać plik, ten z przefiltrowanymi wartościami i przekształcić dane na GRanges

```{r}
vcf_file <- read.vcfR("filtered_vcf_file.vcf")  

chromosomes <- vcf_file@fix[, "CHROM"]  # Kolumna CHROM
positions <- as.integer(vcf_file@fix[, "POS"]) # Kolumna POS (pozycje startowe wariantów)

 valid_rows <- !is.na(positions)
  chromosomes <- chromosomes[valid_rows]
  positions <- positions[valid_rows]

gr_vcf <- GRanges(seqnames = chromosomes, 
                  ranges = IRanges(start = positions, end = positions))
head(gr_vcf)
```

## Ze względu na to, że mogą być wiersze z brakującymi wartościami, usunełam je i po tym stworzono obiekt GRanges, który od razu można podejrzeć. Widać pikerwsze sześć linijek: lokalizacja wariantów - wszystko na chromosmie 22; pozycje tych wariantów - są w bliskiej odległości od siebie; strand - szystkie warianty są na neutralnym strandzie. Brakuje innych metadanych i informacji o długości sekwencji.

# Dodaje brakujące metadane

```{r}
mcols(gr_vcf) <- data.frame(ID = c("var1", "var2", "var3", "var4", "var5", "var6"))
gr_vcf
```

# Następnie trzeba przeprowadzić anotację wariantów, wcześniej należy stworzyć obiekt typu GRanges, który jest używany do przechowywania informacji o regionach genomu

```{r}
gr_vcf <- GRanges(seqnames = chromosomes, ranges = IRanges(start = positions, end = positions))

```

# Można wzbogacić warianty o inne informacje genowe przy wykorzystaniu innych pakietów

```{r}

variants <- rowRanges("filtered_vcf_file.vcf")
annotation_results <- locateVariants(query = variants, subject = txdb, region = AllVariants())

head(annotation_results)
```

# Następnie można wyeksportować wyniki do pliku

```{r}
write.csv(as.data.frame(annotation_results), "annotated_variants.csv", row.names = FALSE)
```

# Kolejnym krokiem jest wybór wariantów z regionów UTR

```{r}
utr5_variants <- locateVariants(query = rowRanges(vcf_data), subject = txdb, region = FiveUTRVariants()) # Warianty w 5'UTR

utr3_variants <- locateVariants(query = rowRanges(vcf_data), subject = txdb, region = ThreeUTRVariants()) # Warianty w 3'UTR

cat("Liczba wariantów w 5'UTR:", num_utr5, "\n")
cat("Liczba wariantów w 3'UTR:", num_utr3, "\n")
```

## Jest znacząca różnica w liczbie wariantów na poszczególnych końcach

# Następnie należy wyodrębnić warianty znajdujące się w regionach międzygenowych

```{r}
intergenic_variants <- locateVariants(query = rowRanges(vcf_data), subject = txdb, region = IntergenicVariants())

num_intergenic <- nrow(intergenic_variants)

cat("Liczba wariantów międzygenowych:", num_intergenic, "\n")
```

## Wynik wynosi 3028 wariantów między genowych

# Podsumowując po przefiltrowaniu pliku nadal zostało dużo wariantów, mimo całkiem wysokiego progu. Z wyższym progiem, liczba wariantów mocno malała, dlatego zostałam przy tym, ponieważ wtedy analiza mogłaby nie wyjść dostatecznie zgodna. W pliku było wiele braków potrzebnych danych, dzięki którym analiza mogła zostać dobrze przeprowadzona, tyle ile mogłam to dodałam brakujących informacji. Jest duża różnica między wariantami 5' i 3' UTR a wariantami międzygenowymi. Może to być spowodowane tym, że regiony 5' i 3' UTR są bardziej narażone na mutacje i są bardziej zróżnicowane niż regiony międzygenowe.
