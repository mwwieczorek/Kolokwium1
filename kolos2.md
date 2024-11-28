---
title: "Kolokwium2"
output: html_document
date: "2024-11-26"
editor_options: 
  markdown: 
    wrap: 72
---

### Pobieranie pakietów potrzebnych do obrbki sekwencji, mapowania genomu oraz instalacja BiocManagera

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE)) 
    install.packages("BiocManager")

#wczytanie BiocManagera
library(BiocManager)
```

### Pobieranie za pomocą Biocmanagera następujących pakietów

Instalacja wszystkich pakietów jednocześnie w celu ułatwienia pracy :)

```{r}
BiocManager::install(c("GenomicFeatures", "AnnotationDbi", 
"ShortRead", "Biostrings","ggplot2", "Rsubread", "Rqc","GenomicAlignments"))
```

### Wczytywanie wszystkich pobranych pakietów

```{r}
library(GenomicFeatures) 
library(AnnotationDbi)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(Rsubread)
library (Rqc)
library(GenomicAlignments)
```

### Zimportowanie pliku FastQ do R i wczytanie odczytów

```{r}
fastq_file <- "C:/Users/Sensumart/Desktop/kolokwium2/ecoli_raw1.fq"

fq_reads <- readFastq(fastq_file)

#wyciągnięcie podstawowych infromacji - liczba odczytów w sekwencji 
length(fq_reads)  

```

**Wynik:309440** –\> tyle jest odczytów w sekwencji z pliku Fastq
ecoli_raw1

### Obliczanie zawartości GC w celu zbadania czy nie występuję wysoka zawartość GC i biasy GC

```{r}
#obliczanie zawrtości GC dla surowego odczytu
gc_content <- letterFrequency(sread(fq_reads), letters = "GC", as.prob = TRUE)
#histogram zawartości GC 
hist(gc_content, breaks = 50, main = "Zawartość GC w oryginalnych odczytach", xlab = "Procent GC")
```

### Długość sekwencji w nukleotydach - aby zobrazować jak się zmienia długość sekwencji przed i po obróce danych do mapowania

```{r}

# Pobranie długości poszczególnych sekwencji
sequence_lengths <- width(sread(fq_reads))

# Oblicz całkowitą długość sekwencji
total_length <- sum(sequence_lengths)

# Wyświetlenie wyniku
cat("Całkowita długość sekwencji w pliku FASTQ:", total_length, "nukleotydów\n")
```

**Wynik:**

```         
Całkowita długość sekwencji w pliku ecoli_raw1: 47223547 nukleotydów
```

### Wykonanie raportu QC dla sekwencji ecoli_raw1 w celu kontroli jakości przed obróbką - w celu znalezienia problemów z jakością odczytów, które mogą wpłynąć na mapowanie sekwencji

```{r}
# Kontrola jakości - generowanie obiektu z wynikami kontroli jakosci - raport QC za pomoca pakietu ShortRead
qa_results <- qa(fastq_file, type = "fastq")
#generowanie raportu QC
report(qa_results, dest = "C:/Users/Sensumart/Desktop/kolokwium2/QA_report1")

```

# Triming&Filtering

## **1.Trimming: przycinanie konkretnych nt o niskiej jakości Phred**

Cel: usunięcie baz o niskiej jakości mogących stanowić przeszkode w
mapowaniu sekwencji do genomu referencyjnego

```{r}
#Przycinanie odczytów o nieskiej jakości na końcach odczytów - zadane parametry funkcji 
#a ---> Phred+33+ kodowanie - wybieram próg Phred30 (?)
#k = 3 tzn.funkcja będzie szukać co najmniej trzech kolejnych baz o jakości poniżej progu Phred30 aby rozpocząć przycinanie 
#halfwidth ustawiam na szerokość 5 baz - funkcja obliczać będzie średnią wartość jakości w oknie o szerokości 5 baz

trimmed_reads <- trimTailw(fq_reads, k = 3, a = "?", halfwidth = 2)

# Liczba odczyhtów
length(fq_reads)        # Przed przycinaniem[1]
length(trimmed_reads)   # Po przycinaniu[2]

#liczba zmodyfikowanych odczytów[3]
sum(width(trimmed_reads) < width(fq_reads))

# procent sekwencji, który nie uległ zmodyfikowaniu[4]
(length(trimmed_reads)/length(fq_reads))* 100

# procent sekwencji, który uległ zmodyfikowaniu[5]
100-((length(trimmed_reads)/length(fq_reads))* 100)


```

**Wynik**i:

[1] 309440

[2] 298001

[3] 193342

[4] 96.30332

[5] 3.696678

## **2.Filtering:**

## **Cel: usuwanie całych odczytów o niskiej jakości lub zbyt krótkich**

```{r}
#60bp - minimalna długość odczytów 
# Filtrowanie odczytów 
filtered_reads <- trimmed_reads[width(trimmed_reads) >= 60]

#Sprawdzanie długośći odczytów po filtracji 
length(trimmed_reads)       # Po przycinaniu[1]
length(filtered_reads)      # Po filtracji[2]


(length(filtered_reads)/length(trimmed_reads))*100 #procent odczytów pozostały po filtracji  [3]

100 - ((length(filtered_reads)/length(trimmed_reads))*100)#procent odczytów odrzuconych podczas filtracji [4] 




```

**Wyniki**:

```         
[1] 298001
[2] 266966
[3] 89.58561
[4] 10.41439
Podczas filtracji odrzucono 10,41 % odczytów 
```

### Kontrola jakości po przycinani i filtrowaniu

Cel: zbdadanie wpływu obróbki sekwencji na jakość odczytów

```{r}
#zapisanie pliku po modyfikacji jako plik fastq
writeFastq(filtered_reads, "C:/Users/Sensumart/Desktop/kolokwium2/ecoli_raw2.fq", "a")

qa_results_2 <- qa("C:/Users/Sensumart/Desktop/kolokwium2/ecoli_raw2.fq", type = "fastq")
#generowanie raportu QC
report(qa_results_2, dest = "C:/Users/Sensumart/Desktop/kolokwium2/QA_report2")



```

### Analiza zawartości GC po modyfikacjach

**Cel: Zbadanie jak obróbka sekwencji wpłynęła na zawartość GC**

```{r}
#obliczanie zawrtości GC dla zmodyfikowanych odczytów
gc_content <- letterFrequency(sread(filtered_reads), letters = "GC", as.prob = TRUE)
#histogram zawartości GC 
hist(gc_content, breaks = 50, main = "Zawartość GC w zmodyfikowanych odczytach", xlab = "Procent GC")
```

### Analiza rozkładu długości odczytów

**Cel: Zbadanie, jak przycinanie wpłynęło na długość odczytów,
generowanie histogramów**

```{r}

# Przed przycinaniem 
hist(width(fq_reads), breaks = 50, main = "Długość odczytów przed przycinaniem", xlab = "Długość (bp)")
   
# Po przycinaniu 
hist(width(filtered_reads), breaks = 50, main = "Długość odczytów po przycinaniu", xlab = "Długość (bp)")
   
   
```

### Wykrywanie i usuwanie sekwencji adapterów

**Cel: Zidentyfikowanie obecności sekwencji adapterów i usunięcie ich z
odczytów.**

```{r}
#Zdefiniowanie sekwencji adaptera 

adapter_seq <- DNAString("AGATCGGAAGAGC")
```

```{r}
# Przycinanie adapterów z odczytów :
trimmed_reads <- trimLRPatterns(
  Lpattern = adapter_seq,
  subject = filtered_reads
)

# Defuniujemy odczyty po przycięciu adapterów:
filtered_reads <- trimmed_reads

```

### Sprawdzanie efektów przycinania

```{r}
# Porównanie długości przed i po przycięciu adapterów
length(filtered_reads) #przed przycieciem adapterów[1]
length(trimmed_reads) #po przycięciu adapterów[2]


# ile odczytów zostało zmodyfikowanych
   sum(width(filtered_reads) < width(trimmed_reads)) #[3]
   
```

**Wyniki**:

[1] 266966

[2] 266966

[3] 0 - liczba zmodyfikowanych odczytów równa 0 sugeruje, że adapterów
nie było.

```{r}

# Pobranie długości poszczególnych sekwencji
sequence_lengths <- width(sread(filtered_reads))

# Oblicz całkowitą długość sekwencji
total_length <- sum(sequence_lengths)

# Wyświetlenie wyniku
cat("Całkowita długość sekwencji w pliku po przycinaniu adapterów:", total_length, "nukleotydów\n")
```

**Wynik:**

```         
Całkowita długość sekwencji w pliku po przycinaniu adapterów: 35300224 nukleotydów
```

### Wpływ obróbki sekwencji na jej długość w nukleotydach

```{r}
47223547 - 35300224
```

**Wynik:** całkowita długość sekwencji surowej po modyfikacjach uległa
skróceniu o 11923323 nukleotydy

### Finalna konrola jakości

**Cel: Zbadanie wpływu usuwania sekwencji adapterów na jakość odczytów**

```{r}
#zapisanie efektu jak nowy plik fastq - po przycinaniu adapterów
writeFastq(filtered_reads, "C:/Users/Sensumart/Desktop/kolokwium2/ecoli_raw3.fq","a")

qa_results3 <- qa("C:/Users/Sensumart/Desktop/kolokwium2/ecoli_raw3.fq", type = "fastq")
report(qa_results3, dest = "C:/Users/Sensumart/Desktop/kolokwium2/QA_report3")

```

## ALIGMENT

```{r}
#importowanie genomu referencyjnefgo 
ref_genome <- readDNAStringSet("C:/Users/Sensumart/Desktop/kolokwium2/ecoli_genome.fna.gz")

#Budowanie indeksu genomu

buildindex(basename = "ecoli_index", reference = "C:/Users/Sensumart/Desktop/kolokwium2/ecoli_genome.fna.gz")


```

**Wszystko zostało zapisane w katalogu:**

```         
Index C:\Users\Sensumart\Downloads\ecoli_index was successfully built.
```

### Mapowanie odczytów do genomu referencyjnego

```{r}
#mapowanie odczytów i utworzenie pliku BAM
align(index = "ecoli_index",
      readfile1 = "C:/Users/Sensumart/Desktop/kolokwium2/ecoli_raw3.fq",
      input_format = "FASTQ",
      output_file = "C:/Users/Sensumart/Desktop/kolokwium2/aligned_sample.BAM")
```

**Wynik:**

Total_reads 266966\
Mapped_reads 266918\
Uniquely_mapped_reads 262013\
Multi_mapping_reads 4905\
Unmapped_reads 48\
Indels 129

**Procent odczytów zmapowanych:**

```{r}
(266918/266966)*100 
```

**Wynik**: 99.98

**Procent odczytów niezmapowanych:**

```{r}
(48/266966)*100
```

**Wynik**: 0.018

### Możliwe przyczyny niezmapowania odczytów:

-   insercje i delecje (określone jako Indels) lub SNP (polimorfizmy
    pojedynczych nukleotydów)

### Analiza wyników mapowania

```{r}
#Import pliku BAM - zmapowanych odczytów do R  
aln <- readGAlignments("C:/Users/Sensumart/Desktop/kolokwium2/aligned_sample.BAM")
```

### Oblicznie pokrycia genomu

```{r}
coverage_data <- coverage(aln)
```

### Wizualizacja pokrycia genomu za pomocą pakietu ggplot2

**Cel: identyfikacja regionó o najwyższym pokryciu**

```{r}

# Konwersja pokrycia do data frame
cov_df <- as.data.frame(coverage_data[[1]])
cov_df$position <- as.numeric(rownames(cov_df))

# Wykres pokrycia
#w celu wizualizacji regionów o najwyżsyzm pokryciu 
cov_df <- as.data.frame(coverage_data[[1]])
cov_df$position <- as.numeric(rownames(cov_df))

ggplot(cov_df[1:25000, ], aes(x = position, y = value)) +
  geom_line(color = "red") +
  labs(title = "Pokrycie genomu E. coli",
       x = "Pozycja w genomie",
       y = "Liczba zmapowanych odczytów")
```

### Pokrycie średnie genomu

```{r}
#Obliczanie pokrycia genomu  
#średnie pokrycie = suma długości wszystkich odczytów/przez długość genomu referencyjnego  
mean(cov_df$value)   
#średnie pokrycie genomu = ile razy dany fragment został odczytany podczas sekwencjonowania
```

**Wynik**: 7.602825 --\> pokrycie genomu poniżej 10 mówi nam,że jest to
bardzo niskie pokrycie, które może prowadzić do pominięcia rzadkich
mutacji i błędów w analizie

### Generowanie zbiorczego raportu dla 3 próbek

**Cel: Analzia wpływu modyfikacji sekwencji na jakość odczytów**

-   **ecoli_raw1** (surowa sekwencja)
-   **ecoli_raw2** (sekwencja po trymowaniu i filtrowaniu)
-   **ecoli_raw3** (sekwencja po przycięciu adaterów)

```{r}
#Generowanie zbiorczego raportu 
fq_files <- list.files(path = "kolokwium2", pattern = "ecoli_raw[1-3].fq", full.names = TRUE)
qa_results4 <- qa(fq_files, type = "fastq") report(qa_results4, dest =
"QA_report_multi")
```
