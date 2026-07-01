# Materiały i metody — analiza wariantów mitochondrialnych na poziomie pojedynczych komórek

## Spis treści

- [Cel analizy](#cel-analizy)
- [Dane wejściowe](#dane-wejściowe)
- [Step 1 — filtrowanie jednoznacznie mapujących się odczytów mitochondrialnych](#step-1--filtrowanie-jednoznacznie-mapujących-się-odczytów-mitochondrialnych)
- [Step 2 — liczenie unikalnych odczytów mtDNA na komórkę](#step-2--liczenie-unikalnych-odczytów-mtdna-na-komórkę)
- [Step 3 — wybór komórek z wystarczającą liczbą odczytów mtDNA](#step-3--wybór-komórek-z-wystarczającą-liczbą-odczytów-mtdna)
- [Step 4 — zliczanie alleli A/C/G/T per komórka i pozycja mtDNA](#step-4--zliczanie-alleli-acgt-per-komórka-i-pozycja-mtdna)
- [Step 5 — empiryczne modelowanie rozkładu jakości baz](#step-5--empiryczne-modelowanie-rozkładu-jakości-baz)
- [Step 6 — wybór progu jakości baz](#step-6--wybór-progu-jakości-baz)
- [Step 7 — filtrowanie pozycji mtDNA po jakości baz](#step-7--filtrowanie-pozycji-mtdna-po-jakości-baz)
- [Step 8 — wywoływanie kandydackich alleli alternatywnych](#step-8--wywoływanie-kandydackich-alleli-alternatywnych)
- [Step 9 — filtrowanie kandydackich wariantów mtDNA](#step-9--filtrowanie-kandydackich-wariantów-mtdna)
- [Step 12 — klasyczne klastrowanie mtDNA](#step-12--klasyczne-klastrowanie-mtdna-metodą-hierarchical-clustering)
- [Step 121 / Step 12B — fuzzy k-medoids](#step-121--step-12b--fuzzy-k-medoids-na-podstawie-dystansu-mtdna)
- [Step 15 — integracja klasycznych klastrów mtDNA z Seurat](#step-15--integracja-klasycznych-klastrów-mtdna-z-obiektem-seurat)
- [Step 151 / Step 15B — integracja fuzzy klastrów mtDNA z Seurat](#step-151--step-15b--integracja-fuzzy-klastrów-mtdna-z-seurat)
- [Podsumowanie przepływu danych](#podsumowanie-przepływu-danych)
- [Najważniejsze założenia analityczne](#najważniejsze-założenia-analityczne)

---

## Cel analizy

Celem pipeline’u było wykrycie informatywnych wariantów mitochondrialnych w danych single-cell oraz wykorzystanie ich do grupowania komórek na podstawie podobieństwa profili heteroplazmii mtDNA. Pipeline wykonywał kolejno: selekcję jednoznacznie mapujących się odczytów mitochondrialnych, zliczanie odczytów mtDNA na komórkę, wybór komórek z wystarczającą liczbą odczytów mitochondrialnych, zliczanie alleli A/C/G/T per komórka i pozycja, empiryczną ocenę jakości baz, filtrowanie pozycji o niskiej jakości, wywoływanie kandydackich alleli alternatywnych, filtrowanie wariantów oraz klastrowanie komórek metodą klasyczną hierarchical clustering albo metodą fuzzy k-medoids.

Pipeline został stworzony na podstawie metody opisanej w publikacji: https://www.sciencedirect.com/science/article/pii/S0092867419300558?via%3Dihub.

Pipeline był uruchamiany w środowisku `~/venvs/scmt_env/bin/activate`, używał programu `samtools` do indeksowania BAM oraz skryptów Python i R z katalogów `BAMs/scripts` i `BAMs/scripts_new`. Domyślne parametry w analizowanym pliku pipeline’u obejmowały m.in. minimalną liczbę unikalnych odczytów mtDNA na komórkę `MIN_MT_READS=200`, minimalną jakość mapowania `MAPQ=30`, nazwę chromosomu mitochondrialnego `CHR_M="chrM"` oraz próg pokrycia do obliczeń dystansu mitochondrialnego `PUBLICATION_COVERAGE_THRESHOLD=9`, co odpowiada regule `depth > 9`, czyli co najmniej 10 odczytów.

---

## Dane wejściowe

Dla każdej próbki pipeline oczekiwał jako głównych danych wejściowych:

1. pliku BAM:

   ```text
   BAMs/all_samples_a1m95/<sample>.bam
   ```

2. pliku z kodami kreskowymi komórek, wyszukiwanego automatycznie w katalogu BAM, np.:

   ```text
   barcodes.<sample>.tsv.gz
   <sample>.barcodes.tsv.gz
   barcodes.tsv.gz
   filtered_feature_bc_matrix/barcodes.tsv.gz
   raw_feature_bc_matrix/barcodes.tsv.gz
   ```

3. tabeli metadanych Seurat zawierającej co najmniej kolumny `sample` i `barcode`:

   ```text
   ~/KSzade/sctanno_ksz/sctanno_ksz_all_cells_with_sample.tsv
   ```

4. obiektu Seurat RDS używanego później do wizualizacji i kontroli jakości klastrów:

   ```text
   ~/KSzade/sctanno_ksz/rzeszow/SCT_annotated_KSz.Rds
   ```

Pipeline najpierw przygotowywał listę barcode’ów zgodnych między plikiem z barcode’ami Cell Ranger i metadanymi Seurat. Z metadanych Seurat wybierane były barcode’y należące do danej próbki, a następnie przecinane z barcode’ami dostępnymi w danych wejściowych. W ten sposób powstawał plik:

```text
<sample>.barcodes.from_seurat.tsv
```

Ten plik był dalej używany jako lista komórek dopuszczonych do zliczania odczytów mtDNA.

---

## Step 1 — filtrowanie jednoznacznie mapujących się odczytów mitochondrialnych

**Skrypt:**

```text
filter_unique_mt_reads.py
```

**Wejście:**

```text
BAMs/all_samples_a1m95/<sample>.bam
```

**Wyjście:**

```text
<sample>.uniq_mt.bam
<sample>.uniq_mt.bam.bai
<sample>.filter_unique_mt_reads.log
```

W tym kroku z pliku BAM wybierane były wyłącznie odczyty spełniające jednocześnie następujące warunki:

```text
read nie jest unmapped
read nie jest secondary
read nie jest supplementary
read.reference_name == chrM
read.mapping_quality >= 30
read ma tag NH
NH == 1
```

Warunek `NH == 1` oznacza, że odczyt został uznany za jednoznacznie mapujący się według alignera. Skrypt zapisuje tylko takie odczyty do nowego pliku BAM. Następnie pipeline indeksuje wynikowy BAM przy użyciu `samtools index`. Log zawiera liczbę wszystkich odczytów widzianych przez skrypt oraz liczbę odczytów zapisanych jako unikalne odczyty mitochondrialne.

---

## Step 2 — liczenie unikalnych odczytów mtDNA na komórkę

**Skrypt:**

```text
count_unique_mt_reads_per_cell.py
```

**Wejście:**

```text
<sample>.uniq_mt.bam
<sample>.barcodes.from_seurat.tsv
```

**Wyjście:**

```text
<sample>.uniq_mt_reads_per_cell.tsv
<sample>.count_unique_mt_reads_per_cell.log
```

Skrypt przechodzi przez odczyty z pliku `uniq_mt.bam` i zlicza liczbę odczytów mitochondrialnych dla każdej komórki z listy barcode’ów. Uwzględniane są tylko odczyty:

```text
mapped
primary alignment
z tagiem CB
z barcode’em obecnym na liście valid barcodes
```

Wynikiem jest tabela:

```text
barcode    unique_mt_reads
```

Każda komórka z listy wejściowej otrzymuje wartość `unique_mt_reads`; jeśli nie miała żadnych odczytów mtDNA, dostaje wartość 0.

---

## Step 3 — wybór komórek z wystarczającą liczbą odczytów mtDNA

**Skrypt:**

```text
select_barcodes_by_mt_reads.py
```

**Wejście:**

```text
<sample>.uniq_mt_reads_per_cell.tsv
```

**Parametr:**

```text
MIN_MT_READS=200
```

**Wyjście:**

```text
<sample>.mt200_barcodes.tsv
<sample>.select_barcodes_by_mt_reads.log
```

W tym kroku wybierane są tylko te komórki, dla których:

```text
unique_mt_reads >= 200
```

Wynikowy plik zawiera jedną kolumnę bez nagłówka — listę barcode’ów komórek przechodzących próg minimalnej liczby unikalnych odczytów mitochondrialnych. Log zapisuje liczbę barcode’ów przed i po filtrowaniu.

---

## Step 4 — zliczanie alleli A/C/G/T per komórka i pozycja mtDNA

**Skrypt:**

```text
per_cell_per_allele_counts1.py
```

**Wejście:**

```text
<sample>.uniq_mt.bam
<sample>.mt200_barcodes.tsv
```

**Wyjście:**

```text
<sample>.per_cell_per_allele_counts1.tsv.gz
<sample>.per_cell_per_allele_counts1.base_quality_distribution.png
<sample>.per_cell_per_allele_counts1.log
```

Skrypt wykonuje pileup po pliku `uniq_mt.bam` i dla każdej pozycji mitochondrialnej zlicza, ile razy w danej komórce obserwowano każdą z zasad:

```text
A, C, G, T
```

Nie stosuje się tutaj filtrowania po jakości bazy, ponieważ `min_base_quality=0`. Jest to ważne, ponieważ rozkład jakości baz ma być oceniony empirycznie w kolejnych krokach, a nie narzucony na tym etapie. Skrypt ignoruje delecje, refskipy, odczyty unmapped, secondary, supplementary, odczyty bez tagu `CB`, barcode’y spoza listy komórek wybranych w Step 3 oraz znaki inne niż A/C/G/T.

Dla każdej kombinacji:

```text
barcode × pozycja mtDNA
```

zapisywane są:

```text
barcode
pos
A
C
G
T
depth
A_mean_bq
C_mean_bq
G_mean_bq
T_mean_bq
```

gdzie `depth = A + C + G + T`, a `*_mean_bq` oznacza średnią jakość baz dla danego allelu w danej komórce i pozycji. Jeśli dany allel nie był obserwowany, wartość jakości zapisywana jest jako `NA`.

Dodatkowo skrypt tworzy wykres globalnego rozkładu jakości baz:

```text
<sample>.per_cell_per_allele_counts1.base_quality_distribution.png
```

---

## Step 5 — empiryczne modelowanie rozkładu jakości baz

**Skrypt:**

```text
gaussian_plot.py
```

**Wejście:**

```text
<sample>.per_cell_per_allele_counts1.tsv.gz
```

**Wyjście:**

```text
<sample>.GMM_meanBQ_density.png
<sample>.GMM_meanBQ_QC.tsv.gz
<sample>.gaussian_plot.log
```

Skrypt agreguje dane z poziomu `komórka × pozycja × allel` do poziomu:

```text
pozycja × allel
```

Dla każdego allelu w każdej pozycji obliczana jest średnia ważona jakość bazy:

```text
mean_bq = sum(count_in_cell × mean_bq_in_cell) / total_count
```

Następnie do rozkładu wartości `mean_bq` dopasowywany jest 3-składnikowy Gaussian Mixture Model. Składnik o najwyższej średniej jest interpretowany jako komponent wysokiej jakości. Dla każdej pozycji i allelu obliczane jest prawdopodobieństwo przynależności do tego komponentu wysokiej jakości, a wynik zapisywany jest w kolumnie:

```text
p_high_conf
```

Jeśli:

```text
p_high_conf > 0.99
```

to pozycja-allel otrzymuje status:

```text
pass_99pct = YES
```

Plik QC ma kolumny:

```text
pos
allele
n_reads
mean_bq
p_high_conf
pass_99pct
```

Dodatkowo tworzony jest wykres gęstości rozkładu `mean_bq` z zaznaczonymi trzema komponentami GMM oraz empirycznym progiem 99% dla komponentu wysokiej jakości.

---

## Step 6 — wybór progu jakości baz

**Skrypt:**

```text
brak osobnego skryptu; logika awk/zcat w pipeline
```

**Wejście:**

```text
<sample>.GMM_meanBQ_QC.tsv.gz
```

**Wyjście:**

```text
<sample>.BQth.txt
<sample>.use_BQth.txt
```

Pipeline odczytuje plik QC z kroku 5 i wybiera minimalną wartość `mean_bq` spośród wierszy, które przeszły próg:

```text
pass_99pct == YES
```

Inaczej mówiąc, wybierany jest najniższy średni BQ, który nadal należy z prawdopodobieństwem >99% do komponentu wysokiej jakości. Następnie wartość ta jest zapisywana jako surowy próg w pliku:

```text
<sample>.BQth.txt
```

oraz zaokrąglana w dół do liczby całkowitej i zapisywana jako:

```text
<sample>.use_BQth.txt
```

Ten całkowity próg jest używany w następnym kroku do filtrowania pozycji o niskiej jakości.

---

## Step 7 — filtrowanie pozycji mtDNA po jakości baz

**Skrypt:**

```text
filter_per_cell_allele_counts_by_bq_plus.py
```

**Wejście:**

```text
<sample>.per_cell_per_allele_counts1.tsv.gz
<sample>.use_BQth.txt
```

**Parametr:**

```text
--low-impact-af 0.002
```

czyli:

```text
1 / 500 = 0.002
```

**Wyjście:**

```text
<sample>.filtered_plus.BQfiltered.per_cell_alleles.tsv.gz
<sample>.filtered_plus.position_allele_BQ_QC.tsv.gz
<sample>.filtered_plus.bad_positions.txt
<sample>.filtered_plus.good_positions.txt
<sample>.filtered_plus.BQ_filter_summary.txt
<sample>.filter_per_cell_allele_counts_by_bq_plus.log
```

Skrypt agreguje dane per pozycja i allel. Dla każdej pozycji i allelu liczone są:

```text
total_count
position_depth
mean_bq
allele_fraction = total_count / position_depth
```

Pozycja jest uznawana za problematyczną i usuwana, jeśli zawiera co najmniej jeden allel spełniający warunek:

```text
mean_bq < bq_threshold
oraz
allele_fraction >= 0.002
```

Zastosowana jest więc reguła wyjątku dla alleli o bardzo małym wpływie na heteroplazmię. Allele o niskiej jakości nie powodują usunięcia pozycji, jeśli ich częstość jest mniejsza niż `1/500`, czyli `0.002`.
Filtrowanie jest konserwatywne na poziomie pozycji: jeśli pozycja zostanie oznaczona jako zła, usuwane są wszystkie wiersze dotyczące tej pozycji ze wszystkich komórek. Wynikowa tabela zachowuje oryginalny format:

```text
barcode
pos
A
C
G
T
depth
A_mean_bq
C_mean_bq
G_mean_bq
T_mean_bq
```

ale zawiera tylko pozycje, które przeszły filtr jakości.

---

## Step 8 — wywoływanie kandydackich alleli alternatywnych

**Skrypt:**

```text
call_candidate_alleles.py
```

**Wejście:**

```text
<sample>.filtered_plus.BQfiltered.per_cell_alleles.tsv.gz
```

**Wyjście:**

```text
<sample>.candidate_alleles.tsv
<sample>.candidate_alleles_by_cell.long.tsv.gz
<sample>.call_candidate_alleles.log
```

Skrypt sumuje zliczenia A/C/G/T po wszystkich komórkach dla każdej pozycji. Allel o największej liczbie odczytów w danej pozycji jest traktowany jako allel „ref-like”, czyli dominujący allel referencyjno-podobny. Wszystkie pozostałe allele z liczbą odczytów większą od zera są traktowane jako potencjalne allele alternatywne.

Dla każdego wariantu tworzona jest nazwa:

```text
<pos>_<ref_like>><alt>
```

np.:

```text
1234_A>G
```

Plik `candidate_alleles.tsv` zawiera podsumowanie wariantu na poziomie populacji komórek:

```text
variant
pos
ref_like
alt
total_depth
total_alt_count
mean_heteroplasmy_population
n_cells_alt_positive
mean_alt_bq_observed
```

gdzie:

```text
mean_heteroplasmy_population = suma alt_count / suma depth
```

Plik `candidate_alleles_by_cell.long.tsv.gz` zawiera długi format danych per komórka i wariant:

```text
barcode
variant
pos
ref_like
alt
depth
alt_count
vaf
alt_mean_bq
```

gdzie:

```text
vaf = alt_count / depth
```

Ten długi format jest później głównym wejściem do filtrowania wariantów i klastrowania komórek.

---

## Step 9 — filtrowanie kandydackich wariantów mtDNA

**Skrypt:**

```text
filter_variants_no_bulk_after_bq_plus.py
```

**Wejście:**

```text
<sample>.candidate_alleles.tsv
<sample>.candidate_alleles_by_cell.long.tsv.gz
```

**Parametry używane przez pipeline:**

```text
--min-heteroplasmy 0.025
--min-cells-alt-positive 1
--min-total-alt-count 1
```

**Wyjście:**

```text
<sample>.filtered_variants.no_bulk.tsv
<sample>.filtered_variants.no_bulk.QC_all_candidates.tsv
<sample>.filtered_variants_by_cell.long.tsv.gz
<sample>.filter_variants_no_bulk_after_bq_plus.log
```

Wariant przechodzi filtr, jeśli spełnia wszystkie warunki:

```text
mean_heteroplasmy_population >= 0.025
n_cells_alt_positive >= 1
total_alt_count >= 1
total_depth >= 1
```

Skrypt tworzy również tabelę diagnostyczną dla wszystkich kandydatów, w której zapisuje osobne kolumny logiczne informujące, czy dany wariant przeszedł każdy z filtrów:

```text
pass_min_heteroplasmy
pass_min_cells_alt_positive
pass_min_total_alt_count
pass_min_total_depth
pass_final
```

Do dalszych etapów przekazywane są tylko warianty z `pass_final == TRUE`. Długi plik per komórka jest filtrowany do zestawu wariantów, które przeszły filtr na poziomie populacji komórek.

---

## Step 10

W analizowanym pipeline’ie nie ma aktywnego kroku oznaczonego jako Step 10.

---

## Step 11 — krok wyłączony

Step 11, opisany historycznie jako tworzenie macierzy heteroplazmii, jest w tej wersji pipeline’u wyłączony. Pipeline wypisuje informację:

```text
STEP 11: disabled in publication-like pipeline version
```

i przechodzi dalej. Dane do klastrowania w Step 12 i Step 121 są pobierane bezpośrednio z pliku:

```text
<sample>.filtered_variants_by_cell.long.tsv.gz
```

---

## Step 12 — klasyczne klastrowanie mtDNA metodą hierarchical clustering

**Skrypt:**

```text
plot_heatmap_and_cluster_publication_mito_distance_auto_k_hclust_REWRITTEN_retained_tiebreak.R
```

**Wejście:**

```text
<sample>.filtered_variants_by_cell.long.tsv.gz
```

**Główne parametry:**

```text
PUBLICATION_COVERAGE_THRESHOLD=9
PUBLICATION_MIN_SHARED_VARIANTS=1
PUBLICATION_NA_DISTANCE_ACTION=one
PUBLICATION_REQUIRE_ALT_SIGNAL_CELL=1
PUBLICATION_MIN_ALT_FOR_CELL=1
AUTO_K_MIN=2
AUTO_K_MAX=30
AUTO_K_MIN_CLUSTER_SIZE=2
```

**Najważniejsze wyjścia:**

```text
<sample>.publication_mito.<DPTH>.auto.publication_mito_distance.D.tsv.gz
<sample>.publication_mito.<DPTH>.auto.publication_mito_relatedness.Kmito.tsv.gz
<sample>.publication_mito.<DPTH>.auto.publication_mito_shared_variants.N.tsv.gz
<sample>.publication_mito.<DPTH>.auto.publication_mito_distance.clustering_ready.D.tsv.gz
<sample>.publication_mito.<DPTH>.auto.auto_k_selection.tsv
<sample>.publication_mito.<DPTH>.auto.best_k.txt
<sample>.clusters.<DPTH>.k<best_k>.tsv
<sample>.publication_mito.<DPTH>.auto.hclust.assignments.k<best_k>.tsv
<sample>.publication_mito.<DPTH>.auto.publication_mito_distance.cluster_sizes.k<best_k>.tsv
<sample>.publication_mito.<DPTH>.auto.publication_mito_distance.heatmap.k<best_k>.pdf
<sample>.best_k.<DPTH>.txt
<sample>.valid_k.<DPTH>.txt
```

Ten krok wykonuje klastrowanie komórek na podstawie profili heteroplazmii mtDNA. Najpierw wybierane są komórki z rzeczywistym sygnałem ALT, jeśli `PUBLICATION_REQUIRE_ALT_SIGNAL_CELL=1`. Komórka jest dopuszczona do klastrowania, jeżeli dla co najmniej jednego wariantu spełnia:

```text
depth > PUBLICATION_COVERAGE_THRESHOLD
alt_count >= PUBLICATION_MIN_ALT_FOR_CELL
vaf > 0
```

Przy domyślnym progu `PUBLICATION_COVERAGE_THRESHOLD=9` oznacza to:

```text
depth >= 10
```

Następnie dla każdej pary komórek obliczany jest dystans mitochondrialny:

```text
D_ij = mean sqrt(abs(VAF_i - VAF_j))
```

Średnia liczona jest tylko po tych wariantach, które mają wystarczające pokrycie w obu komórkach:

```text
depth_i > PUBLICATION_COVERAGE_THRESHOLD
depth_j > PUBLICATION_COVERAGE_THRESHOLD
```

Minimalna liczba wspólnych wariantów wymagana do obliczenia dystansu wynosi domyślnie 1. Pipeline opisuje ten krok jako „publication-like mitochondrial distance clustering”.

Na podstawie macierzy dystansu wykonywane jest klastrowanie hierarchiczne:

```text
hclust(D, method = "average")
```

Następnie testowane są wartości `k` od `AUTO_K_MIN` do `AUTO_K_MAX`, ograniczone przez liczbę komórek. Dla każdej wartości `k` komórki są dzielone na klastry. Klastry mniejsze niż `AUTO_K_MIN_CLUSTER_SIZE` są przenoszone do kategorii:

```text
unassigned
```

Dla pozostałych przypisanych komórek liczona jest średnia silhouette. Najlepsze `k` wybierane jest według reguły:

```text
najwyższa mean_silhouette_assigned
potem najwyższy assigned_fraction
potem najwyższa liczba przypisanych komórek
potem niższa liczba retained clusters
potem najniższe k
```

Jeśli nie ma skończonego wyniku silhouette, skrypt stosuje tryb awaryjny oparty na liczbie zachowanych klastrów i liczbie przypisanych komórek.

Główny plik klastrów dla najlepszego `k` zawiera tylko komórki przypisane do klastrów:

```text
barcode
mt_cluster
```

Komórki przesunięte do `unassigned` są zachowane w pełnej tabeli przypisań, ale nie w głównym pliku klastrów przekazywanym do Step 15.

---

## Step 121 / Step 12B — fuzzy k-medoids na podstawie dystansu mtDNA

**Skrypt:**

```text
plot_heatmap_and_cluster_publication_mito_distance_fuzzy_kmedoids_auto_k_FINAL_retained_cluster_tiebreak.R
```

**Wejście:**

```text
<sample>.filtered_variants_by_cell.long.tsv.gz
```

**Główne parametry:**

```text
PUBLICATION_COVERAGE_THRESHOLD=9
PUBLICATION_MIN_SHARED_VARIANTS=1
PUBLICATION_NA_DISTANCE_ACTION=one
PUBLICATION_REQUIRE_ALT_SIGNAL_CELL=1
PUBLICATION_MIN_ALT_FOR_CELL=1
AUTO_K_MIN=2
AUTO_K_MAX=30
KMEDOIDS_ASSIGNMENT_PROB_THRESHOLD=0.95
KMEDOIDS_MIN_CLUSTER_SIZE=2
KMEDOIDS_MIN_RETAINED_CLUSTERS=1
KMEDOIDS_MIN_ASSIGNED_FRACTION=0
KMEDOIDS_MEMBERSHIP_EXPONENT=2
```

**Najważniejsze wyjścia:**

```text
<sample>.publication_mito.<DPTH>.fuzzy.publication_mito_distance.D.tsv.gz
<sample>.publication_mito.<DPTH>.fuzzy.publication_mito_distance.clustering_ready.D.tsv.gz
<sample>.publication_mito.<DPTH>.fuzzy.auto_k_selection.tsv
<sample>.publication_mito.<DPTH>.fuzzy.best_k.txt
<sample>.clusters_fuzzy.<DPTH>.k<best_k>.tsv
<sample>.publication_mito.<DPTH>.fuzzy.fuzzy_kmedoids.assignments.k<best_k>.tsv
<sample>.publication_mito.<DPTH>.fuzzy.fuzzy_kmedoids.membership.k<best_k>.tsv
<sample>.publication_mito.<DPTH>.fuzzy.fuzzy_kmedoids.cluster_sizes.k<best_k>.tsv
<sample>.publication_mito.<DPTH>.fuzzy.fuzzy_kmedoids.distance_heatmap.k<best_k>.pdf
<sample>.best_k_fuzzy.<DPTH>.txt
<sample>.valid_k_fuzzy.<DPTH>.txt
```

Ten krok używa tej samej tabeli wejściowej i tej samej definicji dystansu mitochondrialnego co Step 12. Różnica polega na metodzie klastrowania. Zamiast twardego cięcia drzewa hierarchicznego używane jest fuzzy k-medoids przez funkcję:

```text
cluster::fanny()
```

z macierzą dystansu jako wejściem i `diss=TRUE`.

Dla każdej testowanej wartości `k` skrypt oblicza macierz membership, czyli prawdopodobieństwo/przynależność każdej komórki do każdego klastra. Komórka jest początkowo przypisywana do klastra tylko wtedy, gdy jej maksymalne membership spełnia:

```text
max_membership >= 0.95
```

Jeśli maksymalne membership jest niższe niż 0.95, komórka pozostaje:

```text
unassigned
```

Następnie klastry z liczbą przypisanych komórek mniejszą niż `KMEDOIDS_MIN_CLUSTER_SIZE` są usuwane z listy klastrów wysokiej pewności, a ich komórki również przenoszone są do `unassigned`.

Najlepsze `k` wybierane jest spośród wyników spełniających warunki:

```text
fanny_ok == TRUE
mean_silhouette_assigned jest skończone
n_retained_clusters >= KMEDOIDS_MIN_RETAINED_CLUSTERS
assigned_fraction >= KMEDOIDS_MIN_ASSIGNED_FRACTION
```

Ranking kandydatów odbywa się według:

```text
najwyższa mean_silhouette_assigned
najwyższy assigned_fraction
najwyższa liczba przypisanych komórek
najniższe k
```

Skrypt zapisuje zarówno główny plik klastrów przypisanych komórek:

```text
barcode
mt_cluster
```

jak i pełną tabelę przypisań, zawierającą m.in.:

```text
barcode
mt_cluster_raw
max_membership
initially_assigned_by_membership
moved_to_unassigned_due_to_small_cluster
mt_cluster
assignment_status
```

oraz pełną macierz membership dla najlepszego `k`.

---

## Step 15 — integracja klasycznych klastrów mtDNA z obiektem Seurat

**Skrypt:**

```text
qc_mt_clusters_on_seurat_safe_terminal_v7_multik.R
```

**Wejście:**

```text
SCT_annotated_KSz.Rds
<sample>.uniq_mt_reads_per_cell.tsv
<sample>.filtered_variants_by_cell.long.tsv.gz
<sample>.step15_cluster_files.<DPTH>.txt
```

Plik `step15_cluster_files` zawiera dwie kolumny:

```text
k
cluster_file
```

Dla Step 15 wskazuje on główny plik klastrów z klasycznego hclust:

```text
<sample>.clusters.<DPTH>.k<best_k>.tsv
```

**Wyjście główne:**

```text
<sample>.mt.best_k.seurat_with_mt_lineage_multik.rds
<sample>.mt.best_k.seurat_only_mt_assigned_cells_union_multik.rds
<sample>.mt.best_k.seurat_metadata_with_mt_lineage_multik.tsv
<sample>.mt.best_k.mt_variant_hover_per_cell.tsv
<sample>.mt.best_k.filtered_variants_by_cell.long.with_seurat_barcodes.tsv
```

**Wyjścia per k:**

```text
<sample>.mt.best_k.k<best_k>.RNA_UMAP_harmony_by_<cluster_col>.pdf
<sample>.mt.best_k.k<best_k>.RNA_UMAP_harmony_<sample>_mt_clusters_highlighted.svg
<sample>.mt.best_k.k<best_k>.RNA_UMAP_harmony_<sample>_mt_clusters_highlighted.interactive.html
<sample>.mt.best_k.k<best_k>.<cluster_col>_vs_sample.tsv
<sample>.mt.best_k.k<best_k>.<cluster_col>_vs_patient.tsv
<sample>.mt.best_k.k<best_k>.<cluster_col>_vs_seurat_clusters.tsv
<sample>.mt.best_k.k<best_k>.<cluster_col>_vs_cluster_annotation.tsv
```

Skrypt ładuje obiekt Seurat, tabelę liczby odczytów mtDNA per komórka, długi plik wariantów mtDNA oraz listę plików klastrów. Barcode’y mtDNA są modyfikowane przez dodanie prefiksu próbki, np.:

```text
<sample>_<barcode>
```

aby odpowiadały nazwom komórek w obiekcie Seurat. Następnie sprawdzane jest dopasowanie barcode’ów między tabelami mtDNA i obiektem Seurat. Tworzone są pliki diagnostyczne z barcode’ami niedopasowanymi między źródłami.

Do metadanych Seurat dodawana jest kolumna:

```text
unique_mt_reads
```

oraz kolumna z klastrami mtDNA. Dla klasycznego Step 15 nazwa kolumny jest oparta na prefiksie:

```text
<sample>.mt_cluster
```

Komórki, które nie występują w pliku klastrów, otrzymują wartość:

```text
unassigned
```

Skrypt wykonuje także kontrolę bezpieczeństwa po metadanych: sprawdza, czy komórki przypisane do klastrów mtDNA należą do oczekiwanej próbki w kolumnie `sample`. Jeśli przypisane komórki pojawią się w nieoczekiwanej próbce, skrypt przerywa działanie, co chroni przed błędnym dopasowaniem barcode’ów.

Na potrzeby interaktywnych wykresów i opisów hover skrypt zachowuje tylko warianty z rzeczywistym sygnałem ALT:

```text
depth > publication_coverage_threshold
vaf > 0
```

i przygotowuje per-komórkowe informacje:

```text
mt_n_variants_detected
mt_variants_hover
mt_top_variant
mt_top_variant_depth
mt_top_variant_alt_count
mt_top_variant_vaf
mt_top_variant_alt_mean_bq
```

Jeśli w obiekcie Seurat dostępna jest redukcja `umap_harmony`, generowane są wykresy UMAP z naniesionymi klastrami mtDNA oraz wersje SVG i HTML z dodatkowymi informacjami o wariantach.

---

## Step 151 / Step 15B — integracja fuzzy klastrów mtDNA z Seurat

**Skrypt:**

```text
qc_mt_clusters_on_seurat_safe_terminal_v7_multik.R
```

**Wejście:**

```text
SCT_annotated_KSz.Rds
<sample>.uniq_mt_reads_per_cell.tsv
<sample>.filtered_variants_by_cell.long.tsv.gz
<sample>.step15B_fuzzy_cluster_files.<DPTH>.txt
```

Plik `step15B_fuzzy_cluster_files` wskazuje wynik fuzzy k-medoids:

```text
<sample>.clusters_fuzzy.<DPTH>.k<fuzzy_best_k>.tsv
```

**Wyjście:**

```text
<sample>.mt.fuzzy_best_k.*
```

Ten krok jest odpowiednikiem Step 15 dla wyników fuzzy k-medoids. Używa tego samego skryptu R, ale innego pliku klastrów i innej nazwy kolumny metadanych:

```text
<sample>.mt_cluster_fuzzy
```

Ważna różnica polega na tym, że plik fuzzy klastrów jest plikiem `assigned-only`. Oznacza to, że zawiera tylko komórki przypisane do klastrów z wysoką pewnością membership oraz po usunięciu małych klastrów. Komórki, które w Step 121 zostały oznaczone jako `unassigned`, nie pojawiają się w głównym pliku klastrów i w Seurat otrzymują wartość:

```text
unassigned
```

Dzięki temu wizualizacja fuzzy pokazuje tylko komórki przypisane z wysoką pewnością, a nie wymusza przypisania wszystkich komórek do klastrów.

---

## Podsumowanie przepływu danych

Cały przepływ można zapisać jako:

```text
BAM
↓
unikalne odczyty chrM, MAPQ >= 30, NH == 1
↓
<sample>.uniq_mt.bam
↓
liczba unikalnych odczytów mtDNA per komórka
↓
<sample>.uniq_mt_reads_per_cell.tsv
↓
komórki z unique_mt_reads >= 200
↓
<sample>.mt200_barcodes.tsv
↓
pileup A/C/G/T per komórka i pozycja
↓
<sample>.per_cell_per_allele_counts1.tsv.gz
↓
GMM na mean base quality
↓
<sample>.GMM_meanBQ_QC.tsv.gz
↓
empiryczny próg BQ
↓
<sample>.use_BQth.txt
↓
filtrowanie pozycji o niskiej jakości
↓
<sample>.filtered_plus.BQfiltered.per_cell_alleles.tsv.gz
↓
kandydackie warianty mtDNA
↓
<sample>.candidate_alleles.tsv
<sample>.candidate_alleles_by_cell.long.tsv.gz
↓
warianty przefiltrowane po heteroplazmii populacyjnej
↓
<sample>.filtered_variants.no_bulk.tsv
<sample>.filtered_variants_by_cell.long.tsv.gz
↓
klastrowanie mtDNA:
  Step 12: hclust average
  Step 121: fuzzy k-medoids
↓
pliki klastrów:
<sample>.clusters.<DPTH>.k<best_k>.tsv
<sample>.clusters_fuzzy.<DPTH>.k<fuzzy_best_k>.tsv
↓
integracja z Seurat i wizualizacja UMAP:
Step 15 / Step 151
```

---

## Najważniejsze założenia analityczne

1. Analiza wariantów mtDNA opiera się wyłącznie na odczytach jednoznacznie mapujących się do chromosomu mitochondrialnego `chrM`, z `MAPQ >= 30` i `NH == 1`.

2. Do dalszej analizy dopuszczane są tylko komórki mające co najmniej 200 unikalnych odczytów mitochondrialnych.

3. Jakość baz nie jest filtrowana na etapie pileupu. Zamiast tego pipeline najpierw zbiera pełny rozkład jakości, dopasowuje model GMM i wyznacza empiryczny próg wysokiej jakości.

4. Pozycje mtDNA są filtrowane konserwatywnie: jeśli istotny allel w danej pozycji ma średnią jakość poniżej progu, cała pozycja jest usuwana. Wyjątek stanowią allele o bardzo małym udziale, `allele_fraction < 1/500`.

5. Warianty kandydackie są definiowane względem allelu dominującego w danej pozycji, oznaczonego jako `ref_like`.

6. Do klastrowania używane są tylko warianty z heteroplazmią populacyjną co najmniej 2,5%.

7. Dystans między komórkami jest liczony jako średnia z `sqrt(abs(VAF_i - VAF_j))` po wariantach pokrytych w obu komórkach.

8. Step 12 wykonuje klasyczne hierarchical clustering metodą average linkage i automatycznie wybiera najlepsze `k`.

9. Step 121 wykonuje fuzzy k-medoids i przypisuje komórki tylko wtedy, gdy maksymalne membership wynosi co najmniej 0.95.

10. Step 15 i Step 151 nie wywołują nowych wariantów, tylko integrują gotowe klastry mtDNA z obiektem Seurat, dodają metadane i generują wizualizacje UMAP oraz tabele QC.
