# Zusammenfassung der NGS-Analyse-Vorlesung und Klausurvorbereitung

## Grundlagen der NGS-Analyse

### 1. Wichtige Konzepte:
- **NGS (Next Generation Sequencing)**: Hochdurchsatz-Sequenzierungstechnologien, die große Mengen an DNA/RNA schnell und kostengünstig sequenzieren
- **OMICs**: Verschiedene "Omics"-Bereiche (Genomik, Transkriptomik, Proteomik etc.) in der biomedizinischen Forschung

### 2. Anwendungsbereiche (Klausuraufgabe 1.1):
- **Tumordiagnostik**: Identifikation von Krebsmutationen
- **Seltene genetische Erkrankungen**: Finden krankheitsverursachender Varianten
- **Pränataldiagnostik**: Untersuchung fetaler DNA
- **Antibiotikaresistenz**: Analyse resistenter Bakterienstämme

### 3. Typische NGS-Pipeline (Klausuraufgabe 1.2):
1. **Library Preparation**: Vorbereitung der DNA-Proben
2. **Sequenzierung**: Erzeugung der Rohdaten (FASTQ-Dateien)
3. **Alignment**: Zuordnung der Reads zum Referenzgenom (BAM-Dateien)
4. **Varianten-Calling**: Identifikation von Abweichungen (VCF-Dateien)
5. **Annotation/Priorisierung**: Bewertung der gefundenen Varianten

## Wichtige Dateiformate (Klausuraufgabe 1.4)

| Format | Beschreibung | Inhalt |
|--------|-------------|--------|
| **FASTQ** | Rohdatenformat | Sequenz + Qualitätswerte pro Read |
| **BAM/SAM** | Alignments | Zuordnung der Reads zum Referenzgenom |
| **VCF** | Varianten | Gefundene Mutationen (SNPs, Indels etc.) |
| **Pileup** | Zwischenformat | Positionelle Übersicht der Alignments |
| **GFF** | Genomannotation | Genpositionen und -funktionen |

## Alignment und Varianten (Klausuraufgabe 2)

### 1. Alignment-Methoden:
- **Seed-and-Extend**: Verwendung kurzer "Seeds" (kmers) zur schnellen Suche, gefolgt von genauerer Ausrichtung (z.B. BWA, Bowtie)
- **Burrows-Wheeler-Transform (BWT)**: Effiziente Indexierung für schnelle Suche (in BWA verwendet)

### 2. Varianten-Typen:
- **Synonyme Varianten**: Veränderte DNA-Sequenz, aber gleiche Aminosäure (keine Proteinveränderung)
- **Nicht-synonyme Varianten**: Veränderte DNA-Sequenz führt zu anderer Aminosäure (Proteinveränderung)
- **Frameshift-Varianten**: Insertion/Deletion, die den Leserahmen verschiebt (meist schwerwiegende Proteinveränderung)

## Tools und Python-Pakete (Klausuraufgaben 1.3 und 3.2)

### Wichtige Tools:
1. **BWA** (Burrows-Wheeler Aligner) - für Alignment
2. **Samtools** - für BAM-Dateimanipulation
3. **GATK** - für Variantencalling
4. **VarScan** - einfacher SNV-Caller
5. **Delly** - für strukturelle Varianten

### Python-Pakete für NGS-Analyse:
- **pysam** - für BAM/VCF-Dateien
- **Biopython** - allgemeine Bioinformatik-Funktionen
- **pandas** - Datenanalyse
- **scikit-learn** - für Machine-Learning-Anwendungen

## Klausurtipps

1. **Pseudocode schreiben können** (mehrere Aufgaben):
   - Variantenerkennung aus BAM-Datei
   - Deletionserkennung
   - CNV-Analyse
   - OLC-Assembly

2. **Methoden verstehen**:
   - de Bruijn Graphen und Euler-Pfade (Assembly)
   - MUMmer-Algorithmus (Comparative Genomics)
   - Phylogenetische Baumkonstruktion

3. **Formeln kennen**:
   - z-Score für CNV-Normalisierung

4. **Visualisierungen beschreiben/zeichnen können**:
   - Paired-End Reads bei strukturellen Varianten
   - de Bruijn Graphen

## Lernstrategie für die Woche

1. **Grundkonzepte verstehen** (Alignment, Variantentypen, Pipeline)
2. **Wichtige Tools und Dateiformate lernen**
3. **Pseudocode-Übungen machen**
4. **Algorithmen (de Bruijn, MUMmer) wiederholen**
5. **Altklausurenfragen durchgehen**

Brauchst du zu einem bestimmten Thema noch mehr Erklärungen oder Beispiele? Ich kann gerne noch detaillierter auf einzelne Aspekte eingehen!




### Lernzusammenfassung: NGS-Analyse mit Python – Lecture 03

#### **1. Einführung in die Variantenerkennung**
- **Ziel**: Identifikation genomischer Varianten (SNVs, InDels) in menschlichen Genomen.
- **Werkzeuge**: 
  - Picard, samtools (PCR-Duplikatentfernung).
  - GATK, ABRA2 (Indel-Realignment).

#### **2. Pileup-Format**
- **Struktur**: Enthält Informationen zu Chromosom, Position, Referenzbase, Tiefe (DP), beobachtete Basen und Qualitätswerte.
- **Beispiel**:
  ```
  seq1 272 T 24 ,.$.....,,.,.,...,,,.,..^+. <<<+;<<<<<<<<<<<=<;<;7<&
  ```
- **Generierung**: 
  ```bash
  samtools mpileup -f $REF sample.bam > sample.pileup
  ```

#### **3. Einfacher SNV-Caller**
- **Kriterien für Varianten**:
  - Mindest-Tiefe (z.B. DP > 10).
  - Mindest-Allelhäufigkeit (z.B. AF > 0.2).
  - Strand-Bias-Filterung.
  - Basenqualität (z.B. Q > 20).
- **Algorithmus**:
  1. Zählen der Referenz- und Alternativbasen.
  2. Berechnung der Allelfrequenz: `AF = alt_count / DP`.
  3. Filterung basierend auf Qualitätsschwellenwerten.

#### **4. Erweiterte Filterung**
- **Basenqualität**: Umwandlung von ASCII zu Phred-Score:
  ```python
  base_quality = ord(qualities[index]) - 33
  ```
- **Qualitätskontrolle**: 
  - Entfernung von Basen mit niedriger Qualität (z.B. Q < 20).
  - Summe der Qualitätswerte für zusätzliche Filterung.

#### **5. VCF-Format (Variant Call Format)**
- **Struktur**: Standardisiertes Format für Variantenspeicherung.
  - Enthält Metadaten (##), Kopfzeile (#CHROM, POS, ID, ...) und Varianteninformationen.
- **Beispiel**:
  ```
  #CHROM POS ID REF ALT QUAL FILTER INFO
  chr1 100 . A G 50 PASS DP=100;AF=0.5
  ```

#### **6. Probabilistische Varianten-Caller**
- **Nachteile einfacher Methoden**:
  - Verlust von Sensitivität durch feste Schwellenwerte.
  - Keine Berücksichtigung von Genotyp-Wahrscheinlichkeiten.
- **Lösungen**: 
  - Bayes'sche Modelle (z.B. BCFtools, GATK).
  - Verwendung von Genotyp-Likelihoods.

#### **7. Praktische Implementierung in Python**
- **Schritte**:
  1. **Argumentparsing**: 
     ```python
     parser.add_argument('-p', '--pileup', required=True, help='Pileup-Datei')
     ```
  2. **Pileup-Parsing**:
     ```python
     with open(infile) as f:
         for line in f:
             chrom, pos, ref, depth, bases, quals = line.split("\t")
     ```
  3. **SNV-Erkennung**:
     ```python
     if alt_count >= min_alt and depth >= min_depth and AF >= min_af:
         print(f"SNV gefunden: {chrom} {pos}")
     ```
  4. **Ausgabe im VCF-Format**.

#### **8. Tools und Ressourcen**
- **Galaxy**: Benutzerfreundliche Oberfläche für Variantenanalyse (z.B. VarScan).
- **VarScan**: Java-basierter Variantencaller für NGS-Daten.
- **BCFtools**: Kommandozeilentool für Variantenaufruf:
  ```bash
  samtools mpileup -uf $REF sample.bam | bcftools call -mv -Oz -o variants.vcf.gz
  ```

#### **9. Python-Grundlagen für die NGS-Analyse**
- **Wichtige Konzepte**:
  - **Listen/Dictionaries**: Für Datenorganisation.
  - **Datei-I/O**: `with open()` für sichere Handhabung.
  - **Klassen**: Modularisierung des Codes (z.B. `SNVCaller`-Klasse).
- **Beispiel**:
  ```python
  class SNVCaller:
      def __init__(self, min_depth, min_af):
          self.min_depth = min_depth
          self.min_af = min_af
      def detect_snv(self, chrom, pos, ref, depth, bases):
          # Implementierung der SNV-Erkennung
  ```

#### **10. Zusammenfassung der Qualitätskriterien**
| **Kriterium**          | **Beispielwert** | **Beschreibung**                          |
|------------------------|------------------|-------------------------------------------|
| Mindest-Tiefe (DP)     | 10               | Minimale Lese-Tiefe für Variantenaufruf.  |
| Mindest-Allelhäufigkeit (AF)| 0.2          | Minimale Häufigkeit des Alternativallels. |
| Basenqualität (Q)      | 20               | Phred-Score für Basenqualität.            |
| Strand-Bias            | Kein Bias        | Vermeidung von einseitiger Abdeckung.     |

#### **11. Weiterführende Themen**
- **Indel-Erkennung**: Ähnlich zu SNVs, aber mit Fokus auf Insertionen/Deletionen.
- **Genotypisierung**: Zuweisung von Genotypen (z.B. 0/0, 0/1, 1/1).
- **Performance-Optimierung**: Nutzung von Multiprocessing für große Datensätze.

Diese Zusammenfassung bietet einen Überblick über die Schlüsselkonzepte der Vorlesung, von grundlegenden Python-Implementierungen bis hin zu fortgeschrittenen Varianten-Calling-Methoden. Für vertiefende Übungen empfiehlt sich die praktische Anwendung mit Beispiel-Datensätzen (z.B. aus dem 1000-Genome-Projekt).




### **Lernzusammenfassung: NGS-Analyse mit Python – Lecture 04**  
*(Relevanz für Klausuraufgaben: Variantenerkennung, VCF-Format, Qualitätsfilterung, Probabilistische Methoden)*  

---

#### **1. Variantenerkennung mit Python**  
**Klausurrelevanz**: Implementierung eines SNV-Callers (Aufgabe: Code-Analyse oder Algorithmus-Design)  
- **Ziel**: Erkennung von SNVs und InDels aus Pileup-Dateien.  
- **Schritte**:  
  1. **Pileup-Parsing**: Extrahiere Chromosom, Position, Referenzbase, Tiefe (DP), Basen und Qualitätswerte.  
  2. **SNV-Detektion**:  
     - Zähle Referenz- (`.,`) und Alternativbasen (`A,a,C,c,G,g,T,t`).  
     - Filtere nach:  
       - Mindest-Tiefe (`-d`, z.B. DP ≥ 10).  
       - Mindest-Allelhäufigkeit (`-af`, z.B. AF ≥ 0.2).  
       - Basenqualität (`-bq`, Phred-Score ≥ 20).  
     ```python
     if alt_count >= min_alt and depth >= min_depth and AF >= min_af:
         print(f"SNV: {chrom}:{pos}, Ref={ref}, Alt={alt_base}")
     ```

---

#### **2. VCF-Format (Variant Call Format)**  
**Klausurrelevanz**: VCF-Header und Einträge erstellen (Aufgabe: Datenformatierung)  
- **Header**: Metadaten mit `##fileformat`, `##INFO`, `##FORMAT`.  
- **Einträge**:  
  - **Spalten**: `CHROM`, `POS`, `ID`, `REF`, `ALT`, `QUAL`, `FILTER`, `INFO`, `FORMAT`, `SAMPLE`.  
  - **Beispiel**:  
    ```python
    vcf_line = f"{chrom}\t{pos}\t.\t{ref}\t{alt_base}\t{qual}\tPASS\tAC={alt_count};AF={AF}\tGT:AD\t0/1:{alt_count}"
    ```
- **Wichtig**: `INFO`-Felder wie `AC` (Allele Count), `AF` (Allele Frequency), `DP` (Depth).

---

#### **3. Qualitätsfilterung**  
**Klausurrelevanz**: Filterstrategien (Aufgabe: Qualitätskriterien bewerten)  
- **Basenqualität**: ASCII zu Phred-Score:  
  ```python
  base_quality = ord(qual_str[index]) - 33  # Sanger-Encoding
  ```  
- **Strand-Bias**: Fisher’s Exact Test (zur Vermeidung von Sequenzierungsartefakten):  
  ```python
  oddsratio, pvalue = stats.fisher_exact([[ref_fwd, ref_rev], [alt_fwd, alt_rev]])
  if pvalue < 0.05:  # Signifikanter Strand-Bias → Filter
  ```

---

#### **4. Probabilistische Varianten-Caller**  
**Klausurrelevanz**: Bayes’sche Methoden vs. Threshold-basierte Ansätze (Aufgabe: Vor-/Nachteile vergleichen)  
- **Problem**: Hard Thresholds (z.B. DP > 10) verlieren Sensitivität bei niedriger Coverage.  
- **Lösung**: Genotyp-Likelihoods berechnen:  
  - **Homozygot Referenz (aa)**: Wahrscheinlichkeit, dass `n-k` Basen Fehler sind.  
  - **Heterozygot (ab)**: Binomialverteilung mit `p=0.5`.  
  - **Bayes’scher Ansatz**: Prior-Wissen (z.B. SNP-Häufigkeiten in Populationen) integrieren.  

---

#### **5. InDel-Erkennung (Hausaufgabe)**  
**Klausurrelevanz**: InDel-Calling in Pileup-Daten (Aufgabe: Algorithmus skizzieren)  
- **Herausforderung**: Erkennung von `+n`/`-n`-Mustern (z.B. `+2AG` für Insertion).  
- **Ansatz**:  
  - Suche nach `+`/`-` in den Basencalls.  
  - Extrahiere Länge und eingefügte/gelöschte Sequenz.  
  - Filtere nach Qualität und Strand-Balance.  

---

#### **6. Tools & Bibliotheken**  
**Klausurrelevanz**: Werkzeugauswahl (Aufgabe: Tool-Vergleich)  
- **BCFtools**: Probabilistischer Caller mit `mpileup` + `call`.  
- **GATK**: Haplotype-Caller für komplexe Varianten.  
- **VarScan**: Einfacher Threshold-basierter Caller.  

---

### **Zusammenfassung der Qualitätskriterien**  
| **Kriterium**          | **Parameter**      | **Klausurrelevanz**                     |  
|------------------------|--------------------|-----------------------------------------|  
| Tiefe (DP)             | `-d 10`            | Mindest-Coverage (Aufgabe: Filterung)   |  
| Allelfrequenz (AF)     | `-af 0.2`          | Heterozygotie-Erkennung                 |  
| Basenqualität (Q)      | `-bq 20`           | Fehlerfilterung (Phred-Score)           |  
| Strand-Bias (SB)       | Fisher-Test (p<0.05)| Artefaktfilterung                       |  

---

### **Wichtige Code-Snippets für die Klausur**  
1. **Pileup-Parsing**:  
   ```python
   with open(pileup) as f:
       for line in f:
           chrom, pos, ref, depth, bases, quals = line.split("\t")
   ```  
2. **VCF-Schreiben**:  
   ```python
   header = "##fileformat=VCFv4.1\n#CHROM\tPOS\tREF\tALT\n"
   with open(out.vcf, "w") as f:
       f.write(header + vcf_line)
   ```  

--- 

**Tipp für die Klausur**:  
- Merken Sie sich die **VCF-Spaltenstruktur** und die **Berechnung des Phred-Scores** (`ord(char) - 33`).  
- Verstehen Sie den Unterschied zwischen **Threshold-basiertem** und **probabilistischem Calling**.  

*Für Lecture 03 siehe vorherige Zusammenfassung (Pileup-Format, SNV-Caller-Implementierung).*


### Zusammenfassung der Lecture 05: NGS-Analyse mit Python

#### **1. Einführung in die Variantenannotation**
- **Problemstellung**: Identifikation der kausalen Variante unter Millionen von SNVs (Single Nucleotide Variants).
- **Priorisierung**: 
  - Genom: ~4M SNVs → Exom: ~80k SNVs → Proteinverändernd: ~15k SNVs → Kausale Variante: 1 SNV.
- **Pathogenitätsklassifikation**: 
  - Populationsfrequenz (z.B. gnomAD), phänotypische Daten (HPO), evolutionäre Konservierung, Spleißdefekte.

#### **2. Variantenannotation mit Ensembl VEP (Variant Effect Predictor)**
- **Funktionen**:
  - Annotation genomischer Regionen (genic, intergenic, coding, intronic).
  - Effekte auf Gen/Protein (synonym, missense, Stop-Gain, Frameshift, Spleißdefekte).
  - Funktionale Auswirkungen (Schadensvorhersage, Domänen, 3D-Struktur).
  - Populationsfrequenz und Krankheitsassoziation.
- **Datenbanken**: 
  - Ensembl, ClinVar, OMIM, COSMIC, gnomAD, SIFT, PolyPhen-2, CADD.

#### **3. Praktische Anwendung von VEP**
- **Web-Tool**: Eingabe von Varianten (HGVS-Notation, Region, ID) zur Annotation.
- **REST-API**: Automatisierte Abfragen in Python (JSON/XML-Format).
  - Beispiel: Abfrage einer Variante auf Chromosom 1 (Position 6524705, Allel T).
- **Python-Code**:
  ```python
  import requests
  server = "https://rest.ensembl.org"
  ext = "/vep/human/region/1:6524705:6524705/T?"
  r = requests.get(server + ext, headers={"Content-Type": "application/json"})
  decoded = r.json()
  print(decoded)
  ```

#### **4. Strukturelle Varianten (SVs) und Copy Number Variants (CNVs)**
- **Typen**: Deletionen, Duplikationen, Inversionen, Translokationen.
- **Detektionsmethoden**:
  - **Read-Count-Analyse**: Abweichungen in der Abdeckung.
  - **Paired-End Mapping (PEM)**: Diskordante Read-Paare (längere/kürzere Insertgrößen).
  - **Split-Reads**: Reads, die über Bruchpunkte hinweg alignen.
  - **Soft-Clipped Reads**: Teilweise Alignments an Bruchpunkten.
- **Pipeline**: Alignment → Clustering diskordanter Reads → Breakpoint-Bestimmung → Merging.

#### **5. Implementierung eines Deletions-Callers in Python**
- **Schritte**:
  1. Insertgrößenverteilung berechnen (Mittelwert + Standardabweichung).
  2. Diskordante Reads identifizieren (Insertgröße > Mittelwert + 3×SD).
  3. Clustering von Reads zur Bestimmung von Deletionsgrenzen.
- **Python-Bibliotheken**:
  - **pysam**: Lesen/Schreiben von BAM-Dateien, Extraktion von Insertgrößen.
  - **matplotlib**: Visualisierung der Insertgrößenverteilung.
  - **statistics**: Berechnung von Mittelwert und Standardabweichung.
- **Beispielcode**:
  ```python
  import pysam
  bam = pysam.AlignmentFile("sample.bam", "rb")
  for read in bam.fetch():
      if read.is_paired and abs(read.template_length) > threshold:
          print(f"Discordant read at {read.reference_name}:{read.reference_start}-{read.reference_end}")
  ```

#### **6. Anwendungsfall: Antibiotikaresistente Bakterien**
- **Beispiel**: Detektion von Deletionen in *Pseudomonas aeruginosa*.
- **Ergebnis**: Identifikation von Genverlusten, die mit Resistenz assoziiert sind.

#### **7. Zusammenfassung**
- **Variantenannotation**: Kritisch für die Priorisierung kausaler Varianten.
- **Strukturelle Varianten**: Erfordern spezielle Methoden (PEM, Split-Reads).
- **Python-Tools**: pysam, requests, matplotlib ermöglichen automatisierte Analysen.

Diese Lecture vermittelt praktische Fähigkeiten zur Analyse von NGS-Daten, von der Annotation einzelner Varianten bis zur Detektion komplexer struktureller Veränderungen.






### Fortsetzung der Lecture 06: Detektion struktureller Varianten (SVs) mit Python

#### **1. Grundlagen struktureller Varianten**
- **Typen genomischer Alterationen**:
  - **Deletionen**: Verlust von DNA-Abschnitten.
  - **Duplikationen**: Vervielfältigung von DNA-Abschnitten.
  - **Inversionen**: Umkehrung der DNA-Sequenz.
  - **Translokationen**: Verlagerung von DNA zwischen Chromosomen.
- **Detektionsmethoden**:
  - **Read-Count-Analyse**: Abweichungen in der Sequenzierabdeckung.
  - **Paired-End Mapping (PEM)**: Analyse diskordanter Read-Paare.
  - **Split-Reads**: Reads, die über Bruchpunkte hinweg alignen.
  - **Soft-Clipped Reads**: Teilweise Alignments an Bruchpunkten.

#### **2. Paired-End Mapping (PEM) für SV-Detektion**
- **Prinzip**:
  - **Konkordante Read-Paare**: Erwartete Insertgröße und Orientierung.
  - **Diskordante Read-Paare**: Abweichende Insertgröße/Orientierung (z. B. durch Deletionen).
- **Insertgrößenverteilung**:
  - Normalverteilung mit Mittelwert (z. B. 300–400 bp) und Standardabweichung.
  - **Schwellenwert für Diskordanz**: Mittelwert + 3× oder 4× Standardabweichung.
- **Beispiel für Deletionen**:
  - Längere Insertgröße → Reads flankieren die Deletion (ein Read upstream, einer downstream).

#### **3. Implementierung eines Deletions-Callers in Python**
- **Schritte**:
  1. **Insertgrößenverteilung berechnen**:
     ```python
     import pysam
     bam = pysam.AlignmentFile("sample.bam", "rb")
     insert_sizes = [read.template_length for read in bam.fetch() if read.is_proper_pair and read.template_length > 0]
     mean, std = statistics.mean(insert_sizes), statistics.stdev(insert_sizes)
     ```
  2. **Diskordante Reads identifizieren**:
     ```python
     discordant_threshold = mean + 4 * std
     discordant_reads = [(read.reference_name, read.reference_end, read.next_reference_start) 
                         for read in bam.fetch() 
                         if read.is_paired and abs(read.template_length) > discordant_threshold]
     ```
  3. **Clustering**:
     - Gruppierung von Reads mit ähnlichen Breakpoints (z. B. maximale Distanz innerhalb von 2× Standardabweichung).
     - **Ausgabe**: Deletionsgrenzen als TSV-Datei.
       ```python
       with open("deletions.tsv", "w") as f:
           f.write("DEL\tchr\tstart\tend\n")
           for chrom, start, end in clusters:
               f.write(f"DEL\t{chrom}\t{start}\t{end}\n")
       ```

#### **4. Tools und Bibliotheken**
- **pysam**:
  - Lesen/Schreiben von BAM-Dateien, Zugriff auf Read-Merkmale (`is_paired`, `template_length`, etc.).
  - Beispiel:
    ```python
    read = next(bam.fetch())
    print(f"Chromosom: {read.reference_name}, Position: {read.reference_start}-{read.reference_end}")
    ```
- **matplotlib/statistics**:
  - Visualisierung der Insertgrößenverteilung und statistische Auswertung.

#### **5. Anwendungsfall: Antibiotikaresistenz in Pseudomonas aeruginosa**
- **Daten**:
  - BAM-Datei (`BK6866.sort.bam`) mit Reads eines resistenten Stammes.
- **Ergebnis**:
  - Detektion von Gen-Deletionen, die mit Resistenz assoziiert sind (z. B. Verlust von Resistenzgenen).

#### **6. Zusammenfassung und Ausblick**
- **Strukturelle Varianten** erfordern spezielle Analysemethoden (PEM, Split-Reads).
- **Python-Tools** (pysam, matplotlib) ermöglichen automatisierte SV-Detektion.
- **Erweiterungen**:
  - Implementierung für Insertionen/Inversionen (Analyse der Orientierung diskordanter Reads).
  - Integration von Split-Read-Analysen für präzisere Breakpoints.

---

### Hinweis zur Gedächtnisklausur
- **Relevante Themen**:
  - Unterschiede zwischen SNVs und SVs.
  - Methoden zur SV-Detektion (PEM, Split-Reads, Read-Count).
  - Bedeutung der Insertgrößenverteilung.
  - Python-Code zur Identifikation diskordanter Reads.
- **Tipp**: Nutzen Sie die praktischen Beispiele (z. B. `deletion_caller.py`) zur Vertiefung!
<img width="644" height="8899" alt="grafik" src="https://github.com/user-attachments/assets/5ad9167b-e10f-4070-8097-76c25076650e" />

