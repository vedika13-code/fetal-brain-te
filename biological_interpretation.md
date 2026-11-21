# Biological interpretation: TE expression in prenatal brain (demo)

Reading notes for the interview. Part 1 is the biology you need. Part 2 interprets each result
from `fetal-brain-te`. Part 3 connects it to the BrainSpan project. Part 4 is likely questions
with short answers.

---

## Part 1. The biology

### What TEs are
- **Transposable elements** are DNA sequences that can move to new positions in the genome.
  About **half of human DNA** is TE-derived, but almost all copies are dead: mutated, truncated,
  unable to move.
- **Class I, retrotransposons (copy-and-paste):** transcribed to RNA, reverse-transcribed to DNA,
  inserted somewhere new. The original stays, so copies accumulate. Nearly all human TEs are this kind.
- **Class II, DNA transposons (cut-and-paste):** move as DNA via a transposase. All inactive in humans.

### The main retrotransposons
| Element | Type | Size | Can it move? |
|---|---|---|---|
| **L1 (LINE-1)** | non-LTR, autonomous | ~6 kb full length | Yes. The only autonomous active family; ~100 intact L1HS copies |
| **Alu** | SINE, non-autonomous | ~300 bp | Yes, using L1's proteins; >1 million copies |
| **SVA** | composite, non-autonomous | ~2 kb | Yes, using L1's proteins; young, hominid-specific |
| **ERVs / LTR elements** | remnants of ancient retroviral infections | variable | Essentially no in humans |
| **MIR, L2** | very old SINE / LINE | – | No; ancient, highly diverged |

### How L1 copies itself
- A full-length L1 has a **5′ UTR containing its promoter**, then **ORF1** (an RNA-binding protein)
  and **ORF2** (endonuclease + reverse transcriptase), then a 3′ UTR with a poly-A tail.
- The L1 RNA and its proteins return to the nucleus. ORF2p **nicks the genomic DNA** and uses the
  nicked end as a primer to **reverse-transcribe the RNA directly at the new site**: target-primed
  reverse transcription (TPRT).
- Reverse transcription often stops early, so most new copies are **5′-truncated**, lack the
  promoter, and are dead on arrival. That is why your L1HS examples in the annotation were only
  ~250 bp.
- The L1 5′ UTR also has an **antisense promoter**, which can drive transcription *out* of the L1
  into neighbouring sequence. Useful when discussing antisense TE reads.

### How cells keep TEs silent
- **DNA methylation** of TE promoters (CpG-rich L1 5′ UTRs).
- **KRAB zinc-finger proteins (KZFPs):** hundreds of them, each recognising specific TE sequences.
  They recruit **KAP1/TRIM28**, which brings **SETDB1** to lay down **H3K9me3** (heterochromatin).
- **HUSH complex** silences young L1s inside active genes.
- **piRNAs** in the germline.

### Age and divergence
- When a TE copy inserts, it matches its "master" sequence exactly. Over time, it accumulates random
  mutations.
- RepeatMasker reports each copy's **divergence from the subfamily consensus**. **Low divergence =
  young.** L1HS and AluY are young; MIR and L2 are very old (>100 million years).
- Young copies are **nearly identical to each other**, so short reads from them map to many places.
  That's the core computational problem.

### TEs in the developing brain
- L1 is **active in neural progenitor cells**: engineered L1 reporters retrotranspose in rat and human
  neural progenitors (Muotri et al. 2005; Coufal et al. 2009). This can make neurons slightly
  genetically different from each other (**somatic mosaicism**).
- **Playfoot et al. 2021 (Genome Research)** analysed TE expression across 16 brain regions from
  8 post-conception weeks to adulthood. Key findings: a distinct **KZFP:TE expression profile at the
  late prenatal to early postnatal transition**, and TE-derived **alternative promoters** driving
  neurogenesis genes. Their analysis was largely at **subfamily** level. The Mitacs project's
  **locus-specific** focus is the step beyond.
- Increased L1 copy number / activity has been **reported** in schizophrenia brain samples
  (e.g. Bundo et al. 2014). Disease links are an active, not settled, area.

### Expression is not mobility
An RNA-seq atlas measures **transcription**, not new insertions. TE-derived reads can come from:
1. a TE using **its own promoter** (what people usually mean by "TE activation");
2. a TE inside an **exon or UTR** of a gene: the TE is just part of the gene's mRNA;
3. a TE in an **intron**, read as part of the gene's **unspliced pre-mRNA** (read-through);
4. a TE acting as an **alternative promoter** for a nearby gene.
Only a handful of full-length, intact L1s could actually retrotranspose. Saying this unprompted is a
strong signal you understand the field.

---

## Part 2. Interpreting your results

### Result 1: a third of "multimappers" were reference artefacts
- After alignment, 7–9% of read pairs were multimapped. Removing alignments to `_alt` contigs
  (alternate haplotype versions of some regions in hg38) made **32–35% of those unique**.
- **Meaning:** some multimapping is a property of the **reference**, not of the biology. A read from a
  region present twice in the reference (main chromosome + `_alt`) is not genuinely ambiguous.
- **Take-home:** check where multimapping comes from before attributing it to repeats.

### Result 2: TE reads by class
| Class | % of TE reads |
|---|---|
| SINE (Alu, MIR) | ~45% |
| LINE (L1, L2) | ~31% |
| DNA transposons | ~12% |
| LTR (ERVs) | ~11% |
| SVA | ~0.3% |

- This roughly follows **how much of each class sits inside genes**, especially in UTRs and introns
  of expressed genes. It **does not** mean SINEs are "more active"; most of these reads are
  host-gene transcripts (see Result 4).
- Both samples look almost identical, which is expected: the bulk TE signal is dominated by
  gene expression, which is broadly similar between the two stages.

### Result 3: young TEs lose the most signal without multimappers
| Subfamily | Kept by unique-only counting |
|---|---|
| L1HS | 40% |
| L1PA2 | 60% |
| AluYa5 | 67% |
| SVA_F | 68% |
| AluYb8 | 77% |
| Overall, all TEs | ~96% |

- **Overall, multimappers add only ~4%**, because most TE reads come from old, diverged copies that
  map uniquely.
- **But for the youngest elements, 30–60% of the signal is in multimapping reads.** L1HS, the only
  family that can still retrotranspose, loses the most.
- **Meaning:** a standard pipeline looks fine on average but silently loses most of the signal from
  exactly the elements a TE atlas cares about. This is the justification for TE-aware tools.
- **Caveat on SVA:** RepeatMasker divergence is inflated for SVA because of its composite,
  repeat-rich structure (VNTR region), so its position on the age axis isn't comparable with L1/Alu.

### Fractional vs EM assignment (how the multimappers were counted)
- I used **fractional assignment**: a read with 10 alignments gives 0.1 to each location.
- **EM** (TEtranscripts, TElocal, Telescope) instead iterates: start evenly, then re-weight each read
  towards locations with more unique-read evidence, until stable. EM concentrates reads on the copies
  that are actually expressed.
- Fractional splitting **over-spreads** reads onto silent copies; at locus level that matters, at
  subfamily level much less (all copies are summed anyway).

### Result 4: where expressed TE loci sit
| Context (loci with ≥5 read pairs) | early_CS20 | mid_14pcw |
|---|---|---|
| Exonic, same strand as gene | 6,307 | 5,893 |
| Intronic, same strand | 2,118 | 2,549 |
| Intergenic | 659 | 626 |
| Intronic, opposite strand | 343 | 325 |
| Exonic, opposite strand | 149 | 170 |

The library is **reverse-stranded** and counting was **strand-specific**: a read counts for a TE only
if it comes from the TE's own strand. That makes the orientation informative:

- **Exonic, same strand (~65%):** TE fragments that became part of genes, mostly in **3′ UTRs**.
  The reads are the **gene's mRNA**. This is why the top subfamilies are old SINEs (MIR, AluJ, AluS)
  and L2: they are the TEs most often embedded in UTRs. Biologically real ("exaptation"), but not TE
  activation.
- **Intronic, same strand (~20%):** most likely **pre-mRNA / intron retention** of the host gene.
  Nuclear RNA contains unspliced transcripts, so reads across intronic TEs are expected.
  Interesting: slightly more of these in the 14 pcw sample. With n = 1 and different regions, I would
  not interpret it.
- **Intronic or exonic, opposite strand (~5%):** **cannot be the host gene's own transcript.** Best
  candidates for **independent TE transcription**, e.g. from a TE promoter (including L1's antisense
  promoter) or from unannotated antisense transcripts.
- **Intergenic (~7%):** no annotated gene to explain them. Candidates for autonomous TE transcription,
  but also possibly **unannotated lncRNAs**, enhancer RNAs, or genes missing from the annotation.

**Take-home:** only a small minority (~12%) of expressed TE loci are good candidates for TE-driven
transcription. Separating them **needs stranded data and context annotation**. A count table alone
would mislead.

### Result 5: autism risk genes (SFARI)
- 1,125 expressed TE loci in 504 SFARI genes; **96% same strand as the gene**, so mostly the autism
  genes' own transcripts. About **40 antisense loci**, in genes such as *PTPRD* and *RBFOX1*.
- **Important confound: gene length.** *PTPRD* and *RBFOX1* are among the longest human genes
  (each well over 1 Mb), and long neuronal genes are full of intronic TEs. Many autism genes are long
  neuronal genes. So **finding many TE loci in SFARI genes is expected from size alone**, and is **not**
  evidence of a TE–autism link.
- To test for real enrichment you'd compare against **length-matched, expression-matched non-SFARI
  genes**.
- What's genuinely interesting: **antisense TE transcription inside long neuronal genes** could, in
  principle, interfere with their expression (antisense transcripts can affect transcription or
  splicing). That's a hypothesis to test, not a finding.

### Early vs mid comparison
- Class composition and top subfamilies are nearly identical. With **n = 1 per stage**, different
  brain regions, and only 4M read pairs, **no developmental conclusion is possible**. The demo
  shows the method, not biology.

---

## Part 3. What this means for the BrainSpan atlas

1. **Keep multimappers and use EM at locus level** (TElocal / Telescope / SQuIRE). Young elements are
   where it matters most.
2. **Report confidence per locus.** Loci supported mainly by multimappers are uncertain; flag them.
3. **Classify every locus by context and strand** to separate TE-driven transcription from host-gene
   read-through. **First question for the lab: is BrainSpan's RNA-seq stranded?** If not, orientation
   can't be used, and read-through is much harder to separate.
4. **Model the developmental design properly:** age as continuous or by developmental period, region
   as a factor, **donor** as a random/batch effect (several regions per brain), and covariates such as
   **RIN, post-mortem interval, sex**.
5. **Reference choices:** remove `_alt` contigs or use an analysis set; consider **T2T-CHM13**, which
   resolves many repeats missing or collapsed in hg38.
6. **Polymorphic insertions:** donors carry non-reference TE insertions; their reads are lost or
   misassigned. Donor genotypes (if available) or tools like MELT/xTea would help.
7. **Validation:** long-read RNA-seq to confirm locus-specific transcripts; chromatin data (H3K9me3,
   H3K4me3 at TE promoters) to support autonomous activity; comparison with Playfoot et al. 2021.

---

## Part 4. Likely questions, short answers

**Why are multimapping reads a problem for TEs?**
Young TE copies are nearly identical, so reads match many locations equally well. Standard pipelines
discard them, which removes most signal from young elements like L1HS (60% in my demo).

**What's the difference between subfamily- and locus-level quantification?**
Subfamily sums all copies (e.g. all L1HS); locus assigns reads to individual copies. Locus is harder
(multimappers must be resolved to one copy) but tells you *which* copy is expressed, which matters for
regulation and for linking to nearby genes or disease loci.

**How does EM work here?**
Start by splitting each multimapper evenly across its locations. Estimate each location's expression,
then re-split every read in proportion to those estimates. Repeat until the estimates stop changing.
Locations with more unique-read support end up with more of the shared reads.

**Does TE expression mean TEs are jumping?**
No. Most TE-derived RNA is TE sequence inside genes' transcripts or pre-mRNA. Only a few intact L1s
can retrotranspose, and expression is necessary but not sufficient.

**How would you tell autonomous TE transcription from read-through?**
Stranded data plus context: antisense-to-gene and intergenic loci are better candidates. Also look
for reads starting at the TE's promoter (5′ end), a drop in coverage at TE boundaries vs the
surrounding intron, and supporting chromatin marks.

**Why did you find so many TEs in autism genes?**
Mostly gene length: autism genes are often very long neuronal genes full of TEs, and 96% of the loci
were on the gene's own strand. A real test needs length- and expression-matched control genes.

**What would you do differently with more resources?**
EM-based locus-level tool on a cluster; all BrainSpan samples with donor and quality covariates;
T2T reference; long-read validation.

**What are the main limitations of your demo?**
n = 1 per stage, different regions; first 4M reads; fractional (not EM) assignment; reference TE
copies only; no statistics.

---

## Key references
- Muotri AR et al. (2005) Somatic mosaicism in neuronal precursor cells mediated by L1
  retrotransposition. *Nature*.
- Coufal NG et al. (2009) L1 retrotransposition in human neural progenitor cells. *Nature*.
- Bundo M et al. (2014) Increased L1 retrotransposition in the neuronal genome in schizophrenia.
  *Neuron*.
- Jin Y et al. (2015) TEtranscripts. *Bioinformatics*.
- Bendall ML et al. (2019) Telescope. *PLoS Computational Biology*.
- Playfoot CJ et al. (2021) Transposable elements and their KZFP controllers are drivers of
  transcriptional innovation in the developing human brain. *Genome Research*.
