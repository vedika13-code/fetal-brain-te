# Samples: HDBR prenatal brain RNA-seq (ArrayExpress E-MTAB-4840), ENA sample accessions
SAMPLES="early_CS20 mid_14pcw"
declare -A ENA=( [early_CS20]=ERS1221416 [mid_14pcw]=ERS1221652 )
MAX_READS=4000000        # reads streamed per sample (first N reads of the run)

DATA=data
REF=ref
WORK=work
RESULTS=results
THREADS=4
