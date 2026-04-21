import sys
import os
import getopt
import pandas as pd
import time
import subprocess
import shlex
from functools import reduce
import re

class LibInfo():
    def __init__(self, libname="", treat="", pair_mode="PE", r1="", r2="", adaptertype=None, qualityscore=""):
        self.libname = libname
        self.fastq_dir_r1 = r1
        self.fastq_dir_r2 = r2
        self.treat = treat
        self.pair_mode = pair_mode
        self.adaptertype = adaptertype
        self.qualityscore = qualityscore

    def outFastqList(self):
        if self.pair_mode == "SE":
            return [self.fastq_dir_r1]
        else:
            return [self.fastq_dir_r1, self.fastq_dir_r2]

class SampleSeries():
    def __init__(self, name="", species="", ref_genome=""):
        self.name = name
        self.species = species
        self.ref_genome = ref_genome
        self.lib_count = 0
        self.lib_list = []

    def parseSampleList(self, sample):
        SampleInfo = pd.DataFrame
        if os.path.exists(sample):
            FileType = os.path.splitext(sample)[-1]
            if FileType in ["", ".txt"]:
                SampleInfo = pd.read_csv(sample, sep="\t", header=0, encoding='utf-8')
            elif FileType in [".xlsx", ".xls"]:
                SampleInfo = pd.read_excel(sample, header=0)
            elif FileType in [".csv"]:
                SampleInfo = pd.read_csv(sample, header=0, encoding='utf-8')
            else:
                sys.stderr.write("Unknown type of file, exiting...")
        else:
            sys.stderr.write("Sample list does not exist, exiting...")

        SampleInfo.columns = SampleInfo.columns.map(lambda x: x.lower())

        for idx, row in SampleInfo.iterrows():
            LIBA = LibInfo()
            if 'lib' in row:
                LIBA.libname = row["lib"]
            if 'adaptertype' in row:
                LIBA.adaptertype = row['adaptertype']
            if 'qualityscore' in row:
                LIBA.qualityscore = row['qualityscore']
            if 'treat' in row:
                LIBA.treat = row['treat']
            if 'r1' in row:
                ReadOneDir = row['r1']
                if ReadOneDir != "":
                    if os.path.exists(ReadOneDir):
                        LIBA.fastq_dir_r1 = ReadOneDir
                    else:
                        raise FileExistsError(ReadOneDir)
            if 'r2' in row:
                ReadTwoDir = row['r2']
                if ReadTwoDir != "":
                    if os.path.exists(ReadTwoDir):
                        if LIBA.fastq_dir_r1 == "":
                            LIBA.fastq_dir_r1 = ReadTwoDir
                            LIBA.pair_mode = "SE"
                        else:
                            LIBA.fastq_dir_r2 = ReadTwoDir
                            LIBA.pair_mode = "PE"
                    else:
                        raise FileExistsError(ReadTwoDir)

            self.lib_list.append(LIBA)
            self.lib_count = self.lib_count + 1

def makeDirectory(target_list):
    for target in target_list:
        if not os.path.exists(target):
            os.makedirs(target)
        else:
            sys.stderr.write(target + " exists.")

def writeLog(target, type, mode, string):
    with open(target + "." + type, mode, encoding='utf-8') as F:
        F.write("==========\nMission time: " + time.asctime(time.localtime()) + "\n")
        F.write(string)

class FastQC():
    command_string = "fastqc -o {outdir} {fastqs}"

    def __init__(self, fastqs=[], outdir=".", lib="FastQC"):
        self.fastqs = fastqs
        self.outdir = outdir
        self.libname = lib
    
    def run(self):
        p = subprocess.Popen(self.command_string.format(outdir=self.outdir, fastqs=" ".join(self.fastqs)), 
                             shell=True, 
                             stderr=subprocess.STDOUT, 
                             stdout=subprocess.PIPE,
                             encoding='utf-8')

        print(time.asctime(time.localtime()) + "\t" + self.command_string.format(outdir=self.outdir, fastqs=" ".join(self.fastqs)))
        if p.stdout:
            writeLog(self.outdir + "/" + self.libname + ".fastqc", "log", "a", p.stdout.read())
        if p.stderr:
            writeLog(self.outdir+ "/" + self.libname + ".fastqc", "err", "a", p.stderr.read())

        if p.wait() != 0:
            raise RuntimeError("FastQC error.")

class MultiQC():
    command_string = "multiqc --interactive -o {outdir} {input}"

    def __init__(self, target, outdir="."):
        self.target = target
        self.outdir = outdir

    def run(self):
        p = subprocess.Popen(self.command_string.format(outdir=self.outdir, input=self.target),
                             shell=True, 
                             stderr=subprocess.STDOUT, 
                             stdout=subprocess.PIPE,
                             encoding='utf-8')

        print(time.asctime(time.localtime()) + "\t" + self.command_string.format(outdir=self.outdir, input=self.target))

        if p.stdout:
            writeLog(self.outdir + "/" + "multiqc", "log", "a", p.stdout.read())
        if p.stderr:
            writeLog(self.outdir+ "/" + "multiqc", "err", "a", p.stderr.read())

        if p.wait() != 0:
            raise RuntimeError("MultiQC error.")

class TrimGalore():
    command_string = "trim_galore --cores {cores} --{QualityScore} {AdapterType} -o {outdir} --dont_gzip {mode} {fastqs}"

    def __init__(self, fastqs=[], outdir=".", mode="PE", cores=10, adaptertype="illumina", adapter=None, adapter2=None, qualityscore="phred33", lib="TrimGalore"):
        self.fastqs = fastqs
        self.outdir = outdir
        self.adapter = adapter
        self.adapter2 = adapter2
        self.mode = mode
        self.cores = cores
        self.adaptertype = adaptertype
        self.qualityscore = qualityscore
        self.libname = lib
    
    def run(self):
        if self.adapter is None and self.adapter2 is None:
            if self.adaptertype and self.adaptertype != "auto":
                self.adaptertype = "--" + self.adaptertype
            else:
                self.adaptertype = ""
        elif self.adapter is not None and self.adapter2 is None:
            self.adaptertype = "-a " + self.adapter
        elif self.adapter is None and self.adapter2 is not None:
            self.adaptertype = "-a " + self.adapter2
        else:
            if self.mode == "PE":
                self.adaptertype = "-a " + self.adapter + " -a2 " + self.adapter2
            else:
                sys.stderr.write("[warnings] adapter2 provided but single end FASTQ provided! Will ignore adapter2...")

        if self.mode == "PE":
            self.mode = "--paired"
        else:
            self.mode = ""
        
        print(time.asctime(time.localtime()) + "\t" + self.command_string.format(outdir=self.outdir, fastqs=" ".join(self.fastqs), 
                                                            AdapterType=self.adaptertype, QualityScore=self.qualityscore, mode=self.mode, cores=self.cores))
        p = subprocess.Popen(self.command_string.format(outdir=self.outdir, fastqs=" ".join(self.fastqs), 
                                                        AdapterType=self.adaptertype, QualityScore=self.qualityscore, mode=self.mode, cores=self.cores), 
                             shell=True, 
                             stderr=subprocess.STDOUT, 
                             stdout=subprocess.PIPE,
                             encoding='utf-8')
        
        if p.stdout:
            writeLog(self.outdir + "/" + self.libname + ".trim_galore", "log", "a", p.stdout.read())
        if p.stderr:
            writeLog(self.outdir+ "/" + self.libname + ".trim_galore", "err", "a", p.stderr.read())

        if p.wait() != 0:
            raise RuntimeError("TrimGalore error.")

class StarAlign():
    command_string = "STAR --twopassMode Basic \
                        {RAMGB} \
                        --genomeDir {ref_genome} \
                        --outSAMunmapped Within \
                        --outFilterType BySJout \
                        --outSAMattributes NH HI AS NM MD MC \
                        --outFilterMultimapNmax {MultiMap} \
                        --outFilterMismatchNmax 999 \
                        --outFilterMismatchNoverReadLmax 0.04 \
                        --alignIntronMin 20 \
                        --alignIntronMax 1000000 \
                        --outWigType bedGraph \
                        --outWigNorm None \
                        {AlignMateGapMax} \
                        --alignSJoverhangMin 8 \
                        --alignSJDBoverhangMin 1 \
                        --sjdbScore 1 \
                        --runThreadN {Thread} \
                        --outSAMtype BAM SortedByCoordinate \
                        --quantMode TranscriptomeSAM GeneCounts \
                        --outSAMheaderHD @HD VN:1.4 SO:coordinate \
                        --outFileNamePrefix {Name} \
                        --readFilesIn {Fastq}"

    def __init__(self, ref_genome="", RAMGB=0, MultiMap=20, Thread=10, Name="STAR", Fastqs=[], mode="PE"):
        self.ref_genome = ref_genome
        self.RAMGB = RAMGB
        self.MultiMap = MultiMap
        self.Thread = Thread
        self.Name = Name
        self.Fastq = Fastqs
        self.pair_mode = mode
    
    def run(self):
        if self.RAMGB > 0:
            self.RAMGB = "--limitBAMsortRAM " + str(self.RAMGB) + "000000000"
        else:
            self.RAMGB = ""
        if self.pair_mode == "PE":
            AlignMateGapMax = "--alignMatesGapMax 1000000"
        else:
            AlignMateGapMax = ""

        Fastq = " ".join(self.Fastq)

        print(time.asctime(time.localtime()) + "\t" + " ".join(shlex.split(self.command_string.format(RAMGB=self.RAMGB,
                                                           ref_genome=self.ref_genome,
                                                           AlignMateGapMax=AlignMateGapMax,
                                                           MultiMap=self.MultiMap,
                                                           Thread=self.Thread,
                                                           Name=self.Name,
                                                           Fastq=Fastq))))
        
        p = subprocess.Popen(" ".join(shlex.split(self.command_string.format(RAMGB=self.RAMGB,
                                                           ref_genome=self.ref_genome,
                                                           AlignMateGapMax=AlignMateGapMax,
                                                           MultiMap=self.MultiMap,
                                                           Thread=self.Thread,
                                                           Name=self.Name,
                                                           Fastq=Fastq))),
                             shell=True,
                             stderr=subprocess.STDOUT,
                             stdout=subprocess.PIPE,
                             encoding='utf-8')

        if p.wait() != 0:
            raise RuntimeError("STAR align error.")

class RsemCount():
    command_string = "rsem-calculate-expression --num-threads {Thread} \
                        --alignments \
                        -q \
                        --bam \
                        {PairMode} \
                        --seed 12345 \
                        --estimate-rspd \
                        --no-bam-output \
                        --forward-prob 0 \
                        {InputBAM} \
                        {RsemRefDir} \
                        {outdir}"
    
    def __init__(self, ref_genome="", Thread=10, inputs="", outdir="", pair_mode="PE"):
        self.ref_genome = ref_genome
        self.Thread = Thread
        self.inputs = inputs
        self.outdir = outdir
        self.pair_mode = pair_mode
    
    def run(self):
        if self.pair_mode == "PE":
            self.pair_mode = "--paired-end"
        else:
            self.pair_mode = ""

        print(time.asctime(time.localtime()) + "\t" +" ".join(shlex.split(self.command_string.format(Thread=self.Thread,
                                                        PairMode=self.pair_mode,
                                                        InputBAM=self.inputs,
                                                        RsemRefDir=self.ref_genome,
                                                        outdir=self.outdir))))
        p = subprocess.call(" ".join(shlex.split(self.command_string.format(Thread=self.Thread,
                                                        PairMode=self.pair_mode,
                                                        InputBAM=self.inputs,
                                                        RsemRefDir=self.ref_genome,
                                                        outdir=self.outdir))),
                            shell=True, stderr=None, stdout=None)

def parseTrimAdapter(inputadapter):
    adapter = None 
    adapter2 = None
    adaptertype = None

    if inputadapter:
        inputadapter = str(inputadapter)
        AdapterSetFlag = 0
        PatternAdapter = re.compile(r'(-a\s*)([atcg]+)', re.I)
        PatternAdapter2 = re.compile(r'(-a2\s*)([atcg]+)', re.I)
        PatternAdapterMatch = PatternAdapter.search(inputadapter)
        PatternAdapterMatch2 = PatternAdapter2.search(inputadapter)
        
        if PatternAdapterMatch:
            adapter = PatternAdapterMatch.group(2)
            AdapterSetFlag = 1
        if PatternAdapterMatch2:
            adapter2 = PatternAdapterMatch2.group(2)
            AdapterSetFlag = 1

        if (AdapterSetFlag == 0) and (inputadapter.lower() in ["illumina", "nextera", "small_rna", "auto"]):
            adaptertype = inputadapter.lower()

    return {"adapter": adapter, "adapter2": adapter2, "adaptertype": adaptertype}

def parseTrimQualityScore(qualityscore):
    if qualityscore in ["phred33", "phred64"]:
        return qualityscore
    else:
        return "phred33"

def get_uniquely_mapped_reads(log_file):
    with open(log_file, 'r') as f:
        for line in f:
            if "Uniquely mapped reads number" in line:
                return int(line.split('|')[-1].strip())
    return None

def calculate_rpm(input_files, uniquely_mapped_reads, output_file):
    command = f"bedtools unionbedg -i {' '.join(input_files)} | awk '{{OFS=\"\\t\"; print $1,$2,$3,$4/{uniquely_mapped_reads}*1000000*100,$5/{uniquely_mapped_reads}*1000000*100}}' > {output_file}"
    subprocess.run(command, shell=True, check=True)

def Usage():
    print("Usage: python RNA-seq.py [options...]\n")
    print("Options:\n")
    print("  -i FILE    Sample information.\n")
    print("  -n STR     Project name.\n")
    print("  -p INT     Core number.\n")
    print("  -r DIR     STAR reference directory.\n")
    print("  -s DIR     RSEM reference directory.\n")

def main():
    Experiment = "RNA-seq"
    SampleDir = ""
    Cores = 1
    RefGenomeDir = ""
    RSEMRefDir = ""

    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:n:p:r:s:", ["input=", "name=", "cores=", "star-reference=", "rsem-reference="])
    except getopt.GetoptError:
        Usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt in ['h', "--help"]:
            Usage()
            sys.exit(2)
        elif opt in ['-i', "--input="]:
            SampleDir = arg
        elif opt in ['-n', "--name="]:
            Experiment = arg
        elif opt in ['-p', "--cores="]:
            Cores = arg
        elif opt in ["-r", "--star-reference="]:
            RefGenomeDir = arg
        elif opt in ["-s", "--rsem-reference="]:
            RSEMRefDir = arg

    if not os.path.exists(SampleDir):
        raise FileExistsError(SampleDir)
    if not os.path.exists(RefGenomeDir):
        raise FileExistsError(RefGenomeDir)
    if not os.path.exists(RSEMRefDir + ".grp"):
        raise FileExistsError(RSEMRefDir)
    
    makeDirectory(["01.QC/raw_fastqc/", "01.QC/raw_multiqc/", "01.QC/trimmed_fastq/", "02.Align", "03.Count", "04.DEG"])

    SampleInfo = SampleSeries()
    SampleInfo.parseSampleList(SampleDir)

    # fastqc for each library
    for i in range(SampleInfo.lib_count):
        FastQC(SampleInfo.lib_list[i].outFastqList(), "01.QC/raw_fastqc", SampleInfo.lib_list[i].libname).run()

    # multi qc
    MultiQC(target="01.QC/raw_fastqc", outdir="01.QC/raw_multiqc").run()

    # trim adapter
    for i in range(SampleInfo.lib_count):
        lib_info = SampleInfo.lib_list[i]
        ParsedAdapter = parseTrimAdapter(lib_info.adaptertype)

        TrimGalore(fastqs=lib_info.outFastqList(), 
                   outdir="01.QC/trimmed_fastq", 
                   mode=lib_info.pair_mode, 
                   cores=Cores, 
                   adaptertype=ParsedAdapter["adaptertype"],
                   adapter=ParsedAdapter["adapter"],
                   adapter2=ParsedAdapter["adapter2"],
                   qualityscore=parseTrimQualityScore(lib_info.qualityscore),
                   lib=lib_info.libname).run()

        SampleInfo.lib_list[i].fastq_dir_r1 = "01.QC/trimmed_fastq/" + lib_info.libname + "_R1_val_1.fq"
        if lib_info.pair_mode == "PE":
            SampleInfo.lib_list[i].fastq_dir_r2 = "01.QC/trimmed_fastq/" + lib_info.libname + "_R2_val_2.fq"

    # trimmed fastqc
    for i in range(SampleInfo.lib_count):
        FastQC(fastqs=SampleInfo.lib_list[i].outFastqList(), outdir="01.QC/raw_fastqc", lib=SampleInfo.lib_list[i].libname + ".trimmed").run()
    MultiQC(target="01.QC/raw_fastqc", outdir="01.QC/raw_multiqc").run()

    # star align
    for i in range(SampleInfo.lib_count):
        lib_info = SampleInfo.lib_list[i]
        StarAlign(ref_genome=RefGenomeDir, Thread=Cores, Name="02.Align/" + lib_info.libname, Fastqs=lib_info.outFastqList()).run()

    # qc star align
    MultiQC(target="02.Align/*Log.final.out", outdir="02.Align/STAR.MultiQC").run()

    # Calculate RPM
    for i in range(SampleInfo.lib_count):
        lib_info = SampleInfo.lib_list[i]
        log_file = f"02.Align/{lib_info.libname}Log.final.out"
        uniquely_mapped_reads = get_uniquely_mapped_reads(log_file)
        if uniquely_mapped_reads:
            input_files = [f"02.Align/{lib_info.libname}Signal.Unique.str2.out.bg", f"02.Align/{lib_info.libname}Signal.Unique.str1.out.bg"]
            output_file = f"02.Align/{lib_info.libname}.Unique.RPM.x100.bdg2"
            calculate_rpm(input_files, uniquely_mapped_reads, output_file)

    # rsem count
    for i in range(SampleInfo.lib_count):
        lib_info = SampleInfo.lib_list[i]
        RsemCount(ref_genome=RSEMRefDir, Thread=Cores, inputs="02.Align/" + lib_info.libname + "Aligned.toTranscriptome.out.bam", pair_mode=lib_info.pair_mode, outdir="03.Count/" + lib_info.libname).run()

    # compile RSEM Count logs
    MultiQC(target="03.Count/*.stat", outdir="03.Count/RSEM.MultiQC").run()

    # generate summary LIB Count table
    SampleTable = pd.DataFrame(columns=["condition", "pairmode", "libname"], )

    GeneTableList = []
    IsoformTableList = []
    for i in range(SampleInfo.lib_count):
        lib_info = SampleInfo.lib_list[i]
        ColumnPrefix = lib_info.libname
        CountTable = pd.read_csv("03.Count/" + lib_info.libname + ".genes.results", sep="\t", header=0)
        CountTable.rename(columns={"length": ColumnPrefix + ".length", 
                                   "effective_length": ColumnPrefix + ".effective_length", 
                                   "expected_count": ColumnPrefix + ".expected_count",
                                   "TPM": ColumnPrefix + ".TPM",
                                   "FPKM": ColumnPrefix + ".FPKM"}, inplace=True)
        GeneTableList.append(CountTable)

        CountTable = pd.read_csv("03.Count/" + lib_info.libname + ".isoforms.results", sep="\t", header=0)
        CountTable.rename(columns={"length": ColumnPrefix + ".length", 
                                   "effective_length": ColumnPrefix + ".effective_length", 
                                   "expected_count": ColumnPrefix + ".expected_count",
                                   "TPM": ColumnPrefix + ".TPM",
                                   "FPKM": ColumnPrefix + ".FPKM"}, inplace=True)
        IsoformTableList.append(CountTable)

        SampleTable.loc[ColumnPrefix] = [lib_info.treat, lib_info.pair_mode, lib_info.libname]

    TotalGeneTable = reduce(lambda left,right: pd.merge(left, right, on=["gene_id", "transcript_id(s)"], how="outer"), GeneTableList).fillna(0)
    TotalIsoformTable = reduce(lambda left,right: pd.merge(left, right, on=["transcript_id", "gene_id"], how="outer"), IsoformTableList).fillna(0)

    pd.DataFrame.to_csv(TotalGeneTable, "03.Count/" + Experiment + ".GeneTotal.Results", sep="\t", index=False)
    pd.DataFrame.to_csv(TotalIsoformTable, "03.Count/" + Experiment + ".IsoformsTotal.Results", sep="\t", index=False)
    pd.DataFrame.to_csv(SampleTable, "03.Count/" + Experiment + ".sample.txt", sep="\t", index=True)

if __name__ == "__main__":
    print("Start processing...")
    main()
    print("End processing.")
