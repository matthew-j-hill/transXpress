import os
import sys
import shutil
import re
import csv
import Bio.SeqIO
import subprocess
from snakemake.utils import min_version
from Bio.Seq import Seq
import pandas as pd

min_version("5.4.1")

configfile: "config.yaml"

# for executable in ["samtools", "bowtie2", "kallisto", "signalp6", "targetp", "Trinity", "blastp", "cmscan", "hmmscan", "fastqc", "seqkit", "R"]:
#    if not shutil.which(executable):
#        sys.stderr.write("Warning: Cannot find %s in your PATH\n" % (executable))

#TRINITY_EXECUTABLE_PATH=shutil.which("Trinity")
#TRINITY_HOME=os.path.dirname(os.path.join(os.path.dirname(TRINITY_EXECUTABLE_PATH), os.readlink(TRINITY_EXECUTABLE_PATH))) ##Have to resolve the symbolic link that conda makes

# https://github.com/griffithlab/rnaseq_tutorial/wiki/Trinity-Assembly-And-Analysis

# These rules don't need to be sent to a cluster
localrules: all, clean, trimmomatic_split, trimmomatic_merge, trinity_butterfly_split, transcriptome_copy, compare_qc_after_trim, fastqc_before_trim_warnings


rule all:
  """
  List of target files of the transxpress pipeline.
  """
  input:
    "samples_trimmed.txt",
    "transcriptome.fasta",
    "transcriptome.pep",
    "transcriptome_stats.txt",
    "ExN50_plot.pdf",
    "transcriptome_annotated.fasta",
    "transcriptome_annotated.pep",
    "transcriptome_TPM_blast.csv",
    "busco",
    "busco_report.txt",
    "fastqc_before_trim",
    "multiqc_before_trim",
    "fastqc_after_trim",
    "multiqc_after_trim",
    "multiqc_before_trim.txt",
    "multiqc_after_trim.txt",
    "WARNING_fastqc_before_trim_overview.txt",
    "FastQC_comparison_after_trim.txt",
    "edgeR_trans"

rule clean:
  """
  Removes files from the transxpress directory when specifically called by the user.
  """

  shell:
    """
    if [ -f samples_trimmed.txt ]; then
      cut -f 2 < samples_trimmed.txt | xargs --no-run-if-empty rm -rf
    fi
    rm -rf trinity_* tmp* log kallisto* transcriptome* pipeliner* annotation* transdecoder* trimmomatic* samples_trimmed* ExN50_plot.pdf multiqc fastqc edgeR_trans
    """


rule fastqc_before_trim:
  """
  Runs fastQC on individual input reads files.
  """
  input:
    samples=config["samples_file"]
  output:
    directory("fastqc_before_trim")
  log:
    "logs/fastqc_before_trim.log"
  conda:
    "/projects/wenglab/testtube/matthew/miniforge3/envs/transxpress-qc"
  params:
    memory="4"
  threads:
    4
  shell:
    """
    mkdir {output} &> {log}
    FILES=$(awk '{{ printf "%s\\n%s\\n", $3,$4}}' {input})
    fastqc -f fastq -t {threads} -o {output} $FILES &>> {log}
    """


rule multiqc_before_trim:
  """
  Creates multiQC report from individual fastQC reports.
  """
  input:
    "fastqc_before_trim"
  output:
    out_dir=directory("multiqc_before_trim"),
    report="multiqc_before_trim.txt"
  log:
    "logs/multiqc_before_trim.log"
  conda:
    "/projects/wenglab/testtube/matthew/miniforge3/envs/transxpress-qc"
  params:
    memory="4"
  threads:
    1
  shell:
    """
    # common error resolved by those two export commands
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    mkdir {output.out_dir} &> {log}
    multiqc -o {output.out_dir} {input} &>> {log}
    cp multiqc_before_trim/multiqc_data/multiqc_fastqc.txt {output.report} &>> {log}
    """

rule fastqc_before_trim_warnings:
  """
  Prints warnings for FastQC modules which produced some warnings or fails of the input data.
  """
  input:
    fastqc_report="multiqc_before_trim.txt"
  output:
    warning_file="WARNING_fastqc_before_trim_overview.txt"
  log:
    "logs/fastqc_before_trim_warnings.log"
  run:
    from fastqc_warnings import mapping, reasons_mapping
    
    with open(output["warning_file"], "w") as out_file:
      df1 = pd.read_csv(input["fastqc_report"], sep='\t')
      df1 = df1.set_index('Sample')

      # filter only columns containing 'pass','warn' or 'fail'
      df1_filtered = df1.loc[:, df1.isin(['pass','warn','fail']).all()]
      df1_filtered = df1_filtered.loc[:, df1_filtered.columns != 'basic_statistics'] # basic statistics module always passes

      n_samples = df1.shape[0]

      for column in df1_filtered.columns:
        if column in mapping.keys():
          out_file.writelines(mapping[column]+'\n')
        else:
          out_file.writelines(column+'\n')
        results = df1_filtered[column].tolist()
        if ('warn' in results) or ('fail' in results):
          n_warn = results.count('warn')
          n_fail = results.count('fail')
          n_warn_perc = round(n_warn/n_samples*100,2)
          n_fail_perc = round(n_fail/n_samples*100,2)
          if n_warn > 0:
            out_file.writelines(f'{n_warn} out of {n_samples} produced "WARN" ({n_warn_perc}%)\n')
          if n_fail > 0:
            out_file.writelines(f'{n_fail} out of {n_samples} produced "FAIL" ({n_fail_perc}%)\n')
          if column in reasons_mapping.keys():
            out_file.writelines('\n' + reasons_mapping[column])
        else:
          out_file.writelines('ok'+'\n')
      out_file.writelines('source: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/\n')
    
    # print to stdout
    with open(output["warning_file"], "r") as out_file:
      for line in out_file.readlines():
        print(line.strip())
      print('See this output in: WARNING_fastqc_before_trim_overview.txt')
    
checkpoint trimmomatic_split:
  """
  Splits file with information about all reads files into separate files
  so they can be processed with trimmomatic in parallel.
  """   
  input:
    samples=config["samples_file"]
  output:
    directory("trimmomatic")
  log:
    "logs/trimmomatic_split.log"
  shell:
    """
    mkdir -p {output} &> {log}
    split -d -l 1 {input} {output}/sample_ &>> {log}
    """


rule trimmomatic_parallel:
  """
  Processes individual read files with trimmomatic.
  """
  input:
    "trimmomatic/sample_{job_index}"
  output:
    "trimmomatic/completed_{job_index}"
  log:
    "logs/trimmomatic_parallel{job_index}.log"
  conda:
    "/projects/wenglab/testtube/matthew/miniforge3/envs/transxpress-trimmomatic"
  params:
    memory="10"
  threads:
    4
  shell:
    """
    # note: trimmomatic can use gzipped files directly
    read SAMPLE REPLICATE F_READS R_READS < {input}
    # If the sample line is empty, ignore it
    if [ -z "$REPLICATE" ]; then
      touch {output}
      exit 0
    fi
    if [ ! -z "$R_READS" ]; then
      trimmomatic PE -threads {threads} $F_READS $R_READS trimmomatic/{wildcards[job_index]}.R1-P.qtrim.fastq.gz trimmomatic/{wildcards[job_index]}.R1-U.qtrim.fastq.gz trimmomatic/{wildcards[job_index]}.R2-P.qtrim.fastq.gz trimmomatic/{wildcards[job_index]}.R2-U.qtrim.fastq.gz {config[trimmomatic_parameters]} &> {log}
      echo $SAMPLE	$REPLICATE	trimmomatic/{wildcards[job_index]}.R1-P.qtrim.fastq.gz	trimmomatic/{wildcards[job_index]}.R2-P.qtrim.fastq.gz > {output} 2>> {log}
    else
      trimmomatic SE -threads {threads} $F_READS trimmomatic/{wildcards[job_index]}.U.qtrim.fastq.gz {config[trimmomatic_parameters]} &> {log}
      echo $SAMPLE      $REPLICATE      trimmomatic/{wildcards[job_index]}.U.qtrim.fastq.gz > {output} 2>> {log}
    fi
    """


def trimmomatic_completed_parallel_jobs(wildcards):
  """
  Returns names of files with information about files processed with trimmomatic.
  """
  parallel_dir = checkpoints.trimmomatic_split.get(**wildcards).output[0]
  job_ids = glob_wildcards(os.path.join(parallel_dir, "sample_{job_index}")).job_index
  completed_ids = expand(os.path.join(parallel_dir,"completed_{job_index}"), job_index=job_ids)
  return completed_ids


rule trimmomatic_merge:
  """
  Creates a file with information about all reads files processed with trimmomatic.
  """
  input:
    trimmomatic_completed_parallel_jobs
  output:
    samples_trimmed="samples_trimmed.txt"
  log:
    "logs/trimmomatic_merge.log"
  shell:
    """
    cat {input} > {output.samples_trimmed} 2> {log}
    """

rule fastqc_after_trim:
  """
  Runs fastQC on individual trimmed reads files.
  """
  input:
    samples="samples_trimmed.txt"
  output:
    directory("fastqc_after_trim")
  log:
    "logs/fastqc_after_trim.log"
  conda:
    "/projects/wenglab/testtube/matthew/miniforge3/envs/transxpress-qc"
  params:
    memory="4"
  threads:
    4
  shell:
    """
    mkdir {output} &> {log}
    FILES=$(awk '{{ printf "%s\\n%s\\n", $3,$4}}' {input})
    fastqc -f fastq -t {threads} -o {output} $FILES &>> {log}
    """


rule multiqc_after_trim:
  """
  Creates multiQC report from individual fastQC reports after the trimming.
  """
  input:
    "fastqc_after_trim"
  output:
    out_directory=directory("multiqc_after_trim"),
    report="multiqc_after_trim.txt"
  log:
    "logs/multiqc_after_trim.log"
  conda:
    "/projects/wenglab/testtube/matthew/miniforge3/envs/transxpress-qc"
  params:
    memory="4"
  threads:
    1
  shell:
    """
    # common error resolved by those two export commands
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    mkdir {output.out_directory} &> {log}
    multiqc -o {output.out_directory} {input} &>> {log}
    cp multiqc_after_trim/multiqc_data/multiqc_fastqc.txt {output.report} &>> {log}
    """

rule compare_qc_after_trim:
  """
  Compares the quality of reads before and after trimming.
  Compares PER BASE SEQUENCE QUALITY, PER SEQUENCE QUALITY SCORES, 
  OVERREPRESENTED SEQUENCES and ADAPTER CONTENT modules of FastQC.
  """
  input:
    before="multiqc_before_trim.txt",
    after="multiqc_after_trim.txt"
  output:
    result="FastQC_comparison_after_trim.txt"
  log:
    "logs/compare_qc_after_trim.log"
  run:

    mapping = {
      "per_base_sequence_quality": "PER BASE SEQUENCE QUALITY",
      "per_sequence_quality_scores": "PER SEQUENCE QUALITY SCORES",
      "overrepresented_sequences": "OVERREPRESENTED SEQUENCES",
      "adapter_content": "ADAPTER CONTENT"
    }

    def compare_module(module_name, before_trim, after_trim, out_file):
      module_name = mapping[module_name]
      out_file.writelines(f"""**{module_name}**
Before trimming:
pass: {str(before_trim.count('pass'))} out of {str(len(before_trim))} ({str(before_trim.count('pass')/len(before_trim)*100)}%)
warn: {str(before_trim.count('warn'))} out of {str(len(before_trim))} ({str(before_trim.count('warn')/len(before_trim)*100)}%)
fail: {str(before_trim.count('fail'))} out of {str(len(before_trim))} ({str(before_trim.count('fail')/len(before_trim)*100)}%)
      
After trimming:
pass: {str(after_trim.count('pass'))} out of {str(len(after_trim))} ({str(after_trim.count('pass')/len(after_trim)*100)}%)
warn: {str(after_trim.count('warn'))} out of {str(len(after_trim))} ({str(after_trim.count('warn')/len(after_trim)*100)}%)
fail: {str(after_trim.count('fail'))} out of {str(len(after_trim))} ({str(after_trim.count('fail')/len(after_trim)*100)}%)
{'-' * 50}
""")

    before = pd.read_csv(input["before"], sep='\t')
    before = before.set_index('Sample')

    after = pd.read_csv(input["after"], sep='\t')
    after = after.set_index('Sample')

    modules = ['per_base_sequence_quality', 'per_sequence_quality_scores', 'overrepresented_sequences', 'adapter_content']

    with open(output["result"], "w") as out_file:
      out_file.write("FASTQC COMPARISON BEFORE AND AFTER TRIMMING\n")
      out_file.write('-' * 50 + '\n')

      for module in modules:
        before_list = before[module].tolist()
        after_list = after[module].tolist()
        compare_module(module, before_list, after_list, out_file)

    # print the results also in the terminal
    with open(output["result"], "r") as out_file:
      for line in out_file.readlines():
        print(line.strip())
      print('See this output in FastQC_comparison_after_trim.txt')


rule trinity_inchworm_chrysalis:
  """
  Runs first two stages of Trinity assembly (Inchworm and Chrysalis).
  """
  input:
    samples="samples_trimmed.txt",
  output:
    "trinity_out_dir/recursive_trinity.cmds"
  log:
    "logs/trinity_inchworm_chrysalis.log"
  conda:
    "/projects/wenglab/testtube/matthew/miniforge3/envs/transxpress-trinityutils"
  params:
    memory="256"
  threads:
    16
  shell:
    """
    Trinity --no_distributed_trinity_exec --max_memory {params.memory}G --CPU {threads} --samples_file {input} {config[trinity_parameters]} {config[strand_specific]} &> {log}
    """


checkpoint trinity_butterfly_split:
  """
  Preparation for last stage of Trinity (Butterfly) parallelization by splitting
  the commands into independent parts.
  """
  input:
    "trinity_out_dir/recursive_trinity.cmds"
  output:
    directory("trinity_out_dir/parallel_jobs")
  log:
    "logs/trinity_split.log"
  params:
    memory="4",
    cpus="2"
  shell:
    """
    mkdir -p {output} &> {log}
    sed 's/--max_memory [^ ]*/--max_memory {params.memory}G/g; s/--CPU [0-9]*/--CPU {params.cpus}/g' {input} | shuf > {output}/shuffled_cmds.tmp 2>> {log}
    split -n l/1000 -e -d {output}/shuffled_cmds.tmp {output}/job_ &>> {log}
    rm {output}/shuffled_cmds.tmp &>> {log}
    """


rule trinity_butterfly_parallel:
  """
  Runs Trinity Butterfly commands (which were split into parts) in parallel.
  """
  input:
    "trinity_out_dir/parallel_jobs/job_{job_index}"
  output:
    "trinity_out_dir/parallel_jobs/completed_{job_index}"
  log:
    "logs/trinity_parallel{job_index}.log"
  conda:
    "/projects/wenglab/testtube/matthew/miniforge3/envs/transxpress-trinityutils"
  params:
    memory="4"
  threads:
    2
  shell:
    """
    set -eo pipefail
    while IFS= read -r cmd; do
        [ -z "$cmd" ] && continue
        echo "=== Running: $cmd" >> {log}
        eval "$cmd" >> {log} 2>&1
    done < {input}
    cp -p {input} {output} >> {log} 2>&1
    """


def trinity_completed_parallel_jobs(wildcards):
  """
  Returns filenames of the files processed in parallel.
  """
  parallel_dir = checkpoints.trinity_butterfly_split.get(**wildcards).output[0]
  job_ids = glob_wildcards(os.path.join(parallel_dir, "job_{job_index}")).job_index
  completed_ids = expand(os.path.join(parallel_dir,"completed_{job_index}"), job_index=job_ids)
  return completed_ids

rule trinity_butterfly_parallel_merge:
  input:
    jobs=trinity_completed_parallel_jobs,
    cmds="trinity_out_dir/recursive_trinity.cmds"
  output:
    cmds_completed="trinity_out_dir/recursive_trinity.cmds.completed"
  log:
    "logs/trinity_butterfly_parallel_merge.log"
  params:
    memory="10"
  threads:
    1
  shell:
    """
    # Can crash if there are too many parallel jobs
    # See https://bitbucket.org/snakemake/snakemake/issues/878/errno-7-argument-list-too-long-path-to-bin 
    cat {input.jobs} > {output.cmds_completed} 2> {log}
    """


rule trinity_final:
  """
  Runs the final Trinity assembly.
  """
  input:
    cmds_completed="trinity_out_dir/recursive_trinity.cmds.completed",
    samples="samples_trimmed.txt"
  output:
    transcriptome="trinity_out_dir.Trinity.fasta", #new output files live in base directory now, not trinity_out_dir.
    gene_trans_map="trinity_out_dir.Trinity.fasta.gene_trans_map" #according to trinity dev, this is a feature, not a bug
  log:
    "logs/trinity_final.log"
  conda:
    "/projects/wenglab/testtube/matthew/miniforge3/envs/transxpress-trinityutils"
  params:
    memory="256"
  threads:
    16
  shell:
    """
    Trinity --max_memory {params.memory}G --CPU {threads} --samples_file {input.samples} {config[trinity_parameters]} {config[strand_specific]} &>> {log}
    """

rule transcriptome_copy:
  """
  Copies the assembled transcriptome to the transxpress directory.
  """
  input:
    transcriptome=rules.trinity_final.output.transcriptome,
    gene_trans_map=rules.trinity_final.output.gene_trans_map,
  output:
    transcriptome="transcriptome.fasta",
    gene_trans_map="transcriptome.gene_trans_map",
    #create copies the same way older Trinity versions did, for parity, in case other dependencies on that path exist?
    redundant_transcriptome="trinity_out_dir/Trinity.fasta",
    redundant_gene_trans_map="trinity_out_dir/Trinity.fasta.gene_trans_map"
  log:
    "logs/transcriptome_copy.log"
  shell:
    """
    cp -p {input.transcriptome} {output.transcriptome} &> {log}
    cp -p {input.gene_trans_map} {output.gene_trans_map} &>> {log}
    #create copies the same way older Trinity versions did, for parity, in case other dependencies on that path exist?
    cp -p {input.transcriptome} {output.redundant_transcriptome} &> {log}
    cp -p {input.gene_trans_map} {output.redundant_gene_trans_map} &>> {log}
    """


rule trinity_stats:
  """
  Runs Trinity script to get statistics about the assembled transcriptome
  (number of transcripts, genes, GC content, EXN50).
  """
  input:
    transcriptome="transcriptome.fasta",
    expression="transcriptome_expression_isoform.tsv"
  output:
    stats="transcriptome_stats.txt",
    exN50="transcriptome_exN50.tsv",
    exN50plot="ExN50_plot.pdf"
  log:
    "logs/trinity_exN50.log"
  conda:
    "/projects/wenglab/testtube/matthew/miniforge3/envs/transxpress-trinityutils"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    TRINITY_HOME=$(python -c 'import os;import shutil;TRINITY_EXECUTABLE_PATH=shutil.which("Trinity");print(os.path.dirname(os.path.join(os.path.dirname(TRINITY_EXECUTABLE_PATH), os.readlink(TRINITY_EXECUTABLE_PATH))))')

    $TRINITY_HOME/util/TrinityStats.pl {input.transcriptome} > {output.stats} 2> {log}
    $TRINITY_HOME/util/misc/contig_ExN50_statistic.pl {input.expression} {input.transcriptome} > {output.exN50} 2>> {log}
    $TRINITY_HOME/util/misc/plot_ExN50_statistic.Rscript {output.exN50} &>> {log}
    """

rule busco:
  """
  Runs BUSCO to assess the completeness of the transcriptome.
  """
  input:
    transcriptome="transcriptome.fasta"
  output:
    out_directory=directory("busco"),
    report="busco_report.txt"
  log:
    "logs/busco.log"
  conda:
    "/projects/wenglab/testtube/matthew/miniforge3/envs/transxpress-busco"
  params:
    memory="10"
  threads:
    4
  shell:
    """
    lineage={config[lineage]} &> {log}
    if [ -z "$lineage"] &>> {log}
    then &>> {log}
      busco -m transcriptome -i {input.transcriptome} -o {output.out_directory} --auto-lineage -c {threads} &>> {log}
    else &>> {log}
      busco -m transcriptome -i {input.transcriptome} -o {output.out_directory} -l $lineage -c {threads} &>> {log}
    fi &>> {log}

    status=$?

    if [ $status -eq 0 ]
    then
      echo "BUSCO run completed successfully" &>> {log}
      cp busco/short_summary*.txt {output.report} &>> {log}
    else
      echo "BUSCO run failed" &>> {log}
      exit 1 &>> {log}
    fi
    """

rule transdecoder_longorfs:
  """
  Runs first stage of Transdecoder extracting the long open reading frames.
  """
  input:
    transcriptome="transcriptome.fasta",
  output:
    orfs="transcriptome.orfs"
  log:
    "logs/transdecoder_longorfs.log"
  conda:
    "/projects/wenglab/testtube/matthew/miniforge3/envs/transxpress-transdecoder"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    rm -rf {input.transcriptome}.transdecoder_dir &> {log}
    TransDecoder.LongOrfs -t {input.transcriptome} --output_dir transdecoder &>> {log} 
    cp -p transdecoder/{input.transcriptome}.transdecoder_dir/longest_orfs.pep {output.orfs} &>> {log}
    """


rule transdecoder_predict:
  """
  Runs second stage of Transdecoder predicting the likely coding regions based
  on Pfam and SwissProt hits.
  """
  input:
    transcriptome="transcriptome.fasta",
    pfam="annotations/pfam_orfs.out",
    blastp="annotations/sprotblastp_orfs.out"
  output:
    "transcriptome.pep"
  log:
    "logs/transdecoder_predict.log"
  conda:
    "/projects/wenglab/testtube/matthew/miniforge3/envs/transxpress-transdecoder"
  params:
    memory="10"
  threads:
    1
  shell:
    """
    TransDecoder.Predict -t {input.transcriptome} --output_dir transdecoder --retain_pfam_hits {input.pfam} --retain_blastp_hits {input.blastp} &> {log}
    cp -p transdecoder/{input.transcriptome}.transdecoder.pep {output} &>> {log}
    """

rule trinity_DE:
  """
  Runs Trinity script to perform differential expression analysis using edgeR.
  """
  input:
    samples=config["samples_file"],
    expression="kallisto.gene.counts.matrix"
  output:
    directory("edgeR_trans")
  log:
    "logs/trinity_DE.log"
  conda:
    "/projects/wenglab/testtube/matthew/miniforge3/envs/transxpress-trinityutils"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    TRINITY_HOME=$(python -c 'import os;import shutil;TRINITY_EXECUTABLE_PATH=shutil.which("Trinity");print(os.path.dirname(os.path.join(os.path.dirname(TRINITY_EXECUTABLE_PATH), os.readlink(TRINITY_EXECUTABLE_PATH))))')


    num_replicates=`awk '{{print $2}}' {input.samples} | sort | uniq | wc -l` &> {log}
    num_samples=`awk '{{print $1}}' {input.samples} | sort | uniq | wc -l` &>> {log}
    num_replicates_minus_samples=$((num_replicates - num_samples)) &>> {log}
    if [ $num_replicates_minus_samples -gt 1 ] &>> {log}
    then &>> {log}
	    $TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix {input.expression} --method edgeR --output {output} --samples_file {input.samples} &>> {log}
    else &>> {log}
	    echo "No biological replicates to run proper differential expression analysis, last-resorting to edgeR with --dispersion 0.1" &>> {log}
      $TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix {input.expression} --method edgeR --output {output} --samples_file {input.samples} --dispersion {config[dispersion]} &>> {log}
    fi &>> {log}
    """


# a more elegant way to do this is:
# seqkit seq -i {input} | seqkit replace -s -p "\*" -r "" | seqkit split -f -s 1 -O {output}
# however, there seems to be a problem in seqkit: https://github.com/shenwei356/seqkit/issues/65
checkpoint fasta_split_fasta:
  """
  Splits the transcriptome.fasta into smaller files so
  they can be processed (annotated) in parallel.
  """
  input:
    "transcriptome.fasta"
  output:
    directory("annotations/chunks_fasta")
  log:
    "logs/fasta_split.log"
  params:
    memory="2"
  threads:
    1
  run:
    os.makedirs(output[0], exist_ok=True)
    with open(input[0], "r") as input_handle:
      output_handle = None
      count = 0
      for record in Bio.SeqIO.parse(input_handle, "fasta"):
        if (count % 5000 == 0):
          if output_handle is not None:
            output_handle.close()
          fileindex = str(int(count / 5000) + 1);
          filename = os.path.join(output[0], fileindex + ".fasta")
          output_handle = open(filename, "w");
        
        # String the description, because some tools (e.g. deeploc) include it in their output
        record.description = ""
        Bio.SeqIO.write(record, output_handle, "fasta")
        count += 1
      if output_handle is not None:
        output_handle.close()
        

checkpoint fasta_split_orfs:
  """
  Splits the transcriptome.orfs into smaller files so
  they can be processed (annotated) in parallel.
  """
  input:
    "transcriptome.orfs"
  output:
    directory("annotations/chunks_orfs")
  log:
    "logs/orfs_split.log"
  params:
    memory="2"
  threads:
    1
  run:
    os.makedirs(output[0], exist_ok=True)
    with open(input[0], "r") as input_handle:
      output_handle = None
      count = 0
      for record in Bio.SeqIO.parse(input_handle, "fasta"):
        if (count % 5000 == 0):
          if output_handle is not None:
            output_handle.close()
          fileindex = str(int(count / 5000) + 1);
          filename = os.path.join(output[0], fileindex + ".orfs")
          output_handle = open(filename, "w");
        
        # String the description, because some tools (e.g. deeploc) include it in their output
        record.description = ""
        Bio.SeqIO.write(record, output_handle, "fasta")
        count += 1
      if output_handle is not None:
        output_handle.close()
        
checkpoint fasta_split_pep:
  """
  Splits the transcriptome.pep into smaller files so
  they can be processed (annotated) in parallel.
  """
  input:
    "transcriptome.pep"
  output:
    directory("annotations/chunks_pep")
  log:
    "logs/pep_split.log"
  params:
    memory="2"
  threads:
    1
  run:
    os.makedirs(output[0], exist_ok=True)
    with open(input[0], "r") as input_handle:
      output_handle = None
      count = 0
      for record in Bio.SeqIO.parse(input_handle, "fasta"):
        if (count % 5000 == 0):
          if output_handle is not None:
            output_handle.close()
          fileindex = str(int(count / 5000) + 1);
          filename = os.path.join(output[0], fileindex + ".pep")
          output_handle = open(filename, "w");
        # Remove predicted stop codons, because some annotation tools do not like them (e.g. InterProScan) 
        #if wildcards["extension"] == "pep":
        record.seq = record.seq.strip("*")
        # String the description, because some tools (e.g. deeploc) include it in their output
        record.description = ""
        Bio.SeqIO.write(record, output_handle, "fasta")
        count += 1
      if output_handle is not None:
        output_handle.close()
        

def parallel_annotation_tasks_fasta(wildcards):
  """
  Returns filenames of files which were annotated in parallel.
  """
  parallel_dir = checkpoints.fasta_split_fasta.get(**wildcards).output[0]
  job_ids = glob_wildcards(os.path.join(parallel_dir, "{index}.fasta" )).index
  completed_files = expand("annotations/{task}/{index}.out",index=job_ids, task=wildcards["task"])
  return completed_files

def parallel_annotation_tasks_orfs(wildcards):
  """
  Returns filenames of files which were annotated in parallel.
  """
  parallel_dir = checkpoints.fasta_split_orfs.get(**wildcards).output[0]
  job_ids = glob_wildcards(os.path.join(parallel_dir, "{index}.orfs")).index
  completed_files = expand("annotations/{task}/{index}.out",index=job_ids, task=wildcards["task"])
  return completed_files

def parallel_annotation_tasks_pep(wildcards):
  """
  Returns filenames of files which were annotated in parallel.
  """
  parallel_dir = checkpoints.fasta_split_pep.get(**wildcards).output[0]
  job_ids = glob_wildcards(os.path.join(parallel_dir, "{index}.pep")).index
  completed_files = expand("annotations/{task}/{index}.out",index=job_ids, task=wildcards["task"])
  return completed_files


rule annotation_merge_fasta:
  """
  Merges files that were annotated in parallel into a single file.
  """
  input:
    parallel_annotation_tasks_fasta
  output:
    "annotations/{task}_fasta.out"
  log:
    "logs/{task}_fasta_merge.log"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    cat {input} > {output} 2> {log}
    """

rule annotation_merge_orfs:
  """
  Merges files that were annotated in parallel into a single file.
  """
  input:
    parallel_annotation_tasks_orfs
  output:
    "annotations/{task}_orfs.out"
  log:
    "logs/{task}_orfs_merge.log"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    cat {input} > {output} 2> {log}
    """


rule annotation_merge_pep:
  """
  Merges files that were annotated in parallel into a single file.
  """
  input:
    parallel_annotation_tasks_pep
  output:
    "annotations/{task}_pep.out"
  log:
    "logs/{task}_pep_merge.log"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    cat {input} > {output} 2> {log}
    """


rule rfam_parallel:
  """
  Runs cmscan on Rfam database on smaller nucleotide files in parallel.
  """
  input:
    fasta="annotations/chunks_fasta/{index}.fasta",
    db="db/Rfam.cm"
  output:
    "annotations/rfam/{index}.out"
  log:
    "logs/rfam_{index}.log"
  conda:
    "/projects/wenglab/testtube/matthew/miniforge3/envs/transxpress-rfam"
  params:
    memory="8"
  threads:
    4
  shell:
    """
    cmscan -E {config[e_value_threshold]} --rfam --cpu {threads} --tblout {output} {input[db]} {input[fasta]} &> {log}
    """


rule pfam_parallel:
  """
  Runs hmmscan on Pfam database on smaller protein files in parallel.
  """
  input:
    fasta="annotations/chunks_orfs/{index}.orfs",
    db="db/Pfam-A.hmm"
  output:
    "annotations/pfam/{index}.out"
  log:
    "logs/pfam_{index}.log"
  conda:
    "/projects/wenglab/testtube/matthew/miniforge3/envs/transxpress-pfam"
  params:
    memory="2"
  threads:
    4
  shell:
    """
    # Transdecoder requires --domtblout output
    hmmscan -E {config[e_value_threshold]} --cpu {threads} --domtblout {output} {input[db]} {input[fasta]} &> {log}
    """


rule sprot_blastp_parallel:
  """
  Runs DIAMOND blastp search on SwissProt database on smaller protein files in parallel.
  """
  input:
    fasta="annotations/chunks_orfs/{index}.orfs",
    db="db/uniprot_sprot.dmnd"
  output:
    "annotations/sprotblastp/{index}.out"
  log:
    "logs/sprotblastp{index}.log"
  conda:
    "/projects/wenglab/testtube/matthew/miniforge3/envs/transxpress-blast"
  params:
    memory="4"
  threads:
    4
  shell:
    """
    diamond blastp --query {input[fasta]} --db {input[db]} --threads {threads} --evalue {config[e_value_threshold]} --max-hsps 1 --max-target-seqs 1 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle --out {output} &> {log}
    """


rule signalp_parallel:
  """
  Runs signalp on smaller protein files in parallel to predict signal peptides.
  """
  input:
    "annotations/chunks_pep/{index}.pep"
  output:
    "annotations/signalp/{index}.out"
  log:
    "logs/signalp_{index}.log"
  conda:
    "/projects/wenglab/testtube/matthew/miniforge3/envs/signalp"
  params:
    memory="8"
  threads:
    1
  shell:
    """
    mkdir -p annotations/signalp &> {log}
    signalp6 --fastafile {input} --organism eukarya --output_dir signalp_{wildcards.index} --format none --mode fast &>> {log}
    mv signalp_{wildcards.index}/prediction_results.txt {output}
    rm -r signalp_{wildcards.index}
    """


rule targetp_parallel:
  """
  Runs TargetP on smaller protein files to predict presence and type of
  targeting peptide.
  """
  input:
    "annotations/chunks_pep/{index}.pep"
  output:
    "annotations/targetp/{index}.out"
  log:
    "logs/targetp_{index}.log"
  conda:
    "/projects/wenglab/testtube/matthew/miniforge3/envs/targetp"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    targetp -fasta {input} -format short -org {config[targetp]} -prefix {wildcards.index} &> {log}
    mv {wildcards.index}_summary.targetp2 {output} &>> {log}
    """


rule kallisto:
  """
  Runs Trinity script to perform transcript expression quantification 
  using Kallisto.
  """
  input:
    samples="samples_trimmed.txt",
    transcriptome="transcriptome.fasta",
    gene_trans_map="transcriptome.gene_trans_map"
  output:
    "transcriptome_expression_isoform.tsv",
    "transcriptome_expression_gene.tsv",
    "kallisto.gene.counts.matrix"
  log:
    "logs/kallisto.log"
  conda:
    "/projects/wenglab/testtube/matthew/miniforge3/envs/transxpress-trinityutils"
  params:
    memory="16" # increased memory from 2 to 8 since it was not sufficient
  threads:
    8
  shell:
    """
    TRINITY_HOME=$(python -c 'import os;import shutil;TRINITY_EXECUTABLE_PATH=shutil.which("Trinity");print(os.path.dirname(os.path.join(os.path.dirname(TRINITY_EXECUTABLE_PATH), os.readlink(TRINITY_EXECUTABLE_PATH))))')

    assembler="{config[assembler]}"
    strand_specific="{config[strand_specific]}"
    $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts {input.transcriptome} {config[strand_specific]} --seqType fq --samples_file {input.samples} --prep_reference --thread_count {threads} --est_method kallisto --gene_trans_map {input.gene_trans_map} &> {log}
    
    $TRINITY_HOME/util/abundance_estimates_to_matrix.pl --est_method kallisto --name_sample_by_basedir --gene_trans_map {input.gene_trans_map} */abundance.tsv &>> {log}
    if [ -f kallisto.isoform.TMM.EXPR.matrix ]; then
      cp -p kallisto.isoform.TMM.EXPR.matrix {output[0]} &>> {log}
    elif [ -f kallisto.isoform.TPM.not_cross_norm ]; then
      cp -p kallisto.isoform.TPM.not_cross_norm {output[0]} &>> {log}
    else
      echo Neither kallisto.isoform.TMM.EXPR.matrix or kallisto.isoform.TPM.not_cross_norm were produced
      exit 1
    fi
    if [ -f kallisto.gene.TMM.EXPR.matrix ]; then
      cp -p kallisto.gene.TMM.EXPR.matrix {output[1]} &>> {log}
    elif [ -f kallisto.gene.TPM.not_cross_norm ]; then
      cp -p kallisto.gene.TPM.not_cross_norm {output[1]} &>> {log}
    else
      echo Neither kallisto.gene.TMM.EXPR.matrix or kallisto.gene.TPM.not_cross_norm were produced
      exit 1
    fi
    """

 
rule download_sprot:
  """
  Downloads and prepares SwissProt database.
  """
  output:
    "db/uniprot_sprot.fasta",
    "db/uniprot_sprot.dmnd"
  log:
    "logs/download_sprot.log"
  conda:
    "/projects/wenglab/testtube/matthew/miniforge3/envs/transxpress-blast"
  params:
    memory="4"
  threads:
    1
  shell:
    """
    wget --directory-prefix db "ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz" &> {log}
    gunzip db/uniprot_sprot.fasta.gz &>> {log}
    diamond makedb --in db/uniprot_sprot.fasta -d db/uniprot_sprot &>> {log}
    """


rule download_pfam:
  """
  Downloads and prepares Pfam database.
  """
  output:
    "db/Pfam-A.hmm"
  log:
    "logs/download_pfam.log"
  conda:
    "/projects/wenglab/testtube/matthew/miniforge3/envs/transxpress-pfam"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    wget --directory-prefix db "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz" &> {log}
    gunzip db/Pfam-A.hmm.gz &>> {log}
    hmmpress db/Pfam-A.hmm &>> {log}
    """


rule download_rfam:
  """
  Downloads and prepares Rfam database.
  """
  output:
    "db/Rfam.cm"
  log:
    "logs/download_rfam.log"
  conda:
    "/projects/wenglab/testtube/matthew/miniforge3/envs/transxpress-rfam"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    wget --directory-prefix db "ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz" &> {log}
    gunzip db/Rfam.cm.gz &>> {log}
    cmpress db/Rfam.cm &>> {log}    
    """


rule annotated_fasta:
  """
  Puts annotations in the headers of transcripts/proteins in the 
  .fasta/.pep transcriptome files.
  """
  input:
    transcriptome="transcriptome.fasta",
    proteome="transcriptome.pep",
    expression="transcriptome_expression_isoform.tsv",
    rfam_results="annotations/rfam_fasta.out",
    blastp_results="annotations/sprotblastp_orfs.out",
    pfam_results="annotations/pfam_orfs.out",
    signalp_results="annotations/signalp_pep.out",
    targetp_results="annotations/targetp_pep.out"
  output:
    transcriptome_annotated="transcriptome_annotated.fasta",
    proteome_annotated="transcriptome_annotated.pep",
    tpm_blast_table="transcriptome_TPM_blast.csv"
  log:
    "logs/annotated_fasta.log"
  params:
    memory="4"
  threads:
    1
  run:
    # Open log file
    with open(log[0], "w") as log_handle:
      ## Annotation map: transcript id -> description
      expression_annotations = {}
      blastp_annotations = {}
      pfam_annotations = {}
      rfam_annotations = {}
      signalp_annotations = {}
      targetp_annotations = {}
  
      ## Load kallisto results
      print ("Loading expression values from", input["expression"], file=log_handle)
      with open(input["expression"]) as input_handle:
        csv_reader = csv.reader(input_handle, delimiter='\t')
        columns = next(csv_reader)
        for row in csv_reader:
          expression_annotations[row[0]] = columns[1] + "=" + str(row[1])
          for i in range(2, len(columns)):
            expression_annotations[row[0]] += " " + columns[i] + "=" + str(row[i])

      ## Load blastp results
      print ("Loading blastp results from", input["blastp_results"], file=log_handle)
      with open(input["blastp_results"]) as input_handle:
        csv_reader = csv.reader(input_handle, delimiter="\t")
        for row in csv_reader:
          if (len(row) < 13): continue
          blastp_annotations[row[0]] = row[12] + " E=" + str(row[10])

      ## make transcript-level blastp lookup (use best/first ORF hit per transcript)
      blastp_by_transcript = {}
      for orf_id, annotation in blastp_annotations.items():
          transcript_id = re.sub(r"\.p[0-9]+$", "", orf_id)
          if transcript_id not in blastp_by_transcript:
              blastp_by_transcript[transcript_id] = annotation

      ## Load pfam results
      print ("Loading pfam predictions from", input["pfam_results"], file=log_handle)
      with open(input["pfam_results"]) as input_handle:
        for line in input_handle:
          if (line.startswith("#")): continue
          row = re.split(" +", line, 22)
          if (len(row) < 23): continue
          if row[3] not in pfam_annotations:
            pfam_annotations[row[3]] = row[1] + " " + row[22] + "E=" + str(row[6])

      ## Load rfam results
      print ("Loading rfam predictions from", input["rfam_results"], file=log_handle)
      with open(input["rfam_results"]) as input_handle:
        for line in input_handle:
          if (line.startswith("#")): continue
          row = re.split(" +", line, 17)
          if (len(row) < 18): continue
          rfam_annotations[row[2]] = row[1] + " " + row[17] + "E=" + str(row[15])
  
      ## Load signalp results
      print ("Loading signalp predictions from", input["signalp_results"], file=log_handle)
      with open(input["signalp_results"]) as input_handle:
        csv_reader = csv.reader(input_handle, delimiter="\t")
        for row in csv_reader:
          if (len(row) < 2): continue
          prediction = str(row[1])
          if prediction == "OTHER":
            signalp_annotations[row[0]] = "OTHER (No SP)"
          elif prediction == "SP":
            signalp_annotations[row[0]] = "standard secretory signal peptides (Sec/SPI)"
          elif prediction == "LIPO":
            signalp_annotations[row[0]] = "lipoprotein signal peptides (Sec/SPII)"
          elif prediction == "TAT":
            signalp_annotations[row[0]] = "Tat signal peptides (Tat/SPI)"
          elif prediction == "TATLIPO":
            signalp_annotations[row[0]] = "Tat lipoprotein signal peptides (Tat/SPII)"
          elif prediction == "PILIN":
            signalp_annotations[row[0]] = "Pilin and pilin-like signal peptides (Sec/SPIII)"
          else:
            signalp_annotations[row[0]] = str(row[1])

      ## Load targetp results
      translation =    {"SP": "Signal peptide",
                        "mTP": "Mitochondrial transit peptide",
                        "cTP": "chloroplast transit peptide",
                        "luTP": "thylakoidal lumen composite transit peptide",
                        "noTP": "no targeting peptide" }
      print("Loading targetp predictions from", input["targetp_results"], file=log_handle)
      with open(input["targetp_results"]) as input_handle:
        for line in input_handle:
          if (line.startswith("#")): continue
          row = line.split()
          if (len(row) >= 8):
            targetp_annotations[row[0]] = translation[row[1]] + ", " + " ".join(row[7:])
          elif (len(row) >= 2):
            targetp_annotations[row[0]] = translation[row[1]]
      
      ## Do the work
      print ("Annotating FASTA file", input["transcriptome"], "to", output["transcriptome_annotated"], file=log_handle)
      with open(input["transcriptome"], "r") as input_fasta_handle, open(output["transcriptome_annotated"], "w") as output_fasta_handle:
        for record in Bio.SeqIO.parse(input_fasta_handle, "fasta"):
          transcript_id = record.id
          record.description = "TPM: " + expression_annotations.get(transcript_id)
          if transcript_id in blastp_by_transcript:
            record.description += "; blastp: " + blastp_by_transcript.get(transcript_id)
          if transcript_id in rfam_annotations:
            record.description += "; rfam: " + rfam_annotations.get(transcript_id)
          Bio.SeqIO.write(record, output_fasta_handle, "fasta")
      
      print ("Annotating FASTA file", input["proteome"], "to", output["proteome_annotated"], file=log_handle)
      with open(input["proteome"], "r") as input_fasta_handle, open(output["proteome_annotated"], "w") as output_fasta_handle:
        for record in Bio.SeqIO.parse(input_fasta_handle, "fasta"):
          transcript_id = re.sub("\\.p[0-9]+$", "", record.id)
          record.description = "transdecoder: " + re.search("ORF type:([^,]+,score=[^,]+)", record.description).group(1)
          if transcript_id in expression_annotations:
            record.description += "; TPM: " + expression_annotations.get(transcript_id)
          if record.id in blastp_annotations:
            record.description += "; blastp: " + blastp_annotations.get(record.id)
          if record.id in pfam_annotations:
            record.description += "; pfam: " + pfam_annotations.get(record.id)
          if transcript_id in rfam_annotations:
            record.description += "; rfam: " + rfam_annotations.get(transcript_id)
          if record.id in signalp_annotations:
            record.description += "; signalp: " + signalp_annotations.get(record.id)
          if record.id in targetp_annotations:
            record.description += "; targetp: " + targetp_annotations.get(record.id)
          # Add sequence ID prefix from configuration
          if config["annotated_fasta_prefix"]:
            record.id = config["annotated_fasta_prefix"] + "|" + record.id
          Bio.SeqIO.write(record, output_fasta_handle, "fasta")

      print ("Generating transcriptome_TPM_blast.csv table", file=log_handle)
      with open(input["expression"], "r") as input_csv_handle, open(output["tpm_blast_table"], "w") as output_csv_handle:
        csv_reader = csv.reader(input_csv_handle, delimiter="\t")
        csv_writer = csv.writer(output_csv_handle, delimiter=",")
        csv_columns = next(csv_reader)
        csv_columns[0] = "transcript"
        for i in range(1, len(csv_columns)):
          csv_columns[i] = "TPM(" + csv_columns[i] + ")"
        csv_columns.append("blastp")
        csv_writer.writerow(csv_columns)
        for row in csv_reader:
          row.append(blastp_by_transcript.get(row[0], ""))
          csv_writer.writerow(row)


