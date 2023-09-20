def get_signalp_splits(wildcards):
    import os
    checkpoint_output = checkpoints.split_fasta.get(**wildcards).output[0]
    return expand('./split_files/{i}.fasta',
           i=glob_wildcards(os.path.join(checkpoint_output, '{i}.fasta')).i)



#configfile: "config.yaml"



rule trim_reads:
#   Description: trims the provided raw reads
#   todo: switch trimmomatic for fastp?
#   todo: implement autodetection of compositional bias trimming?
#   todo: do we provide the adapter list? switching to fastp would provide automatic adapter identification
    input:
        r1 = config['R1'],
        r2 = config['R2'],
        adapters = config['adapters'],
    output:
        r1 = "trimmed_reads/" + config['R1'].split("/")[-1],
        r1_unpaired = "trimmed_reads/unpaired." + config['R1'].split("/")[-1],
        r2 = "trimmed_reads/" + config['R2'].split("/")[-1],
        r2_unpaired = "trimmed_reads/unpaired." + config['R2'].split("/")[-1],
    threads: config['threads']
    shell:
        """
        mkdir -p trimmed_reads
        trimmomatic PE -threads {threads} {input.r1} {input.r2} {output.r1} {output.r1_unpaired} {output.r2} {output.r2_unpaired} ILLUMINACLIP:{input.adapters}:2:40:15 LEADING:15 TRAILING:15 MINLEN:25 SLIDINGWINDOW:4:15"
        """

rule assemble_transcriptome:
#   Description: Assembles a transcriptome if it is not provided. Uses Trinity 
#   In this case sequencing reads MUST be provided in the config.
    input:
        r1 = rules.trim_reads.output.r1,
        r2 = rules.trim_reads.output.r2,
    output:
        assembly = "trinity_out_dir/Trinity.fasta"
    params:
        memory = str(config['memory']) + "G"
    threads: config['threads']
    shell:
        """
        Trinity --seqType fq --left {rules.trim_reads.output.r1} --right {rules.trim_reads.output.r2} --CPU {threads} --max_memory {params.memory} --KMER_SIZE 31
        """
    

rule build_contaminants_database:
#   Description: builds blast database for the removal of contaminants   
#   todo: make this optional like in the original code
    input:
        fasta_db = config['contaminants']
    output:
        blast_db = config['contaminants'] + ".nin"
    shell:
        """
        makeblastdb -dbtype nucl -in {input.fasta_db}
        """

rule blast_on_contaminants:
#   Description: performs the actual blast of the contigs against the contaminants database
    input:
        blast_db = rules.build_contaminants_database.output.blast_db,
        blast_db_alias =  config['contaminants'],
        contigs = config['transcriptome'] if 'transcriptome' in config else rules.assemble_transcriptome.output.assembly,
    output:
        blast_result = config['basename'] + ".blastsnuc.out"
    params:
        evalue = config['contamination_evalue']
    threads: config['threads']
    shell:
        """
        blastn -db {input.blast_db_alias} -query {input.contigs} -out {output.blast_result} -outfmt 6 -evalue {params.evalue} -max_target_seqs 1 -num_threads {threads}
        """

rule filter_contaminants:
#   Description: performs the actual filtering
    input: 
        blast_result = rules.blast_on_contaminants.output.blast_result,
        contigs = config['transcriptome'] if 'transcriptome' in config else rules.assemble_transcriptome.output.assembly,
    output:
        filtered_contigs = config['basename'] + ".filtered.fasta"
    run:
        from Bio import SeqIO
        records = []
        infile = open(input.blast_result, 'r')
        for line in infile:
            line = line.rstrip()
            if line[0] != '#':
                blast = line.split()                                        
                records.append(blast[0]) # we recover the ID of the significan hits
        infile.close()
        recordIter = SeqIO.parse(open(input.contigs), "fasta")
        with open(output.filtered_contigs, "w") as handle:
            for rec in recordIter:
                if rec.id not in records:
                    SeqIO.write(rec, handle, "fasta")

rule detect_orfs:
#   Description: finds complete orfs within the input nucleotide sequences. 
#   i'm testing this with orfipy instead of orffinder to leverage multithreading
    input:
        nucleotide_sequences = rules.filter_contaminants.output.filtered_contigs
    output:
        aa_sequences = config['basename'] + ".faa"
    params:
        minlen = "66" if "minlen" not in config else config['minlen'],
        maxlen = "30000000" if "maxlen" not in config else config['maxlen']
    threads: config['threads']
    shell:
        """
        orfipy --procs {threads} --start ATG --pep {output.aa_sequences} --min {params.minlen} --max {params.maxlen} {input.nucleotide_sequences} --outdir .
        """

rule cluster_peptides:
#   Description: runs cd-hit on predicted peptide to remove excess redundancy
    input:
        aa_sequences = rules.detect_orfs.output.aa_sequences
    output:
        filtered_aa_sequences = config['basename'] + ".filtered.faa" #todo might want to change the basename to a param also in the other cd-hit rule if we decide on keeping it
    params:
        threshold = config['clustering_threshold'], #todo : we might want to use separate thresholds if we're going to run cd-hit on transcripts and peptides
        memory = str(int(config['memory'])*1000),
        basename = config['basename'] + ".filtered"
    threads: config['threads']
    shell:
        """
        cd-hit -i {input.aa_sequences} -o {output.filtered_aa_sequences} -c {params.threshold} -M {params.memory} -T {threads} 
        """

rule trim_peptides:
#   Description: this rule trims all peptides to only the first 50 aminoacids, as they are the only useful part for signalp. This step improves load time.
    input:
        aa_sequences = rules.cluster_peptides.output.filtered_aa_sequences
    output:
        trimmed_sequences = config['basename'] + ".trimmed.faa",
    threads: 
        config['threads']
    run:
        from Bio import SeqIO
        import subprocess
        with open(output.trimmed_sequences, "w") as outfile:
            for seq in SeqIO.parse(input.aa_sequences, "fasta"):
                outfile.write(f">{seq.id}\n{seq.seq[:31]}\n")


checkpoint split_fasta:
    input:
        fasta_file = rules.trim_peptides.output.trimmed_sequences
    output:
        split_dir = directory("split_files"),
        marker = "split_fasta.done"
    params:
        chunk_size = 9000 # using 9000 instead of 50000 for usability in normal desktop/laptop pcs. May be user defined.
    run:
        from Bio import SeqIO
        import os
        def batch_iterator(iterator, batch_size):
            """Returns lists of length batch_size.

            This is a generator function, and it returns lists of the
            entries from the supplied iterator.  Each list will have
            batch_size entries, although the final list may be shorter.
            
            src: https://biopython.org/wiki/Split_large_file
            """
            entry = True  # Make sure we loop once
            while entry:
                batch = []
                while len(batch) < batch_size:
                    try:
                        entry = next(iterator)
                    except StopIteration:
                        entry = None
                    if entry is None:
                        break
                    batch.append(entry)
                if batch:
                    yield batch
        # Open the large fasta file and use batch_iterator to split the file into batches of params.chunk_size sequences.
        os.makedirs(output.split_dir, exist_ok=True)
        record_iter = SeqIO.parse(open(input.fasta_file), "fasta")
        for i, batch in enumerate(batch_iterator(record_iter, params.chunk_size)):
            # Write the current batch to a split fasta file.
            output_file = f"{output.split_dir}/{i + 1}.fasta"
            handle = open(output_file, "w")
            SeqIO.write(batch, handle, "fasta")
            handle.close()
        with open(output.marker, "w") as f: # touches the file. Might come useful for storing information later
            f.write('')

import glob

def get_fasta_splits(wildcards):
    chkpt_output = checkpoints.split_fasta.get(**wildcards).output[0]
    return glob.glob(f"{chkpt_output}/*.fasta")

rule run_signalp:
    input: 
        fasta_file = "split_files/{file_id}.fasta",
        marker = "split_fasta.done"
    output:
        outfile = "split_sigp/{file_id}_summary.signalp5"
    params:
        prefix = "split_sigp/{file_id}"
    shell:
        """
        mkdir -p split_sigp
        signalp -batch 5000 -fasta {input.fasta_file} -org euk -format short -verbose -prefix {params.prefix}
        """

rule filter_signalp_outputs:
#   Description: this rule filters the output of the multiple signalp runs and extracts only those with a probability of signal peptide greater than a threshold. Only one file should be produced from the multiple signalp files. Two outputs are expected: a filtered (or not?) table with the signalp results and a filtered fasta of only those peptides with a signal
    input:
        files = expand("{file_id}", file_id = glob.glob("split_sigp/*.signalp5"))
    output:
        "filtered_sigp.tsv"
    params:
        threshold = 0.8 # todo: user defined
    shell:
        """
        cat {input.files} | sed '/^#/d' | awk '$3 > {params.threshold}' > {output}
        """

rule extract_secreted_peptides:
    input:
        signalp_result = rules.filter_signalp_outputs.output,
        fasta_file = rules.cluster_peptides.output.filtered_aa_sequences
    output:
        secreted_peptides = "secreted_peptides.fasta",
        non_secreted_peptides = "non_secreted_peptides.fasta"
    run:
        from Bio import SeqIO
        with open(str(input.signalp_result)) as infile:
            records = []
            for line in infile:
                records.append(line.rstrip().split("\t")[0])
        with open(output.non_secreted_peptides, "w") as n_outfile:
            with open(output.secreted_peptides, "w") as outfile:
                for seq in SeqIO.parse(input.fasta_file, "fasta"):
                    if seq.id in records:
                        SeqIO.write(seq, outfile, "fasta")
                    else:
                        SeqIO.write(seq, n_outfile, "fasta")



rule run_phobius: #todo: remember to inform the user about the installation procedure. I added a dependency in the conda env with a convenient installation script
#   Description: runs phobius
    input: 
        rules.extract_secreted_peptides.output.secreted_peptides
    output:
        table = "phobius_predictions.tsv"
    shell:
        """
        phobius -short {input} | sed 's/\s\+/\t/g' | awk '$2 == 0' > {output.table}
        """

rule extract_non_TM_peptides:
#   Description: extracts non-TM peptides from the phobius output
    input:
        phobius_result = rules.run_phobius.output.table,
        fasta_file = rules.extract_secreted_peptides.output.secreted_peptides
    output:
        non_TM_peptides = "non_TM_peptides.fasta"
    run:
        from Bio import SeqIO
        with open(str(input.phobius_result)) as infile:
            records = []
            for line in infile:
                records.append(line.rstrip().split("\t")[0])
        with open(output.non_TM_peptides, "w") as outfile:
            for seq in SeqIO.parse(input.fasta_file, "fasta"):
                if seq.id in records:
                    SeqIO.write(seq, outfile, "fasta")


rule build_toxin_blast_db: #todo: do we switch to diamond for max speed?
#   Description: builds a blast database for toxin prediction
    input:
        db = config['toxin_db']
    output:
        outfile = config['toxin_db'] + ".nin",
    shell:
        """
        makeblastdb -dbtype prot -in {input.db} 
        """

rule blast_on_toxins:
#   Description: runs blastp against the toxin database. The query are the peptides without any signal sequence. The rule runs blast and extracts the fasta at the same time. might be split in two rules for easier management.
    input:
        fasta_file = rules.extract_non_TM_peptides.output.non_TM_peptides,
        db_file = rules.build_toxin_blast_db.output.outfile,
        blast_db_alias = config['toxin_db'],
    output:
        blast_result = "toxin_blast_results.tsv",
        hits_fasta = config['basename'] + '_toxins_by_similarity.fasta'
    params:
        evalue = config['toxins_evalue'] if 'toxins_evalue' in config else "1e-10"
    threads: 
        config['threads']
    run:
        import subprocess
        from Bio import SeqIO
        command_line = "blastp -query {input.fasta_file} -evalue {params.evalue} -max_target_seqs 1 -threads {threads} -db {input.blast_db_alias} -outfmt 6 -out {output.blast_result}"
        subprocess.run(command_line, shell=True)
        with open(str(output.blast_result)) as infile:
            records = []
            for line in infile:
                records.append(line.rstrip().split("\t")[0])
        with open(output.hits_fasta, "w") as outfile:
            for seq in SeqIO.parse(input.fasta_file, "fasta"):
                if seq.id in records:
                    SeqIO.write(seq, outfile, "fasta")


#todo: finish implementation of this. 
rule retrieve_candidate_toxins:
#   Description: this rule just creates a fasta from the positive hits in the toxin similarity and structure searches. 
    input:
        structure_based = rules.extract_non_TM_peptides.output.non_TM_peptides,
        similarity_based = rules.blast_on_toxins.output.hits_fasta
    output:
        config["basename"] + "_candidate_toxins.fasta"
    shell:
        """
        cat {input.structure_based} {input.similarity_based} > {output}
        """

rule download_pfam:
#   Description: downloads pfam database. I'd like to leave it this way as the database is fairly small and will be downloaded in parallel with slow steps. 
    output:
        pfam_db = "Pfam-A.hmm"
    shell:
        """
        wget -c https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz && gunzip Pfam-A.hmm.gz
        """

rule run_hmmer:
#   Description: runs hmmer against the pfam database. 
    input:
        pfam_db = rules.download_pfam.output.pfam_db
    output:
        tblout = config['basename'] + ".tblout",
        domtblout = config['basename'] + ".domtblout"
    params:
        evalue = config['pfam_evalue'] if 'pfam_evalue' in config else "1e-5"
    threads: 
        config['threads']
    shell:
        """
        hmmsearch --cut_ga --cpu {threads} --domtblout {output.domtblout} --tblout {output.tblout} {input.pfam_db} {input.fasta_file} 
        """

rule parse_hmmsearch_output:
#   Description: parses and aggregates the hmmer output, uses the domtblout file
    input: 
        domtblout = rules.run_hmmer.output.domtblout
    output:
        filtered_table = config['basename'] + ".domtblout.tsv"
    run:
        df_domtblout = pandas.read_csv("{input.domtblout}", comment="#", delim_whitespace=True, names=["target name","accession_t","tlen","query name","accession_Q","qlen","E-value","score_1","bias_1","#","of","c-Evalue","i-Evalue","score_2","bias_2","from_1","to_1","from_2","to_2","from_3","to_3","acc","description of target"])
        aggregated_domains = df_domtblout.groupby('target name')['query name'].apply(list).reset_index()
        aggregated_domains['query name'] = aggregated_domains['query name'].apply(lambda x: list(set(x)))
        aggregated_domains.to_csv("{output.filtered_table}", sep="\t", index=False)

rule run_wolfpsort:
#   Description: runs wolfpsort on secreted peptides inferred by signalp 
    input:
        rules.extract_secreted_peptides.output.secreted_peptides
    output:
        config['basename'] + "_secreted_wolfpsort_prediction.tsv"
    params:
        wps_path = config['wolfPsort_path'],
        awk = "awk '{print $1\"\t\"$2}'"
    shell:
        """
        {params.wps_path} animal < {input} | {params.awk} > {output}
        """




### conditional rules


if config.get("quant"): # only run this rule if the quant option is active
    rule run_salmon:
    #   Description: run quantification on the entire transcriptome
        input: 
            transcriptome = config['transcriptome'],
            r1 = rules.trim_reads.output.r1,
            r2 = rules.trim_reads.output.r2
        output:
            quant_dir = directory(config['basename'] + "_quant"),
            quantification = config['basename'] + "_quant/quant.sf"
        threads:
            config['threads']
        params:
            index = "{input.transcriptome}.idx"
        shell:
            """
            salmon index -t {input.transcriptome} -i {params.index} -p {threads}
            salmon quant -i {params.index} -l A -1 {input.r1} -2 {input.r2} --validateMappings -o {output.quant_dir}
            """


if config.get("blast_uniprot"):
    rule download_uniprot:
    #   Description: downloads the uniprot fasta
        output:
            db_dir = directory("databases"),
            database = "databases/uniprot.fasta.gz"
        shell:
            """
            curl https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz -o {output.database}
            """
    
    rule make_uniprot_blast_database:
    #   Description: builds a blast database from the uniprot fasta
        input:
            fasta_file = rules.download_uniprot.output.database
        output:
            db_file = "databases/uniprot_blast_db"
        shell:
            """
            makeblastdb -in {input.fasta_file} -dbtype prot -out {output.db_file}
            """

    rule blast_on_uniprot:
    #   Description: run blast against the uniprot database, return only the best hit
        input:
            fasta_file = rules.cluster_peptides.output.filtered_aa_sequences,
            db_file = rules.make_uniprot_blast_database.output.db_file
        output:
            blast_result = config['basename'] + "_uniprot_blast_results.tsv"
        params:
            evalue = config['uniprot_evalue'] if 'uniprot_evalue' in config else "1e-10"
        threads: 
            config['threads']
        shell:
            """
            blastp -query {input.fasta_file} -evalue {params.evalue} -max_target_seqs 1 -threads {threads} -db {input.db_file} -outfmt 6 -out {output.blast_result}
            """


# TODO: follow this comment for the rule that will wraps everything up and create the final table. -> Also, in my opinion these peptides should be marked with a warning flag in the output, specifying which issue affects them (e.g. “this peptide lacks a signal peptide”, “this peptide contains a transmembrane domain”, etc.)


#todo: try to run signalp during the split rule to avoid problems. issue: if the process is interrupted abnormally during the run the rule is almost certain to misbehave and rerun the whole thing

#todo: rule run_signalp: # requires some testing. 
#todo: test with stripped sequences. This means that all sequences are preprocessed to be cut to a fixed length that would contain a signal peptide (like 50 or so). this might save memory and improve time. Moreover, we could try deduplicating these cut sequences and rereplicate afterwards to avoid predicting the same signal over and over. 


rule all: #todo: there is a bug that makes this rule run even if no split file is done. might solve by adding a checkpoint file at the end of the split rule
    input: 
        rules.extract_secreted_peptides.output,
        rules.blast_on_toxins.output,
        rules.run_wolfpsort.output,
        rules.parse_hmmsearch_output.output,
#        split_done = "split_fasta.done",
#        signalp_files = expand("split_sigp/{f}_summary.signalp5", f=[i.split("/")[1].split(".")[0] for i in  glob.glob("split_files/*.fasta")])
