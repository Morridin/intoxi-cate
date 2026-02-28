from lib import config, assemble_transcriptome, cluster_peptides, hmmer, blast_on_toxins, blast_on_uniprot, \
    signalp, retrieve_candidate_toxins, build_output_table, run_salmon, detect_by_structure

if __name__ == "__main__":
    # --- Validate Config ----------------------------------------------------------
    has_reads = config.get("R1") is not None
    has_transcripts = config.get("transcriptome") is not None
    has_proteome = config.get("proteome_fasta") is not None
    quant = config.get("quant", True)

    # --- 1. Proteome OR 2. Reads (+/- transcriptome)  -----------------------------
    if has_proteome and (has_reads or has_transcripts):
        raise ValueError("Choose between providing only a proteome or a transcriptome/read")
    if not (has_reads or has_transcripts or has_proteome):
        raise ValueError(
            "No entry point provided in config.yaml. Please provide either a proteome, a transcriptome or reads.")

    if has_proteome:
        if quant:
            raise ValueError("quant must be false in config.yaml if you provide proteome_fasta.")
        if config.get("R1") is not None or config.get("R2") is not None:
            raise ValueError("You must not fill in R1/R2 in config.yaml when proteome_fasta is used.")

    transcriptome = assemble_transcriptome()
    clustered_peptides = cluster_peptides(transcriptome)

    toxins_blast_result = blast_on_toxins(clustered_peptides)

    signalp_result = detect_by_structure(clustered_peptides).set_index("ID")

    toxin_candidates = retrieve_candidate_toxins(clustered_peptides, toxins_blast_result, signalp_result)

    if quant:
        salmon_result = run_salmon()
    else:
        salmon_result = None

    hmmer_result = hmmer(toxin_candidates)

    if config.get("swissprot", False):
        uniprot_blast_result = blast_on_uniprot(toxin_candidates)
    else:
        uniprot_blast_result = None

    print(
        f"\n\n"
        f"# ============================================================================ #\n"
        f"#                              Pipeline complete                               #\n"
        f"# ============================================================================ #\n"
        f"The final pipeline output can be found under {
        build_output_table(toxin_candidates, hmmer_result, toxins_blast_result.reset_index(), signalp_result, uniprot_blast_result.reset_index(), salmon_result)
        }"
    )
