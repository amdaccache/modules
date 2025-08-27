process PARABRICKS_STARFUSION {
    tag "$meta.id"
    label 'process_high'
    label 'process_gpu'
    stageInMode 'copy'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.5.1-1"

    input:
    tuple val(meta), path(chimeric_junction)
    tuple val(meta2), path(genome_lib_dir)

    output:
    // tuple val(meta), path("*.bam")                  , emit: bam              , optional:true
    // tuple val(meta), path("*.bai")                  , emit: bai              , optional:true
    // tuple val(meta), path("${prefix}_genome_lib_build_dir"), emit: reference
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.15.1' // WARN: This is the actual version of the STAR-FUSION, but version information of tool is not updated and prints '1.15.0'
    
    def num_gpus = task.accelerator ? "--num-gpus $task.accelerator.request" : ''
    """
    pbrun \\
        --chimeric-junction ${chimeric_junction} \\
        --genome-lib-dir ${genome_lib_dir} \\
        --output-dir $prefix \\
        $num_gpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunzip: \$(echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//')
        hmmer: \$(echo \$(hmmpress -h | grep HMMER | sed 's/# HMMER //' | sed 's/ .*//' 2>&1))
        STAR-Fusion: $VERSION
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.15.1' // WARN: This is the actual version of the STAR-FUSION, but version information of tool is not updated and prints '1.15.0'
    """
    mkdir -p ${prefix}_genome_lib_build_dir

    touch ${prefix}_genome_lib_build_dir/AnnotFilterRule.pm
    echo | gzip > ${prefix}_genome_lib_build_dir/blast_pairs.dat.gz
    touch ${prefix}_genome_lib_build_dir/blast_pairs.idx

    mkdir -p ${prefix}_genome_lib_build_dir/__chkpts
    touch ${prefix}_genome_lib_build_dir/__chkpts/annotfiltrule_cp.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/blast_pairs.idx.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/cp_gene_blast_pairs.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/cp_pfam_dat.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/cp_ref_annot_cdna.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/fusion_annot_lib.cp.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/_fusion_annot_lib.idx.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/index_pfam_hits.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/index_ref_annot_cdna.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/makeblastdb.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/mm2_genome_idx.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/mm2.splice_bed.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/_prot_info_db.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/ref_annot.cdsplus.dfam_masked.fa.cp.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/ref_annot.cdsplus.dfam_masked.fa.idx.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/ref_annot.gtf.gene_spans.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/ref_annot.gtf.mini.sortu.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/ref_annot.gtf.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/ref_genome_fai.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/ref_genome.fa.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/trans.blast.dat.cp.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/trans.blast.dat.index.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/validate_ctat_genome_lib.ok

    echo | gzip > ${prefix}_genome_lib_build_dir/fusion_annot_lib.gz
    touch ${prefix}_genome_lib_build_dir/fusion_annot_lib.idx
    touch ${prefix}_genome_lib_build_dir/pfam_domains.dbm
    echo | gzip > ${prefix}_genome_lib_build_dir/PFAM.domtblout.dat.gz

    touch ${prefix}_genome_lib_build_dir/ref_annot.cdna.fa
    touch ${prefix}_genome_lib_build_dir/ref_annot.cdna.fa.idx
    touch ${prefix}_genome_lib_build_dir/ref_annot.cds
    touch ${prefix}_genome_lib_build_dir/ref_annot.cdsplus.fa
    touch ${prefix}_genome_lib_build_dir/ref_annot.cdsplus.fa.idx
    touch ${prefix}_genome_lib_build_dir/ref_annot.gtf
    touch ${prefix}_genome_lib_build_dir/ref_annot.gtf.gene_spans
    touch ${prefix}_genome_lib_build_dir/ref_annot.gtf.mini.sortu
    touch ${prefix}_genome_lib_build_dir/ref_annot.gtf.mm2.splice.bed
    touch ${prefix}_genome_lib_build_dir/ref_annot.pep
    touch ${prefix}_genome_lib_build_dir/ref_annot.prot_info.dbm

    touch ${prefix}_genome_lib_build_dir/ref_genome.fa
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.fai
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.mm2
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.ndb
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.nhr
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.nin
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.njs
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.not
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.nsq
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.ntf
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.nto

    mkdir -p ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/build.ok
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/chrLength.txt
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/chrNameLength.txt
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/chrName.txt
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/chrStart.txt
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/exonGeTrInfo.tab
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/exonInfo.tab
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/geneInfo.tab
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/Genome
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/genomeParameters.txt
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/Log.out
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/SA
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/SAindex
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/sjdbInfo.txt
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/sjdbList.fromGTF.out.tab
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/sjdbList.out.tab
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/transcriptInfo.tab

    touch ${prefix}_genome_lib_build_dir/trans.blast.align_coords.align_coords.dat
    touch ${prefix}_genome_lib_build_dir/trans.blast.align_coords.align_coords.dbm
    echo | gzip > ${prefix}_genome_lib_build_dir/trans.blast.dat.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunzip: \$(echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//')
        hmmer: \$(echo \$(hmmpress -h | grep HMMER | sed 's/# HMMER //' | sed 's/ .*//' 2>&1))
        STAR-Fusion: $VERSION
    END_VERSIONS
    """

}
