nextflow.enable.dsl=2

include {FASTQC} from './modules/fastqc'
include {TRIM} from './modules/trimmomatic'
include {BOWTIE2_BUILD} from './modules/bowtie2_build'
include {BOWTIE2_ALIGN} from './modules/bowtie2_align'
include {SAMTOOLS_SORT} from './modules/samtools_sort'
include {SAMTOOLS_IDX} from './modules/samtools_idx'
include {SAMTOOLS_FLAGSTAT} from './modules/samtools_flagstat'
include {BAMCOVERAGE} from './modules/deeptools_bamcoverage'
include {MULTIQC} from './modules/multiqc'
include {MULTIBWSUMMARY} from './modules/deeptools_multibwsummary'
include {PLOTCORRELATION} from './modules/deeptools_plotcorrelation'
include {TAGDIR} from './modules/homer_maketagdir'
include {FINDPEAKS} from './modules/homer_findpeaks'
include {POS2BED} from './modules/homer_pos2bed'
include {BEDTOOLS_INTERSECT} from './modules/bedtools_intersect'
include {BEDTOOLS_REMOVE} from './modules/bedtools_remove'
include {ANNOTATE} from './modules/homer_annotatepeaks'
include {COMPUTEMATRIX} from './modules/deeptools_computematrix'
include {PLOTPROFILE} from './modules/deeptools_plotprofile'
include {FIND_MOTIFS_GENOME} from './modules/homer_findmotifsgenome'

workflow {

    
    //Here we construct the initial channels we need
    
    read_ch = Channel.fromPath(params.samplesheet)
    | splitCsv( header: true )
    | map{ row -> tuple(row.name, file(row.path)) }
    


    FASTQC(read_ch)
    TRIM(read_ch, params.adapter_fa)

    BOWTIE2_BUILD(params.genome) 
    BOWTIE2_ALIGN(TRIM.out.trimmed, BOWTIE2_BUILD.out.index)


    SAMTOOLS_SORT(BOWTIE2_ALIGN.out.bam)
    SAMTOOLS_FLAGSTAT(SAMTOOLS_SORT.out.sorted)

    BAMCOVERAGE(SAMTOOLS_SORT.out.sorted)

    multiqc_ch = FASTQC.out.zip.map{ it[1] }.flatten() \
    .mix(TRIM.out.log.map{ it[1] }) \
    .mix(SAMTOOLS_FLAGSTAT.out.flagstat.map{ it[1] })
    .collect()
    
    MULTIQC(multiqc_ch)

    bigwig_files_ch = BAMCOVERAGE.out.bw.map { it[1] }.collect()
    sample_labels_ch = BAMCOVERAGE.out.bw.map { it[0] }.collect()

    MULTIBWSUMMARY(bigwig_files_ch, sample_labels_ch)
 

    PLOTCORRELATION(MULTIBWSUMMARY.out.npz)

    TAGDIR(BOWTIE2_ALIGN.out)

    TAGDIR.out
    .map { sample_id, tagdir ->
        def parts = sample_id.split('_')
        def type = parts[0]          // "IP" or "INPUT"
        def rep  = parts[1]          // "rep1"
        tuple(rep, type, tagdir)
    }
    .groupTuple(by: 0)
    .map { rep, types, tagdirs ->
        // After groupTuple, you get separate lists for each grouped element
        def ip_idx = types.findIndexOf { it == 'IP' }
        def input_idx = types.findIndexOf { it == 'INPUT' }
        tuple(rep, tagdirs[ip_idx], tagdirs[input_idx])
    }
    .set { peakcalling_ch }

    FINDPEAKS(peakcalling_ch)

    POS2BED(FINDPEAKS.out)

    POS2BED.out
    .map { sample_id, bed_file -> bed_file }  // Extract just the bed files
    .collect()  // Gather all bed files into a list
    .map { beds -> tuple(['rep1', 'rep2'], beds) }  // Create tuple for process
    .set { intersect_ch }

    BEDTOOLS_INTERSECT(intersect_ch)

    BEDTOOLS_REMOVE(BEDTOOLS_INTERSECT.out, params.blacklist)

    ANNOTATE(BEDTOOLS_REMOVE.out, params.genome, params.gtf)

    BAMCOVERAGE.out.bw
    .filter { sample_id, bw -> sample_id.startsWith('IP') }
    .set { ip_bw_ch }

    COMPUTEMATRIX(ip_bw_ch, params.ucsc_genes)

    PLOTPROFILE(COMPUTEMATRIX.out)

    FIND_MOTIFS_GENOME(BEDTOOLS_REMOVE.out, params.genome)


}