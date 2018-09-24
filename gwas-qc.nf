#!/usr/bin/env nextflow

bedBimFamPublishDir = (params.bedbimfam == null ? "$params.output/bed-bim-fam" : params.bedbimfam)

def toPrefixTuple = { file -> tuple(file.name.take(file.name.lastIndexOf('.')), file.toRealPath()) }
def flatGroupTuple = { it -> tuple(it[0], *it[1].sort()) }

mapPeds = Channel.fromFilePairs("${params.mapped}/*.{map,ped}", flat: true)

process convertBedBimFam {
    publishDir bedBimFamPublishDir

    input:
    set fileSetId, 'file.map', 'file.ped' from mapPeds

    output:
    set fileSetId, "${fileSetId}.bed", "${fileSetId}.bim", "${fileSetId}.fam" into convertedBedBimFams

    when:
    params.mapped != null

    shell:
    '''
    plink --file file --make-bed --out !{fileSetId}
    '''
}

Channel
    .fromPath("$bedBimFamPublishDir/*.{bed,bim,fam}")
    .map(toPrefixTuple)
    .groupTuple()
    .map(flatGroupTuple) 
    .concat(convertedBedBimFams)
    .unique { it[0] }
    .into { bedBimFams; debugBedBimFams }

debugBedBimFams.subscribe {
    println it
}

process checkSexDiscordance {
    // *.hh, *.sexcheck, *.log, *.sexprobs, *.sexcheck.fail
    publishDir "$params.output/sex-discordance"

    input:
    set fileSetId, "file.bed", "file.bim", "file.fam" from bedBimFams

    output:
    set fileSetId, "${fileSetId}.sexcheck.fail", "${fileSetId}.sexcheck" into sexDiscordanceFails

    shell:
    '''
    plink --bfile file --check-sex --out !{fileSetId}
    grep PROBLEM !{fileSetId}.sexcheck > !{fileSetId}.sexprobs
    awk '{ print $1, $2 }' !{fileSetId}.sexprobs >> !{fileSetId}.sexcheck.fail
    '''
}

//     input:
//     val meta from mapPeds
//     file 'ref.fasta' from references
//     file "${meta.id}.vcf" from vcfs

//     output:
//     set val(meta), "decomposed.normalized.vcf" into decompoedNormalizedVCFs

//     when:
//     params.map-ped-in != null

//     shell:
//         "plink --file $file_set_name --make-bed --out $file_set_name",
//         "mv *.bed *.bim *.fam /pfs/out/$id",
//     '''
//     zless !{meta.id}.vcf |
//         sed 's/ID=AD,Number=./ID=AD,Number=R/' | 
//         vt decompose -s - |
//         vt normalize -r ref.fasta - > decomposed.normalized.vcf
//     '''

// }

// outsToPick=["annotated.vcf.gz", "snpEff_genes.txt", "snpEff_summary.html"]


// annotatedVCFs.subscribe {
//     println it
// }

// workflow.onComplete {
//     println "${workflow.scriptId}, ${workflow.commitId}, ${workflow.sessionId}"
//     println "${workflow.projectDir}, ${workflow.runName}, ${workflow.workDir}"
// }
