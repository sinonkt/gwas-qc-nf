#!/usr/bin/env nextflow

def toPrefixTuple = { file -> tuple(file.name.take(file.name.lastIndexOf('.')), file.toRealPath()) }
def flatGroupTuple = { it -> tuple(it[0], *it[1].sort()) }

mapPeds = Channel.fromFilePairs("${params.mapped}/*.{map,ped}", flat: true)
bedBimFamPublishDir = (params.bedbimfam == null ? "$params.output/bed-bim-fam" : params.bedbimfam)
sexDiscordancePublishDir = "$params.output/sex-discordance"
imissHetPublishDir = "$params.output/imiss-het"

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
    .into { 
        bedBimFamsForSexDiscordance;
        bedBimFamsForImissHet;
        debugBedBimFams
    }

debugBedBimFams.subscribe {
    println it
}

process checkSexDiscordance {
    // *.hh, *.sexcheck, *.log, *.sexprobs, *.sexcheck.fail
    publishDir sexDiscordancePublishDir

    input:
    set fileSetId, "file.bed", "file.bim", "file.fam" from bedBimFamsForSexDiscordance

    output:
    set fileSetId, "${fileSetId}.sexcheck.fail", "${fileSetId}.sexcheck" into sexDiscordanceFails

    shell:
    '''
    plink --bfile file --check-sex --out !{fileSetId}
    grep PROBLEM !{fileSetId}.sexcheck > !{fileSetId}.sexprobs
    awk '{ print $1, $2 }' !{fileSetId}.sexprobs >> !{fileSetId}.sexcheck.fail
    '''
}
// IndividualMissingAndHeterozygosityRate
process imissHet {
    
    publishDir imissHetPublishDir

    input:
    set fileSetId, "file.bed", "file.bim", "file.fam" from bedBimFamsForImissHet

    output:
    set fileSetId, "${fileSetId}.imiss", "${fileSetId}.lmiss", "${fileSetId}.het" into imissHets

    shell:
    '''
    plink --bfile file --missing --out !{fileSetId}
    plink --bfile file --het --out !{fileSetId}
    '''
}

imissHets.into { imissHetsForOutlying; imissHetsForGraph }

process imissHetOutlying {

    input:
    set fileSetId, "file.imiss", "file.lmiss", "file.het" from imissHetsForOutlying

    output:
    set fileSetId, "${fileSetId}.imisshet.fail", "${fileSetId}.gfr.fail", "${fileSetId}.het.fail" into imissHetsFails

    when:
    params.gfr != null && params.hetrsdinterval != null

    shell:
    '''
    awk -v gfr=!{params.gfr} 'NR > 1 && $6 > gfr { printf \"%s\t%s\\n\", $1, $2 }' file.imiss > !{fileSetId}.gfr.fail
    mu=$(awk 'NR > 1 {acc_h_rate+=(($5 - $3)/$5)} END {printf (acc_h_rate/(NR-1))}' file.het)
    sd=$(awk -v mu=$mu 'NR > 1 {h_rate=(($5 - $3)/$5); acc_squared_diff+=(h_rate-mu)^2} END {printf sqrt(acc_squared_diff/(NR-1))}' file.het)
    awk -v mu=$mu -v sd=$sd -v sd_interval=!{params.hetrsdinterval} 'NR > 1 { h_rate=(($5 - $3)/$5) } NR > 1 && (h_rate < (mu - sd_interval * sd) || h_rate > (mu + sd_interval * sd)) { printf \"%s\t%s\\n\", $1, $2 }' file.het > !{fileSetId}.het.fail
    cat *.fail | sort -k1 | uniq > !{fileSetId}.imisshet.fail
    '''
}

process imissHetGraph {

    input:
    set fileSetId, "file.imiss", "file.lmiss", "file.het" from imissHetsForGraph

    output:
    stdout graphOut

    shell:
    '''
    echo 'test'
    '''
}

imissHetsFails.subscribe {
    println it
}

//   "script=\"/usr/local/src/imiss-vs-het.Rscript\"",
//   "while IFS=, read -r cid gfr sd_interval; do",
//   "args=( $script '-i' $file_set_name )",
//   "([ \"$sd_interval\" != \"NULL\" ]) && args+=( '-s' $sd_interval )",
//   "([ \"$gfr\" != \"NULL\" ]) && args+=( '-f' $gfr )",
//   "args+=( '-o' \"/pfs/out/${id}/${cid}.pdf\" )",
//   "echo \"${args[@]}\"",
//   "Rscript \"${args[@]}\"",

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
