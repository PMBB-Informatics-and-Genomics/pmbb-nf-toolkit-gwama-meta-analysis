
workflow {
    pheno_list = ['T2D', 'AAA']
    pheno = Channel.fromList(pheno_list)

    analyses = ['AFR_EUR': ['AFR', 'EUR'],
                'ALL': ['AFR', 'EUR', 'EAS', 'AMR'],
                'Leave_AFR_Out': ['EUR', 'EAS', 'AMR']]

    println(pheno_list)
    println(analyses)
    println('\n\n')
    
    analyses_keys = analyses.keySet().collect()
    analysis = Channel.fromList(analyses_keys)

    pheno_analysis = pheno.combine(analysis)

    make_infile_output = make_input_file_list(pheno_analysis)
    gwama_output = call_gwama(make_infile_output)

    plots = plot_meta_results(gwama_output)
    filtered_results = filter_sumstats(gwama_output)
    make_summary_table(filtered_results.collect())

}

process make_input_file_list {
    publishDir "${launchDir}/Meta/${analysis}/"

    input:
        tuple val(pheno), val(analysis)
    output:
        tuple val(pheno), val(analysis), path("${pheno}.sumstats.in")
    stub:
        """
        touch ${pheno}.sumstats.in
        """
}

process call_gwama {
    publishDir "${launchDir}/Meta/${analysis}"

    input:
        tuple val(pheno), val(analysis), path(sumstats_infile)
    output:
        tuple val(pheno), val(analysis), path("${pheno}.meta.txt")
    stub:
        """
        touch ${pheno}.meta.txt
        """
}

process plot_meta_results {
    publishDir "${launchDir}/Plots/"

    input:
        tuple val(pheno), val(analysis), path("${pheno}.meta.txt")
    output:
        path("{manhattan,qq}.${analysis}.${pheno}.meta.png")
    stub:
        """
        touch manhattan.${analysis}.${pheno}.meta.png
        touch qq.${analysis}.${pheno}.meta.png
        """
}

process filter_sumstats {
    publishDir "${launchDir}/Meta/Suggestive/"

    input:
        tuple val(pheno), val(analysis), path("${pheno}.meta.txt")
    output:
        tuple val(pheno), val(analysis), path("${analysis}.${pheno}.meta_filtered.txt")
    stub:
        """
        touch ${analysis}.${pheno}.meta_filtered.txt
        """
}

process make_summary_table {
    publishDir "${launchDir}/Summary/"

    input:
        val all_filtered_sumstats
    output:
        'meta_all_suggestive.csv'
    stub:
        """
        touch meta_all_suggestive.csv
        """
}