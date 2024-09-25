nextflow.enable.dsl = 2

// log info
log.info """\
NEXTFLOW - DSL2 - GWAMA Meta - P I P E L I N E
============================================
run_as            :  ${workflow.commandLine}
run_location      :  ${launchDir}
started_at        :  ${workflow.start}
container         :  ${workflow.containerEngine}:${workflow.container}

Phenotype Lists
============================================
binary phenos     :  ${params.bin_pheno_list}
quant. phenos     :  ${params.quant_pheno_list}

Input Information
============================================
sumstats_suffix   :  ${params.sumstats_suffix}
input_col_names   :
${params.input_col_names.collect().join('\n').replace('=', '\t=\t') }

GWAMA Parameters
============================================
gwama_path        :  ${params.gwama_path}
min EAF for meta  :  ${params.min_meta_input_EAF}

Post-Processing Parameters
============================================
annotate          :  ${params.annotate}
p thresh          :  ${params.p_cutoff_summarize}

Meta-Analysis Cohort Scheme
============================================
${params.analyses.collect().join('\n').replace('=', '\t:\t')}

"""
.stripIndent()

workflow {
    cohort_pheno_sumstats = GWAMA_SETUP()
    gwama_meta_sumstats = GWAMA_META(cohort_pheno_sumstats)
}

include { BIOFILTER_POSITIONS } from './biofilter_wrapper.nf'

import groovyx.gpars.dataflow.DataflowBroadcast

DataflowBroadcast get_pheno_channel(Map params) {
    all_phenos = []
    all_phenos.addAll(params.bin_pheno_list)
    all_phenos.addAll(params.quant_pheno_list)
    pheno = Channel.fromList(all_phenos)
    return pheno
}

workflow GWAMA_SETUP {
    main:
        /*
        This is the code chunk for starting with GWAMA
        If you have summary statistics from somewhere else,
        you can pass them directly
        */

        // Make phenotype Channel
        pheno = get_pheno_channel(params)

        // Make cohort Channel by iterating over analysis groups from config
        // to get unique cohorts
        all_unique_cohorts = []
        params.analyses.each { analysis, cohort_list -> all_unique_cohorts.addAll(cohort_list) }
        all_unique_cohorts = all_unique_cohorts.unique()

        // Make cohort x pheno combinations
        cohort_pheno = Channel.fromList(all_unique_cohorts).combine(pheno)

        // Map cohort x pheno combinations to file paths
        cohort_pheno_sumstats = cohort_pheno.map { c, p ->
            new Tuple(c, p, "${launchDir}/${c}/Sumstats/${c}.${p}${params.sumstats_suffix}")
        }
    emit:
        cohort_pheno_sumstats
}

workflow GWAMA_META {
    take:
        cohort_pheno_sumstats // channel of tuples with 3 items

    main:
        // Get all of the phenotypes in one list
        pheno = get_pheno_channel(params)

        // Iterate over input channel to see which cohort/phenotype files exist
        sumstats_keep = cohort_pheno_sumstats.filter { c, p, sumstats -> file(sumstats).exists() }
        sumstats_dropped = cohort_pheno_sumstats.filter { c, p, sumstats -> !file(sumstats).exists() } \
            .map { c, p, sumstats -> new Tuple(c, p) }
        write_dropped_file(sumstats_dropped.collect(flat: false))

        // Pass each individual cohort/phenotype summary stats file to be munged
        gwama_munge_script = "${launchDir}/scripts/munge_sumstats_for_gwama.py"
        munge_output = munge_sumstats_file(sumstats_keep, gwama_munge_script)

        // Get a List of Tuples of ALL (cohort, analysis) combinations from params
        // Examples: (AFR_ALL, AFR_EUR), (EUR_ALL, AFR_EUR)
        // Examples: (AFR_M, ALL_M), (AMR_M, ALL_M), (EUR_M, ALL_M)
        all_analyses_cohorts = []
        params.analyses.each { analysis, cohort_list ->
            cohort_list.each { cohort ->
                all_analyses_cohorts.add(new Tuple(cohort, analysis))
            }
        }

        /*
        Goal: get a channel to join with the munge output channel.
            create Channel of (cohort, analysis) pairs,
            combine (cohort, analysis) pairs with phenotype, -> (cohort, analysis, pheno)
            reorder Tuples to be (cohort, pheno, analysis),
            group tuples by cohort and pheno to get (cohort, pheno, analysis_list)
        */
        all_analyses_cohorts = Channel.fromList(all_analyses_cohorts).combine(pheno).map {
            cohort, analysis, pheno ->
            new Tuple(cohort, pheno, analysis)
            }.groupTuple(by: [0, 1])

        /*
        Goal: get a Channel with (cohort, pheno, analysis, sumstats file)
            join (cohort, pheno, analysis_list) with (cohort, pheno, sumstats_file)
            to get (cohort, pheno, analysis_list, sumstats_file)
            transpose to get new Channel with (cohort, pheno, analysis, file path)
        */
        //
        analysis_cohort_files = all_analyses_cohorts.join(munge_output, by: [0, 1]).transpose()

        /*
        Take (cohort, pheno, analysis, file path),
            group by pheno and analysis to get (cohort_list, pheno, analysis, file_list),
            reorder to get new Channel with (analysis, pheno, cohort_list, file_list)
        */
        analysis_inputs = analysis_cohort_files.groupTuple(by: [1, 2]).map { cl, p, a, fl -> new Tuple(a, p, cl, fl) }

        // Call process to make GWAMA input file list
        make_infile_output = make_infile(analysis_inputs)

        /*
        Join (analysis, pheno, infile) output Channel with
            (analysis, pheno, cohort_list, file_list) to get
            new Channel with (analysis, pheno, infile, cohort_list, file_list)
        */
        gwama_input = make_infile_output.join(analysis_inputs, by: [0, 1])

        // Call GWAMA tool
        gwama_output = call_gwama(gwama_input)
        gwama_meta_output = gwama_output.map { analysis, pheno, meta, gc, snps -> new Tuple(analysis, pheno, meta, snps) }
        gwama_gc_output = gwama_output.map { analysis, pheno, meta, gc, snps -> new Tuple(analysis, pheno, gc) }

        // Add chromosome and position columns back into summary stats
        gwama_meta_sumstats = add_chr_pos_to_meta_sumstats(gwama_meta_output)

        // Plots and report post-processing
        plotting_script = "${launchDir}/scripts/plot_meta_results_manhattan.py"

        // Filter summary stats for summary tables
        (filtered_results, meta_analysis_Ns) = filter_sumstats(gwama_meta_sumstats)
        filtered_sumstats_list = filtered_results.map { analysis, pheno, sumstats -> sumstats }.collect()
        sample_sizes_list = meta_analysis_Ns.map { analysis, pheno, info -> info }.collect()
        sample_size_table = make_analysis_size_table(sample_sizes_list)

        if (params['annotate']) {
            // Use the biofilter mini-workflow to get RSIDs and nearest genes
            biofilter_input = make_biofilter_positions_input(filtered_sumstats_list)
            bf_input_channel = Channel.of('gwama_meta').combine(biofilter_input)
            biofilter_annots = BIOFILTER_POSITIONS(bf_input_channel)

            plots = plot_meta_results_with_annot(gwama_meta_sumstats.combine(biofilter_annots), plotting_script)
            make_summary_table_with_annot(filtered_sumstats_list, biofilter_annots)
        }
        else {
            // Skip biofilter, go straight to plots and tables
            plots = plot_meta_results(gwama_meta_sumstats, plotting_script)
            make_summary_table(filtered_sumstats_list)
        }
    emit:
        gwama_meta_sumstats // three-part tuple of (analysis, phenotype, sumstats file path)
        sample_size_table // csv with analysis sample sizes
}

// If any of the expected summary stats files didn't exist, write a log
process write_dropped_file {
    publishDir "${launchDir}"

    input:
        val cohort_pheno_list
    output:
        path('sumstats_dropped.txt')
    shell:
        """
        echo "${cohort_pheno_list}" \
          | sed 's|], \\[|\\n|g' \
          | sed 's|\\[\\[||g' \
          | sed 's|]]||g' \
          | sort \
          > sumstats_dropped.txt
        """
}

// Munge the summary stats files to format the input for GWAMA
process munge_sumstats_file {
    publishDir "${launchDir}/Meta/MSS/"
    memory '18GB'
    input:
        tuple val(cohort), val(pheno), path(sumstats)
        path(gwama_munge_script)
    output:
        tuple val(cohort), val(pheno), path("${cohort}.${pheno}.munged.txt")
    shell:
        """
        echo "${params.input_col_names.collect().join('\n')}" > colnames.txt
        ${params.my_python} ${gwama_munge_script} \
          --sumstats ${sumstats} \
          --maf ${params.min_meta_input_EAF} \
          --pheno ${pheno} \
          --cohort ${cohort} \
          --colnames colnames.txt \
          --traitType ${params.quant_pheno_list.contains(pheno) ? 'quant' : 'bin'}

        """
    stub:
        """
        touch ${cohort}.${pheno}.munged.txt
        """
}

// Make a list of sumstats file paths that GWAMA will use as input
process make_infile {
    publishDir "${launchDir}/Meta/${analysis}/"

    input:
        tuple val(analysis), val(pheno), val(cohort_list), path(input_files)
    output:
        tuple val(analysis), val(pheno), path("${pheno}.sumstats.in")
    shell:
        """
        echo "${input_files.join('\n')}" > ${pheno}.sumstats.in
        """
}

// Call the GWAMA tool
process call_gwama {
    publishDir "${launchDir}/Meta/${analysis}"
    memory '25GB'

    input:
        tuple val(analysis), val(pheno), path(sumstats_infile), val(cohort_list), path(sumstats_file_list)
    output:
        tuple val(analysis), val(pheno), path("${pheno}.meta.out"), \
            path("${pheno}.meta.gc.out"), path("${pheno}.meta.snp_coords.txt")
    shell:
        """
        # First save chromosome, position, and variant ID
        cut -f1-3 *.munged.txt | sort -n | uniq > ${pheno}.meta.snp_coords.txt

        ${params.gwama_path} \
          --filelist ${sumstats_infile} \
          --output ${pheno}.meta \
          --genomic_control \
          ${params.quant_pheno_list.contains(pheno) ? '-qt' : ''}
        """
    stub:
        """
        touch ${pheno}.meta.out
        touch ${pheno}.meta.gc.out
        touch ${pheno}.meta.snp_coords.txt
        """
}

// Add CHR and POS columns back to the summary stats
process add_chr_pos_to_meta_sumstats {
    publishDir "${launchDir}/Meta/Sumstats/"

    input:
        tuple val(analysis), val(pheno), path(sumstats_file), path(snp_coords)
    output:
        tuple val(analysis), val(pheno), path("${analysis}.${pheno}.meta.gz")
    script:
        """
        #! ${params.my_python}
        import pandas as pd

        coords = pd.read_table('${snp_coords}', index_col='MARKERNAME')
        results = pd.read_table('${sumstats_file}', index_col='rs_number')

        df = pd.concat([coords, results], axis=1)
        df.index.name = 'variant_id'
        df.to_csv('${analysis}.${pheno}.meta.gz', sep='\\t')
        """
    stub:
        """
        touch ${analysis}.${pheno}.meta.gz
        """
}

// Filter the sumstats for summary tables
process filter_sumstats {
    publishDir "${launchDir}/Meta/Suggestive/"

    input:
        tuple val(analysis), val(pheno), path(sumstats_file)
    output:
        tuple val(analysis), val(pheno), path("${analysis}.${pheno}.meta_filtered.txt")
        tuple val(analysis), val(pheno), path("${analysis}.${pheno}.meta_N.txt")
    script:
        """
        #! ${params.my_python}

        import pandas as pd
        import sys

        print('${pheno}', '${analysis}')

        output_file = '${analysis}.${pheno}.meta_filtered.txt'
        output_file_N = '${analysis}.${pheno}.meta_N.txt'
        max_meta_N = 0

        chunks = []
        for chunk in pd.read_table('${sumstats_file}', chunksize=4E5):
            chunk_max_N = chunk['n_samples'].max() if not pd.isnull(chunk['n_samples'].max()) else 0
            max_meta_N = max(max_meta_N, chunk_max_N)
            chunk = chunk[chunk['p-value'] <= ${params.p_cutoff_summarize}]
            chunks.append(chunk)

        df = pd.concat(chunks)
        df['PHENO'] = '${pheno}'
        df['ANALYSIS'] = '${analysis}'
        df.to_csv(output_file, sep='\\t', index=False, na_rep='NA')

        N_df = pd.Series(index=['PHENO', 'ANALYSIS', 'N_Samples'])
        N_df['PHENO'] = '${pheno}'
        N_df['ANALYSIS'] = '${analysis}'
        N_df['N_Samples'] = max_meta_N
        N_df.index.name = 'col_name'
        N_df.name = 'value'
        N_df.to_csv(output_file_N, sep='\\t')
        """
    stub:
        """
        touch ${analysis}.${pheno}.meta_filtered.txt
        touch ${analysis}.${pheno}.meta_N.txt
        """
}

// Create a file to pass to biofilter
process make_biofilter_positions_input {
    publishDir "${launchDir}/Annotations/"

    input:
        path(filtered_sumstats)
    output:
        path('gwama_meta_biofilter_input_positions.txt')
    script:
        """
        #! ${params.my_python}
        import pandas as pd

        dfs = []
        input_list = '${filtered_sumstats.join(' ')}'.split()
        for f in input_list:
            dfs.append(pd.read_table(f))
        all = pd.concat(dfs)
        keep_cols = ['CHR', 'variant_id', 'POS']
        all[keep_cols].to_csv('gwama_meta_biofilter_input_positions.txt', header=False, index=False, sep=' ')
        """
    stub:
        '''
        touch gwama_meta_biofilter_input_positions.txt
        '''
}

// Make Manhattan and QQ plots
process plot_meta_results {
    publishDir "${launchDir}/Plots/"

    input:
        tuple val(analysis), val(pheno), path(sumstats)
        path(plotting_script)
    output:
        path "${analysis}.${pheno}.{manhattan.png,qq.png,qq.csv}"
    shell:
        """
        ${params.my_python} ${plotting_script} \
          --pheno ${pheno} \
          --analysis ${analysis} \
          --sumstats ${sumstats}
        """
    stub:
        """
        touch ${analysis}.${pheno}.manhattan.png
        touch ${analysis}.${pheno}.qq.png
        touch ${analysis}.${pheno}.qq.csv
        """
}

// Make Manhattan and QQ plots with nearest gene annotations
process plot_meta_results_with_annot {
    publishDir "${launchDir}/Plots/"
    memory '25GB'

    input:
        tuple val(analysis), val(pheno), path(sumstats), val(data_nickname), path(biofilter_annots)
        path plotting_script
    output:
        path "${analysis}.${pheno}.{manhattan.png,qq.png,qq.csv}"
    shell:
        """
        ${params.my_python} ${plotting_script} \
          --pheno ${pheno} \
          --analysis ${analysis} \
          --sumstats ${sumstats} \
          --annot ${biofilter_annots}
        """
    stub:
        """
        touch ${analysis}.${pheno}.manhattan.png
        touch ${analysis}.${pheno}.qq.png
        touch ${analysis}.${pheno}.qq.csv
        """
}

// Make top hits summary table
process make_summary_table {
    publishDir "${launchDir}/Summary/"

    input:
        path all_filtered_sumstats
    output:
        path('gwama_meta_all_suggestive.csv')
    script:
        """
        #! ${params.my_python}

        import pandas as pd
        dfs = []
        input_list = '${all_filtered_sumstats.join(' ')}'.split()
        output = "gwama_meta_all_suggestive.csv"

        for f in input_list:
            dfs.append(pd.read_table(f))

        pd.concat(dfs).sort_values(by=['CHR', 'POS']).to_csv(output, index=False)
        """
    stub:
        '''
        touch gwama_meta_all_suggestive.csv
        '''
}

// Make top hits summary table with RSIDs and nearest genes
process make_summary_table_with_annot {
    publishDir "${launchDir}/Summary/"

    input:
        path all_filtered_sumstats
        tuple val(data_nickname), path(biofilter_annots)
    output:
        path('gwama_meta_all_suggestive.csv')
    script:
        """
        #! ${params.my_python}

        import pandas as pd
        dfs = []
        input_list = '${all_filtered_sumstats.join(' ')}'.split()
        output = "gwama_meta_all_suggestive.csv"

        for f in input_list:
            dfs.append(pd.read_table(f))

        all_meta = pd.concat(dfs).sort_values(by=['CHR', 'POS'])

        annot_df = pd.read_csv('${biofilter_annots}', index_col='Var_ID')
        all_meta[['Gene', 'RSID']] = annot_df.loc[all_meta['variant_id'], ['Gene', 'RSID']].values

        all_meta.to_csv(output, index=False)
        """
    stub:
        '''
        touch gwama_meta_all_suggestive.csv
        '''
}


// Combine the N Samples Files to get a Table of Sample Size
process make_analysis_size_table {
    publishDir "${launchDir}/Summary/"

    input:
        path(all_sample_sizes)
    output:
        path('gwama_meta_analysis_sizes.csv')
    script:
        """
        #! ${params.my_python}

        import pandas as pd
        dfs = []
        input_list = '${all_sample_sizes.join(' ')}'.split()
        output = "gwama_meta_analysis_sizes.csv"

        for f in input_list:
            dfs.append(pd.read_table(f, index_col='col_name', dtype=str))

        df = pd.concat(dfs, axis=1).transpose()
        df.columns = ['PHENO', 'ANALYSIS', 'N_Samples']
        df = df[['ANALYSIS', 'PHENO', 'N_Samples']]
        df['N_Samples'] = df['N_Samples'].str.replace('.0', '')
        df['N_Samples'] = df['N_Samples'].mask(df['N_Samples'].isin(['0', '-9']))
        print(df)

        df.sort_values(by=['ANALYSIS', 'PHENO']).to_csv(output, index=False, na_rep='NA')
        """
    stub:
        '''
        touch gwama_meta_analysis_sizes.csv
        '''

}
