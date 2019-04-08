/*
Identify immunogenic drivers
*/
WITH
selected_tumor_pairs AS
(
	SELECT
		tumor_pair_barcode,
		ps.case_barcode,
		idh_codel_subtype,
		tumor_barcode_a,
		tumor_barcode_b--,
		--row_number() OVER (PARTITION BY case_barcode ORDER BY surgical_interval_mo DESC, portion_a ASC, portion_b ASC, substring(tumor_pair_barcode from 27 for 3) ASC) AS priority
	FROM analysis.tumor_pairs ps
	LEFT JOIN analysis.blocklist b1 ON b1.aliquot_barcode = ps.tumor_barcode_a
	LEFT JOIN analysis.blocklist b2 ON b2.aliquot_barcode = ps.tumor_barcode_b
	LEFT JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
	WHERE
	--	comparison_type = 'longitudinal' AND
		sample_type_b <> 'M1' AND
		b1.coverage_exclusion = 'allow' AND b2.coverage_exclusion = 'allow' 
),
selected_aliquots AS
(
	SELECT tumor_barcode_a AS aliquot_barcode, case_barcode FROM selected_tumor_pairs --WHERE priority = 1
	UNION
	SELECT tumor_barcode_b AS aliquot_barcode, case_barcode FROM selected_tumor_pairs --WHERE priority = 1
),
selected_genes AS
(
	SELECT DISTINCT sn.gene_symbol, chrom, pos, alt, sn.variant_classification, variant_classification_priority, hgvs_p
	FROM analysis.snvs sn
	INNER JOIN analysis.dnds_fraction_sel_cv ds ON ds.gene_symbol = sn.gene_symbol AND (ds.qglobal_cv < 0.05 OR ds.gene_symbol IN ('TERT'))
	LEFT JOIN analysis.variant_classifications vc ON sn.variant_classification = vc.variant_classification
	WHERE
		(sn.gene_symbol NOT IN ('TERT','IDH1') AND variant_classification_priority IS NOT NULL) OR
		(sn.gene_symbol = 'TERT' AND sn.variant_classification = '5''Flank' AND lower(sn.pos) IN (1295228,1295250)) OR
		(sn.gene_symbol = 'IDH1' AND sn.hgvs_p IN ('p.R132C','p.R132G','p.R132H','p.R132S'))
),
selected_genes_samples AS
(
	SELECT aliquot_barcode, case_barcode, gene_symbol, chrom, lower(pos) AS start_pos, upper(pos)-1 AS end_pos, alt
	FROM selected_aliquots, selected_genes
),
gene_sample_coverage AS
(
	SELECT
		gene_symbol,
		sg.aliquot_barcode,
		sg.case_barcode,
		round(median(alt_count + ref_count)) AS median_cov
	FROM analysis.full_genotypes fgt
	INNER JOIN selected_genes_samples sg ON fgt.aliquot_barcode = sg.aliquot_barcode AND fgt.chrom = sg.chrom AND fgt.start = sg.start_pos AND fgt.end = sg.end_pos AND fgt.alt = sg.alt
	GROUP BY 1,2,3
),
gene_pair_coverage AS
(
	SELECT
		stp.tumor_pair_barcode,
		stp.case_barcode,
		c1.gene_symbol,
		c1.median_cov AS cov_a,
		c2.median_cov AS cov_b
	FROM selected_tumor_pairs stp
	LEFT JOIN gene_sample_coverage c1 ON c1.aliquot_barcode = stp.tumor_barcode_a
	LEFT JOIN gene_sample_coverage c2 ON c2.aliquot_barcode = stp.tumor_barcode_b AND c1.gene_symbol = c2.gene_symbol
),
neoag_mutation AS
(
	SELECT DISTINCT chrom, pos, alt, mutation AS aa_change
	FROM analysis.pvacseq_fraction
),
variants_by_case_and_gene AS
(
	SELECT
		gtc.gene_symbol,
		gtc.case_barcode,
		gtc.tumor_pair_barcode,
		gtc.tumor_barcode_a,
		gtc.tumor_barcode_b,
		gtc.chrom,
		gtc.pos,
		gtc.alt,
		gtc.variant_classification,
		sg.hgvs_p,
		pv.mutation AS aa_change,
		pv.netmhcpan_mt_score,
		--mutect2_call_a AS selected_call_a,
		--mutect2_call_b AS selected_call_b,
		--lag(mutect2_call_a) OVER w = mutect2_call_a OR lag(mutect2_call_b) OVER w = mutect2_call_b AS is_same_variant,
		--vaf_corrected_call_a AS selected_call_a,
		--vaf_corrected_call_b AS selected_call_b,
		--lag(vaf_corrected_call_a) OVER w = vaf_corrected_call_a OR lag(vaf_corrected_call_b) OVER w = vaf_corrected_call_b AS is_same_variant,
		(alt_count_a > 0) AS selected_call_a,
		(alt_count_b > 0) AS selected_call_b,
		lag((alt_count_a > 0)) OVER w = (alt_count_a > 0) OR lag((alt_count_b > 0)) OVER w = (alt_count_b > 0) AS is_same_variant,
		row_number() OVER w AS priority
	FROM analysis.master_genotype_comparison gtc
	INNER JOIN selected_tumor_pairs stp ON stp.tumor_pair_barcode = gtc.tumor_pair_barcode 
	INNER JOIN selected_genes sg ON sg.chrom = gtc.chrom AND sg.pos = gtc.pos AND sg.alt = gtc.alt
	LEFT JOIN analysis.pvacseq_fraction pv ON pv.tumor_pair_barcode = gtc.tumor_pair_barcode AND pv.chrom = gtc.chrom AND pv.pos = gtc.pos AND pv.alt = gtc.alt
	WHERE 
		(mutect2_call_a OR mutect2_call_b) AND 
		(ref_count_a+alt_count_a) >= 15 AND
		(ref_count_b+alt_count_b) >= 15
	WINDOW w AS (PARTITION BY gtc.gene_symbol, gtc.tumor_pair_barcode ORDER BY mutect2_call_a::integer + mutect2_call_b::integer DESC, variant_classification_priority, (ref_count_a+alt_count_a) + (ref_count_b+alt_count_b) DESC, netmhcpan_mt_score ASC)
),
variants_agg AS
(
	SELECT
		gene_symbol, case_barcode, tumor_pair_barcode, tumor_barcode_a, tumor_barcode_b, chrom, pos, alt, variant_classification, hgvs_p, aa_change, netmhcpan_mt_score,
		aa_change IS NOT NULL AS is_immunogenic,
		(CASE
		 WHEN selected_call_b AND selected_call_a THEN 'S'
		 WHEN selected_call_a AND NOT selected_call_b THEN 'P'
		 WHEN selected_call_b AND NOT selected_call_a THEN 'R' END) fraction
	FROM variants_by_case_and_gene
	WHERE priority = 1
)
SELECT va.* FROM variants_agg va
INNER JOIN analysis.silver_set ss ON ss.tumor_pair_barcode = va.tumor_pair_barcode