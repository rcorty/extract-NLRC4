version 1.0

workflow NLRC4_processing {
    input {
        Array[File] vcf_files
        Int batch_size = 1000
        String chr = "chr2"
        Int start_pos = 32224000
        Int end_pos = 32257000
    }

    # Calculate number of batches needed
    Int num_files = length(vcf_files)
    Int num_batches = ceil(num_files / batch_size)

    # Create array of batch indices
    scatter (i in range(num_batches)) {
        # Calculate start and end indices for this batch
        Int start_idx = i * batch_size
        Int end_idx = min(start_idx + batch_size, num_files)
        Array[File] batch = slice(vcf_files, start_idx, end_idx)

        call process_vcf_batch {
            input:
                vcf_files = batch,
                chromosome = chr,
                start = start_pos,
                end = end_pos
        }
    }
}

task process_vcf_batch {
    input {
        Array[File] vcf_files
        String chromosome
        Int start
        Int end
    }

    command {
        for vcf in ${sep=' ' vcf_files}; do
            bcftools view -r ${chromosome}:${start}-${end} $vcf | \
            bcftools view -f PASS | \
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%GT\t%AD\t%AF\t%DP\t%F1R2\t%F2R1\n' \
            >> batch_output.csv
        done
    }

    output {
        File filtered_data = "batch_output.csv"
    }

    runtime {
        docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
        memory: "4 GB"
        cpu: "1"
    }
}
