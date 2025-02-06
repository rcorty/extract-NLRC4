version 1.0

workflow NLRC4_processing {
    input {
        Array[String] vcf_files
        Array[Int] batch_numbers  # We'll provide this as range(0, 30) for 1000 files per batch
        String chr = "chr2"
        Int start_pos = 32224000
        Int end_pos = 32257000
        Int files_per_batch = 1000
    }

    scatter (batch_num in batch_numbers) {
        call process_batch {
            input:
                all_files = vcf_files,
                batch_number = batch_num,
                files_per_batch = files_per_batch,
                chromosome = chr,
                start = start_pos,
                end = end_pos
        }
    }

    call merge_outputs {
        input:
            input_files = process_batch.filtered_data
    }
}

task process_batch {
    input {
        Array[String] all_files
        Int batch_number
        Int files_per_batch
        String chromosome
        Int start
        Int end
    }

    command <
        start_idx=$((~{batch_number} * ~{files_per_batch}))
        files=(~{sep=' ' all_files})
        
        # Process each file in the batch
        for ((i=0; i<~{files_per_batch}; i++)); do
            idx=$((start_idx + i))
            if [ $idx -lt ${#files[@]} ]; then
                vcf=${files[$idx]}
                bcftools view -r ~{chromosome}:~{start}-~{end} $vcf | \
                bcftools view -f PASS | \
                bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%GT\t%AD\t%AF\t%DP\t%F1R2\t%F2R1\n' \
                >> batch_output.csv
            fi
        done
    >>>

    output {
        File filtered_data = "batch_output.csv"
    }

    runtime {
        docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
        memory: "8 GB"  # Increased for batch processing
        cpu: "2"        # Increased for batch processing
        disks: "local-disk 50 SSD"  # Increased for batch processing
    }
}

task merge_outputs {
    input {
        Array[File] input_files
    }

    command {
        cat ${sep=' ' input_files} > merged_output.csv
    }

    output {
        File merged_data = "merged_output.csv"
    }

    runtime {
        docker: "ubuntu:latest"
        memory: "4 GB"
        cpu: "1"
        disks: "local-disk 50 SSD"
    }
}
