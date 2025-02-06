version 1.0

workflow NLRC4_processing {
    input {
        Array[String] vcf_files
        Array[Int] batch_numbers
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

    command {
        BATCH_START=$(( ${batch_number} * ${files_per_batch} ))
        
        # Create a file list for this batch
        printf "%s\n" "${sep='\n' all_files}" > all_files.txt
        sed -n "$((BATCH_START + 1)),$((BATCH_START + ${files_per_batch}))p" all_files.txt > batch_files.txt
        
        # Process each file in the batch
        while read vcf; do
            bcftools view -r ${chromosome}:${start}-${end} $vcf | \
            bcftools view -f PASS | \
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%GT\t%AD\t%AF\t%DP\t%F1R2\t%F2R1\n' \
            >> batch_output.csv
        done < batch_files.txt
    }

    output {
        File filtered_data = "batch_output.csv"
    }

    runtime {
        docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
        memory: "8 GB"
        cpu: "2"
        disks: "local-disk 50 SSD"
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
