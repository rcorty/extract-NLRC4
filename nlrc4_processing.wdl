version 1.0

workflow NLRC4_processing {
    input {
        File? vcf_file_list
        
        # Optional parameters with defaults
        Int files_per_batch = 1000
        String chr = "chr2"
        Int start_pos = 32224000
        Int end_pos = 32257000
    }
    
    # Handle input flexibility
    Array[String] vcf_files = read_lines(select_first([vcf_file_list]))
    
    # Calculate number of batches needed
    Int num_files = length(vcf_files)
    Int num_batches = ceil(num_files / files_per_batch)
    Array[Int] batch_numbers = range(num_batches)

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

    output {
        File merged_output = merge_outputs.merged_data
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
    
    Int batch_start = batch_number * files_per_batch
    Int total_files = length(all_files)
    Int actual_batch_size = if (batch_start + files_per_batch > total_files) then (total_files - batch_start) else files_per_batch
    
    command <
        set -e  # Exit on error
        set -o pipefail  # Exit if any command in a pipe fails
        
        # Create a file list for this batch
        printf "%s\n" "~{sep='\n' all_files}" > all_files.txt
        sed -n "$((~{batch_start} + 1)),$((~{batch_start} + ~{actual_batch_size}))p" all_files.txt > batch_files.txt
        
        # Add header to output file
        echo -e "CHROM\tPOS\tREF\tALT\tGT\tAD\tAF\tDP\tF1R2\tF2R1" > batch_output.csv
        
        # Process files with error handling and logging
        echo "Processing batch ~{batch_number} (files ~{
