# A wrapper script for analysis of relabeling incidence wthin the analyzed reads 
# When running the script, pass '1' to the script to execute the analysis. Run with empty argument to check that all input/output parameters are correct and the script does not contain any execution errors. 

params_1="$PWD/config_files/params_1.sh"
params_2="$PWD/config_files/samples_table.sh"
params_3="$PWD/config_files/summarize_params.txt"

###########################################
# Gather pipeline parameter data
function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }

. $params_1
assert_
. $params_2
assert_

genes_analyzed="$2"
oligoD="$3"
id_len="$4"

echo "genes_analyzed = $genes_analyzed
relabeling oligo = $oligoD
sample identifier length = $id_len 
#########
"
###########################################
# Iterate over analyzed data
# If BC3_min_cutoff different from 1 was used in the first part of the pipeline, please modify the line below with relevant value(s)
for BC3_min_cutoff in $BC3_min; do
    
    # Check that input data exists and create the output directory
    if [ ! -d "$params_dir_out_1/tables.BC3cutoff$BC3_min_cutoff" ]; then
        exit "Error: not all parameters/files exit"
    fi

    outdir="$params_dir_out_1/tables.BC3cutoff$BC3_min_cutoff.relabeling"
    if [ ! -d $outdir ]; then 
        mkdir $outdir; assert_; 
    fi
    
    runtime=$(date +"%y-%m-%dT%H%M")
    logdir="$outdir/logs_$runtime"
    
    if [ ! -d $logdir ]; then 
        mkdir $logdir; assert_; 
    fi
    
    # Iterate over analyzed sample treatment (Cont/Exp)
    for i in ${!title[@]}; do
        
        # Gather input data parameters 
        title1="${title[i]}"_"idx$i"
        refs1=$(echo "${sort_ref[i]}"";""${sort_refs[i]}" | tr ";" "\n")
        
        size_f1="${size_f[$i]}"; assert_
        read_seq="${read_seq[$i]}"; assert_
        read_pos="${read_pos[$i]}"; assert_    
        
        echo "size_f = $size_f1"
        echo "read_seq = $read_seq"
        echo "read_pos = $read_pos"
        echo "#########" 
        
        # Iterate over sorted read mutation data and run consensus jobs
        k=0
        for r1 in $refs1; do
            let k=$k+1
            echo "k = "$k
            r2="$r1"
            
            if [ $k -eq 1 ]; then r2="$r2.others"; fi # First item in loop refers to "others" (ambigous origin reads)
            
            #  Skip data that do not match user defined gene of origin
            flag=0
            for gene in $(echo $genes_analyzed | tr "," "\n"); do
                if [ $r1 == $gene ]; then
                    let flag+=1
                fi  
            done
            
            if [ $flag == 0 ]; then
                continue
            fi            
            
            # Delete three lines below to also analyze "others" files. By default they are skipped, since they usually do not contain enough data for meaningful analysis
            if [[ $r2 == *"others"* ]]; then 
                continue
            fi           
            
            # Analyzed file and reference names
            bam_name="$title1.$r2.bwa.sorted.bam"
            barcode_name="$title1.assembled.filtered.fastq.trimmed.wrongId.barcodes"
            ref1="$params_dir_reference"/"$r1".fa
            
            echo "bam_name = $bam_name"
            echo "ref1 = $ref1"
            
            # Iterate over reads having specific BC3s
            for sam_BX_tag in ''; do
                # Define BC3 tag variable added to the output file name
                echo "sam_BX_tag = \"$sam_BX_tag\""
                tag1=""
                if [ "$sam_BX_tag" != "" ]; then tag1="tag-$sam_BX_tag."; fi
                
                # Define I/O variables
                barcode_table="$params_dir_out_1/filtered/$barcode_name"
                mut_table="$params_dir_out_1/tables.BC3cutoff$BC3_min_cutoff/$bam_name.""$tag1""mutationFrequncyPerBarcode.txt"
                table_out1="$outdir/$bam_name.""$tag1""relabeling.txt"
                log_out1="$logdir/$bam_name.""$tag1""relabeling"
                
                echo "barcode_table = $barcode_table"
                echo "mut_table = $mut_table"
                echo "table_out1 = $table_out1"
                echo "log file = $log_out1"
                echo "#########" 
                
                # Check that all required parameters are OK before executing the jobs
                if [ -f "$barcode_table" ] && [ -f "$mut_table" ] && [ -d "$params_dir_out_1/tables.BC3cutoff$BC3_min_cutoff" ]; then
                    if [ ! -z "$log_out1" ] && [ ! -z "$outdir" ] && [ ! -z "$table_out1" ] && [ ! -z "$BC3_min_cutoff" ]; then
                        echo 'All files and parameters are found'
                        
                        # Run jobs to analyze relabeling incidence
                        if [ $1 -eq 1 ]; then
                            sleep 1

                            sbatch -o "$log_out1.out" -e "$log_out1.err" \
                                   -p hive1d,hive7d,hiveunlim,queen \
                                   --wrap "module load R/3.6.0-foss-2019a; Rscript relabeling_dist.R \"$barcode_table\" \"$table_out1\" \"$oligoD\" \"$size_f1\" \"$read_seq\" \"$read_pos\" \"$id_len\" \"$mut_table\""
                            
                        fi

                    else
                        echo "Error: some of the required parameters are empty"
                        exit
                    fi
                else
                    echo "Error: not all I/O data exists"
                    exit
                fi
            done
            echo '----------------'
        done
    done
done

conda deactivate
