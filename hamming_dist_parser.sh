# A wrapper script for analyzing hamming distance of read families' BC5s grouped by family size
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

###########################################
# Iterate over analyzed data
# If BC3_min_cutoff different from 1 was used in the first part of the pipeline, please modify the line below with relevant value(s)
for BC3_min_cutoff in $BC3_min; do

    # Create the output directory
    outdir="$params_dir_out_1/tables_consensus.BC3cutoff$BC3_min_cutoff.hamming"
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
        
        echo "BC3_min_cutoff = $BC3_min_cutoff"
        
        # Iterate over sorted read data and run consensus jobs
        k=0
        for r1 in $refs1; do
            let k=$k+1
            echo "k = "$k
            r2="$r1"
            if [ $k -eq 1 ]; then r2="$r2.others"; fi # First item in loop refers to "others" (ambigous origin reads)
            
            # Delete three lines below to also analyze "others" files. By default they are skipped, since they usually do not contain enough data for meaningful analysis
            if [[ $r2 == *"others"* ]]; then 
                continue
            fi            
            
            # Analyzed file and reference names            
            bam_name="$title1.$r2.bwa.sorted.bam"
            ref1="$params_dir_reference"/"$r1".fa
            
            echo "bam_name = $bam_name"
            echo "ref1 = $ref1"
            
            # Iterate over reads having specific BC3s
            for sam_BX_tag in ''; do
                echo "sam_BX_tag = \"$sam_BX_tag\""
                tag1=""
                if [ "$sam_BX_tag" != "" ]; then tag1="tag-$sam_BX_tag."; fi
                
                # Define I/O variables
                table_in1="$params_dir_out_1/tables_consensus.BC3cutoff$BC3_min_cutoff/$bam_name.""$tag1""consensus.txt"
                table_out1="$params_dir_out_1/tables_consensus.BC3cutoff$BC3_min_cutoff.hamming/$bam_name.""$tag1""hamming.txt"
                log_out1="$logdir/$bam_name.""$tag1""hmm"
                                
                echo "table_in1 = $table_in1"
                echo "table_out1 = $table_out1"
                echo "log_out1 = $log_out1"
                
                # Check input file presence
                if [ ! -f $table_in1 ]; then 
                    echo "Input $table_in1 is missing!"; assert_;
                    continue 
                fi
                
                # Check that all required parameters are OK before executing the jobs
                if [ -d "$params_dir_out_1/tables_consensus.BC3cutoff$BC3_min_cutoff" ]; then
                    if [ ! -z "$table_in1" ] && [ ! -z "$table_out1" ] && [ ! -z "$log_out1" ] && [ ! -z "$BC3_min_cutoff" ]; then
                        echo 'All files and parameters are found'
                         
                        # Run jobs to analyze hamming distance between read families' BC5s grouped by family size
                        if [ $1 -eq 1 ]; then
                           sleep 1
                           cores=20
                                
                           sbatch -N1 -n1 -o "$log_out1.out" -e "$log_out1.err" --mem=25000 \
                                  -p hive1d,hive7d,hiveunlim,queen \
                                  --wrap "module load R/3.6.0-foss-2019a; Rscript hamming_dist_parser.R $table_in1 $table_out1"                                  
                        fi             
                    
                    else
                        echo "Error: some of the required parameters are empty"
                        exit
                    fi
                else
                    echo "Error: not all I/O data exists"
                fi
            done
            echo '----------------'
        done
	done
done

conda deactivate
