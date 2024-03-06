# A wrapper script for collecting error rate data from summary files produced by "consensus-count_summary" script
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
    outdir="$params_dir_out_1/tables_consensus.BC3cutoff$BC3_min_cutoff.err_rate"
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
        
        # Iterate over sorted read data
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
                # Define BC3 tag variable added to the output file name
                echo "sam_BX_tag = \"$sam_BX_tag\""
                tag1=""
                if [ "$sam_BX_tag" != "" ]; then tag1="tag-$sam_BX_tag."; fi
                
                # Define I/O variables and check that input directory exists
                table_in1_dir="$params_dir_out_1/tables_consensus.BC3cutoff$BC3_min_cutoff.summary/$bam_name""$tag1"".cutoffs-update"
                log_out1="$logdir/$bam_name.""$tag1""err"
                
                if [ ! -d $table_in1_dir ]; then 
                    echo "One or more $bam_name.$tag1 dirs are missing!"; assert_;
                    continue 
                fi
                
                echo "table_in1_dir = $table_in1_dir"
                echo "log_out1 = $log_out1"
                
                # Check that all required parameters are OK before executing the jobs
                if [ ! -z "$table_in1_dir" ] && [ ! -z "$BC3_min_cutoff" ]; then
                    echo 'All parameters are found'
                    
                    # Run jobs to collect error rate data
                    if [ $1 -eq 1 ]; then
                       sleep 1
                            
                       sbatch -o "$log_out1.out" -e "$log_out1.err" --mem=28000 \
                              -p hive1d,hive7d,hiveunlim,queen \
                              --wrap "module load R/3.6.0-foss-2019a; Rscript err_rate_collate.R $table_in1_dir $outdir"
                    fi
                
                else
                    echo "Error: some of the required parameters are empty"
                    exit
                fi
            done
            echo '----------------'
        done
	done
done

conda deactivate
