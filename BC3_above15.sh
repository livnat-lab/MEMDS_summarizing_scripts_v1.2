# A wrapper script for reporting amount of read families per sample associated with more than 15 different 3' barcodes
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
    
    # Check that barcode distribution data exists and create the output directory
    outdir="$params_dir_out_1/tables.BC3cutoff$BC3_min_cutoff.BC_dist"
    if [ ! -d $outdir ]; then
        echo "Error: missing BC_dist folder in $params_dir_out_1"
        exit
    fi

    runtime=$(date +"%y-%m-%dT%H%M")
    logdir="$outdir/logs_$runtime"
    if [ ! -d $logdir ]; then 
        mkdir $logdir; assert_;
    fi
    
    # Extract name of the analyzed sample (for output table)
    parts=($(echo "$params_dir_out_1" | tr "/" "\n"))
    sample=${parts[${#parts[@]}-2]}
    
    # Define I/O variables
    input_dir=$outdir
    table_out1="$input_dir/$sample""_BC3_above15.txt"
    log_out1="$logdir/$sample""_BC3_above15_log"
                
    echo "input = $input_dir"
    echo "table_out1 = $table_out1"
    echo "sample = $sample"
    echo "log_out = $log_out1"
                
    # Check that all required parameters are OK before executing the jobs
    if [ ! -z "$input_dir" ] && [ ! -z "$table_out1" ] && [ ! -z "$BC3_min_cutoff" ]; then
        echo 'All files and parameters are defined'
                    
        # Run jobs to analyze presence of read families with more than 15 different BC3 groups
        if [ $1 -eq 1 ]; then
            sbatch -o "$log_out1.out" -e "$log_out1.err" \
                   -p hive1d,hive7d,hiveunlim,queen \
                    --wrap "module load R/3.6.0-foss-2019a; Rscript BC3_above15.R $input_dir $table_out1"
                        
        fi

     else
        echo "Error: some of the required parameters are empty"
        exit
     fi

done

conda deactivate
