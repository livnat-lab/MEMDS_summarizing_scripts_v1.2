# An interactive wrapper for running mutation analysis scripts across multiple samples

######## Set title ########
fgRed=$(tput setaf 1)
txBold=$(tput bold)
txReset=$(tput sgr0)

termwidth="$(tput cols)"
mssg="Welcome to the MEMDS analysis pipeline"
let padding=($termwidth - ${#mssg} - 3)/2

printf "\n%${padding}s ${fgRed}${txBold} $mssg %${padding}s ${txReset}\n"

######## Gather path infromation ########
echo "Please provide path(s) to the MEMDS pipeline script folder. E.g: My_computer/My_experiment/scripts"
printf "Multiple paths can be provided, separated by semi-colon (;)\n\n"

read path_str

paths=$(echo "$path_str" | tr ";" "\n")

printf "\n${fgRed}${txBold}These paths were provided:${txReset}\n"

for path in $paths; do
    echo "$path"
done

printf "\nAre these paths correct? Y/N: \n"
while read input; do    
        case $input in
            [yY])
                echo 'Continuing'
                break
                ;;
            [nN])
                echo "Incorrect paths encountered. Exiting"
                exit
                ;;
             *)
                printf "Please input Y or N\n" >&2
        esac
done

######## List possible steps ########
explanation="
${fgRed}${txBold}Choose one of the following steps:${txReset}\n
1. Analyze barcode distribution\n
2. Parse Hamming distance of 5' barcodes\n
3. Summarize and analyze consensus counts for each cut-off set\n
4. Check nucleotide frequencies within primary barcodes (5')\n
5. Check relabeling incidence in read families\n\n"

printf "$explanation"

######## Execute the steps ########
for path in $paths; do

    ########
    printf "\n^v^v^v^v^v^v^v^v^v^\n"
    printf "${txBold}Working on $path${txReset}\n"
    
    echo "Please choose step to run [1-5]: "
    read step
    
    while true; do
        if [[ ! $step =~ ^[0-9]+$ ]] || [ $step -gt 5 ]; then
            echo "Please choose ${txBold}valid${txReset} step to run [1-5]: "
            read step
            continue
        else 
            break 
        fi
    done
    
    while true; do
    printf "\nIs the step correct? Y/N: \n"	
        
        read input
        case $input in
            [yY])
                echo "Executing $step"
                break
                ;;
            [nN])
                echo "Please choose step to run [1-5]: "
                read step
                
                while true; do
                    if [[ ! $step =~ ^[0-9]+$ ]] || [ $step -gt 5 ]; then
                        echo "Please choose ${txBold}valid${txReset} step to run [1-5]: "
                        read step
                        continue
                    else 
                        break 
                    fi
                done
                    
                continue
                ;;
             *)
                printf "Please input Y or N" >&2
        esac
    done
    ########
    
    case $step in
       1)
         printf "\n${fgRed}${txBold}Choose one of the following sub-steps:${txReset}\n
1. Analyze barcode distribution, based on mutation tables\n
2. Calculate percent of families with 15+ 3' barcodes\n
3. Check eveness of 3' barcode group size distribution\n
4. Check how many 3' barcode groups have more than one read in them\n
5. Report frequencies of N largest 3' barcode groups\n\n"
         printf "Please choose sub-step to run [1-5]: "
         
         read substep
        
         while true; do
            if [[ ! $substep =~ ^[0-9]+$ ]] || [ $substep -gt 5 ]; then
                echo "Please choose ${txBold}valid${txReset} sub-step to run [1-5]: "
                read substep
                continue
            else 
                break 
            fi
         done
         
         case $substep in
            1)
            cd $path
            bash BC_dist.sh 1
            ;;
            2)
            cd $path
            bash BC3_above15.sh 1
            ;;            
            3)
            cd $path
            bash BC_evenness.sh 1
            ;;
            4)
            cd $path
            bash BC3_size_above1.sh 1
            ;;
            5)
            cd $path
            
            printf "${txBold}How many BC3 groups should be reported? ${txReset}\n"
            read col_report
                        
            bash BC3_freq.sh 1 $col_report                                    
         esac            
         ;;
       2)
         cd $path
         bash hamming_dist_parser.sh 1
         ;;
       3)
         printf "\n${fgRed}${txBold}Choose one of the following sub-steps:${txReset}\n
1. Summarise consensus counts per cut-off criteria set\n
2. Create a table of counted substitutions per cut-off criteria set\n
3. Create a table of expected substitution counts inside ROI per cut-off criteria set\n
4. Create a table of expected error and mutation rates per cut-off criteria set\n
5. Recalculate expected error and mutation rates using custom parameters\n
6. Plot substitution count per read position, as a percent of all substitutions\n
7. Create a table of plasmid to WT ratios per cut-off criteria set\n\n"
         printf "Please choose sub-step to run [1-7]: "
         
         read substep
   
         while true; do
            if [[ ! $substep =~ ^[0-9]+$ ]] || [ $substep -gt 7 ]; then
                echo "Please choose ${txBold}valid${txReset} sub-step to run [1-7]: "
                read substep
                continue
            else 
                break 
            fi
         done
         
         case $substep in
             1)
             cd $path
             bash consensus-count_summary.sh 1
             ;;
             2)
             cd $path
            
             printf "${txBold}Please list gene names, comma-separated, as listed in parameter table, to which apply the script: ${txReset}\n"
             read genes_analyzed
             
             printf "${txBold}Please list positions to exclude from counts, comma-separated (e.g.: 39,52,...) or \"NONE\": ${txReset}\n"
             read excl_pos
             
             printf "${txBold}Please list substitutions to exclude from counts, comma-separated (e.g.: CA,GA,...) or \"NONE\": ${txReset}\n"
             read excl_sub
             
             printf "${txBold}Please list variants to exclude from counts, comma-separated (e.g.: 39CT,47GA,...) or \"NONE\": ${txReset}\n"
             read excl_var         	
             
             bash substitution_collate.sh 1 $genes_analyzed $excl_pos $excl_sub $excl_var
             ;;
             3)
             cd $path
            
             printf "${txBold}Please list gene names, comma-separated, as listed in parameter table, to which apply the script: ${txReset}\n"
             read genes_analyzed            
            
             printf "${txBold}Please list positions to exclude from counts, comma-separated (e.g.: 39,52,...) or \"NONE\": ${txReset}\n"
             read excl_pos
            
             printf "${txBold}Please list substitutions to exclude from counts, comma-separated (e.g.: CA,GA,...) or \"NONE\": ${txReset}\n"
             read excl_sub
            
             printf "${txBold}Please list variants to exclude from counts, comma-separated (e.g.: 39CT,47GA,...) or \"NONE\": ${txReset}\n"
             read excl_var           
            
             bash expected_mut.sh 1 $genes_analyzed $excl_pos $excl_sub $excl_var
             ;;         	
             4)
             cd $path
             bash err_rate_collate.sh 1
             ;;
             5)
             cd $path
             
             printf "${txBold}Please list gene names, comma-separated, as listed in parameter table, to which apply the script: ${txReset}\n"
             read genes_analyzed                        
            
             printf "${txBold}Please list positions to exclude from counts, comma-separated (e.g.: 39,52,...) or \"NONE\": ${txReset}\n"
             read excl_pos
            
             printf "${txBold}Please list substitutions to exclude from counts, comma-separated (e.g.: CA,GA,...) or \"NONE\": ${txReset}\n"
             read excl_sub
            
             printf "${txBold}Please list variants to exclude from counts, comma-separated (e.g.: 39CT,47GA,...) or \"NONE\": ${txReset}\n"
             read excl_var           
            
             bash re_calc_error_rate.sh 1 $genes_analyzed $excl_pos $excl_sub $excl_var
             ;;
             6)
             cd $path
             
             printf "${txBold}Please list gene names, comma-separated, as listed in parameter table, to which apply the script: ${txReset}\n"
             read genes_analyzed                        
             
             printf "${txBold}Please list positions to exclude from counts, comma-separated (e.g.: 39,52,...) or \"NONE\": ${txReset}\n"
             read excl_pos
             
             printf "${txBold}Please list substitutions to exclude from counts, comma-separated (e.g.: CA,GA,...) or \"NONE\": ${txReset}\n"
             read excl_sub
            
             printf "${txBold}Please list variants to exclude from counts, comma-separated (e.g.: 39CT,47GA,...) or \"NONE\": ${txReset}\n"
             read excl_var           
            
             bash plot_subs.sh 1 $genes_analyzed $excl_pos $excl_sub $excl_var
             ;;
             7)
             cd $path
             
             printf "${txBold}Please list gene names, comma-separated, as listed in parameter table, to which apply the script: ${txReset}\n"
             read genes_analyzed                        
             
             printf "${txBold}Please list plasmid profiles to report, with multiple profiles separated by semicolon(e.g.: 39CT,46C-,47C-;39CT,48T-,51-A): ${txReset}\n"
             read plasmid             
             bash plasmid-WT_ratio.sh 1 $genes_analyzed $plasmid                     	
         esac
         ;;  
       4)
         cd $path
         
         printf "${txBold}Please list gene names, comma-separated, as listed in parameter table, to which apply the script: ${txReset}\n"
         read genes_analyzed          
        
         printf "${txBold}Please list length of sample identifier sequence in BC5. If \"factors table\" lists multiple options, supply a list of lengths separated by \"|\": ${txReset}\n"
         read id_len          
         
         bash BC5_nucleotide_dist.sh 1 $genes_analyzed $id_len
         ;;
       5)
         cd $path
         
         printf "${txBold}Please list gene names, comma-separated, as listed in parameter table, to which apply the script: ${txReset}\n"
         read genes_analyzed          
         
         printf "${txBold}Please provide relabeling oligo sequence: ${txReset}\n"
         read oligoD
        
         printf "${txBold}Please list length of sample identifier sequence in BC5. If \"factors table\" lists multiple options, supply a list of lengths separated by \"|\": ${txReset}\n"
         read id_len               
        
         bash relabeling_dist.sh 1 $genes_analyzed $oligoD $id_len                                                 		
    esac
done
