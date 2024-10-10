#!/usr/bin/bash

# Function to display help
display_help() {
    echo "Usage: $0 -b BAMFILE -e BEDFILE -o OUTDIR"
    echo
    echo "   -b BAMFILE    Path to the input BAM file"
    echo "   -e BEDFILE    Path to the input BED file"
    echo "   -o OUTDIR     Directory to save the output files"
    echo
    echo "Options:"
    echo "   --help        Display this help message and exit"
    exit 0
}

# Parse command-line options
while getopts ":b:e:o:h" opt; do
    case ${opt} in
        b )
            bamfile=$OPTARG
            ;;
        e )
            bedfile=$OPTARG
            ;;
        o )
            outdir=$OPTARG
            ;;
        h )
            display_help
            ;;
        \? )
            echo "Invalid option: -$OPTARG" >&2
            display_help
            ;;
        : )
            echo "Option -$OPTARG requires an argument." >&2
            display_help
            ;;
    esac
done
shift $((OPTIND -1))

# Ensure all required arguments are provided
if [ -z "$bamfile" ] || [ -z "$bedfile" ] || [ -z "$outdir" ]; then
    echo "Error: Missing required arguments"
    display_help
fi

/usr/bin/bedtools bamtobed -i "${bamfile}" > "${outdir}/temp.bed"
/usr/bin/bedtools coverage -d -a "${bedfile}" -b "${outdir}/temp.bed" > "${outdir}/temp.uniform.bed"
/home/pipelines/Consensus_pipeline_with_espresso/scripts/coverage_plots/coverage_plot.py "${outdir}/temp.uniform.bed" "${outdir}/coverage_plot.pdf"
rm "${outdir}/temp.bed" "${outdir}/temp.uniform.bed"

