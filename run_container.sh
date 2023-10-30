#!/bin/sh

_impprove_container_runtime() {
    if [ -n "$IMPPROVE_CONTAINER_RUNTIME" ]; then
        if command -v "$IMPPROVE_CONTAINER_RUNTIME" 1>/dev/null; then
            printf '%s' "$IMPPROVE_CONTAINER_RUNTIME"
        else
          # printf '\e[31m%s\e[m not found on the system path!\n' "$IMPPROVE_CONTAINER_RUNTIME"
          unset IMPPROVE_CONTAINER_RUNTIME
          _impprove_container_runtime
        fi
    elif command -v podman 1>/dev/null; then
        printf "podman"
    elif command -v apptainer 1>/dev/null; then
        printf "apptainer"
    elif command -v singularity 1>/dev/null; then
        printf "singularity"
    elif command -v docker 1>/dev/null; then
        printf "docker"
    else
      printf '\e[31mError:\e[m No container runtime found!\n' >&2
      exit 1
    fi
}

_impprove_container_bindarg() {
    if [ "$1" = "podman" ] || [ "$1" = "docker" ]; then
        printf '%s' "-v"
    elif  [ "$1" = "apptainer" ] || [ "$1" = "singularity" ]; then
        printf '%s' "-B"
    fi
}

_impprove_check_dirs() {
    for dir in "$@"; do
        if [ ! -d "$(dirname "$0")/$dir" ]; then
            printf " \e[31;1m!\e[m \e[2mThe \e[36m/%s\e[0;2m directory is missing\e[m\n" "$dir" >&2
            exit 1
        fi
    done
}

impprove__filter() {
    if [ $# -eq 0 ] || [ "$1" = "--help" ] || [ "$1" = "-h" ]; then
        printf '
    %s filter -i INFILE.VCF -o OUTFILE.VCF [-e] [--info INFO]

This command runs INFILE.VCF through bcftools then slivar, writing the output to
OUTFILE.vcf. The input VCF must use the CRGh38/Hg38 genome assembly.

When using Ensembl identifiers, supply the \e[34m-e\e[m flag.

It is equivalent to running:
    bcftools csq --samples - --ncsq 40 --gff-annot $annot_data_file --local-csq --fasta-ref $genome_data $vcf_input_file --output-type u |
    slivar expr --vcf - -g $gnomad_data_file --info ''INFO.impactful && INFO.gnomad_popmax_af < 0.01 && variant.QUAL >= 20 && variant.ALT[0] != "*"'' -o $outfile

The \e[34m--info\e[m passed to slivar can be customised by providing your own value.
' "$(basename "$0")"
        exit 0
    fi
    _impprove_check_dirs "data"
    input_file=""
    output_file=""
    sinfo='INFO.impactful && INFO.gnomad_popmax_af < 0.01 && variant.QUAL >= 20 && variant.ALT[0] != "*"'
    genomefile="Hg38p14.fa"
    annotfile="Homo_sapiens.GRCh38.108.gff3.gz"
    gnomadfile="gnomad.hg38.v2.zip"
    while [ "$#" -gt 0 ]; do
        case "$1" in
            -i)
                shift
                input_file="$1"
                ;;
            -o)
                shift
                output_file="$1"
                ;;
            -e)
                genomefile="Hg38p14e.fa"
                ;;
            --info)
                shift
                sinfo="$1"
                ;;
            *)
                printf "Flag \e[31m%s\e[m unrecognised\n" "$1"
                ;;
        esac
        shift
    done
    # Check for mandatory options
    if [ -z "$input_file" ] || [ -z "$output_file" ]; then
        printf "\e[31mError:\e[m Mandatory options -i INFILE.VCF and -o OUTFILE.VCF are required.\n"
        exit 1
    fi
    if [ ! -f "$input_file" ]; then
        printf '\e[31mError:\e[m \e[36m%s\e[m is not a file\n' "$input_file"
        exit 1
    fi
    tempdir=$(mktemp -d)
    cp "$input_file" "$tempdir"
    basefolder="$(dirname "$(realpath -s "$0")")"
    runtime="$(_impprove_container_runtime)"
    bind="$(_impprove_container_bindarg "$runtime")"
    logfile="$tempdir/log.txt"
    didlogproblem=false
    set -- "$bind" "$basefolder/data:/data" "$bind" "$tempdir:/workdir"
    printf '[\e[33m%s\e[m] ' "$runtime" >&2
    if [ "$runtime" = "podman" ] || [ "$runtime" = "docker" ]; then
        printf "Running bcftools... "
        printf ':: BCFTOOLS ::\n' > "$logfile"
        $runtime run "$@" --entrypoint="bcftools" localhost/apply-impprove \
            csq --samples - \
            --ncsq 40 --gff-annot "/data/$annotfile" --local-csq \
            --fasta-ref "/data/$genomefile" "/workdir/$(basename "$input_file")" \
            --output-type u -o "/workdir/intermediate.bcf" \
            2>>"$logfile"
        exitcode=$?
        if [ $exitcode -ne 0 ]; then
            printf '\e[31mexit code %s\e[m\n  See the logfile %s for more information\n' "$exitcode" "$logfile" >&2
            didlogproblem=true
        fi
        printf "slivar... "
        printf '\n\n:: SLIVAR ::\n' >> "$logfile"
        $runtime run "$@" --entrypoint="slivar" localhost/apply-impprove \
            expr --vcf "/workdir/intermediate.bcf" \
            -g "/data/$gnomadfile" --info "$sinfo" \
            -o "/workdir/outfile.vcf" \
            2>>"$logfile"
        exitcode=$?
        if [ $exitcode -ne 0 ]; then
            printf '\e[31mexit code %s\e[m\n  See the logfile %s for more information\n' "$exitcode" "$logfile" >&2
            didlogproblem=true
        fi
    elif [ "$runtime" = "apptainer" ] || [ "$runtime" = "singularity" ]; then
        printf "Running bcftools... "
        printf ':: BCFTOOLS ::\n' > "$logfile"
        $runtime exec "$@" "$basefolder/apply-impprove.sif" \
            bcftools csq --samples - \
            --ncsq 40 --gff-annot "/data/$annotfile" --local-csq \
            --fasta-ref "/data/$genomefile" "/workdir/$(basename "$input_file")" \
            --output-type u -o "/workdir/intermediate.bcf" \
            2>>"$logfile"
        exitcode=$?
        if [ $exitcode -ne 0 ]; then
            printf '\e[31mexit code %s\e[m\n  See the logfile %s for more information\n' "$exitcode" "$logfile" >&2
            didlogproblem=true
        fi
        printf "slivar... "
        printf '\n\n:: SLIVAR ::\n' >> "$logfile"
        $runtime exec "$@" "$basefolder/apply-impprove.sif" \
            slivar expr --vcf "/workdir/intermediate.bcf" \
            -g "/data/$gnomadfile" --info "$sinfo" \
            -o "/workdir/outfile.vcf" \
            2>>"$logfile"
        exitcode=$?
        if [ $exitcode -ne 0 ]; then
            printf '\e[31mexit code %s\e[m\n  See the logfile %s for more information\n' "$exitcode" "$logfile" >&2
            didlogproblem=true
        fi
    else
      printf ' \e[1;31m!\e[m Container runtime \e[31m%s\e[m is not supported\n' "$runtime" >&2
    fi
    if ! $didlogproblem; then
        # TODO could replace this with a call to `sed` inside the container
        sed -i 's;##bcftools/csq;##bcftools_csq;' "$tempdir/outfile.vcf"
        mv "$tempdir/outfile.vcf" "$output_file"
        rm -rf "$tempdir"
        printf "\e[32mdone\e[m\n" >&2
    fi
}

impprove__predict() {
    if [ $# -eq 0 ] || [ "$1" = "--help" ] || [ "$1" = "-h" ]; then
        printf '
    %s predict -i INFILE.VCF -o OUTFOLDER [-m MODELS] HPOS...

Apply IMPPROVE models to INFILE.VCF, writing results to a folder within OUTFOLDER.

Restrict the models used with \e[34m-m\e[m.
Specify the HPO terms to use with \e[34mHPOS...\e[m.

The \e[34m-i\e[m argument can also point to a folder of VCF files.
' "$(basename "$0")"
        exit 0
    fi
    _impprove_check_dirs "data" "models"
    input_file=""
    output_folder=""
    models_arg=""
    hpo_args=""
    while [ "$#" -gt 0 ]; do
        case "$1" in
            -i)
                shift
                input_file="$1"
                ;;
            -o)
                shift
                output_folder="$1"
                ;;
            -m)
                shift
                models_arg="$1"
                ;;
            *)
                hpo_args="$hpo_args $1"
                ;;
        esac
        shift
    done
    # Check for mandatory options
    if [ -z "$input_file" ] || [ -z "$output_folder" ]; then
        printf "\e[31mError:\e[m Mandatory options -i INFILE.VCF and -o OUTFILE.VCF are required.\n" >&2
        exit 1
    fi
    if [ ! -e "$input_file" ]; then
        printf '\e[31mError:\e[m \e[36m%s\e[m is not a file or folder\n' "$input_file" >&2
        exit 1
    fi
    if [ ! -d "$output_folder" ]; then
        mkdir -p "$output_folder"
    fi
    tempdir=$(mktemp -d)
    if [ -d "$input_file" ]; then
        input_folder="$input_file"
    else
        input_folder="$tempdir"
        cp "$input_file" "$tempdir"
    fi
    basefolder="$(dirname "$(realpath -s "$0")")"
    runtime="$(_impprove_container_runtime)"
    bind="$(_impprove_container_bindarg "$runtime")"
    set -- "$bind" "$basefolder/models:/models" \
        "$bind" "$basefolder/data:/data" \
        "$bind" "$input_folder:/vcfs" \
        "$bind" "$output_folder:/predictions"
    printf '[\e[33m%s\e[m] ' "$runtime" >&2
    printf "Predicting...\n" >&2
    if [ "$runtime" = "podman" ] || [ "$runtime" = "docker" ]; then
        if [ -z "$models_arg" ]; then
           $runtime run "$@" "localhost/apply-impprove" $hpo_args
        else
           $runtime run "$@" "localhost/apply-impprove" \
               -m "$models_arg" $hpo_args
        fi
    elif [ "$runtime" = "apptainer" ] || [ "$runtime" = "singularity" ]; then
        if [ -z "$models_arg" ]; then
           $runtime run "$@" "$basefolder/apply-impprove.sif" $hpo_args
        else
           $runtime run "$@" "$basefolder/apply-impprove.sif" \
               -m "$models_arg" $hpo_args
        fi
    else
      printf ' \e[1;31m!\e[m Container runtime \e[31m%s\e[m is not supported\n' "$runtime" >&2
    fi
    rm -rf "$tempdir"
    printf "\e[32mDone\e[m\n" >&2
}

impprove() {
  cmdname="$1"; shift
  if type "impprove__$cmdname" >/dev/null 2>&1; then
    "impprove__$cmdname" "$@"
  else
    command impprove "$cmdname" "$@" # call the *real* impprove command
  fi
}

# if the functions above are sourced into an interactive interpreter, the user can
# just call "impprove predict" or "impprove filter" with no further code needed.

# if invoked as a script rather than sourced, call function named on argv via the below;
# note that this must be the first operation other than a function definition
# for $_ to successfully distinguish between sourcing and invocation:
# [ "$_" != "$0" ] && return

# make sure we actually *did* get passed a valid subcommand name
if type "impprove__$1" >/dev/null 2>&1; then
  # invoke that subcommand, passing arguments through
  impprove "$@" # same as "$1" "$2" "$3" ... for full argument list
elif [ -z "$1" ]; then
    printf '
    %s filter|predict [args...]

Either filter a VCF file in preperation for prediction, or create predictions
from a VCF file. See the \e[34mfilter\e[m and \e[34mpredict\e[m subcommands'' help
for more information.

This scripts acts as a small shim over the container, for convenience.
Run the container itself with --help to see more information on it.
' "$(basename "$0")"
    exit 0
else
  printf 'IMPPROVE subcommand \e[31m%s\e[m not recognized. Try \e[34mfilter\e[m or \e[34mpredict\e[m.\n' "$1" >&2
  exit 1
fi
