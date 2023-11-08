#!/bin/bash

escval() {
    printf \'
    unescaped=$1
    while :
    do
        case $unescaped in
        *\'*)
            printf %s "${unescaped%%\'*}""'\''"
            unescaped=${unescaped#*\'}
            ;;
        *)
            printf %s "$unescaped"
            break
        esac
    done
    printf \'
}

contains() {
    target=$1
    shift
    for pattern in "$@"; do
        case $target in
            *$pattern*) return 0;;
        esac
    done
    return 1
}

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

_impprove_container_exec() {
    runtime=""
    bind=""
    donebinds=false
    doneexec=false
    for arg in "$@"; do
        if [ -z "$runtime" ]; then
            set --
            runtime="$arg"
            bind="$(_impprove_container_bindarg "$runtime")"
        elif ! $donebinds && [ "$arg" = "--" ]; then
            donebinds=true
        elif $donebinds && ! $doneexec; then
            doneexec=true
            if [ "$runtime" = "podman" ] || [ "$runtime" = "docker" ]; then
                set -- "$@" '--entrypoint' "$arg" 'localhost/apply-impprove'
            elif [ "$runtime" = "apptainer" ] || [ "$runtime" = "singularity" ]; then
                set -- "$@" 'exec' 'apply-impprove.sif' "$arg"
            else
              :
            fi
        elif ! $donebinds; then
            set -- "$@" "$bind" "$arg"
        else
          set -- "$@" "$arg"
        fi
    done
    printf '\n\nruntime: %s\n\n' "$runtime"
    if [ "$runtime" = "podman" ] || [ "$runtime" = "docker" ]; then
        echo "$runtime run $*"
        "$runtime" run "$@"
    elif [ "$runtime" = "apptainer" ] || [ "$runtime" = "singularity" ]; then
        echo "$runtime $*"
        "$runtime" "$@"
    else
      printf ' \e[1;31m!\e[m Container runtime \e[31m%s\e[m is not supported\n' "$runtime" >&2
      return 17
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
    slivar expr --vcf - -g $gnomad_data_file --info '"'"'INFO.impactful && INFO.gnomad_popmax_af < 0.01 && variant.QUAL >= 20 && variant.ALT[0] != "*"'"'"' -o $outfile

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
    if [ -z "$input_file" ]; then
        printf "\e[31mError:\e[m Mandatory options \e[33m-i INFILE.VCF\e[m and -o OUTFILE.VCF are required.\n"
        exit 1
    elif [ -z "$output_file" ]; then
        printf "\e[31mError:\e[m Mandatory options -i INFILE.VCF and \e[33m-o OUTFILE.VCF\e[m are required.\n"
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
    logfile="$tempdir/log.txt"
    printf '[\e[33m%s\e[m] ' "$runtime" >&2
    printf "Checking... " >&2
    checkref=false
    reference=$(_impprove_container_exec "$runtime" "$tempdir:/workdir" -- \
        bcftools view --header-only "/workdir/$(basename "$input_file")" 2>/dev/null \
        | grep '##reference')
    exitcode=$?
    if [ $exitcode -ne 0 ]; then
        printf "\e[33m(could not locate reference header in VCF)\e[m " >&2
    elif contains "$reference" "hg38" "Hg38" "HG38"; then
        :
    elif contains "$reference" "grch38" "Grch38" "GRCh38" "GRCH38"; then
        :
    elif contains "$reference" "hg19" "Hg19" "HG19"; then
        printf "\e[33m(seems to be using Hg19, not Hg38)\e[m " >&2
        checkref=true
    elif contains "$reference" "grch37" "Grch37" "GRCh37" "GRCH37"; then
        printf "\e[33m(seems to be using GRCh37, not GRCh38)\e[m " >&2
        checkref=true
    else
        printf '\e[33m(unrecognised reference genome: %s)\e[m ' "${reference##\#\#reference=file://}" >&2
        checkref=true
    fi
    if $checkref; then
        num_records=$(_impprove_container_exec "$runtime" \
            "$basefolder/data:/data" "$tempdir:/workdir" -- \
            bcftools stats "/workdir/$(basename "$input_file")" 2>&1 |
                   grep "number of records:")
        num_records="${num_records##*number of records:	}"
        num_mismatch=$(_impprove_container_exec "$runtime" \
            "$basefolder/data:/data" "$tempdir:/workdir" -- \
            bcftools norm --check-ref w --fasta-ref /data/Hg38p14.fa \
            "/workdir/$(basename "$input_file")" 2>&1 | grep -c "REF_MISMATCH")
        if [ "$(("$num_mismatch" + "$num_mismatch"))" -gt "$num_records" ]; then
            printf "\n \e[1;31m!\e[m More than half of the variants do not match \
the Hg38 reference genome (%s out of %s). It is likely the input VCF is built against \
an older reference genome. If the VCF uses Hg19/GRCh37 you can use the \
\e[36mliftover\e[m subcommand to lift it to Hg38 like so:\n\n  \
%s liftover -i %s -o %s.hg38.vcf\n\n\
Then try again.\n" "$num_mismatch" "$num_records" "$(basename "$0")" "$input_file" "${input_file%.vcf}"
            return 18
        fi
    fi
    printf "Running bcftools... " >&2
    printf ':: BCFTOOLS ::\n' > "$logfile"
    _impprove_container_exec "$runtime" \
        "$basefolder/data:/data" "$tempdir:/workdir" -- \
        bcftools csq --samples - \
        --ncsq 40 --gff-annot "/data/$annotfile" --local-csq \
        --fasta-ref "/data/$genomefile" "/workdir/$(basename "$input_file")" \
        --output-type u -o "/workdir/intermediate.bcf" \
        2>>"$logfile"
    exitcode=$?
    printf '\nexit code = %s' "$exitcode"
    if [ $exitcode -ne 0 ]; then
        printf '\e[31mexit code %s\e[m\n  See the logfile %s for more information\n' "$exitcode" "$logfile" >&2
        return 1
    fi
    printf "slivar... " >&2
    printf '\n\n:: SLIVAR ::\n' >> "$logfile"
    _impprove_container_exec \
        "$basefolder/data:/data" "$tempdir:/workdir" -- \
        slivar expr --vcf "/workdir/intermediate.bcf" \
        -g "/data/$gnomadfile" --info "$sinfo" \
        -o "/workdir/outfile.vcf" \
        2>>"$logfile"
    exitcode=$?
    if [ $exitcode -ne 0 ]; then
        printf '\e[31mexit code %s\e[m\n  See the logfile %s for more information\n' "$exitcode" "$logfile" >&2
        return 1
    fi
    # TODO could replace this with a call to `sed` inside the container
    sed -i 's;##bcftools/csq;##bcftools_csq;' "$tempdir/outfile.vcf"
    mv "$tempdir/outfile.vcf" "$output_file"
    rm -rf "$tempdir"
    printf "\e[32mdone\e[m\n" >&2
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
                hpo_args="$hpo_args $(escval "$1")"
                ;;
        esac
        shift
    done
    # Check for mandatory options
    if [ -z "$input_file" ]; then
        printf "\e[31mError:\e[m Mandatory options \e[33m-i INFILE.VCF\e[m and -o OUTFOLDER are required.\n" >&2
        exit 1
    elif [ -z "$output_folder" ]; then
        printf "\e[31mError:\e[m Mandatory options -i INFILE.VCF and \e[33m-o OUTFOLDER\e[m are required.\n" >&2
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

impprove__liftover() {
    if [ $# -eq 0 ] || [ "$1" = "--help" ] || [ "$1" = "-h" ]; then
        printf '
    %s liftover [-f FROM -t TO] -i INFILE.VCF -o OUTFILE.VCF [-r REJFILE.VCF]

Lift \e[36mINFILE.VCF\e[m from genome build \e[36mFROM\e[m (Hg19 by default) to \e[36mTO\e[m (Hg38 by default).
The lifted result is written to \e[36mOUTFILE.VCF\e[m, and rejected variants (that could not
be lifted over) to OUTFILE.VCF.rej if \e[36mREJFILE.VCF\e[m is not specified.
' "$(basename "$0")"
        exit 0
    fi
    _impprove_check_dirs "data/liftover"
    hg_from="19"
    hg_to="38"
    input_file=""
    output_file=""
    rejects_file=""
    while [ "$#" -gt 0 ]; do
        case "$1" in
            -f)
                shift
                hg_from="$1"
                ;;
            -t)
                shift
                hg_to="$1"
                ;;
            -i)
                shift
                input_file="$1"
                ;;
            -o)
                shift
                output_file="$1"
                ;;
            -r)
                shift
                rejects_file="$1"
                ;;
            *)
                printf "Flag \e[31m%s\e[m unrecognised\n" "$1"
                ;;
        esac
        shift
    done
    # Check for mandatory options
    if [ -z "$input_file" ]; then
        printf "\e[31mError:\e[m Mandatory options \e[33m-i INFILE.VCF\e[m and -o OUTFILE.VCF are required.\n" >&2
        exit 1
    elif [ -z "$output_file" ]; then
        printf "\e[31mError:\e[m Mandatory options -i INFILE.VCF and \e[33m-o OUTFILE.VCF\e[m are required.\n" >&2
        exit 1
    fi
    if [ ! -e "$input_file" ]; then
        printf '\e[31mError:\e[m \e[36m%s\e[m is not a file or folder\n' "$input_file" >&2
        exit 1
    fi
    if [ -z "$rejects_file" ]; then
        rejects_file="$output_file.rej"
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
    printf "Lifting... " >&2
    printf ':: Liftover::\n' > "$logfile"
    _impprove_container_exec "$runtime" \
        "$basefolder/data:/data" "$tempdir:/workdir" -- \
        "/liftover/liftover.jl" "$hg_from" "$hg_to" \
        "/workdir/$(basename "$input_file")" "/workdir/lifted.vcf" \
        "/workdir/rej.vcf" >>"$logfile"
    exitcode=$?
    if [ $exitcode -ne 0 ]; then
        printf '\e[31mexit code %s\e[m\n  See the logfile %s for more information\n' "$exitcode" "$logfile" >&2
        didlogproblem=true
    fi
    if ! $didlogproblem; then
        mv "$tempdir/lifted.vcf" "$output_file"
        [ -s "$tempdir/rej.vcf" ] &&
            mv "$tempdir/rej.vcf" "$rejects_file"
        rm -rf "$tempdir"
        printf "\e[32mdone\e[m\n" >&2
    fi
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
    %s liftover|filter|predict [args...]

Either liftover an Hg19 VCF to Hg38, filter a VCF file in preparation for
prediction, or create predictions from a VCF file. See the \e[34mliftover\e[m,
\e[34mfilter\e[m, and \e[34mpredict\e[m subcommand'"'"'s help for more
information.

This scripts acts as a small shim over the container, for convenience.
Run the container itself with --help to see more information on it.
' "$(basename "$0")"
    exit 0
else
  printf 'IMPPROVE subcommand \e[31m%s\e[m not recognized. Try \e[34mfilter\e[m or \e[34mpredict\e[m.\n' "$1" >&2
  exit 1
fi
