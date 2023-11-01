#!/usr/bin/env sh

help () {
    printf '
  \e[1mapply-impprove -m MODELS [HPOs...]\e[m

For each VCF file in \e[94m/vcfs\e[m, predict pathogenic variants
using a baseline score and expression data, on a per-\e[32mHPO\e[m basis..

Output files will be created in \e[94m/predictions\e[m.

\e[1mSwitches\e[m
  \e[34m-m, --models\e[m  The model sets to use (default is everything)

\e[1mFlags\e[m
 • \e[34m--readme\e[m  \e[2mPrint the project README.org\e[m
 • \e[34m--sample\e[m  \e[2mPrint the sample VCF file\e[m

\e[3mSample run\e[m
  \e[2m┆\e[m \e[95mpodman run \\\e[m
  \e[2m┆\e[m \e[95m  -v ~/Downloads/impprove_models:/models \\\e[m
  \e[2m┆\e[m \e[95m  -v ~/Downloads/impprove_data:/data \\\e[m
  \e[2m┆\e[m \e[95m  -v ~/Documents/impprove_vcfs:/vcfs \\\e[m
  \e[2m┆\e[m \e[95m  -v ~/Documents/impprove_predictions:/predictions \\\e[m
  \e[2m┆\e[m \e[95m  localhost/apply-impprove\e[m

See the README (\e[34m--readme\e[m) for information on how the input VCF should
be processed.

The data files needed must be placed in \e[94m/data\e[m.

For bulk processing, a \e[35mhpo.map\e[0m file can be placed under \e[94m/vcfs\e[m, of the format:
  \e[2m┆\e[m \e[95mVCF_FILENAME1 HPO1, HPO2...\e[m
  \e[2m┆\e[m \e[95mVCF_FILENAME2 HPOS...\e[m
  \e[2m┆\e[m \e[95m...\e[m

For each \e[1mHPO\e[m term used, a matching model must exist in \e[94m/models\e[m.

\e[1mSupported model baseline scores\e[m
 • \e[34mCADD\e[m  \e[2mCombined Annotation Dependent Depletion\e[m
 • \e[34mAM\e[m    \e[2mAlphaMissense\e[m

\e[1mSupported model expression data\e[m
 • \e[34mGTEx\e[m             \e[2mGenotype-Tissue Expression bulk\e[m
 • \e[34mTabSapTissue\e[m     \e[2mTabula Sapiens Tissue pseudobulk\e[m
 • \e[34mTabSapTissueM\e[m    \e[2mTabula Sapiens Tissue mean pseudobulk\e[m
 • \e[34mTabSapTissueD\e[m    \e[2mTabula Sapiens Tissue, dispersion\e[m
 • \e[34mTabSapCelltype\e[m   \e[2mTabula Sapiens Celltype pseudobulk\e[m
 • \e[34mTabSapCelltypeM\e[m  \e[2mTabula Sapiens Celltype mean pseudobulk\e[m
 • \e[34mTabSapCelltypeD\e[m  \e[2mTabula Sapiens Celltype dispersion\e[m
';
};

while getopts ":h" option; do
    case $option in
        h) help
           exit 0;;
        *) ;;
    esac
done

if [ "$1" = "--readme" ]; then
    exec cat "$(dirname "$0")/README"
elif [ "$1" = "--sample" ]; then
    exec cat "$(dirname "$0")/sample.vcf"
fi

[ -d "/data" ] ||
    printf '\e[33;1m[ Warning:\e[m The \e[94m/data\e[0m directory is not mounted\n'
[ -d "/models" ] ||
    printf '\e[33;1m[ Warning:\e[m The \e[94m/models\e[0m directory is not mounted\n'
[ -d "/vcfs" ] ||
    printf '\e[33;1m[ Warning:\e[m The \e[94m/vcfs\e[0m directory is not mounted\n'
[ -d "/predictions" ] ||
    printf '\e[33;1m[ Warning:\e[m The \e[94m/predictions\e[0m directory is not mounted\n'

export JULIA_DEPOT_PATH="$(mktemp -d --suffix=".4apptainer"):/usr/local/julia/local/share/julia:/usr/local/julia/share/julia"
exec julia -J"$(dirname "$0")/sysimage.so" --project="$(dirname "$0")" "$(dirname "$0")/apply.jl" "$@"
