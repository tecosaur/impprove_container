€€[38;5;61m                          â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”
                           IMPPROVE CONTAINER

                                Timothy
                          â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”


                               2023-10-18


€€[0m


€€[1m€€[38;5;67m1 €€[0m€€[1m€€[38;5;67mOverview€€[0m

  This container applies IMPPROVE-trained random forest models to VCF files, then
  saves and aggregates the predictions.

  The input VCFs €€[3mmust€€[0m use the Hg38 reference genome, but this is not (currently)
  validated.

  To filter out common variants, use of €€[1m€€[4m€€[38;5;61m€€[1m€€[4m€€[38;5;61m[slivar]€€[0m€€[0m is recommended. This tool has
  previously been tested with €€[1m€€[4m€€[38;5;61m€€[1m€€[4m€€[38;5;61m[the recommended rare-disease slivar filter]€€[0m€€[0m:

  bcftools csq --samples - --ncsq 40 --gff-annot $€€[38;5;160mannot_data_file€€[0m --local-csq --fasta-ref $€€[38;5;160mgenome_data€€[0m $€€[38;5;160mvcf_input_file€€[0m --output-type u |
  slivar expr --vcf - -g $€€[38;5;160mgnomad_data_file€€[0m --info €€[38;5;64m'INFO.impactful && INFO.gnomad_popmax_af < 0.01 && variant.QUAL >= 20 && variant.ALT[0] != "*"'€€[0m -o $€€[38;5;160moutfile€€[0m

  In the command above, €€[38;5;64m$vcf_input_file€€[0m is a placeholder for the VCF file being
  processed. €€[38;5;64m$annot_data_file€€[0m is the Ensembl Hg38 annotation file
  (€€[38;5;64mHomo_sapiens.GRCh38.108.gff3.gz€€[0m), and €€[38;5;64m$gnomad_data_file€€[0m is either an
  Ensembl-style (€€[38;5;64mHg38p14e.fa€€[0m) or UCSC-style (€€[38;5;64mHg38p14.fa€€[0m) Hg38 FASTA file. All
  three reference files can be found in the distributed data directory.


€€[1m€€[4m€€[38;5;61m[slivar]€€[0m <https://github.com/brentp/slivar>

€€[1m€€[4m€€[38;5;61m[the recommended rare-disease slivar filter]€€[0m
<https://github.com/brentp/slivar/wiki/rare-disease>


€€[1m€€[38;5;67m2 €€[0m€€[1m€€[38;5;67mRunning bcftools and slivar€€[0m

  €€[38;5;64mbcftools€€[0m and €€[38;5;64mslivar€€[0m can both be downloaded from the internet, but are also baked
  into this container.


€€[1m€€[38;5;67m3 €€[0m€€[1m€€[38;5;67mRunning the container image€€[0m

  This container can be run with any docker-archive compatible software (e.g.
  docker, podman, apptainer). Here I will use €€[38;5;64mpodman€€[0m.

  The image can be loaded with the €€[38;5;64mload€€[0m command, then run.

  gunzip apply-impprove.tar.gz
  podman load --input apply-impprove.tar
  podman run localhost/apply-impprove --help

  It is also possible to create an Apptainer/Singularity container from this
  container image tarball, see
  €€[1m€€[4m€€[38;5;61m€€[1m€€[4m€€[38;5;61m<https://apptainer.org/docs/user/latest/docker_and_oci.html#containers-in-docker-archive-files>€€[0m€€[0m
  for more information.

  apptainer build apply-impprove.sif docker-archive:apply-impprove.tar


€€[1m€€[38;5;67m4 €€[0m€€[1m€€[38;5;67mHelp€€[0m

  To see the help of the main entrypoint, simply run the container with the €€[38;5;64m--help€€[0m
  flag. As mentioned in the help information, this README can be retrieved by
  running the container with the argument €€[38;5;64m--readme€€[0m.


€€[1m€€[38;5;67m5 €€[0m€€[1m€€[38;5;67mSetup€€[0m

  In order to apply IMPPROVE random forest models, this container needs two sets
  of data:
  €€[38;5;61mâƒ €€[0mIMPPROVE models to apply, and
  €€[38;5;61mâƒ €€[0mA small collection of datasets needed to prepare the input

  The â€œsmall collectionâ€ of datasets consists of ~1 GB of tabular files found at
  €€[38;5;64mtkiHOHPC2111:/data/scratch/timothy/minimal_impprove_data/€€[0m, and a €€[38;5;64mcadd.hdf5€€[0m file
  (see €€[38;5;64mtkiHOHPC2111:/data/scratch/timothy/cadd_small.hdf5/€€[0m). The €€[38;5;64mcadd.hdf5€€[0m should
  be put in (not symlinked) the €€[38;5;64mminimal_impprove_data€€[0m folder. This can then be
  mounted as the €€[38;5;64m/data€€[0m volume.

  podman run €€[38;5;64m\€€[0m
      -v ~/Downloads/impprove_models/some_particular_set:/models €€[38;5;64m\€€[0m
      -v ~/Downloads/minimal_impprove_data:/data €€[38;5;64m\€€[0m
      -v ~/Documents/impprove_run/predictions:/predictions €€[38;5;64m\€€[0m
      -v ~/Documents/impprove_run/vcfs:/vcfs €€[38;5;64m\€€[0m
      localhost/apply-impprove

  Each €€[38;5;64m-v localpath:bindpath€€[0m argument creates a volume seen by the container at
  â€œbindpathâ€ from â€œlocalpathâ€. With apptainer €€[38;5;64m-B localpath:bindpath€€[0m can be used
  instead.

  The €€[38;5;64m/models€€[0m volume should be a directory of per-HPO models, and look something
  like this:

  /models
  â”œâ”€â”€ 2.jlso
  â”œâ”€â”€ 3.jlso
  â”œâ”€â”€ 8.jlso
  â‹®
  â””â”€â”€ 3000050.jlso

  A selection of pre-trained models are available under
  €€[38;5;64mtkiHOHPC2111:/data/scratch/timothy/impprove_pretrained_models€€[0m.


€€[1m€€[38;5;67m6 €€[0m€€[1m€€[38;5;67mOutput€€[0m

  For each input VCF, an output folder is created within the €€[38;5;64m/predictions€€[0m volume.
  Each folder will have a structure something like this:

  /predictions/input_file_name.vcf
  â”œâ”€â”€ all-predictions.csv
  â”œâ”€â”€ all-preds.png
  â”œâ”€â”€ hpo-descriptions.csv
  â”œâ”€â”€ model-scores.csv
  â”œâ”€â”€ top-50-preds.pngâ””
  â”œâ”€â”€ variable-importances.csv
  â””â”€â”€ wmean-good-prediction.csv

  €€[38;5;64mall-predictions.csv€€[0m will contain information on each variant of the VCF,
  starting with the Ensemble Gene ID and the gene name, and then after giving the
  variant itself (in terms of the chromosome, location, ref and alt) each
  subsequent â€œHP:XXXXXXXâ€ column gives the predictions from the model build for
  that HPO term.

  Each model used has an associated self-evaluation score from the bootstrapped
  test/train process, these scores are given in €€[38;5;64mmodel-scores.csv€€[0m. The variable
  importances for each model are found in €€[38;5;64mvariable-importances.csv€€[0m. Each model is
  labelled according to its HPO id, descriptions for the HPO ids used can be found
  in €€[38;5;64mhpo-descriptions.csv€€[0m.

  The score-weighted average prediction across models with a score of at least 0.4
  is given in €€[38;5;64mwmean-good-prediction.csv€€[0m.

  €€[38;5;64mtop-50-preds.png€€[0m is a heatmap of all predictions for the top 50 variants
  overall, while €€[38;5;64mall-preds.png€€[0m is the heatmap for everything.


€€[1m€€[38;5;67m7 €€[0m€€[1m€€[38;5;67mSample€€[0m

  A small sample VCF is bundled with this container for testing, and can be
  fetched by running the container with the argument €€[38;5;64m--sample€€[0m. It contains three
  known pathogenic variants, and 5000 benign variants.
