conda_version=$(conda --version | awk '{print $2}')
conda_major=$(echo $conda_version | awk -F. '{print $1}')
conda_minor=$(echo $conda_version | awk -F. '{print $2}')

if [[ $conda_major -lt 23 ]] || [[ $conda_major -lt 24 && $conda_minor -lt 10 ]]; then
    CONDA_FRONTEND_DEFAULT=mamba
else
    CONDA_FRONTEND_DEFAULT=conda
fi

CONDA_FRONTEND="${CONDA_FRONTEND:=${CONDA_FRONTEND_DEFAULT}}"
echo "Conda version ${conda_version}: using conda frontend $CONDA_FRONTEND"

if [[ -z "$@" ]]; then
    echo "No snakemake options supplied. To run the tests, at least set the number of threads (e.g. -j 1)"
    echo "Example usage: "
    echo
    echo "./runtest.sh -j 1"
    echo
    echo "UPDATE: Starting from conda version 23.10.0, by default conda uses the libmamba-solver."
    echo "If you have an older conda version, consider running"
    echo
    echo "CONDA_FRONTEND=mamba ./runtest.sh -j 1"
    echo
    exit
fi
# check if we are running CI; if not, initialize databases
if [[ -z "$CI" ]]; then
    if [[ ! -e ".initdb" ]]; then
        echo "This looks like the first test run... Installing bioconda packages..."
        snakemake --conda-frontend $CONDA_FRONTEND --use-conda --show-failed-logs -j 1 --conda-cleanup-pkgs cache --conda-create-envs-only -s ../workflow/Snakefile

        source $(dirname $(dirname $CONDA_EXE))/etc/profile.d/conda.sh

        ##############################
        # Krakenuniq database
        ##############################
        echo Building krakenuniq data
        env=$(grep krakenuniq .snakemake/conda/*yaml | awk '{print $1}' | sed -e "s/.yaml://g")
        conda activate $env
        krakenuniq-build --db resources/KrakenUniq_DB --kmer-len 21 --minimizer-len 11 --jellyfish-bin $(pwd)/$env/bin/jellyfish
        conda deactivate

        ##############################
        # Krona taxonomy
        ##############################
        echo Building krona taxonomy
        env=$(grep krona .snakemake/conda/*yaml | awk '{print $1}' | sed -e "s/.yaml://g" | head -1)
        conda activate $env
        cd $env/opt/krona
        ./updateTaxonomy.sh taxonomy
        cd -
        conda deactivate

        ##############################
        # Adjust malt max memory usage
        ##############################
        echo Adjusting malt max memory usage
        env=$(grep hops .snakemake/conda/*yaml | awk '{print $1}' | sed -e "s/.yaml://g" | head -1)
        conda activate $env
        version=$(conda list malt --json | grep version | sed -e "s/\"//g" | awk '{print $2}')
        cd $env/opt/malt-$version
        sed -i -e "s/-Xmx64G/-Xmx3G/" malt-build.vmoptions
        sed -i -e "s/-Xmx64G/-Xmx3G/" malt-run.vmoptions
        cd -
        conda deactivate

        touch .initdb
    fi
fi


echo Running workflow excluding bowtie...
echo snakemake --conda-frontend $CONDA_FRONTEND --use-conda --show-failed-logs --conda-cleanup-pkgs cache -s ../workflow/Snakefile $@
snakemake --conda-frontend $CONDA_FRONTEND --use-conda --show-failed-logs --conda-cleanup-pkgs cache -s ../workflow/Snakefile --configfile config/config.bowtie.yaml $@

if [ $? -ne 0 ]; then
    echo ERROR: Workflow test failed!
    exit
fi


echo Running workflow...
echo snakemake --conda-frontend $CONDA_FRONTEND --use-conda --show-failed-logs --conda-cleanup-pkgs cache -s ../workflow/Snakefile $@
snakemake --conda-frontend $CONDA_FRONTEND --use-conda --show-failed-logs --conda-cleanup-pkgs cache -s ../workflow/Snakefile $@

if [ $? -ne 0 ]; then
    echo ERROR: Workflow test failed!
    exit
fi

echo Generating report...
echo snakemake --conda-frontend $CONDA_FRONTEND -s ../workflow/Snakefile --report --report-stylesheet ../workflow/report/custom.css
snakemake --conda-frontend $CONDA_FRONTEND -s ../workflow/Snakefile --report --report-stylesheet ../workflow/report/custom.css

if [ $? -ne 0 ]; then
    echo ERROR: Report test failed!
    exit
fi

# Add test for missing taxid
echo Running test for missing taxid...
echo 'echo -e "632\n42862" >> results/KRAKENUNIQ/foo/taxID.species'
echo -e "632\n42862" >> results/KRAKENUNIQ/foo/taxID.species
echo snakemake --conda-frontend $CONDA_FRONTEND --use-conda --show-failed-logs --conda-cleanup-pkgs cache -s ../workflow/Snakefile $@
snakemake --conda-frontend $CONDA_FRONTEND --use-conda --show-failed-logs --conda-cleanup-pkgs cache -s ../workflow/Snakefile $@
rm results/KRAKENUNIQ/foo/taxID.species
