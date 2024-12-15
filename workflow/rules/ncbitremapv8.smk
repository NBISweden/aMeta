rule NCBIMapTre:
    """Download ncbi.map and ncbi.tre from https://github.com/husonlab/megan-ce/tree/master/src/megan/resources/files"""
    output:
        tre=os.path.join(config["ncbi_db"], "ncbi.tre"),
        map=os.path.join(config["ncbi_db"], "ncbi.map"),
    input:
        tre=storage(
            "https://github.com/husonlab/megan-ce/raw/master/src/megan/resources/files/ncbi.tre",
        ),
        map=storage(
            "https://github.com/husonlab/megan-ce/raw/master/src/megan/resources/files/ncbi.map",
        ),
    log:
        "logs/NCBI/ncbi.log",
    threads: 1
    shell:
        "mv {input.tre} {output.tre};"
        "mv {input.map} {output.map};"
