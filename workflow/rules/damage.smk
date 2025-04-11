rule MapDamage:
    output:
        dir=directory("results/MAPDAMAGE/{sample}"),
    input:
        species_tax_id="results/KRAKENUNIQ/{sample}/taxID.species",
        bam="results/BOWTIE2/{sample}/AlignedToBowtie2DB.bam",
        bai="results/BOWTIE2/{sample}/AlignedToBowtie2DB.bam.bai",
    params:
        bowtie2_seqid2taxid=config["bowtie2_seqid2taxid_db"],
        BOWTIE2_DB=config["bowtie2_db"],
        options=config["options"].get("MapDamage", ""),
    threads: 10
    log:
        "logs/MAPDAMAGE/{sample}.log",
    conda:
        "../envs/mapdamage.yaml"
    envmodules:
        *config["envmodules"]["mapdamage"],
    benchmark:
        "benchmarks/MAPDAMAGE/{sample}.benchmark.txt"
    message:
        "MapDamage: RUNNING MAPDAMAGE ON SPECIES IDENTIFIED BY KRAKENUNIQ IN SAMPLE {input.bam}"
    shell:
        "mkdir {output.dir}; "
        "if [ -s {input.species_tax_id} ]; then "
        'cat {input.species_tax_id} | parallel -j {threads} "grep -w {{}} {params.bowtie2_seqid2taxid} | cut -f1 > {output.dir}/{{}}.seq_ids" ; '
        "for i in $(cat {input.species_tax_id}); do if [[ $(wc -l < {output.dir}/${{i}}.seq_ids) -ge 0 ]]; then xargs --arg-file={output.dir}/${{i}}.seq_ids samtools view -bh {input.bam} --write-index -@ {threads} -o {output.dir}/${{i}}.tax.bam; fi; done >> {log} 2>&1; "
        "find {output.dir} -name '*.tax.bam' | parallel -j {threads} \"mapDamage {params.options} -i {{}} -r {params.BOWTIE2_DB} --merge-reference-sequences -d {output.dir}/mapDamage_{{}}\" >> {log} 2>&1 || true; "
        "for filename in {output.dir}/*.tax.bam; do newname=`echo $filename | sed 's/tax\.//g'`; mv $filename $newname; done >> {log} 2>&1; "
        "mv {output.dir}/mapDamage_{output.dir}/* {output.dir} >> {log} 2>&1; "
        "rm -r {output.dir}/mapDamage_results >> {log} 2>&1; "
        "else echo NO MICROBES TO AUTHENTICATE > {output.dir}/README.txt; fi"
