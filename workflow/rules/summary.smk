rule Plot_Authentication_Score:
    output:
        heatmap="results/overview_heatmap_scores.pdf",    
    input:
        scores=expand("results/AUTHENTICATION/.{sample}_done",sample=SAMPLES)
    message:
        "Plot_Authentication_Score: PLOTTING HEATMAP OF AUTHENTICATION SCORES"
    params:
        exe=WORKFLOW_DIR / "scripts/plot_score.R",
    log:
        "logs/PLOT_AUTHENTICATION_SCORE/plot_authentication_score.log",
    conda:
        "../envs/r.yaml"
    envmodules:
        *config["envmodules"]["r"],
    shell:
        "Rscript {params.exe} results/AUTHENTICATION $(dirname {output.heatmap}) &> {log}"

