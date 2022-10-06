---
name: Workflow Error
about: I got atlas running but then it encountered this error. For other errors see
  the other template.
title: Error in rule
labels: ''
assignees: ''

---

- [ ] I checked and didn't found a related issue,e.g. while typing the title
- [ ] ** I got an error in the following rule(s):** `   `

<!-- The error messages are highlighted in red in the stdout. Otherwise you can get all the errors occurred in the workflow by running the following command in the atlas working directory:
`grep  -A 8 rule $(ls -t .snakemake/log/* | head -1)`
  -->



- [ ] I checked the log files indicated indicated in the error message (and the cluster logs if submitted to a cluster)

<!-- 
The error message should look something like the rule below 
```

Error in rule run_das_tool:
    jobid: 91
    output: sample2/binning/DASTool/sample2_DASTool_summary.txt, sample2/binning/DASTool/cluster_attribution.tsv
    log: sample2/logs/binning/DASTool.log (check log file(s) for error message)
    conda-env: /data/users/ocetiner/ABERRANT_data/test3/databases/conda_envs/83d50c467b0444c29b1510fd5b6ba829
    shell:
         DAS_Tool --outputbasename sample2/binning/DASTool/sample2  --bins 
        (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)
    cluster_jobid: 7010304

Error executing rule run_das_tool on cluster (jobid: 91, external: 7010304). For error details see the cluster log and the log files of the involved rule(s).

```
Do what is says **For error details see the cluster log and the log files of the involved rule(s).**

cluster logs are usually in `cluster_log/*external_id*` or so.

-->

Here is the relevant log output:

```



```


** Atlas version**

**Additional context**
Add any other context about the problem here.
