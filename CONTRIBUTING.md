# Submit a bug report

*First look in the all open and closed issues for the error message you got.*

As atlas is based on snakemake, check if you find a similar bug already discussed for other snakemake-workflows.

If you don't find any help. Submit a issue:

- specify the system you are working on: linux, cluster, shared filesystem?
- copy the error message
- add the log file of the rule, which produced the error.
- If you run the atlas on a cluster join also the log file of the cluster.

I hope we can help you...

# Contribute to the metagenome-atlas code

## Prerequisites

- know the basic about git and gitHub
- know how snakemake works, otherwise check the [tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html)

## Setup

You can ask the maintainers to be added to the repository and work from a *branch* of the main atlas repository or you can work from a fork of the atlas repository.

Follow the [steps](https://github.com/metagenome-atlas/atlas#install-the-development-version-from-github) to set up the development version of atlas. This allows you to work with the code you have in the git repository.

## Test the code

### Locally

Idelly you should have some test prpject on your local machine.
When you created a new rule and you want to test the output of this rule `my_target.tsv` you can do this by running:

``` atlas run None my_target.tsv ```

### Continuous integration

When you make a pull request to the master branch. Each change in your code gets checked by continuous integration (CI). The tests should make sure that your modification don't break any other use of atlas. However due to the requeirements needed during the execution of atlas, it is not possible to test all functionalities via CI. If you add functionalities to atlas, they should also be tested. Have a look at the scripts in `.test`.
