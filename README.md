# Basestack consensus

This pipeline is an implementation of the ARTIC network pipeline that is deployed in a Docker container that can be easily installed using the Basestack platform (https://github.com/jhuapl-bio/Basestack).

## Normal usage

Install Basestack as described here: https://github.com/jhuapl-bio/Basestack and install the 'Consensus' module.

## Development usage

The pipeline can be directly run, for development or other command line use, with the following commands:

1. Download the latest version of the pipeline with:

```sh
docker pull jhuaplbio/basestack_consensus:latest
```

2. Connect a MinION sequencing run folder to a container of that image:

```sh
docker container run \
  --rm \
  -it \
  -v /path/to/local/MinION/data:/opt/basestack_consensus/sequencing_runs/example-run \
  jhuaplbio/basestack_consensus:latest \
  bash
```

3. Within the container, run the pipeline on the MinION run folder:
```sh
artic-module1-barcode-demux.sh -i /opt/basestack_consensus/sequencing_runs/example-run
```


# misc

To run a dev process on the image on a mounted path, do: 

`docker container run -it --rm -v $ROOT_DATA_PATH:/opt/run-data -v $PWD:/opt/basestack_consensus/code/ncov --name basestack_consensus jhuaplbio/basestack_consensus`

where `$ROOT_DATA_PATH` is the path to the run directory. The pipeline will output to `artic-pipeline` in that folder