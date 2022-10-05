# Pharmacogenetics

## INTRODUCTION

The **PharmProfiler** aims to provide a convenient way to extract pharmacogenetic data from your patient of choice.
It includes a method to create a database based on the PharmGKB website, and methods to use this database to create a report for your individual.

## PRE-REQUISITES
* You will require a [PharmGKB](https://www.pharmgkb.org/) account.
* Python 2.7 with the following packages:
  * BioPython >= 0.6.2
  * tqdm >= 4.19.4
  * requests_oauthlib >= 0.8.0
  * Flask
  * Clustalo
  * FastTree

## BUILDING THE DATABASE

The following command can be used to get general input help:

```bash
python pharmprofiler.py -h
```

To create all necessary tables in the proper order, execute the following command:

```bash
python pharmprofiler.py -b all
```

This will prompt your PharmGKB e-mail address and password multiple times.
You can also (re)create individual tables replacing '-t all' with various forms of -t YOURTABLE'.
(options include: pairs, genes, variants, haplotypes, guidelines, annotations and more.)
If you do this route, keep in mind that the following order is required (all valid options for the '-b' command).

_pairs > drugs > drugmatches > drugvars > genes > haplotypes > hapvars > etcvars > indels > annotations > guidelines > guideoptions_

This means that the table before the > is necessary to build the table after the >.

## PRE-PROCESSING

In order to process patient data, a GVCF file is required.
This file can be generated through the [GATK](https://software.broadinstitute.org/gatk/) pipeline.
To execute this query, you will require a BED file.
To create this file, please first complete database creation
and then execute the following command.

```bash
python pharmprofiler.py -b bed
```

It will create a file that GATK can use to determine the gene locations to search in.
Please execute the GATK VariantCaller with the options --emitRefConfidence 'BP_RESOLUTION' and -L BED_FILE_LOCATION to create the GVCF file needed for further steps.

## IMPORTING GVCF

Once you have your GVCF file and your **PharmProfiler** database has been created, you can proceed with importing your individual's data.
If you want to do everything in one row, including generating a report, execute the following:

```bash
python pharmprofiler.py -p GVCF_LOCATION -a all
```

You can also execute steps separately, from importing to haplotyping to finding matching annotations.

```bash
python pharmprofiler.py -p GVCF_LOCATION -a YOUR_ACTION
```

If you do this route, keep in mind that the following order is required (all valid options for the '-a' command).

_import > haplotype > score > genotype > guidelines > annotate > report_

This means that the action before the > is necessary to perform the action after the >.

## REPORT generation
Using the command:

```
python pharmprofiler.py -p GVCF_LOCATION -a report
```

After executing all other patient options as listed above will net you a Latex file with a report in it. You will require TexMaker or another latex renderer to convert this to PDF. Other formats will be added in the future. ( ;) )

------
### Copyright notice
Copyright 2016 Joanna Wolthuis & Joep de Ligt
