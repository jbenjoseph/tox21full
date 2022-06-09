Tox21 Full
~~~~~~~~~~

This is a Tox21-like dataset created from the raw NIH assay data. The Tox21 dataset we know and love only includes 12 assays. This one includes 64 different assays!

Simply by using this dataset you can train machine learning models with metrics you never thought was possible.


Downloads NIH raw assay data and creates a clean CSV.GZ file ready for import into pandas:

::

    tox21full ~/Downloads/tox21full.csv.gz


You can also create it as a parquet file (more efficent):

:: 

    tox21full --format parquet ~/Downloads/tox21full.parquet