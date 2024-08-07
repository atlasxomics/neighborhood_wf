import subprocess

from enum import Enum
from typing import List

from latch import custom_task, workflow
from latch.types import (
    LatchAuthor,
    LatchDir,
    LatchMetadata,
    LatchParameter,
    LatchRule
)

from wf.upload_to_registry import upload_to_registry, Run


class Genome(Enum):
    mm10 = 'mm10'
    hg38 = 'hg38'
    rnor6 = 'rnor6'


@custom_task(cpu=8, memory=384, storage_gib=500)
def neighborhood_task(
    runs: List[Run],
    project_name: str,
    genome: Genome
) -> LatchDir:

    _r_cmd = ['Rscript', '/root/wf/neighborhood.R', project_name, genome.value]

    runs = [
        (
            f'{run.run_id},'
            f'{run.seurat_obj.local_path},'
            f'{run.condition},'
            f'{run.positions_file.local_path},'
            f'{run.spatial_dir.local_path},'
        )
        for run in runs
    ]

    _r_cmd.extend(runs)
    subprocess.run(_r_cmd)

    _mv_cmd = ['mv', '/root/neighborhood/figures', f'/root/{project_name}']

    subprocess.run(_mv_cmd)

    return LatchDir(
        f'/root/{project_name}',
        f'latch:///Neighborhood_Analysis_outputs/{project_name}'
    )


metadata = LatchMetadata(
    display_name='neighborhood Analysis',
    author=LatchAuthor(
        name='AtlasXomics, Inc.',
        email='noorisotudeh@gmail.com',
        github='https://github.com/atlasxomics/neighborhood_wf',
    ),
    repository='https://github.com/atlasxomics/neighborhood_wf',
    license='MIT',
    parameters={
        'runs': LatchParameter(
            display_name='runs',
            description='List of runs to be analyzed; each run must contain a \
                        run_id and seurat rds file; optional: condition, \
                        tissue position file for filtering on/off tissue, \
                        spatial folder for SpatialDimPlot. Note that multiple \
                        Coditions must be separted by -, for example: \
                        Non_responder-post_treatment-One_month.',
            batch_table_column=True,
            samplesheet=True
        ),
        'project_name': LatchParameter(
            display_name='project name',
            description='Name of output directory in \
                        Neighborhood_Analysis_outputs',
            batch_table_column=True,
            rules=[
                LatchRule(
                    regex="^[^/].*",
                    message="project name cannot start with a '/'"
                )
            ]
        ),
        'genome': LatchParameter(
            display_name='genome',
            description='Reference genome to be used for geneAnnotation and \
                        genomeAnnotation',
            batch_table_column=True,
        ),
        'run_table_id': LatchParameter(
            display_name='Registry Table ID',
            description='The runs will be updated in Registry with its \
                        corresponding outs.'
        ),
        'project_table_id': LatchParameter(
            display_name='The ID of the SOWs Registry table',
            description='The outs will be inserted into the SOW \
                        table for the corresponding runs.'
        )
    },
    tags=[],
)


@workflow(metadata)
def neighborhoodanalysis_workflow(
    runs: List[Run],
    genome: Genome,
    project_name: str,
    run_table_id: str = '761',
    project_table_id: str = '779'
) -> LatchDir:
    '''Workflow for neighborhood analysis.

    # neighborhood Analysis

    **neighborhood Analysis** is a [latch.bio](https://latch.bio/) workflow for
    generating R objects and data for downstream analysis of epigenomic
    [Neighborhood Analysis](https://www.cell.com/cell-metabolism/fulltext/S1550-4131(21)00363-6?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS1550413121003636%3Fshowall%3Dtrue).

    > This workflow does not support combined feature count matrices with >
    2^31 - 1 elements (~100,000 cells); we are working to add support for large
    matrices in a future version.

    ## Inputs
    workflow takes the following parameters,
    * Seurat object file: A seurat object created by [create ArchRProject](https://wiki.latch.bio/workflows/overview) workflow per sample
    * [tissue_positions_list.csv](https://docs.atlasxomics.com/projects/AtlasXbrowser/en/latest/SpatialFolder.html): A comma-separated file in which each row contains a unique barcode, an indicator for whether the tixel is 'on-tissue' (1, 0), and a row/column index
    * [Spatial folder](https://docs.atlasxomics.com/projects/AtlasXbrowser/en/latest/SpatialFolder.html): A directory containing tissue images and experiment metadata
    * Run ID: An identifier for the run
    * Condition (_optional_):  An experimental Condition descriptor
    (ie. 'control', 'diseased')
    Individual runs are batched in a Project with the following global
    parameters,
    * Project Name: A name for the output folder
    * Genome: A reference genome to be used for alignment
    * Upload to SLIMS _(ATX-internal only)_: A T/F toggle for whether to push
    QC results to LIMS

    ## Running the workflow (_in progress_)

    ## Outputs (_in progress_)
    ## Next Steps
    ## Support
    Questions? Comments?  Contact support@atlasxomics.com or post in
    AtlasXomics [Discord](https://discord.com/channels/1004748539827597413/1005222888384770108).

    '''
    project = neighborhood_task(
        runs=runs,
        project_name=project_name,
        genome=genome
    )

    upload_to_registry(
        runs=runs,
        project=project,
        run_table_id=run_table_id,
        project_table_id=project_table_id
    )

    return project
