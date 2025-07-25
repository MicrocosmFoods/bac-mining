
from dataclasses import dataclass
import typing
import typing_extensions

from flytekit.core.annotation import FlyteAnnotation

from latch.types.metadata import NextflowParameter
from latch.types.file import LatchFile
from latch.types.directory import LatchDir, LatchOutputDir

# Import these into your `__init__.py` file:
#
# from .parameters import generated_parameters

generated_parameters = {
    'input_genomes': NextflowParameter(
        type=LatchDir,
        default=None,
        section_title='Input/output options',
        description='Input directory of genomes in fasta format ending in .fa',
    ),
    'genome_list': NextflowParameter(
        type=LatchFile,
        default=None,
        section_title=None,
        description='List of subset genomes to run the pipeline through that are in the input_genomes directory',
    ),
    'outdir': NextflowParameter(
        type=typing_extensions.Annotated[LatchDir, FlyteAnnotation({'output': True})],
        default=None,
        section_title=None,
        description='The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.',
    ),
    'antismash_db': NextflowParameter(
        type=LatchDir,
        default=None,
        section_title='Databases',
        description='Path to directory of pre-downloaded antismash databases',
    ),
    'kofam_db': NextflowParameter(
        type=LatchDir,
        default=None,
        section_title=None,
        description='Path to directory of Kofam KEGG HMM database',
    ),
    'threads': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title='Other',
        description=None,
    ),
    'functional_annotation': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Optionally run kofamscan functional annotation, default set to false. To run annotation set to true.',
    ),
    'smorfinder_mode': NextflowParameter(
        type=typing.Optional[str],
        default='pre_called',
        section_title=None,
        description='Mode for running smorfinder. "single" runs smorfinder directly on genomes, "pre_called" uses pre-called genes from pyrodigal.',
    ),
}

