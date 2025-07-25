import sys
from dataclasses import dataclass
from enum import Enum
import os
import subprocess
import requests
import shutil
from pathlib import Path
import typing
import typing_extensions

from latch.resources.workflow import workflow
from latch.resources.tasks import nextflow_runtime_task, custom_task
from latch.types.file import LatchFile
from latch.types.directory import LatchDir, LatchOutputDir
from latch.ldata.path import LPath
from latch.executions import report_nextflow_used_storage
from latch_cli.nextflow.workflow import get_flag
from latch_cli.nextflow.utils import _get_execution_name
from latch_cli.utils import urljoins
from latch.types import metadata
from flytekit.core.annotation import FlyteAnnotation

from latch_cli.services.register.utils import import_module_by_path

meta = Path("latch_metadata") / "__init__.py"
import_module_by_path(meta)
import latch_metadata

@custom_task(cpu=0.25, memory=0.5, storage_gib=1)
def initialize() -> str:
    token = os.environ.get("FLYTE_INTERNAL_EXECUTION_ID")
    if token is None:
        raise RuntimeError("failed to get execution token")

    headers = {"Authorization": f"Latch-Execution-Token {token}"}

    print("Provisioning shared storage volume... ", end="")
    resp = requests.post(
        "http://nf-dispatcher-service.flyte.svc.cluster.local/provision-storage-ofs",
        headers=headers,
        json={
            "storage_expiration_hours": 15,
            "version": 2,
            "fs_size_tb": 45,
        },
    )
    resp.raise_for_status()
    print("Done.")

    return resp.json()["name"]






@nextflow_runtime_task(cpu=16, memory=50, storage_gib=200)
def nextflow_runtime(pvc_name: str, input_genomes: LatchDir, genome_list: typing.Optional[LatchFile], outdir: typing_extensions.Annotated[LatchDir, FlyteAnnotation({'output': True})], antismash_db: LatchDir, kofam_db: LatchDir, threads: typing.Optional[str], functional_annotation: typing.Optional[bool], smorfinder_mode: typing.Optional[str]) -> None:
    shared_dir = Path("/nf-workdir")

    exec_name = _get_execution_name()
    if exec_name is None:
        print("Failed to get execution name.")
        exec_name = "unknown"

    latch_log_dir = urljoins("latch:///your_log_dir/nf_bacmagmining", exec_name)
    print(f"Log directory: {latch_log_dir}")




    ignore_list = [
        "latch",
        ".latch",
        ".git",
        "nextflow",
        ".nextflow",
        "work",
        "results",
        "miniconda",
        "anaconda3",
        "mambaforge",
    ]

    shutil.copytree(
        Path("/root"),
        shared_dir,
        ignore=lambda src, names: ignore_list,
        ignore_dangling_symlinks=True,
        dirs_exist_ok=True,
    )

    profile_list = ['docker']
    if False:
        profile_list.extend([p.value for p in execution_profiles])

    if len(profile_list) == 0:
        profile_list.append("standard")

    profiles = ','.join(profile_list)

    antismash_db_dir = Path(antismash_db.local_path)
    antismash_shared_dir = shared_dir / antismash_db_dir.name
    shutil.move(str(antismash_db_dir), str(antismash_shared_dir))
 
    kofam_db_dir = Path(kofam_db.local_path)
    kofam_shared_dir = shared_dir / kofam_db_dir.name
    shutil.move(str(kofam_db_dir), str(kofam_shared_dir))

    cmd = [
        "/root/nextflow",
        "run",
        str(shared_dir / "main.nf"),
        "-work-dir",
        str(shared_dir),
        "-profile",
        profiles,
        "-c",
        "latch.config",
        "-resume",
                *get_flag('input_genomes', input_genomes),
                *get_flag('genome_list', genome_list),
                *get_flag('outdir', outdir),
                *get_flag('antismash_db', antismash_shared_dir),
                *get_flag('kofam_db', kofam_shared_dir),
                *get_flag('threads', threads),
                *get_flag('functional_annotation', functional_annotation),
                *get_flag('smorfinder_mode', smorfinder_mode)
    ]

    print("Launching Nextflow Runtime")
    print(' '.join(cmd))
    print(flush=True)

    failed = False
    try:
        env = {
            **os.environ,
            "NXF_ANSI_LOG": "false",
            "NXF_HOME": "/root/.nextflow",
            "NXF_OPTS": "-Xms9600M -Xmx38400M -XX:ActiveProcessorCount=16",
            "NXF_DISABLE_CHECK_LATEST": "true",
            "NXF_ENABLE_VIRTUAL_THREADS": "false",
            "NXF_ENABLE_FS_SYNC": "true",
        }

        if False:
            env["LATCH_LOG_DIR"] = latch_log_dir

        subprocess.run(
            cmd,
            env=env,
            check=True,
            cwd=str(shared_dir),
        )
    except subprocess.CalledProcessError:
        failed = True
    finally:
        print()

        nextflow_log = shared_dir / ".nextflow.log"
        if nextflow_log.exists():
            remote = LPath(urljoins(latch_log_dir, "nextflow.log"))
            print(f"Uploading .nextflow.log to {remote.path}")
            remote.upload_from(nextflow_log)

        print("Computing size of workdir... ", end="")
        try:
            result = subprocess.run(
                ['du', '-sb', str(shared_dir)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=5 * 60
            )

            size = int(result.stdout.split()[0])
            report_nextflow_used_storage(size)
            print(f"Done. Workdir size: {size / 1024 / 1024 / 1024: .2f} GiB")
        except subprocess.TimeoutExpired:
            print("Failed to compute storage size: Operation timed out after 5 minutes.")
        except subprocess.CalledProcessError as e:
            print(f"Failed to compute storage size: {e.stderr}")
        except Exception as e:
            print(f"Failed to compute storage size: {e}")

    if failed:
        sys.exit(1)


@workflow(metadata._nextflow_metadata)
def nf_bacmagmining(input_genomes: LatchDir, genome_list: typing.Optional[LatchFile], outdir: typing_extensions.Annotated[LatchDir, FlyteAnnotation({'output': True})], antismash_db: LatchDir, kofam_db: LatchDir, threads: typing.Optional[str], functional_annotation: typing.Optional[bool], smorfinder_mode: typing.Optional[str]) -> None:
    """
    bacMAGmining

    Sample Description
    """

    pvc_name: str = initialize()
    nextflow_runtime(pvc_name=pvc_name, input_genomes=input_genomes, genome_list=genome_list, outdir=outdir, antismash_db=antismash_db, kofam_db=kofam_db, threads=threads, functional_annotation=functional_annotation, smorfinder_mode=smorfinder_mode)

