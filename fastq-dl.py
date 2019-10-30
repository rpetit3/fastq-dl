#! /usr/bin/env python3
"""
usage: fastq-dl [-h] [--aspera STRING] [--aspera_key STRING]
                [--aspera_speed STRING] [--is_study] [--is_experiment]
                [--is_run] [--group_by_experiment] [--group_by_sample]
                [--outdir OUTPUT_DIR] [--max_attempts INT] [--cpus INT]
                [--ftp_only] [--silent] [--debug] [--version]
                ACCESSION {sra,SRA,ena,ENA}

fastq-dl - Download FASTQs from ENA or SRA

optional arguments:
  -h, --help            show this help message and exit

Required Options:

  ACCESSION             ENA/SRA accession to query. (Study, Experiment, or Run
                        accession)
  {sra,SRA,ena,ENA}     Specify which provider (ENA or SRA) to use. Accepted
                        Values: ENA SRA

Aspera Connect Options:
  --aspera STRING       Path to the Aspera Connect tool "ascp" (Default:
                        "which ascp")
  --aspera_key STRING   Path to Aspera Connect private key, if not given,
                        guess based on ascp path
  --aspera_speed STRING
                        Speed at which Aspera Connect will download. (Default:
                        100M)

Query Related Options:
  --is_study            Query is a Study.
  --is_experiment       Query is an Experiment.
  --is_run              Query is a Run.
  --group_by_experiment
                        Group Runs by experiment accession.
  --group_by_sample     Group Runs by sample accession.

Helpful Options:
  --outdir OUTPUT_DIR   Directory to output downloads to. (Default: ./)
  --max_attempts INT    Maximum number of download attempts (Default: 10)
  --cpus INT            Total cpus used for downloading from SRA (Default: 1)
  --ftp_only            FTP only downloads.
  --silent              Only critical errors will be printed.
  --debug               Skip downloads, print what will be downloaded.
  --version             show program's version number and exit

Example:
fastq-dl SRX477044 SRA
"""
import logging
import os
import time

PROGRAM = "fastq-dl"
VERSION = "1.0.2"
STDOUT = 11
STDERR = 12
logging.addLevelName(STDOUT, "STDOUT")
logging.addLevelName(STDERR, "STDERR")


ENA_URL = ('https://www.ebi.ac.uk/ena/data/warehouse/search?result=read_run&'
           'display=report')
FIELDS = [
    'study_accession', 'secondary_study_accession', 'sample_accession',
    'secondary_sample_accession', 'experiment_accession', 'run_accession',
    'submission_accession', 'tax_id', 'scientific_name',
    'instrument_platform', 'instrument_model', 'library_name',
    'library_layout', 'nominal_length', 'library_strategy',
    'library_source', 'library_selection', 'read_count',
    'base_count', 'center_name', 'first_public', 'last_updated',
    'experiment_title', 'study_title', 'study_alias', 'experiment_alias',
    'run_alias', 'fastq_bytes', 'fastq_md5', 'fastq_ftp', 'fastq_aspera',
    'fastq_galaxy', 'submitted_bytes', 'submitted_md5', 'submitted_ftp',
    'submitted_aspera', 'submitted_galaxy', 'submitted_format',
    'sra_bytes', 'sra_md5', 'sra_ftp', 'sra_aspera', 'sra_galaxy',
    'cram_index_ftp', 'cram_index_aspera', 'cram_index_galaxy',
    'sample_alias', 'broker_name', 'sample_title', 'nominal_sdev',
    'first_created'
]

def set_log_level(error, debug):
    """Set the output log level."""
    return logging.ERROR if error else logging.DEBUG if debug else logging.INFO


def get_log_level():
    """Return logging level name."""
    return logging.getLevelName(logging.getLogger().getEffectiveLevel())


def execute(cmd, directory=os.getcwd(), capture_stdout=False, stdout_file=None,
            stderr_file=None, max_attempts=1, is_sra=False):
    """A simple wrapper around executor."""
    from executor import ExternalCommand, ExternalCommandFailed
    attempt = 0
    while attempt < max_attempts:
        attempt += 1
        try:
            command = ExternalCommand(
                cmd, directory=directory, capture=True, capture_stderr=True,
                stdout_file=stdout_file, stderr_file=stderr_file
            )

            command.start()
            if get_log_level() == 'DEBUG':
                logging.log(STDOUT, command.decoded_stdout)
                logging.log(STDERR, command.decoded_stderr)

            if capture_stdout:
                return command.decoded_stdout
            else:
                return command.returncode
        except ExternalCommandFailed as error:
            logging.error(f'"{cmd}" return exit code {command.returncode}')
            

            if is_sra and command.returncode == 3:
                # The FASTQ isn't on SRA for some reason, try to download from ENA
                error_msg = command.decoded_stderr.split("\n")[0]
                logging.error(error_msg)
                return 'SRA_NOT_FOUND'

            if attempt < max_attempts:
                logging.error(f'Retry execution ({attempt} of {max_attempts})')
                time.sleep(10)
            else:
                raise error


def sra_download(accession, outdir, cpus=1, max_attempts=10):
    """Download FASTQs from SRA using fasterq-dump."""
    fastqs = {'r1':'', 'r2': '', 'single_end': True}
    se = f'{outdir}/{accession}.fastq.gz'
    pe = f'{outdir}/{accession}_2.fastq.gz'

    if not os.path.exists(se) and not os.path.exists(pe):
        execute(f'mkdir -p {outdir}')
        outcome = execute(f'fasterq-dump {accession} --split-files --threads {cpus}',
                            max_attempts=max_attempts, directory=outdir, is_sra=True)
        if outcome == "SRA_NOT_FOUND":
            return outcome
        else:
            execute(f'pigz -p {cpus} -n --fast *.fastq', directory=outdir)

    if os.path.exists(f'{outdir}/{accession}_2.fastq.gz'):
        # Paired end
        fastqs['r1'] = f'{outdir}/{accession}_1.fastq.gz'
        fastqs['r2'] = f'{outdir}/{accession}_2.fastq.gz'
        fastqs['single_end'] = False
    else:
        fastqs['r1'] = f'{outdir}/{accession}.fastq.gz'

    return fastqs


def ena_download(run, outdir, aspera=None, max_attempts=10, ftp_only=False):
    fastqs = {'r1':'', 'r2': '', 'single_end': True}
    fasp = run['fastq_aspera'].split(';')
    ftp = run['fastq_ftp'].split(';')
    md5 = run['fastq_md5'].split(';')
    for i in range(len(fasp)):
        is_r2 = False
        # If run is paired only include *_1.fastq and *_2.fastq, rarely a
        # run can have 3 files.
        # Example:ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR114/007/ERR1143237
        if run['library_layout'] == 'PAIRED':
            if fasp[i].endswith('_2.fastq.gz'):
                # Example: ERR1143237_2.fastq.gz
                is_r2 = True
            elif fasp[i].endswith('_1.fastq.gz'):
                # Example: ERR1143237_1.fastq.gz
                pass
            else:
                # Example: ERR1143237.fastq.gz
                # Not apart of the paired end read, so skip this file. Or,
                # its the only fastq file, and its not a paired
                obs_fq = os.path.basename(fasp[i])
                exp_fq = f'{run["run_accession"]}.fastq.gz'
                if (len(fasp) != 1 and obs_fq != exp_fq):
                    continue

        # Download Run
        if md5[i]:
            fastq = download_ena_fastq(
                fasp[i], ftp[i], outdir, md5[i], aspera,
                max_attempts=max_attempts, ftp_only=ftp_only
            )

            if is_r2:
                fastqs['r2'] = fastq
                fastqs['single_end'] = False
            else:
                fastqs['r1'] = fastq

    return fastqs


def md5sum(fastq):
    """Return the MD5SUM of an input file."""
    if os.path.exists(fastq):
        stdout = execute(f'md5sum {fastq}', capture_stdout=True)
        if stdout:
            md5sum, filename = stdout.split()
            return md5sum
        else:
            return None
    else:
        return None


def download_ena_fastq(fasp, ftp, outdir, md5, aspera, max_attempts=10, ftp_only=False):
    """Download FASTQs from ENA using Apera Connect or FTP."""
    success = False
    attempt = 0
    fastq = f'{outdir}/{os.path.basename(fasp)}'

    if not os.path.exists(fastq):
        execute(f'mkdir -p {outdir}')

        while not success:
            if ftp_only:
                logging.info(f'\t\tFTP download attempt {attempt + 1}')
                execute(f'wget --quiet -O {fastq} {ftp}', max_attempts=max_attempts)
            else:
                logging.info(f'\t\tAspera Connect download attempt {attempt + 1}')
                execute((f'{aspera["ascp"]} -QT -l {aspera["speed"]} -P33001 '
                         f'-i {aspera["private_key"]} era-fasp@{fasp} ./'),
                    directory=outdir,
                    max_attempts=max_attempts
                )

            fastq_md5 = md5sum(fastq)
            if fastq_md5 != md5:
                logging.log(STDOUT, f'MD5s, Observed: {fastq_md5}, Expected: {md5}')
                attempt += 1
                if os.path.exists(fastq):
                    os.remove(fastq)
                if attempt > max_attempts:
                    if not ftp_only:
                        ftp_only = True
                        attempt = 0
                    else:
                        logging.error(
                            f'Download failed after {max_attempts} attempts. '
                            'Please try again later or manually from SRA/ENA.'
                        )
                        sys.exit(1)
                time.sleep(10)
            else:
                success = True

    return fastq


def merge_runs(runs, output):
    """Merge runs from an experiment."""
    if len(runs) > 1:
        cat_cmd = ['cat']
        rm_cmd = ['rm']
        for run in runs:
            cat_cmd.append(run)
            rm_cmd.append(run)
        execute(' '.join(cat_cmd), stdout=output)
        execute(' '.join(rm_cmd))
    else:
        execute(f'mv {runs[0]} {output}')


def get_run_info(experiment):
    """Retreive a list of unprocessed samples avalible from ENA."""
    import requests
    url = f'{ENA_URL}&query="{query}"&fields={",".join(FIELDS)}'
    r = requests.get(url)
    if r.status_code == requests.codes.ok:
        data = []
        col_names = None
        for line in r.text.split('\n'):
            cols = line.rstrip().split('\t')
            if line:
                if col_names:
                    data.append(dict(zip(col_names, cols)))
                else:
                    col_names = cols
        return data
    else:
        return False


def write_json(data, output):
    """Write input data structure to a json file."""
    import json
    with open(output, 'w') as fh:
        json.dump(data, fh, indent=4, sort_keys=True)


def parse_query(query, is_study, is_experiment, is_run):
    "Parse user query, to determine search field value."
    if is_study:
        return f'study_accession={query}'
    elif is_experiment:
        return f'experiment_accession={query}'
    elif args.is_run:
        return f'run_accession={query}'
    else:
        # Try to guess...
        if query[1:3] == 'RR':
            return f'run_accession={query}'
        elif query[1:3] == 'RX':
            return f'experiment_accession={query}'
        else:
            return f'study_accession={query}'


def check_aspera(ascp, private_key, speed):
    "Verify Aspera Connect is available, not if it works."
    error_message = None
    if not os.path.exists(ascp):
        error_message = f'cannot access "{ascp}": No such file or directory'
    else:
        if private_key:
            # User provided path to private key
            if not os.path.exists(private_key):
                error_message = f'cannot access "{private_key}": No such file or directory'
        else:
            # Try to guess private key path, based on ascp path
            key_path = os.path.dirname(ascp).replace('/bin', '/etc')
            private_key = f'{key_path}/asperaweb_id_dsa.openssh'
            if not os.path.exists(private_key):
                error_message = f'cannot access "{private_key}": No such file or directory'

    if error_message:
        logging.error(f'Aspera Related Error: {error_message}')
        sys.exit(1)
    else:
        return {'ascp': ascp, 'private_key': private_key, 'speed': speed}


if __name__ == '__main__':
    import argparse as ap
    import sys

    parser = ap.ArgumentParser(
        prog=PROGRAM,
        conflict_handler='resolve',
        description=(f'{PROGRAM} (v{VERSION}) - Download FASTQs from ENA or SRA')
    )
    group1 = parser.add_argument_group('Required Options', '')
    group1.add_argument('query', metavar="ACCESSION", type=str,
                        help=('ENA/SRA accession to query. (Study, Experiment, or '
                              'Run accession)'))
    group1.add_argument(
        'provider', choices=['sra', 'SRA', 'ena', 'ENA'],
        help='Specify which provider (ENA or SRA) to use. Accepted Values: ENA SRA'
    )

    group2 = parser.add_argument_group('Aspera Connect Options')
    group2.add_argument('--aspera', metavar="STRING", type=str,
                        help='Path to the Aspera Connect tool "ascp" (Default: "which ascp")')
    group2.add_argument(
        '--aspera_key', metavar="STRING", type=str,
        help='Path to Aspera Connect private key, if not given, guess based on ascp path'
    )
    group2.add_argument('--aspera_speed', metavar="STRING", type=str, default="100M",
                        help='Speed at which Aspera Connect will download. (Default: 100M)')

    group3 = parser.add_argument_group('Query Related Options')
    group3.add_argument('--is_study', action='store_true', help='Query is a Study.')
    group3.add_argument('--is_experiment', action='store_true', help='Query is an Experiment.')
    group3.add_argument('--is_run', action='store_true', help='Query is a Run.')
    group3.add_argument('--group_by_experiment', action='store_true',
                        help='Group Runs by experiment accession.')
    group3.add_argument('--group_by_sample', action='store_true',
                        help='Group Runs by sample accession.')

    group4 = parser.add_argument_group('Helpful Options')
    group4.add_argument('--outdir', metavar="OUTPUT_DIR", type=str, default='./',
                        help=('Directory to output downloads to. (Default: ./)'))
    group4.add_argument('--prefix', metavar="PREFIX", type=str, default='fastq',
                        help=('Prefix to use for naming log files (Default: fastq)'))
    group4.add_argument('--max_attempts', metavar="INT", type=int, default=10,
                        help='Maximum number of download attempts (Default: 10)')
    group4.add_argument('--cpus', metavar="INT", type=int, default=1,
                        help='Total cpus used for downloading from SRA (Default: 1)')
    group4.add_argument('--ftp_only', action='store_true', help='FTP only downloads.')
    group4.add_argument(
        '--sra_only', action='store_true', 
        help='Do not attempt to fall back on ENA if SRA download does not work (e.g. missing FASTQ).'
    )
    group4.add_argument('--silent', action='store_true',
                        help='Only critical errors will be printed.')
    group4.add_argument('--verbose', action='store_true',
                        help='Print debug related text.')
    group4.add_argument('--debug', action='store_true',
                        help='Skip downloads, print what will be downloaded.')
    group4.add_argument('--version', action='version', version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    # Setup logs
    FORMAT = '%(asctime)s:%(name)s:%(levelname)s - %(message)s'
    logging.basicConfig(format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S',)
    logging.getLogger().setLevel(set_log_level(args.silent, args.verbose))

    aspera = check_aspera(args.aspera, args.aspera_key, args.aspera_speed) if args.aspera else None
    if not aspera:
        if args.provider == 'ENA':
            logging.info("Aspera Connect not available, using FTP for ENA downloads")
        args.ftp_only = True

    outdir = os.getcwd() if args.outdir == './' else f'{args.outdir}'
    query = parse_query(args.query, args.is_study, args.is_experiment, args.is_run)

    # Start Download Process
    ena_data = get_run_info(query)
    logging.info(f'Query: {args.query}')
    logging.info(f'Archive: {args.provider}')
    logging.info(f'Total Runs To Download: {len(ena_data)}')
    runs = {} if args.group_by_experiment or args.group_by_sample else None
    for i, run in enumerate(ena_data):
        logging.info(f'\tWorking on run {run["run_accession"]}...')
        fastqs = None
        if args.provider.lower() == 'ena':
            fastqs = ena_download(run, outdir, aspera=aspera,
                                   max_attempts=args.max_attempts,
                                   ftp_only=args.ftp_only)
        else:
            fastqs = sra_download(run["run_accession"], outdir, cpus=args.cpus,
                                  max_attempts=args.max_attempts)
            if fastqs == "SRA_NOT_FOUND":
                if args.sra_only:
                    logging.error(f'\t{run["run_accession"]} not found on SRA')
                    ena_data[i]["error"] = 'SRA_NOT_FOUND'
                    fastqs = None
                else:
                    # Retry download from ENA
                    logging.info(f'\t{run["run_accession"]} not found on SRA, retrying from ENA')
                    fastqs = ena_download(run, outdir, aspera=aspera,
                                          max_attempts=args.max_attempts,
                                          ftp_only=args.ftp_only)


        # Add the download results
        if fastqs:
            if args.group_by_experiment or args.group_by_sample:
                name = run["sample_accession"]
                if args.group_by_experiment:
                    name = run["experiment_accession"]

                if name not in runs:
                    runs[name] = {'r1': [], 'r2': []}

                if fastqs['single_end']:
                    runs[name]['r1'].append(fastqs['r1'])
                else:
                    runs[name]['r1'].append(fastqs['r1'])
                    runs[name]['r2'].append(fastqs['r2'])

    # If applicable, merge runs
    if runs and not args.debug:
        for name, vals in runs.items():
            if len(vals['r1']) and len(vals['r2']):
                # Not all runs labled as paired are actually paired.
                if len(vals['r1']) == len(vals['r2']):
                    logging.info(f'\tMerging paired end runs to {name}...')
                    merge_runs(vals['r1'], f'{outdir}/{name}_R1.fastq.gz')
                    merge_runs(vals['r2'], f'{outdir}/{name}_R2.fastq.gz')
                else:
                    logging.info('\tMerging single end runs to experiment...')
                    merge_runs(vals['r1'], f'{outdir}/{name}.fastq.gz')
            else:
                logging.info('\tMerging single end runs to experiment...')
                merge_runs(vals['r1'], f'{outdir}/{name}.fastq.gz')
        write_json(runs, f'{outdir}/{args.prefix}-run-mergers.json')
    write_json(ena_data, f'{outdir}/{args.prefix}-run-info.json')
