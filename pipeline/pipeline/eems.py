import os
import numpy as np
from load import unwrap_america
import pandas as pd
from utils.plink import run_plink
from random import randint


TMP_PLINK = 'tmp_plink'


def run_all(args):
    if not args.dry:
        run_eems(args.eems_snps, ini_file=args.proj, n_demes=args.nDemes)
        create_plots(plot_scripts_folder=args.eems + '/runeems_snps/plot/',
                     out_path=args.output_folder, in_path=args.proj,
                     tmp_script=args.tmp_folder + '/plot.r',
                     wrap=args.wrap)


def create_plots(plot_scripts_folder, out_path, in_path, tmp_script, wrap):
    plot_script = "%s/sed-plots.R" % plot_scripts_folder
    tpl = out_path, in_path, plot_scripts_folder, plot_script, tmp_script
    s = 'sed -e "s!MCMCPATH!%s!g; s!PLOTPATH!%s!g; s!EEMSWD!%s!g;" %s > %s' \
        % (tpl)
    os.system(s)

    if wrap:
        coords_file = '%s.coord' % in_path
    data = np.loadtxt(coords_file)
    data[:, 1] = unwrap_america(data[:, 1])
    np.savetxt(coords_file, data, fmt="%2.2f")
    os.system('Rscript %s' % tmp_script)


def run_eems(eems_exe, ini_file, n_demes, dry=False):
    s = ""
    for nd in n_demes:
        tpl = eems_exe, ini_file, nd
        s += "%s --params %s_%s.ini &\n" % tpl
    s += 'wait'
    os.system(s)


def write_all_files(args, meta_data, polygon):
    if args.diffs is None:
        create_eems_files(args, meta_data=meta_data, polygon=polygon,
                          bed2diffs=args.bed2diffs,
                          bedfile=TMP_PLINK,
                          eems_input_name=args.proj,
                          eems_output_name=args.output_folder,
                          n_runs=args.n_runs,
                          **args.eems_args)
    else:
        create_eems_files(args, meta_data=meta_data, polygon=polygon,
                          bedfile=TMP_PLINK,
                          bed2diffs=args.bed2diffs,
                          diffs_file=args.diffs,
                          eems_input_name=args.proj,
                          eems_output_name=args.output_folder,
                          n_runs=args.n_runs,
                          **args.eems_args)


def filter_data(meta_data, bedfile, plink="plink"):
    """filter_data
    filters bed file to only keep the individuals in meta_data, uses plink
    
    the data is read from bedfile, and written to the file in the TMP_PLINK
    global variable
    
    Parameters
    ----------
    meta_data : pd.DataFrame
        pandas data frame with individuals to keep
    bedfile : path
        the bedfile to be filtered
    plink : path
        the plink executable to be used
    
    """
    include_name = '%s.incl' % TMP_PLINK

    fam = pd.read_table("%s.fam" % bedfile, header=None,
                        skipinitialspace=True, sep=" ")
    fam.columns = ['FAM', 'IND', 'a', 'b', 'c', 'd']

    extract_data = meta_data.merge(fam, on='IND', how='inner')
    extract_data.to_csv(include_name, sep=' ',
                        columns=('FAM', 'IND'),
                        header=None, index=None)

    meta_data = extract_data[['IND', 'POP', 'LAT', 'LONG', 'FAM']]
    meta_data = meta_data[np.logical_not(np.isnan(meta_data['LAT']))]
    meta_data = meta_data[np.logical_not(np.isnan(meta_data['LONG']))]

    flags = dict()
    flags['make-bed'] = ''
    flags['allow-extra-chr'] = ''
    flags['bfile'] = bedfile
    flags['out'] = TMP_PLINK
    flags['keep'] = include_name
    flags['indiv-sort'] = 'f %s' % include_name

    run_plink(plink, flags)

    return meta_data


def make_grid(grid, gridfile, outname):
    scriptpath = '/data/eems-project/eems/pipeline/pipeline/maps/extract_points.r'
    grid2 = '/data/eems-project/eems/pipeline/pipeline/maps/grid%s' % grid
    outer = '%s.outer' % outname
    sample = '%s.coord' % outname
    tpl = (scriptpath, grid2, outer, sample, gridfile)
    s = "Rscript %s %s %s %s %s" % tpl
    print s
    os.system(s)


def create_eems_files(args, meta_data=None, polygon=None,
                      bed2diffs="bed2diffs",
                      bedfile=None, diffs_file=None,
                      eems_input_name="eems_input",
                      eems_output_name="eems_output",
                      n_runs=1,
                      **kwargs):
    """create_eems_files

    creates the eems files
        - diff file with pairwise differences (using bed2diffs)
        - polygon file with the location to run run on
        - sample file with sampling locations
        - ini file with all remaining parameters

    if there is more than one argument passed to nDemes, a different ini file
    is created for each of them.
    
    Parameters
    ----------
    bedfile : path
        Path to the bed files to be converted
    meta_data : pd.DataFrame
        the location data for all individuals
    polygon : list of (str, str)
        the polygon to be written
    bed2diffs : path
        the bed2diffs executable to be used
    
    """
    if diffs_file is not None:
        filter_diffs_file(diffs_file, meta_data.IND, eems_input_name)
    elif bedfile is not None:
        create_diffs_file(bedfile=bedfile, bed2diffs=bed2diffs,
                          outname=eems_input_name)
    else:
        raise ValueError("either bedfile or diffs_file needs to be specified")

    create_polygon_file(polygon, eems_input_name)
    if args.sd is None:
        create_sample_file(meta_data, eems_input_name, 
                           order_file=eems_input_name)
    else:
        print meta_data.columns
        create_sample_file_jitter(meta_data, eems_input_name,
                                  order_file=eems_input_name,
                                  n_runs=int(args.n_runs))

    if args.grid != 0 and args.sd is None:
        kwargs['gridpath'] = "%s_%s" % (eems_input_name, args.grid)
        print "GRID AT", kwargs['gridpath']
        make_grid(args.grid, kwargs['gridpath'], eems_input_name)

    if args.submit_script:
        ss_name = "%s-submit.sh" % args.proj
        ss = open(ss_name, 'w')
        ss.write(submit_script_header())
    if args.run_script:
        rs_name = "%s-run.sh" % args.proj
        rs = open(rs_name, 'w')

    nSites = sum(1 for line in open("%s.bim" % bedfile))
    nDemes = kwargs['nDemes']
    del kwargs['nDemes']
    for nd in nDemes:
        for i in xrange(int(n_runs)):
            if args.sd is not None:
                ini_name = "%s_%s_run%s" % (eems_input_name, nd, i)
                out_name = "%s/%srun%s" % (eems_output_name, nd, i)
                datapath = "%s%d" % (eems_input_name, i)
                os.link("%s.diffs" % eems_input_name, "%s.diffs" % datapath)
                os.link("%s.outer" % eems_input_name, "%s.outer" % datapath)
                os.link("%s.order" % eems_input_name, "%s.order" % datapath)


                if args.grid != 0 :
                    kwargs['gridpath'] = "%s_%s" % (datapath, args.grid)
                    print "GRID AT", kwargs['gridpath']
                    make_grid(args.grid, kwargs['gridpath'], datapath)

            elif n_runs == 1:
                ini_name = "%s_%s" % (eems_input_name, nd)
                out_name = "%s/%s" % (eems_output_name, nd)
                datapath = eems_input_name
            else:
                ini_name = "%s_%s_run%s" % (eems_input_name, nd, i)
                out_name = "%s/%s_run%s" % (eems_output_name, nd, i)
                datapath = eems_input_name

            create_ini_file(ini_name, datapath=datapath,
                            mcmcpath=out_name,
                            meta_data=meta_data,
                            nSites=nSites, n_demes=nd,
                            **kwargs)

            if args.run_script:
                rs.write("%s --params %s.ini &\n" % (args.eems_snps, ini_name))
            if args.submit_script:
                ss.write('submit %s.ini\n' % ini_name)

    if args.run_script:
        rs.close()
    if args.submit_script:
        ss.close()
    kwargs['nDemes'] = nDemes


def submit_script_header():
    s = """
#!/bin/bash
EEMS="$HOME/eems/runeems_snps/src/runeems_snps"
mkdir -p logs/
submit()
{
        inifile=`basename $1`
        extension="${inifile##*.}"
        fileid="${inifile%.*}"

        echo "#!/bin/bash" > submit.sh
        echo "#$ -N $fileid" >> submit.sh
        echo "#$ -l h_vmem=4g" >> submit.sh
        echo "#$ -l s_rt=168:00:00" >> submit.sh
        echo "#$ -l h_rt=168:00:00" >> submit.sh
        echo "#$ -cwd" >> submit.sh
        echo "#$ -V" >> submit.sh
        echo "#$ -e logs/\$JOB_NAME_\$JOB_ID.err" >> submit.sh
        echo "#$ -o logs/\$JOB_NAME_\$JOB_ID.out" >> submit.sh
        echo "$EEMS --params $1 --seed \$JOB_ID" >> submit.sh

        qsub submit.sh
}

"""

    return s



def create_diffs_file(bedfile, bed2diffs, outname, nthreads=4):
    """create a file with pairwise differences using the bed2diffs executable
        
    Parameters
    ----------
    bedfile : path
        the bed file to run bed2diffs on
    bed2diffs : path
        path to the bed2diffs executable
    outname : str
        output file name
    
    """
    tpl = bed2diffs, nthreads, bedfile
    s = "%s --nthreads %d --bfile %s" % (tpl)
    s += " && mv %s.order %s.order " % (bedfile, outname)
    s += " && mv %s.diffs %s.diffs " % (bedfile, outname)
    os.system(s)


def create_ini_file(ini_name, mcmcpath, datapath, meta_data,
                    n_demes, run_id=1,
                    **kwargs):
    """create_ini_file
    
    Parameters
    ----------
    ini_name : path
        the name of the ini file to be created
    mcmcpath : path
        path for the mcmc output of eems
    datapath : path
        path where the input data is loaded from
    meta_data : pd.DataFrame
        data structure with sample meta data
    n_demes : int
        deme parameter, will be added to folder name
    
    """
    kwargs['nSites'] = kwargs.get('nSites', 10000)
    kwargs['nIndiv'] = len(meta_data)
    kwargs['diploid'] = kwargs.get('diploid', 'True')
    kwargs['numMCMCIter'] = kwargs.get('numMCMCIter', 20000)
    kwargs['numBurnIter'] = kwargs.get('numBurnIter', 10000)
    kwargs['numThinIter'] = kwargs.get('numThinIter', 99)
    kwargs['datapath'] = datapath
    kwargs['mcmcpath'] = mcmcpath
    kwargs['seed'] = randint(0, 10000000)

    with open("%s.ini" % ini_name, 'w') as f:
        for k, v in kwargs.iteritems():
            f.write("%s = %s\n" % (k, v))
        f.write("%s = %s\n" % ('nDemes', n_demes))


def filter_diffs_file(diffs_file, individuals, outname):
    """Filters an existing diffs an order file combo to only retain the
    individuals in `individuals`
    
    Parameters
    ----------
    diffs_file : path
        path to a diffs/order file (created by bed2diffs)
    individuals : list of strings
        the individuals to be retained
    outname : str
        output file name
    
    """
    diff_matrix = np.loadtxt("%s.diffs" % diffs_file)
    sample_order = pd.read_table("%s.order" % diffs_file, header=None, sep=" ")
    individuals = pd.DataFrame(individuals)
    individuals.columns = [0]

    print individuals[0]

    for i in individuals[0]:
        print i
        print np.where(i == sample_order[1])[0][0]

    to_keep = [np.where(i == sample_order[1])[0][0] for i in individuals[0]]
    new_diff_matrix = diff_matrix[to_keep][:, to_keep]
    individuals = individuals.merge(sample_order, how="left",
                                    left_on=0, right_on=1)

    np.savetxt("%s.diffs" % outname, new_diff_matrix, fmt="%f")
    individuals.to_csv("%s.order" % outname, header=False, index=False,
                       columns=('0_y', '0_x'), sep=" ")


def create_sample_file(meta_data, outname, order_file=None):
    """create_sample_file
    
    Parameters
    ----------
    meta_data : pd.DataFrame
        data frame of the meta data
    outname : path
        the file of the outputed sample file (without extension)
    order_file : path
        a file with the sample ordering (without extension .order).
        This file is created by bed2diffs. If present, the sample
        coordinates are aligned with this file
    """
    out = "%s.coord" % outname
    if order_file is not None:
        sample_order = pd.read_table("%s.order" % order_file, header=None,
                                     sep=" ")
        meta_data = sample_order.merge(meta_data, how="left", left_on=1,
                                       right_on='IND')

    meta_data.to_csv(out, sep=" ", header=False, index=False, float_format="%2.2f",
                     columns=('LONG', 'LAT'))


def create_sample_file_jitter(meta_data, outname, order_file=None,
                              n_runs=1):
    """create_sample_file_jitter

    creates a set of sample files with some normal jitter added
    
    Parameters
    ----------
    meta_data : pd.DataFrame
        data frame of the meta data
    outname : path
        the file of the outputed sample file (without extension)
    order_file : path
        a file with the sample ordering (without extension .order).
        This file is created by bed2diffs. If present, the sample
        coordinates are aligned with this file
    n_runs : int
        the number of input files created
    """
    if order_file is not None:
        sample_order = pd.read_table("%s.order" % order_file, header=None,
                                     sep=" ")
        meta_data = sample_order.merge(meta_data, how="left", left_on=1,
                                       right_on='IND')


    out = "%s.sample" % (outname)
    meta_data.to_csv(out, sep=" ", header=True, index=False,
                     columns=('IND', 'POP', 'LONG', 'LAT', 'SD'))
    for i in range(n_runs):
        out = "%s%d.coord" % (outname, i)
        a = np.random.normal(meta_data['LONG'], meta_data['SD'])
        b = np.random.normal(meta_data['LAT'], meta_data['SD'])
        a = ["%2.2f" % i for i in a]
        b = ["%2.2f" % i for i in b]
        temp_data = pd.DataFrame(zip(a, b), columns=('LONG', 'LAT'))

        temp_data.to_csv(out, sep=" ", header=False, index=False,
                         columns=('LONG', 'LAT'))


def create_polygon_file(polygon, outname):
    """create_polygon_file
    writes a polygon to outname, for input in eems
    
    Parameters
    ----------
    polygon : shapely Polygon
        the polygon to be written
        
    outname : str
        file name of the polygon
    
    Returns
    -------
    """
    points = list(polygon.exterior.coords)
    np.savetxt("%s.outer" % outname, points, fmt="%f")
