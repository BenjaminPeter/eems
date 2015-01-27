from parameters import Parameters
from load import load_location_file, load_sample_file
from polygon import get_eems_area
from eems import write_all_files, run_all, filter_data


def run(params):
    location_data = load_location_file(params.loc, params.location_header)
    sample_data = load_sample_file(params.sample, params.sample_header)

    meta_data = sample_data.merge(location_data)
    print "loaded meta data"
    polygon, meta_data = get_eems_area(params, meta_data)
    print "got area"
    filter_data(meta_data, params.bed)
    print "filtered data"
    write_all_files(params, meta_data, polygon)
    run_all(params)


def run_cli():
    """run_cli

    runs the pipeline from the command line interface.
    see -h flag for arguments
    """
    params = Parameters.from_command_line()
    run(params)
