# ###################################
# #                                 #
# #    SETI_CATALOG_FUNCTIONS.PY    #
# #                                 #
# ###################################

# Functions to support Optical SETI searches using HARPS data.
# HARPS: https://www.eso.org/sci/facilities/lasilla/instruments/harps.html
#
# These functions operate on entire catalogs of stars (CSV text files).
#
# They assume you have an account at the ESO archive,
# and have configured astroquery for automatic password entry.
# See http://archive.eso.org/eso/eso_archive_main.html
# https://astroquery.readthedocs.io/en/latest/eso/eso.html#automatic-password

# Authors: Benjamin Fields, Jason Goodman
# Repo: https://github.com/spacetraveller42/optical_seti
#

eso_login = "Spaceboy42"
star_list = "OSETI_targets.txt"
results_list = "onelineperhit.txt"
predownload_list = "multistarwidefield.txt"

from astroquery.eso import Eso
from astropy.io import fits
import astropy.config
import optical_seti_functions
from pathlib import Path
eso_cache_path = Path(astropy.config.get_cache_dir()) / "astroquery" / "Eso"

# parse_star_list: Extract important data from input file
#
# Inputs:
#   withstarlist: filename of input list of stars to analyze.
#   predownloader_format: If true, assumes the star list is in
#               the format created by predownloader().  If false,
#               assumes the star list is in the OSETI_targets.txt
#               format.
# Outputs:
#   target_list: List of stars to analyze
#   spectral_types, temperatures, distances: star data
#   harps_objects: HARPS's name for the stars
#   data_files: names of predownloaded files (blank if predownloader_format=False)

def parse_star_list(withstarlist=star_list,predownloader_format=False):
    if (predownloader_format):
        target_list = [x.split(',')[0] for x in open(withstarlist).readlines()]
        spectral_types = [x.split(',')[1] for x in open(withstarlist).readlines()]
        data_files = [x.split(',')[2] for x in open(withstarlist).readlines()] 
        harps_objects = [x.split(',')[3] for x in open(withstarlist).readlines()]
        temperatures = [x.split(',')[4] for x in open(withstarlist).readlines()]
        distances = [x.split(',')[5].strip() for x in open(withstarlist).readlines()]
    else: # OSETI_targets.txt format
        target_list = [x.split('\t')[0] for x in open(withstarlist).readlines()]
        spectral_types = [x.split('\t')[3] for x in open(withstarlist).readlines()] 
        data_files = [""] * len(target_list)  # empty
        harps_objects = target_list
        temperatures = [x.split('\t')[4] for x in open(withstarlist).readlines()] 
        distances = [x.split('\t')[8] for x in open(withstarlist).readlines()]
    return target_list,spectral_types,data_files,harps_objects,temperatures,distances

# download_one_obs: Download one observation file from HARPS
# 
# Inputs:
#   star: Queryable name of star
#   eso: astroquery.eso object (eso.login must already be completed)
# Outputs:
#   file: name of downloaded file
#   harps_object: HARPS's name for the object

def download_one_obs(star,eso):
    tbl = eso.query_surveys('HARPS', target= star,box=0.1)
    harps_object = tbl['Object'][1:2][0]
    arcfile = tbl['ARCFILE'][1:2][0]
    cached_file = eso_cache_path /  (arcfile.replace(":","_")+".fits")
    if (cached_file.exists()):
        print("using cached file " + str(cached_file))
        file = cached_file
    else:
        print("downloading " + str(arcfile))
        file = eso.retrieve_data(arcfile)
    return file,harps_object

# do_search: Search for spikes in all stars in a catalog.
#
# Inputs:
#   withlogin: ESO archive login username.
#   withstarlist: filename of input list of stars to analyze.
#   withresults: filename of output CSV list of identified spikes.
#   categorize: If true, writes separate files categorizing
#               the files by number of spikes found
#               (no spikes, <4 spikes, 4+ spikes)
#   predownloader_format: If true, assumes the star list is in
#               the format created by predownloader().  If false,
#               assumes the star list is in the OSETI_targets.txt
#               format.

def do_search(withlogin=eso_login,withstarlist=star_list,withresults=results_list,categorize=False,predownloader_format=False):
    eso = Eso()
    eso.login(username=withlogin)

    # Get key data from the input catalog
    target_list,spectral_types,data_files,harps_objects,temperatures,distances = parse_star_list(withstarlist,predownloader_format)

    # Create output file(s)
    header_string = "STAR,SPECTRAL TYPE,START INDEX,END INDEX,OBSERVATION FILE,START WAVELENGTH,END WAVELENGTH,HARPS OBJECT,TEMPERATURE,DISTANCE\n"
    if (categorize):  # Create two files, one for stars with few spikes (possible lasers), one for stars with many (possible stellar emission lines)
        output_fewspikes = open("results_fewspikes.txt","w")
        output_fewspikes.write(header_string) 
        output_manyspikes = open("results_manyspikes.txt","w")
        output_manyspikes.write(header_string) 
    else:
        output = open(withresults, "w")
        output.write(header_string)

    # Loop over all stars in the targets list
    for (star, harps_object, spectral_type, temperature, file, distance) in zip(target_list[1:], harps_objects[1:], spectral_types[1:], temperatures[1:], data_files[1:], distances[1:]):
        print(star)
        if not file: # If our star list doesn't have a predownloaded filename, query ESO archive to get one.
            file,harps_object = download_one_obs(star,eso)
        file = Path(file)
        print(file.name)
        # Read in file data
        [wave,arr1] = optical_seti_functions.read_harps_file(file)
        if len(arr1)==0:
            continue
        # Perform SETI spike analysis
        hits_start, hits_end, count  = optical_seti_functions.seti_spike_analyzer(arr1, min_count = 4, max_count = 60, threshold_multiplier = 3.5, window_size = 101)
        # For each hit we found...
        for (start, end) in zip(hits_start, hits_end):
            wavelength_start = wave[start]
            wavelength_end = wave[end]
            if (categorize): # Pick which file to write out to
                if len(hits_start) <= 3:
                    output_dest = output_fewspikes
                else:
                    output_dest = output_manyspikes
            else:
                output_dest = output
            output_dest.write("{},{},{},{},{},{},{},{},{},{}\n".format(star, spectral_type, start, end, file.name, wavelength_start, wavelength_end, harps_object, temperature,distance))
    if (categorize):
        output_fewspikes.close()
        output_manyspikes.close()
    else:
        output.close()

# predownloader: Predownload one spectrum from each of the stars in the catalog.
#
# Inputs:
#   withlogin: ESO archive login username.
#   withstarlist: filename of input list of stars to analyze.
#   withresults: filename of output CSV list of filenames downloaded.

def predownload(withlogin=eso_login,withstarlist=star_list,withresults=predownload_list):
    eso = Eso()
    eso.login(username=withlogin)
    
    # Get key data from the input catalog
    target_list,spectral_types,data_files,harps_objects,temperatures,distances = parse_star_list(withstarlist,predownloader_format=False)

    # Create output file
    output = open(withresults, "w")
    output.write("STAR,SPECTRAL TYPE,OBSERVATION FILE,HARPS OBJECT,TEMPERATURE,DISTANCES\n") 
    for (star, spectral_type, temperature, distance) in zip(target_list[1:], spectral_types[1:], temperatures[1:], distances[1:]):
        print(star)
        file,harps_object = download_one_obs(star,eso)
        output.write("{},{},{},{},{},{}\n".format(star,spectral_type,Path(file).name,harps_object,temperature, distance))
    output.close()

# bulk_predownloader: Predownload one or more spectra from each of the stars in the catalog.
#   Group up files to download to reduce the number of requests we make to ESO servers.
#
# Inputs:
#   withlogin: ESO archive login username.
#   withstarlist: filename of input list of stars to analyze.
#   withresults: filename of output CSV list of filenames downloaded.
#   obs_per_star: number of obs to retrieve for each star.

def bulk_predownload(withlogin=eso_login,withstarlist=star_list,withresults=predownload_list,obs_per_star=1):
    eso = Eso()
    eso.login(username=withlogin)

    download_queue = []         # list of files to download
    download_queue_size = 40    # maximum number of files to queue up before starting download

    # Get key data from the input catalog
    target_list,spectral_types,data_files,harps_objects,temperatures,distances = parse_star_list(withstarlist,predownloader_format=False)

    # Create output file
    output = open(withresults, "w")
    output.write("STAR,SPECTRAL TYPE,OBSERVATION FILE,HARPS OBJECT,TEMPERATURE,DISTANCES\n") 
    # Loop over all stars in the input catalog
    for (star, spectral_type, temperature, distance) in zip(target_list[1:], spectral_types[1:], temperatures[1:], distances[1:]):
        print(star)
        tbl = eso.query_surveys('HARPS', target= star,box=0.1)   # Find all obs for that star
        if not (tbl):                                            # Go to next star if here are no obs
            continue
        if len(tbl)<2:                                           # Go to next star if there's only one obs
            continue
        harps_object = tbl['Object'][1:2][0]
        arcfiles = tbl['ARCFILE'][1:(obs_per_star+1)]            # Grab the first obs_per_star observations
        for i in range(len(arcfiles)):
            cached_file = eso_cache_path /  (arcfiles[i].replace(":","_")+".fits")  # figure out name of local file
            if not (cached_file.exists()):                       # If we don't have the file in cache
                download_queue.append(arcfiles[i])               # Add it to the download queue and record it in output file
            output.write("{},{},{},{},{},{}\n".format(star,spectral_type,cached_file.name,harps_object,temperature, distance))
        if (len(download_queue) > download_queue_size):          # If the download queue is full
            try:
                eso.retrieve_data(download_queue)
            except:                                             # Maybe login expired, relogin and try again
                eso.login(username=withlogin)
                eso.retrieve_data(download_queue)
            download_queue = []                                 # Clear download queue
    if (len(download_queue) > 0):
        eso.retrieve_data(download_queue)
    output.close()