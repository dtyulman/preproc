#!/usr/bin/env python

import argparse
import ConfigParser
from datetime import timedelta
import dicom  
import glob
import multiprocessing as mp
import nipype.interfaces.io as nio
import nipype.pipeline.engine as pe
import nipype.interfaces.fsl as fsl
import nipype.interfaces.spm as spm
import nipype.interfaces.dcm2nii as d2n
import nipype.interfaces.utility as util
import os
import shutil
import sys
from time import time, ctime
## other imports done inside functions due to NiPype engine constraints:
#from ast import literal_eval
#import math
#import nibabel
#from neuro.apps.xnat import Xnat

## --- Uncomment for debug mode --- #
#from nipype import config, logging
#config.enable_debug_mode()
#logging.update_logging(config)

def fatal(message, exit_status=1, *args):
    print "[FATAL ERROR]:", message, args
    print "[--EXITING--]"
    sys.exit(exit_status)

def info(message, *args):
    if args != ():
        print "[INFO]:", message, args
    else:
        print "[INFO]:", message

def normpath(path):
    return os.path.abspath(os.path.expandvars(os.path.expandvars(path)))
    
def get_options():
    #parse command line inputs   
    parser = argparse.ArgumentParser(description="Enhanced reimplementation of the procfast preprocessing pipeline (http://cnlwiki.pbworks.com/w/page/13165336/Data%20Preprocessing:%20procfast) using NiPype. See tokeefe.neuroinfo.org/wiki/Preprocessing for detailed usage information.")    
    dicoms = parser.add_mutually_exclusive_group(required=True)
    dicoms.add_argument('-xnat', '--xnat_alias', help="Alias of XNAT server used in .xnat_auth") 
    dicoms.add_argument('-d', '--dicom_dir', help="Directory on your hard disk where the DICOM files are located.")
#    dicoms.add_argument('-url', '--xnat_url', help="URL of the XNAT server from which to download the data.") 
##    parser.add_argument('-u', '--username', help="Username to log into XNAT. (required if using xnat_url)") 
##    parser.add_argument('-p', '--password', help="Password to log into XNAT. (required if using xnat_url)")   
    parser.add_argument('-id', required=True, help="Subject ID")
    parser.add_argument('-bold', nargs='+', type=int, required=True, help="Specify functional scan numbers (BOLD) separated by spaces.")
    parser.add_argument('-anat', nargs=1, type=int, required=True, help="Specify the anatomical scan number.")     
    parser.add_argument('-o', '-out', '--output_dir', default=None, help="Directory to place processed files. (Default is './<subject_id>')")
    parser.add_argument('-w', '-work', '--working_dir', help="Directory to place all intermediate files (Default is output_dir)")
    parser.add_argument('-T1', '-T1target', '--anat_template', help="Specify T1 target image for normalizing ANAT (MPRAGE). (Default is MNIT1S)")
    parser.add_argument('-T2', '-T2target', '--bold_template', help="Specify T2 target image for normalizing BOLD. (Default is SPM8EPI)")
    parser.add_argument('-skip', type=int, help="Specify number of frames to skip. (Default is 4)")
    parser.add_argument('-tr', type=int, help="Specify sampling rate in seconds. (Will try to read from file if not specified)")
    parser.add_argument('-useconfig', '-cfg', help="Specify a custom configuration file to use.")
    parser.add_argument('-smooth', type=float, help="Specify (in mm) the FWHM (Full Width Half Max) of a Gaussian kernel to smooth preprocessed data.")
    parser.add_argument('-clean', action='store_true', help="Remove working directory after completing preprocessing.")
    parser.add_argument('-slicetime', '-faln', choices=['Y','N'], help="Choose whether to perform slice time correction (a.k.a. frame alignment). (Default is yes)")
    parser.add_argument('-spm', '--spmversion', type=int, choices=[2, 8], help="Choose which version of SPM to use with slice time correction and normalization. (Default is SPM8)")    
    parser.add_argument('-g', '--writegraph', choices=['orig','flat', 'exec', 'hierarchical'], help="Create a .dot file and a .png image of the workflow graph. Requires dot. See NiPype website for more info. (Default is hierarchical)")    
    parser.add_argument('-l', '--legacy', action='store_true', help="Use the calculation for the .rdat file from the original fsl_preprocess.sh (Default is False)")
    args = parser.parse_args()
     
    #initialize options with defaults config
    config = ConfigParser.SafeConfigParser()
    config.optionxform=str
    options = {}
    with open(os.path.join(os.path.dirname(__file__), '.preproc.ini')) as fp:
        config.readfp(fp)
        for sect in config.sections():
            for name, value in config.items(sect):
                options[name]=value
    
    #update options from custom config            
    if args.useconfig != None: #config file was specified by user
        info("Using custom config file")
        try: 
            with open(os.path.expanduser(args.useconfig)) as fp:
                config.readfp(fp)
        except: fatal('Could not find specified configuration file:', args.useconfig)
        for sect in config.sections():
            for name, value in config.items(sect):
                if name not in options.keys(): #check for errors in user config file
                    fatal('Error in custom configuration file: %d is not a valid entry' %name)
                options[name]=value
    
    #update options from command line args         
    for name, value in vars(args).iteritems():
        if name not in options.keys() or value != None:
            options[name] = value
        
    #clean up options dictionary (typecasting, normalize pathnames, etc)
    options['bold'].sort()
    options['code_dir'] = os.path.dirname(__file__) #for internal use
    for name, value in options.iteritems(): 
        if type(value) == str:           
            if value.lower() == 'none' or value.lower() == 'null':
                options[name] = None
            elif value.lower() == 'true':
                options[name] = True
            elif value.lower() == 'false':
                options[name] = False
            else:
                try: options[name] = int(value)
                except: 
                    try: options[name] = float(value)
                    except: pass                          
    if options['output_dir'] == None:
        options['output_dir'] = './%s' %options['id']                               
    if options['working_dir'] == None:
        options['working_dir'] = options['output_dir']
    paths_list = ['output_dir', 'working_dir', 'spm8_path','spm2_path','spm_params']
    for path in paths_list:
        options[path] = normpath(options[path]) 
    
    #set T1 and T2 targets for normalization (choose which template to use)   
    anat_template = options['anat_template'].strip().upper()
    bold_template = options['bold_template'].strip().upper()
    if anat_template in options:
        options['anat_template'] = normpath(options[anat_template])
    if bold_template in options:
        options['bold_template'] = normpath(options[bold_template])
             
    #some basic user input error checking
    if options['spmversion'] == 2 and options['bold_template'].endswith('nii'):
        fatal("Cannot use template %s with SPM2. See online documentation for more information" %options['bold_template'])
    elif options['spmversion'] == 8 and options['bold_template'].endswith('mnc'):
        fatal("Cannot use template %s with SPM8. See online documentation for more information" %options['bold_template'])
    
    return options     

def downloadDicoms(options): #download DICOM files from XNAT and save to hard disk
    info("Opening connection to XNAT...")
    from neuro.apps.xnat import Xnat#, Xnat_15
    if options['xnat_alias']:
        xnat = Xnat.create(options['xnat_alias'])
#    elif options['xnat_url']:       
##        if options['username'] == None or options['password'] == None:
##            fatal("You must specify XNAT username and password if using XNAT URL") 
##        xnat = Xnat_15({"url": options['xnat_url'], "username": options['username'], "password": options['password']})
#        xnat = Xnat.create(options['xnat_url'])
    else:
        fatal("You must specify the location of the DICOM files. Use the -h or --help flag for more information.")
    arcget = xnat.createArcGet()
    arcget.addScanNumbers(options['anat']+options['bold'])
    arcget.setSessionID(options['id'])
    arcget.setOutputDirectory(os.path.join(options['working_dir'], options['id'], 'preproc', 'dicoms'))
    info("Downloading DICOMs from XNAT to: %s" %os.path.join(options['working_dir'], options['id'], 'preproc', 'dicoms'))
    options['dicom_dir'] = arcget.execute()
   
# --- Function to build a workflow for DICOM to NIFTI conversion        
def make_converter_workflow(options):
    # --- Node (dcm2nii): DICOM to NIFTI conversion    
    example_dicom_files = []
    options['max_scan_len'] = 0
    scan_files = options['dicoms_dict']
    for series_num,value in scan_files.iteritems():
        example_dicom_files.append(value['example_file'])
        if value['n_files'] > options['max_scan_len']:
            options['max_scan_len'] = value['n_files']
    
    dcm2nii = pe.Node(interface=d2n.Dcm2nii(), name='dcm2nii')
    dcm2nii.inputs.args = '-e y -f n -d n -p n -a n' #file output naming convention: 's{seriesnumber}a{acqnumber}.nii'
    dcm2nii.inputs.source_names = example_dicom_files
    dcm2nii.inputs.reorient_and_crop = False
    dcm2nii.inputs.reorient = False
    
    # --- Node: (orient): force the given orientation on the scans
    orient = pe.MapNode(interface=fsl.SwapDimensions(), name='orient', iterfield='in_file')
    orient.inputs.new_dims = ('RL', 'PA', 'IS')
    orient.inputs.output_type = 'NIFTI'
    
    # --- Node: (niftisort): given a list of NIFTI files, separate them into anatomical and BOLD scans. 
    niftisorter = pe.Node(interface=util.Function(input_names=['nifti_list', 'anat_number'], 
                                                  output_names=['anat', 'bold'],
                                                  function=nifti_sort), 
                          name='niftisorter')
    niftisorter.inputs.anat_number = options['anat'][0]
    
    # --- Node (convertsink): collect the relevant outputs from the converter workflow
    convertsink = pe.Node(interface=nio.DataSink(), name='convertsink')
    convertsink.inputs.base_directory = os.path.join(options['output_dir'], options['id'])
    convertsink.inputs.parameterization = False
    
    # --- Workflow (converter): build the workflow and make required node connections
    converter = pe.Workflow(name='converter')
    converter.base_dir = os.path.join(options['working_dir'], options['id'], 'preproc')
    
    
    converter.connect(dcm2nii, 'converted_files', orient, 'in_file')
    converter.connect([(orient, convertsink, [(('out_file', get_scan, options['anat'][0]), 'anat')])
                     ])
    for bold in options['bold']:
        converter.connect([(orient, convertsink, [(('out_file', get_scan, bold), 'bold.'+str(bold))])
                         ])
        
    return converter

# --- Converter workflow helper functions:
def getDicoms(dicom_dir): #get DICOM files from hard disk
    print dicom_dir
    scan_files = {};
    info("Walking directory, %s looking for DICOM files" %dicom_dir)
    for root, dirs, files in os.walk(dicom_dir): #walk directory tree, beginning at dicom_dir
        for f in files:
            full_file = os.path.join(root, f)
            try:
                d = dicom.read_file(full_file, stop_before_pixels=True)
            except:
                info("   -> not DICOM '" + full_file + "'")
                continue        
            if(d.SeriesNumber not in scan_files): #store the first file we find of all scans (by series number)
                scan_files[d.SeriesNumber] = {'n_files': 1, 'example_file': full_file}
            else:
                scan_files[d.SeriesNumber]['n_files'] += 1
    return scan_files

def get_scan(list, scan_number):
    import re #regular expressions
    pattern = re.compile(r's0*%da\d{3,}(_newdims)*\.nii(.gz)*$' %scan_number) 
    for nifti in list:
        if re.search(pattern, nifti):
            return nifti
        
def nifti_sort(nifti_list, anat_number):
    import re #regular expressions
    anat = []
    bold = []
    pattern = re.compile(r's0*%da\d{3,}(_newdims)*\.nii(.gz)*$' %anat_number) 
    for nifti in nifti_list:
        if re.search(pattern, nifti):
            anat.append(nifti)
        else:
            bold.append(nifti)
    if len(anat) != 1:
        import sys
        print "[FATAL ERROR]: Incorrect amount of anatomical scans:", anat  
        print "[--EXITING--]"
        sys.exit(1)      
    return anat[0], bold


# --- Function to build anatomical preprocessing workflow
def make_anat_workflow(scan_number,  options):        
    # --- Node (inputanat): used to input appropriate scans into other nodes
    inputanat = pe.Node(interface=util.IdentityInterface(fields=['anat']), 
                        name='inputanat') 
    scan = glob.glob(os.path.join(options['output_dir'], options['id'], 'anat', 's*'+str(scan_number)+'*'))    
    inputanat.inputs.anat = scan[0]
        
    # --- Node (bet): Anatomical brain extraction
    bet = pe.Node(interface=fsl.BET(), name='bet')
    bet.inputs.out_file = 'bhighres.nii'
    bet.inputs.output_type = 'NIFTI'
    bet.inputs.vertical_gradient = options['vertical_gradient']
    
    # --- Node (flirt): Normalize skull-stripped anatomical to target
    bhighres2target = pe.Node(interface=fsl.FLIRT(), name='bhighres2target')
    bhighres2target.inputs.reference = options['anat_template']
    bhighres2target.inputs.out_file = 'bhighres2target.nii'
    bhighres2target.inputs.output_type = 'NIFTI'
    bhighres2target.inputs.out_matrix_file = 'bhighres2target.mat'
    bhighres2target.inputs.cost = 'corratio'
    bhighres2target.inputs.dof = 12
    bhighres2target.inputs.searchr_x = [-90, 90]
    bhighres2target.inputs.searchr_y = [-90, 90]
    bhighres2target.inputs.searchr_z = [-90, 90]
    bhighres2target.inputs.interp = 'trilinear'
    
    # --- Node (convertxfm): invert the transformation matrix
    convert_xfm = pe.Node(interface=fsl.ConvertXFM(), name='convert_xfm')
    convert_xfm.inputs.invert_xfm = True
    convert_xfm.inputs.out_file = 'target2bhighres.mat' 
    convert_xfm.inputs.output_type = 'NIFTI'
    
    # --- Node (flirt): Normalize anatomical to target using affine matrix from skull-stripped normalization
    highres2target = pe.Node(interface=fsl.FLIRT(), name='highres2target')
    highres2target.inputs.reference = options['anat_template']
    highres2target.inputs.out_file = 'highres2target.nii'
    highres2target.inputs.output_type = 'NIFTI'
    highres2target.inputs.apply_xfm = True
    
    # --- Node (flirt)
    highres1112target = pe.Node(interface=fsl.FLIRT(), name='highres1112target')
    highres1112target.inputs.reference = options['anat_template']
    highres1112target.inputs.out_file = 'highres1112target.nii'
    highres1112target.inputs.output_type = 'NIFTI'
    highres1112target.inputs.force_scaling = True
    highres1112target.inputs.args = '-applyisoxfm 1'
    
    # --- Node (anatsink): collect the relevant outputs from the anatomical workflow
    anatsink = pe.Node(interface=nio.DataSink(), name='anatsink')
    anatsink.inputs.base_directory = os.path.join(options['output_dir'], options['id'])
    
    # --- Workflow (anatwork): build the workflow and make required node connections
    anatwork = pe.Workflow(name='anatwork')
    anatwork.base_dir = os.path.join(options['working_dir'], options['id'], 'preproc')
    anatwork.connect([(inputanat, bet, [('anat', 'in_file')]), #anatomical normalization nodes
                      (bet, bhighres2target, [('out_file', 'in_file')]),
                      (bhighres2target, convert_xfm, [('out_matrix_file', 'in_file')]),
                      (inputanat, highres2target, [('anat', 'in_file')]),
                      (bhighres2target, highres2target, [('out_matrix_file', 'in_matrix_file')]),
                      (inputanat, highres1112target, [('anat', 'in_file')]),
                      (bhighres2target, highres1112target, [('out_matrix_file', 'in_matrix_file')])
                    ])
    anatwork.connect([(bet, anatsink, [('out_file', 'qc.@bet')]), #datasink nodes
                      (bhighres2target, anatsink, [('out_file', 'anat.@bhighres2target')]),
                      (highres2target, anatsink, [('out_file', 'anat.@highres2target')]),
                      (highres1112target, anatsink, [('out_file', 'anat.@highres1112target')])
                    ]) 
    
    return anatwork

# --- Function to build the BOLD preprocessing workflow                    
def make_bold_workflow(scan_number, workflow_name, options):                    
    # --- Initialize the BOLD workflow    
    boldwork = pe.Workflow(name=workflow_name)
    boldwork.base_dir = os.path.join(options['working_dir'], options['id'], 'preproc')    
    
    # --- Node (boldsink): collect the relevant outputs from the BOLD workflow
    boldsink = pe.Node(interface=nio.DataSink(), name='boldsink')
    boldsink.inputs.parameterization = False
    boldsink.inputs.base_directory = os.path.join(options['output_dir'], options['id'], 'bold', str(scan_number))
    boldsink_folder = '' 
    
#-------------------------------------------------------------------------------#
    
    # --- Node (frameskip): skip the first (default=4) frames of each of the BOLD files
    frameskip = pe.Node(interface=fsl.ExtractROI(), name='frameskip')
    frameskip.inputs.t_min = options['skip']
    frameskip.inputs.t_size = options['max_scan_len']
    frameskip.inputs.output_type = 'NIFTI'        
    scan = glob.glob(os.path.join(options['output_dir'], options['id'], 'bold', str(scan_number), 's*'+str(scan_number)+'*'))[0]      
    frameskip.inputs.in_file = scan      
    
    boldsink_folder = boldsink_folder + "skip"    
    boldwork.connect(frameskip, 'roi_file', boldsink, boldsink_folder)
    
 
#----------------------------SLICE-TIME CORRECTION------------------------------#  

    # --- If-else: choose whether to build a placeholder node or a set of nodes to perform slicetime correction. 
    if options['slicetime'] == 'N':
        # --- Node: fixes the input/output name to match that of slice time correction so that
        #           the workflow connections remain although slice time correction is not performed        
        slicetime = pe.Node(interface=util.Function(input_names=['input'],
                                                    output_names=['timecorrected_files'],
                                                    function=slicetime_placeholder), name='slicetime_placeholder')      
        boldwork.connect(frameskip, 'roi_file', slicetime, 'input')                              
    
    else:                   
        slicetime_params = pe.Node(interface=util.Function(input_names=['file', 'options'],
                                                           output_names=['nslices', 'tr', 'ta', 'sliceorder', 'refslice'],
                                                           function=slicetime_parameters), 
                                   name='slicetime_parameters')
        slicetime_params.inputs.options = options
        boldwork.connect(frameskip, 'roi_file', slicetime_params, 'file')    
        
        # --- If-else: choose whether to use SPM2 or SPM8 to perform slice-time correction    
        if options['spmversion'] == 2:                  
            slicetime = pe.Node(interface=util.Function(input_names=['file', 'nslices', 'tr', 'ta', 'refslice', 'options'],
                                                       output_names=['timecorrected_files'],
                                                       function=spm2_slicetimecorrect), 
                                name='spm2_slicetime') 
            slicetime.inputs.options = options
            boldsink_folder = boldsink_folder + '_spm2slicetime'
            boldwork.connect([(frameskip, slicetime, [('roi_file', 'file')]),
                              (slicetime_params, slicetime, [('nslices', 'nslices')]),
                              (slicetime_params, slicetime, [('tr', 'tr')]),
                              (slicetime_params, slicetime, [('ta', 'ta')]),
                              (slicetime_params, slicetime, [('refslice', 'refslice')]),
                            ])
            boldwork.connect([(slicetime, boldsink, [('timecorrected_files', boldsink_folder)]),
                            ])
            
        else:  
            print "[INFO]: Using SPM8 to perform slice-timing correction."     
            # --- Node (slicetime): slice-time correction, a.k.a frame alignment 
            slicetime = pe.Node(interface=spm.SliceTiming(), name='slicetime')
            slicetime.inputs.paths = options['spm8_path']
            boldsink_folder = boldsink_folder + '_spm8slicetime'
            # --- Make the connections in the BOLD workflow for slice time correction and all helper nodes
            boldwork.connect([(slicetime_params, slicetime, [('nslices', 'num_slices')]),
                              (slicetime_params, slicetime, [('tr', 'time_repetition')]),
                              (slicetime_params, slicetime, [('ta', 'time_acquisition')]),
                              (slicetime_params, slicetime, [('sliceorder', 'slice_order')]),
                              (slicetime_params, slicetime, [('refslice', 'ref_slice')]),
                              (frameskip, slicetime, [('roi_file', 'in_files')])        
                            ])
            boldwork.connect([ (slicetime, boldsink, [('timecorrected_files', boldsink_folder)]) ])
     
#-----------------------------MOTION CORRECTION-------------------------------#
    
    # --- Node (mcflirt): performs motion correction relative to the refvol frame of the scan
    mcflirt = pe.Node(interface=fsl.MCFLIRT(), name='mcflirt')
    mcflirt.inputs.save_plots = True
    mcflirt.inputs.save_rms = True 
    mcflirt.inputs.save_mats = True
#    mcflirt.inputs.stats_imgs = True
    mcflirt.inputs.output_type = 'NIFTI'
    boldsink_folder = boldsink_folder + '_mc'
    boldwork.connect(slicetime, 'timecorrected_files', mcflirt, 'in_file')
    boldwork.connect([(mcflirt, boldsink, [('mat_file', boldsink_folder+'.MATS')]),
                      (mcflirt, boldsink, [('par_file', boldsink_folder+'.@2')]),
                      (mcflirt, boldsink, [('rms_files', boldsink_folder+'.@3')]),
                      (mcflirt, boldsink, [('mean_img', boldsink_folder+'.@4')]),
                      (mcflirt, boldsink, [('std_img', boldsink_folder+'.@5')]),
                      (mcflirt, boldsink, [('variance_img', boldsink_folder+'.@6')]),
                      (mcflirt, boldsink, [('out_file', boldsink_folder+'.@7')])
                    ])
    # --- If-else: choose whether to use the existing refvol or to save one 
    if scan_number == options['bold'][0]:
        # --- Node (refvol): get the 10th volume of the BOLD scan
        refvol = pe.Node(interface=fsl.ExtractROI(), name='refvol')
        refvol.inputs.t_min = options['ref_vol']
        refvol.inputs.t_size = 1
        refvol.inputs.output_type = 'NIFTI'
        refvol.inputs.roi_file = 'refvol.nii'
        boldwork.connect(slicetime, 'timecorrected_files', refvol, 'in_file')
        boldwork.connect(refvol, 'roi_file', boldsink, '@')      
    else:
        refvol = pe.Node(interface=util.IdentityInterface(fields=['roi_file']), name='refvol')
        boldwork.add_nodes([refvol])
        preproc.connect(globals()['bold'+str(options['bold'][0])], 'refvol.roi_file', boldwork, 'refvol.roi_file')
    
    boldwork.connect(refvol, 'roi_file', mcflirt, 'ref_file') 
       
#----------------------------NORMALIZATION---------------------------------#
   
    # --- If-else: choose whether to use SPM2 or SPM8 to perform BOLD normalization
    if options['spmversion'] == 2: 
        spm_normalise = pe.Node(interface=util.Function(input_names=['apply_to_files', 'source', 'options'], 
                                                   output_names=['normalized_files', 'normalized_source', 'normalization_parameters'], 
                                                   function=spm2_norm), name='spm2_normalise') 
        spm_normalise.inputs.options = options
        boldsink_folder = boldsink_folder + "_spm2norm"
    else: 
        print "[INFO]: Using SPM8 to perform BOLD normalization."     
        # --- Node (spm normalize): Determine normalization parameters
        spm_normalise = pe.Node(interface=spm.Normalize(), name='spm8_normalise')
        spm_normalise.inputs.template = options['bold_template']
        spm_normalise.inputs.source_image_smoothing = 8
        spm_normalise.inputs.template_image_smoothing = 0
        spm_normalise.inputs.affine_regularization_type = 'mni'
        spm_normalise.inputs.DCT_period_cutoff = 25 
        spm_normalise.inputs.nonlinear_iterations = 16
        spm_normalise.inputs.nonlinear_regularization = 1
        spm_normalise.inputs.write_bounding_box = [[-90, -126, -72],[90, 90, 108]]
        spm_normalise.inputs.write_interp= 7
        spm_normalise.inputs.write_voxel_sizes = [2.0, 2.0, 2.0]
        spm_normalise.inputs.write_wrap = [0, 1, 0]
        spm_normalise.inputs.write_preserve = False
        spm_normalise.inputs.paths = [options['spm8_path']] 
        boldsink_folder = boldsink_folder + "_spm8norm"
        
    boldwork.connect([(refvol, spm_normalise, [('roi_file', 'source')]),
                      (mcflirt, spm_normalise, [('out_file', 'apply_to_files')])
                    ]) 
        
    # --- Node (merge): utility node to merge inputs into a single list
    merge = pe.Node(interface=util.Merge(2), name='merge')
    
    # --- Node (fslmaths): replace NaNs (improper numbers) with 0s
    nan2zero = pe.MapNode(interface=fsl.maths.MathsCommand(), name='nan2zero', iterfield='in_file')
    nan2zero.inputs.args = '-nan'
    nan2zero.inputs.output_type = 'NIFTI'  
    
    # --- Node (unmerge): take the list of BOLDs and refvol and separate them into fields
    unmerge = pe.Node(interface=util.Function(input_names=['list'],
                                              output_names=['refvol', 'bold'],
                                              function=unMerge), 
                     name='unmerge')
    

    boldwork.connect([(spm_normalise, merge, [('normalized_files', 'in1'),
                                              ('normalized_source', 'in2')]),
                      (merge, nan2zero, [('out', 'in_file')]),
                      (nan2zero, unmerge, [('out_file', 'list')]),
                    ])
    boldwork.connect([(unmerge, boldsink, [('bold', boldsink_folder)]),
                      (spm_normalise, boldsink, [('normalization_parameters', boldsink_folder+'.@1')]),
                    ])

#-----------------------------SMOOTHING-------------------------------------#
       
    if options['smooth'] != None:
        # --- Node (smooth): smooth the preprocessed data
        smooth = pe.Node(interface=fsl.Smooth(), name='smooth') 
        smooth.inputs.fwhm = options['smooth']
        smooth.inputs.output_type = 'NIFTI'  
        boldsink_folder = boldsink_folder + "_smooth" + str(options['smooth'])
        boldwork.connect(smooth, 'smoothed_file', boldsink, boldsink_folder)
        boldwork.connect(unmerge, 'bold', smooth, 'in_file')
    
#-------------------------create a file for fcMRI_preproc------------------#
    movement = pe.Node(interface=util.Function(input_names=['par_file', 'options'],
                                                           output_names=['dat', 'ddat', 'rdat'],
                                                           function=par_to_dat), 
                                   name='movement')
    movement.inputs.options = options
    boldwork.connect([(mcflirt, movement, [('par_file', 'par_file')]),
                      (movement, boldsink, [('dat', 'movement.@dat'),
                                            ('ddat', 'movement.@ddat'),
                                            ('rdat', 'movement.@rdat')
                                            ])
                    ])
  
    return boldwork

# --- BOLD workflow helper functions
def slicetime_placeholder(input): 
    return input

def slicetime_parameters(file, options):
    #this needs to be all one function due to NiPype engine constraints
            
    #set nslices
    if isinstance(options['num_slices'], int):
        nslices = options['num_slices']
    else:
        import nibabel
        nslices = nibabel.load(file).get_shape()[2]
        # TODO: In fsl_proprocess.sh sliceTimingDim=3 is hardcoded by default and is never changed although  
        # there is code for the cases if sliceTimingDim == 2 or sliceTimingDim == 1 (does fslswapdim)
    
    #set tr
    tr = options['tr']
                        
    #set ta
    ta = tr - tr / float(nslices)
    
    #set sliceorder         
    if options['sliceorder'].find('ascend') != -1:
        sliceorder = range(1, nslices + 1, 1)
    elif options['sliceorder'].find('descend') != -1:
        sliceorder = range(nslices, 0, -1)
    elif options['sliceorder'].find('inter') != -1 and options['sliceorder'].find('middle') != -1 and options['sliceorder'].find('top') != -1:
        sliceorder = [(nslices + 1 - k) / 2 + ((nslices - k) % 2 * (nslices - 1) / 2) + 1 for k in range(1, nslices + 1)]
    elif options['sliceorder'].find('inter') != -1 and options['sliceorder'].find('bottom') != -1 and options['sliceorder'].find('up') != -1:
        sliceorder = range(1, nslices + 1, 2) + range(2, nslices + 1, 2)
    elif options['sliceorder'].find('inter') != -1 and options['sliceorder'].find('top') != -1 and options['sliceorder'].find('down') != -1: 
        sliceorder = range(nslices, 0, -2) + range(nslices - 1, 0, -2) 
    else:
        try:                        
            from ast import literal_eval
            sliceorder = literal_eval(options['sliceorder'])
            options['sliceorder'] = "custom"
        except:
            print "[WARNING]: %s is not a valid entry for slice order. Defaulting to interleaved (bottom -> up)." %options['sliceorder']
            options['sliceorder'] = 'interleaved (bottom -> up)'
            sliceorder = range(1, nslices + 1, 2) + range(2, nslices + 1, 2) #default to interleaved (bottom -> up)

    #set refslice
    if isinstance(options['ref_slice'], int) and options['ref_slice'] > nslices:
        refslice = options['ref_slice']
    else:       
        import math
        refslice = sliceorder[int(math.floor(nslices / 2)) - 1] #middle slice
        print "[INFO]: Using default reference slice: %d" % refslice 
        #NOTE: not actually middle slice but this is what preproc does
          
    return nslices, tr, ta, sliceorder, refslice

def spm2_slicetimecorrect(file, nslices, tr, ta, refslice, options):
    import os
    import subprocess 
    import shutil    
    print "[INFO]: Using SPM2 to perform slice-timing correction."
    shutil.copy(file, os.getcwd())
    os.rename(os.path.join(os.getcwd(), os.path.basename(file)), 'bold.nii')
    os.symlink(os.path.join(options['code_dir'], 'spm2_slicetime.m'), 'spm2_slicetime.m')
    os.listdir(os.getcwd())
    matlabcmd = "matlab -nojvm -nodisplay -nodesktop -nosplash -r \"spm2_slicetime('%s', %d, %d, %d, '%s', '%s')\"" \
                % (options['sliceorder'], refslice, tr, nslices, options['spm_params'], options['spm2_path'])
    subprocess.call(matlabcmd, shell=True)   
    return os.path.abspath('bold_faln.nii.gz')   

def spm2_norm(apply_to_files, source, options):
    import os
    import subprocess 
    import shutil    
    print "[INFO]: Using SPM2 to normalize BOLD scan."
    shutil.copy(apply_to_files, os.getcwd())
    shutil.copy(source, os.getcwd())
    os.rename(os.path.join(os.getcwd(), os.path.basename(apply_to_files)), 'bold.nii')
    os.rename(os.path.join(os.getcwd(), os.path.basename(source)), 'refvol.nii')
    os.symlink(os.path.join(options['code_dir'], 'spm2_norm.m'), 'spm2_norm.m')  
    matlabcmd = "matlab -nojvm -nodisplay -nodesktop -nosplash -r \"spm2_norm('%s', '%s', '%s')\"" \
                         %(options['bold_template'],options['spm_params'],options['spm2_path'])
    subprocess.call(matlabcmd, shell=True)   
    return os.path.abspath('bold_atl.nii.gz'), os.path.abspath('wrefvol.nii.gz'), os.path.abspath('refvol-001_sn.mat')

def unMerge(list):
    import re
    for file in list: 
        if re.search(r'.*wrefvol.*\.nii(.gz)*$', file):
            refvol=file
        else: 
            bold=file
    return refvol, bold

def par_to_dat(par_file, options):
    from numpy import array
    import os
    
    name = os.path.basename(par_file)
    
    par = []
    with open(par_file) as f:
        for line in f:
            row = [float(x) for x in line.split()]
            row = [row[3], row[4], row[5], row[0], row[1], row[2], 1.]
            par.append(row)
    
    out = open("%s.dat" %name,'wb')
    for i in range(len(par)):
        y = ["%9.6f" %a for a in par[i]]
        out.write("%d %s" %(i+1," ".join(y)))
        out.write("\n")
    out.close()
           
    out = open("%s.ddat" %name,'wb')
    prev = par[0]
    for i in range(len(par)):
        y = []
        for j in range(len(par[i])):
            y.append("%9.6f" %(par[i][j] - prev[j]))
        if 0 == i:
            y[-1] = "%9.6f" %(1.0)
        out.write("%d %s" %(i+1," ".join(y)))
        out.write("\n")
        prev = par[i]
    out.close()
    
    if options['legacy']:
        out = open("%s.rdat" %name,'wb')
        data = par
        data_average = [0.0 for x in par[0]]
        for n in range(len(data)):
            for j in range(len(data[0])):
                data_average[j] = (data[n][j] + data_average[j]*(n+1)) / (n+2)
        for n in range(len(data)):
            y = []
            for j in range(len(data[0])):
                y.append("%9.6f" %(data[n][j] - data_average[j]))
            y[-1] = "%9.6f" %(1.0)
            out.write("%d %s" %(n+1," ".join(y)))
            out.write("\n")
        out.close()    
    else:
        out = open("%s.rdat" %name,'wb')
        par = array(par)
        par_min_mean = par - par.mean()
        for i in range(len(par)):
            y = ["%9.6f" %a for a in par[i]]
            out.write("%d %s" %(i+1," ".join(y)))
            out.write("\n")
        out.close()
    
    return os.path.abspath('%s.dat' %name), os.path.abspath('%s.ddat' %name), os.path.abspath('%s.rdat' %name)
    
if __name__ == '__main__':           
    starttime = time()
    info("Start time: %s" %str(ctime(starttime)))
   
    # --- Create options dictionary from config files and parsed command line inputs
    options = get_options()
    
    # --- Get DICOM files
    #if local dicom_dir unspecified, get data from XNAT server and create a dicom_dir
    if options['dicom_dir'] == None: 
        downloadDicoms(options)
    
    #get the first DICOM file from each scan
    options['dicoms_dict'] = getDicoms(options['dicom_dir'])
        
    #get the TR from the DICOM file of one of the BOLD scans
    if not isinstance(options['tr'], int):
        for value in options['dicoms_dict'].itervalues():
            dcm = dicom.read_file(value['example_file'])
            if dcm.data_element('SeriesNumber').value in options['bold']:
                options['tr'] = int(dcm.data_element('RepetitionTime').value)/1000
                break
         
    # --- Create the main workflow which will encompass the other workflows   
    preproc = pe.Workflow(name='preproc')
    preproc.base_dir = os.path.join(options['working_dir'], options['id'])   
    
    # --- Make and run the workflow to convert DICOMs to NIFTIs
    convertflow = make_converter_workflow(options)
    preproc.add_nodes([convertflow])
    convertflow.run()
                      
    # --- Make and run the bold workflow to get the reference volume from the first BOLD scan
    vars()['bold'+str(options['bold'][0])] = make_bold_workflow(options['bold'][0], 
                                                                'bold'+str(options['bold'][0])+'_refvol', 
                                                                options)   
    preproc.add_nodes([vars()['bold'+str(options['bold'][0])]])
    vars()['bold'+str(options['bold'][0])].run()
      
    # --- Make the bold workflows for the rest of the BOLD scans using the acquired reference volume 
    for bold in options['bold'][1:]:
        boldflow = make_bold_workflow(bold, 'bold'+str(bold), options)
    
    # --- Make the workflow to preprocess the anatomical scan
    anatflow = make_anat_workflow(options['anat'][0], options)
    preproc.add_nodes([anatflow])
    
    # --- Create workflow graph
    try:
        info("Creating workflow graph...")
        preproc.write_graph(simple_form=False, graph2use=options['writegraph'])
    except:
        info("Could not write graph. Proceeding to execute workflow.")
    
    # --- Run all nodes in the workflow (that have not yet been run)
    try:
        info("Attempting to run the workflow in parallel...")
        preproc.run(plugin='MultiProc', plugin_args={'n_procs': mp.cpu_count()})
    except:
        info("Could not run the workflow in parallel. Running serially.")
        preproc.run()
         
    # --- Optionally clean up
    if options['clean'] == True:
        info("Finished running the workflow. Cleaning up...")
        shutil.rmtree(os.path.join(options['working_dir'], options['id'], 'preproc'))
    else: 
        info("Not doing cleanup...")      
    
    info("---DONE---")
    
    runtime = time() - starttime  
    info("Finish time: %s" %str(ctime(time())))
    info("Runtime: %s" %str(timedelta(seconds=runtime)))
    

    
