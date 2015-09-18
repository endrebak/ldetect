#!/usr/bin/env python3

'''
Created on Aug 7, 2014

@author: tberisa
'''

import sys
import getopt
import decimal
import pickle
import datetime

import ldetect.baselib.flat_file_consts as cnst
import ldetect.baselib.flat_file as flat
import ldetect.baselib.read_data as rd
import ldetect.baselib.filters as filt

import ldetect.pipeline_elements.E03_matrix_to_vector as matrix_to_vector

import commanderline.commander_line as cl

# Script start

# def main(argv=None):
#     '''
#     -d dataset (orig_data or synth_data)
#     -c chromosome name (chrN, where N in {1,2,...,22} for orig_data)
#     -b (optional) default: first available locus in chromosome ..."begin"
#     -e (optional) default: last available locus in chromosome ..."end"
#     -i (optional) no* / yes ..."image"
#     -o (optional) diag* / vert ..."orientation"
#     -r (optional) sum* / avg ..."reduce"

#     * - default
#     '''
    
#     options = {}
#     options['dataset'] = {}
#     options['dataset']['val'] = None
#     options['dataset']['flag'] = '-d'
#     options['chr_name'] = {}
#     options['chr_name']['val'] = None
#     options['chr_name']['flag'] = '-c'
#     options['begin'] = {}
#     options['begin']['val'] = -1
#     options['begin']['flag'] = '-b'
#     options['end'] = {}
#     options['end']['val'] = -1
#     options['end']['flag'] = '-e'
#     options['img'] = {}
#     options['img']['val'] = 'no'
#     options['img']['flag'] = '-i'
#     options['img']['valid'] = ('yes', 'no')
#     options['orient'] = {}
#     options['orient']['val'] = 'diag'
#     options['orient']['flag'] = '-o'
#     options['orient']['valid'] = ('diag', 'vert')
#     options['red'] = {}
#     options['red']['val'] = 'sum'
#     options['red']['flag'] = '-r'
#     options['red']['valid'] = ('sum', 'avg')

#     # Mandatory options
#     # mandatory_opts = {'-d', '-c', '-b', '-e'}

#     # Credit for arg parsing to Guido van Rossum
#     if argv is None:
#         argv = sys.argv

#     flat.print_log_msg('Run parameters: ' + repr(argv))

#     try:
# #         opts, args = getopt.getopt(argv[1:], 'b:c:d:e:i:h:n', ['ini=','help'])
#         opts, args = getopt.getopt(argv[1:], 'd:c:b:e:i:o:r:h', ['help'])
#     except getopt.error as msg:
#         print(msg)
#         return print_opt_arg_error()

#     # process options

#     for o, a in opts:
#         if o in ('-h', '--help'):
#             print(main.__doc__)
#             return(0)
#         elif o == options['dataset']['flag']:
#             options['dataset']['val'] = a
#         elif o == options['chr_name']['flag']:
#             options['chr_name']['val'] = a
#         elif o == options['begin']['flag']:
#             options['begin']['val'] = int(a)
#         elif o == options['end']['flag']:
#             options['end']['val'] = int(a)
#         elif o == options['img']['flag']:
#             options['img']['val'] = a
#         elif o == options['orient']['flag']:
#             options['orient']['val'] = a
#         elif o == options['red']['flag']:
#             options['red']['val'] = a
# #         elif o == '--ini':
# #             ini_fname = a
#         else:
#             print('Error: Unkown option '+o)
#             return print_opt_arg_error()

#     # process arguments
#     for arg in args:
#         #process(arg) # process() is defined elsewhere 
#         print("Error: Unkown argument: "+arg)
#         return print_opt_arg_error() 

#     # verification
#     for arg in options:
#         if options[arg]['val'] is None:
#             raise Exception('Missing flag: '+options[arg]['flag'])

#         if 'valid' in options[arg]:
#             if options[arg]['val'] not in options[arg]['valid']:
#                 raise Exception('Unkown value: '+options[arg]['val']+' for flag: '+options[arg]['flag'])

#     # for o in mandatory_opts:
#     #     if o not in user_opts:
#     #         print('Error: Missing mandatory option: '+o)
#     #         return print_opt_arg_error()

#     pipeline_lean(options['dataset']['val'], options['chr_name']['val'], options['begin']['val'], options['end']['val'], options['img']['val'], options['orient']['val'], options['red']['val'])

#     print('Done.')

def multiple_pipelines(dataset, name, begin_list, end_list, img='no', orient='diag', red='sum', delta=1000000):
    begin_list = cl.parse_list_from_cmd_line(begin_list)
    end_list = cl.parse_list_from_cmd_line(end_list)

    run = 0

    for begin, end in zip(begin_list, end_list):
        # begin = min(i1,i2)
        # end = max(i1,i2) 
        # pipeline(dataset, name, begin-int(delta), end+int(delta), img, orient, red, begin, end)
        pipeline(dataset, name, begin-int(delta), begin+int(delta), img, orient, red, begin, comment=str(run)+'-begin')
        pipeline(dataset, name, end-int(delta), end+int(delta), img, orient, red, end, comment=str(run)+'-end') 

        run+=1

def pipeline(dataset, name, begin=-1, end=-1, img='no', orient='diag', red='sum', snp=None, comment=''):
    '''
    pipeline(dataset, name, begin=-1, end=-1, img='no', orient='diag', red='sum', snp=None, comment='')

    snp1 and snp2 are loci of two SNPs that need to be converted into ordinal numbers representing row/col in image of matrix
    '''

    analysis = matrix_to_vector.MatrixAnalysis(name, cnst.const[dataset], begin, end)

    print(analysis.snp_first)
    print(analysis.snp_last)

    if(img=='yes'):
        generate_img = True
    elif(img=='no'):
        generate_img = False
    else:
        raise Exception('Error: Unknown argument: '+img)

    if(orient=='vert'):
        analysis.calc_vert(not generate_img) 
    elif(orient=='diag'):
        analysis.calc_diag(not generate_img) 
    else:
        raise Exception('Error: Unknown argument: '+orient)

    if(red=='avg'):
        avg = True
        raise Exception('Average used, but its output is not always consistent - especially for diag!')
    elif(red=='sum'):
        avg = False
    else:
        raise Exception('Error: Unknown argument: '+red)

    t = datetime.datetime.now() 
    t_formatted = t.strftime('%Y_%m_%d_%H_%M_%S')

    out_fname = 'vector-'+dataset+'-'+name+'-'+comment+'-'+str(analysis.snp_first)+'-'+str(analysis.snp_last)+'-'+orient+'-'+red+'-img_'+img+'-'+t_formatted

    analysis.write_output_to_file(out_fname+'.txt.gz', cnst.const['out_delim'], avg)

    if generate_img:
        # flat.print_log_msg('x_values: '+repr(x_values))
        if snp is not None:
            analysis.generate_img('img-'+out_fname+cnst.const['img_out_ext'], snp)
        else:
            analysis.generate_img('img-'+out_fname+cnst.const['img_out_ext'])



    flat.print_log_msg('Done') 

def pipeline_lean(dataset, name, begin=-1, end=-1, img='no', orient='diag', red='sum'):
    '''
    pipeline_lean(dataset, name, begin=-1, end=-1, img='no', orient='diag', red='sum')
    '''
    
    analysis = matrix_to_vector.MatrixAnalysis(name, cnst.const[dataset], begin, end)

    print(analysis.snp_first)
    print(analysis.snp_last)

    t = datetime.datetime.now() 
    t_formatted = t.strftime('%Y_%m_%d_%H_%M_%S')

    out_fname = 'vector-'+dataset+'-'+name+'-'+str(analysis.snp_first)+'-'+str(analysis.snp_last)+'-'+orient+'-'+red+'-img_'+img+'-'+t_formatted
    out_fname += '.txt.gz'
    flat.print_log_msg('out_fname: '+out_fname)

    if(img=='yes'):
        generate_img = True
    elif(img=='no'):
        generate_img = False
    else:
        raise Exception('Error: Unknown argument: '+img)

    if(orient=='vert'):
        analysis.calc_vert(not generate_img) 
    elif(orient=='diag'):
        analysis.calc_diag_lean(out_fname, cnst.const['out_delim'], not generate_img) 
    else:
        raise Exception('Error: Unknown argument: '+orient)

    if(red=='avg'):
        avg = True
        raise Exception('Average used, but its output is not always consistent - especially for diag!')
    elif(red=='sum'):
        avg = False
    else:
        raise Exception('Error: Unknown argument: '+red)

    # Output is done step-by-step
    # analysis.write_output_to_file(out_fname+'.txt.gz', cnst.const['out_delim'], avg)

    if generate_img:
        analysis.generate_img(out_fname+cnst.const['img_out_ext'])

    flat.print_log_msg('Done')

def print_opt_arg_error():
    print('For help use --help')
    return(2)

# if __name__ == '__main__':
#     main()

if __name__ == '__main__':
    cl.commander_line((pipeline_lean, pipeline, multiple_pipelines)) 
