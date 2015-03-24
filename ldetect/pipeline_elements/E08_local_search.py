#!/usr/bin/env python3.4

'''
Created on Aug 6, 2014

@author: tberisa
'''

import ldetect.baselib.flat_file_consts as cnst
import ldetect.baselib.flat_file as flat
import ldetect.baselib.binary_search as binsrch

import sys
import math

import decimal

class LocalSearch:
    def __init__(self, name, start_search, stop_search, initial_breakpoint_index, breakpoints, total_sum, total_N, input_config):
        decimal.getcontext().prec=50
        
        self.name = name
        self.start_search = start_search
        self.stop_search = stop_search
        self.initial_breakpoint_index = initial_breakpoint_index
        self.breakpoints = breakpoints
        self.total_sum = total_sum
        self.total_N = total_N
        self.input_config = input_config
        
        self.matrix = {}
        self.locus_list = []
        self.locus_list_deleted = []
        
        self.precomputed = {}
        self.precomputed['locus_list'] = [] # keep the ordering of loci -> allow for efficient iterating
        self.precomputed['data'] = {} # allow ~O(1) access to each element by it's locus
        
        self.dynamic_delete = True
        self.init_complete = False
        self.search_complete = False
        
        if start_search >= stop_search:
            raise Exception('Error: start_search >= stop_search')
        
        if initial_breakpoint_index>=len(breakpoints) or initial_breakpoint_index<0:
            raise Exception('Error: initial_breakpoint_index index out of bounds')
        
        if breakpoints[initial_breakpoint_index] >= stop_search:
            raise Exception('Error: breakpoint >= stop_search')
        
        if breakpoints[initial_breakpoint_index] <= start_search:
            raise Exception('Error: breakpoint <= start_search')
        
        # tmp_partitions = flat.read_partitions(self.name, self.input_config)
        tmp_partitions = flat.get_final_partitions(self.input_config, self.name, start_search, stop_search)
        
        
        if start_search < tmp_partitions[0][0] or start_search > tmp_partitions[len(tmp_partitions)-1][1]:
            raise Exception('Error: start_search is out of bounds')
        
        if stop_search < tmp_partitions[0][0] or stop_search > tmp_partitions[len(tmp_partitions)-1][1]:
            raise Exception('Error: stop_search is out of bounds')
        
        if initial_breakpoint_index > 0:
            if start_search < breakpoints[initial_breakpoint_index-1]:
                raise Exception('Error: start_search cannot be further than a neighboring breakpoint')
        else:
            pass # this is just to emphasize that this has been thought through and covered. It's taken care of when testing for start_search < tmp_partitions[0][0]
        
        if initial_breakpoint_index < (len(breakpoints)-1):
            if stop_search > breakpoints[initial_breakpoint_index+1]:
                raise Exception('Error: stop_search cannot be further than a neighboring breakpoint')
        else:
            pass # this is just to emphasize that this has been thought through and covered. It's taken care of when testing for stop_search > tmp_partitions[len(tmp_partitions)-1][1]
        
        # work out snp_first, snp_last - watch out if it's the first or last breakpoint
        
        # # snp_first defines where to start reading data
        # if initial_breakpoint_index > 0:
        #     self.snp_first = breakpoints[initial_breakpoint_index-1]
        # else:
        #     self.snp_first = tmp_partitions[0][0] # this gets the first SNP in the chromosome (setting it just to 1 would cause flat.relevant_subpartitions() and consequently flat.get_final_partitions() to fail)
  
        # The previous (above) was not taking into account start_search, but just assumed where search started!
        self.snp_first = start_search

        flat.print_log_msg('snp_first: '+repr(self.snp_first))

        # snp_last defined where to stop reading data
        self.snp_last = stop_search

        flat.print_log_msg('snp_last: '+repr(self.snp_last))
        
        # This is the upper bound for the search space (upper border)
        if initial_breakpoint_index+1 < len(breakpoints):
            self.snp_top = breakpoints[initial_breakpoint_index+1]
        else:
            self.snp_top = tmp_partitions[len(tmp_partitions)-1][1]
        
        flat.print_log_msg('snp_top: '+repr(self.snp_top))

        # This is the bottom bound for the search space (bottom border)
        if initial_breakpoint_index-1 >= 0:
            self.snp_bottom = breakpoints[initial_breakpoint_index-1]
        else:
            self.snp_bottom = tmp_partitions[0][0]

        flat.print_log_msg('snp_bottom: '+repr(self.snp_bottom))

        # flat.print_log_msg('In local search: ')
        # flat.print_log_msg(repr(self.snp_first)+' '+repr(self.snp_last)+' '+repr(self.snp_top))

        # Data must be read until snp_top!
        self.partitions = flat.get_final_partitions(self.input_config, self.name, self.snp_bottom, self.snp_top)

        # flat.print_log_msg('self.partitions: ')
        # flat.print_log_msg(repr(self.partitions))
        
        self.start_locus = -1
        self.start_locus_index = -1
        self.end_locus = -1
        self.end_locus_index = -1
    
    def init_search(self, lean=True):
        if lean:
            return self.init_search_lean()
        else:
            return self.init_search_full() 

    # Based on calc_vert:
    def init_search_full(self):
        # flat.print_log_msg('Removing existing matrix output file')
        # try:
        #     os.remove(cnst.const['out_matrix_delim'])
        # except OSError:
        #     pass
    
        if not self.dynamic_delete:
            raise Exception('Error: dynamic_delete should be True for local search!') 
    
        flat.print_log_msg('Start local search init') 

        # pre-read all relevant partitions at beginning!
        last_p_num = -1
        for p_num_init in range(0, len(self.partitions)-1):
            if self.snp_bottom >= self.partitions[p_num_init+1][0]:
                flat.print_log_msg('Pre-reading partition: '+str(self.partitions[p_num_init])) 
                flat.read_partition_into_matrix(self.partitions, p_num_init, self.matrix, self.locus_list, self.name, self.input_config, self.snp_bottom, self.snp_top)
                last_p_num = p_num_init
            else:
                break

        curr_locus = -1
        # for p_num, p in enumerate(self.partitions):
        for p_num in range(last_p_num+1, len(self.partitions)):
            p = self.partitions[p_num]

            flat.print_log_msg('Reading partition: '+str(p))
            # Data must be read until snp_top!
            flat.read_partition_into_matrix(self.partitions, p_num, self.matrix, self.locus_list, self.name, self.input_config, self.snp_bottom, self.snp_top)
    
            # Determine first locus
            if curr_locus<0: # Either first partition or not found in first partition
                # curr_locus = -1 # <- this should have been set to -1 before entering the main for loop
                if len(self.locus_list)>0:
                    # Find first locus >= snp_bottom
                    for i, locus in enumerate(self.locus_list):
                        if locus >= self.snp_bottom:
                            curr_locus = locus
                            start_locus = locus
                            curr_locus_index = i
                            start_locus_index = i
                            break
                else:
                    raise Exception('Error: locus_list seems to be empty') 
            # else:
            #   if len(self.locus_list)>0:
            #       curr_locus = self.locus_list[0]
            #       curr_locus_index = 0
            #   else:
            #       raise Exception('Error: locus_list seems to be empty')
            else:
                try:
                    curr_locus_index = self.locus_list.index(curr_locus)
                    # curr_locus is carried from prev iteration, but index has changed since part of matrix (and locus_list) has been deleted
                except ValueError:
                    if len(self.locus_list)>0:
                        curr_locus = self.locus_list[0]
                        curr_locus_index = 0
                    else:
                        raise Exception('Error: locus_list seems to be empty')

            if curr_locus<0:
                flat.print_log_msg('Warning: curr_locus not found! Continuing to next partition.')
                flat.print_log_msg('Comment: This is possibly due to snp_bottom being very close to end of partition.')
                flat.print_log_msg('Details: ')
                flat.print_log_msg('Partition: '+repr(p))
                flat.print_log_msg('snp_bottom: '+repr(self.snp_bottom))
                flat.print_log_msg('curr_locus: '+repr(curr_locus)) 
                continue #continue to next partition
                # raise Exception('Error: curr_locus not found!')   
            
            if p_num+1 < len(self.partitions):
                end_locus = self.partitions[p_num+1][0]
                end_locus_index = -1
            else:
                # end_locus = self.partitions[p_num][1]

                # Find last locus <= snp_last
                end_locus_found = False
                for i in reversed(range(0, len(self.locus_list))):
                # for locus in reversed(locus_list):
                    if self.locus_list[i] <= self.snp_last:
                        end_locus = self.locus_list[i]
                        end_locus_index = i
                        end_locus_found = True
                        break

                if not end_locus_found:
                    end_locus_index = 0
                    end_locus = self.locus_list[end_locus_index]
            
            # flat.print_log_msg('self.locus_list control output: '+repr(self.locus_list))

            flat.print_log_msg('Running precompute for partition: '+str(p))

            flat.print_log_msg('start_locus: '+repr(start_locus)+' end_locus: '+repr(end_locus)+' end_locus_index '+repr(end_locus_index))
            # This will not include the very last SNP of the complete range, but that shouldn't be too important since the end of the range shouldn't be a defining location for LD
            while curr_locus <= end_locus:                     
                self.add_locus_to_precomputed(curr_locus) # We want snp_bottom to be added here always (for later use). Same thing for snp_top
                
                # flat.print_log_msg('curr_locus: '+repr(curr_locus)+' end_locus: '+repr(end_locus))

                if (curr_locus > self.snp_first or self.initial_breakpoint_index == 0) and (curr_locus <= self.snp_last): # Do not include snp_first in the calculation unless the very first block is being taken into account. Do not calculate anything above snp_last, just insert dummies
                    for key, el in self.matrix[curr_locus]['data'].items():
                        # don't take into account anything over snp_top
                        if key <= self.snp_top:                        
                            corr_coeff = self.matrix[curr_locus]['data'][key]['shrink'] / math.sqrt( self.matrix[curr_locus]['data'][curr_locus]['shrink'] * self.matrix[key]['data'][key]['shrink'] )
                            
    #                         if curr_locus != key: # Don't include diagonal! ...although not that important.
                            self.add_val_to_precomputed(decimal.Decimal(corr_coeff**2), curr_locus, key) # If the diagonal is included, it doesn't matter because later we add and subtract is exactly once when adding and subra
    #                         else:
    #                             self.add_val_to_precomputed(decimal.Decimal(0), curr_locus, key)
                else:
                    self.add_val_to_precomputed(decimal.Decimal(0), curr_locus, curr_locus) # Dummy value for snp_first! ...in order to be consistent for some other future use of these data structures
                    
                if curr_locus_index+1 < len(self.locus_list):
                    curr_locus_index+=1
                    curr_locus = self.locus_list[curr_locus_index]
                else:
                    flat.print_log_msg('curr_locus_index out of bounds') # The possibility of this happening is only at the end of the range [usually chromosome] (end of last partition)
                    break
    
            # flat.delete_loci_smaller_than_and_output_matrix_to_file(end_locus, self.matrix, locus_list, locus_list_deleted, cnst.const['out_matrix_filename'])
            if self.dynamic_delete:
                flat.print_log_msg('Deleting loci not required any more')
                flat.delete_loci_smaller_than(end_locus, self.matrix, self.locus_list, self.locus_list_deleted)
    
        self.start_locus = start_locus
        self.start_locus_index = start_locus_index
        self.end_locus = end_locus
        self.end_locus_index = end_locus_index
    
        self.init_complete = True

    # Based on calc_vert:
    def init_search_lean(self):
        # flat.print_log_msg('Removing existing matrix output file')
        # try:
        #     os.remove(cnst.const['out_matrix_delim'])
        # except OSError:
        #     pass
    
        if not self.dynamic_delete:
            raise Exception('Error: dynamic_delete should be True for local search!') 
    
        flat.print_log_msg('Start local search init') 

        # pre-read all relevant partitions at beginning!
        last_p_num = -1
        for p_num_init in range(0, len(self.partitions)-1):
            if self.snp_bottom >= self.partitions[p_num_init+1][0]:
                flat.print_log_msg('Pre-reading partition: '+str(self.partitions[p_num_init])) 
                flat.read_partition_into_matrix_lean(self.partitions, p_num_init, self.matrix, self.locus_list, self.name, self.input_config, self.snp_bottom, self.snp_top)
                last_p_num = p_num_init
            else:
                break

        curr_locus = -1
        # for p_num, p in enumerate(self.partitions):
        for p_num in range(last_p_num+1, len(self.partitions)):
            p = self.partitions[p_num]

            flat.print_log_msg('Reading partition: '+str(p))
            # Data must be read until snp_top!
            flat.read_partition_into_matrix_lean(self.partitions, p_num, self.matrix, self.locus_list, self.name, self.input_config, self.snp_bottom, self.snp_top)
    
            # Determine first locus
            if curr_locus<0: # Either first partition or not found in first partition
                # curr_locus = -1 # <- this should have been set to -1 before entering the main for loop
                if len(self.locus_list)>0:
                    # Find first locus >= snp_bottom
                    for i, locus in enumerate(self.locus_list):
                        if locus >= self.snp_bottom:
                            curr_locus = locus
                            start_locus = locus
                            curr_locus_index = i
                            start_locus_index = i
                            break
                else:
                    raise Exception('Error: locus_list seems to be empty') 
            # else:
            #   if len(self.locus_list)>0:
            #       curr_locus = self.locus_list[0]
            #       curr_locus_index = 0
            #   else:
            #       raise Exception('Error: locus_list seems to be empty')
            else:
                try:
                    curr_locus_index = self.locus_list.index(curr_locus)
                    # curr_locus is carried from prev iteration, but index has changed since part of matrix (and locus_list) has been deleted
                except ValueError:
                    if len(self.locus_list)>0:
                        curr_locus = self.locus_list[0]
                        curr_locus_index = 0
                    else:
                        raise Exception('Error: locus_list seems to be empty')

            if curr_locus<0:
                flat.print_log_msg('Warning: curr_locus not found! Continuing to next partition.')
                flat.print_log_msg('Comment: This is possibly due to snp_bottom being very close to end of partition.')
                flat.print_log_msg('Details: ')
                flat.print_log_msg('Partition: '+repr(p))
                flat.print_log_msg('snp_bottom: '+repr(self.snp_bottom))
                flat.print_log_msg('curr_locus: '+repr(curr_locus)) 
                continue #continue to next partition
                # raise Exception('Error: curr_locus not found!')   
            
            if p_num+1 < len(self.partitions):
                end_locus = self.partitions[p_num+1][0]
                end_locus_index = -1
            else:
                # end_locus = self.partitions[p_num][1]

                # Find last locus <= snp_last
                end_locus_found = False
                for i in reversed(range(0, len(self.locus_list))):
                # for locus in reversed(locus_list):
                    if self.locus_list[i] <= self.snp_last:
                        end_locus = self.locus_list[i]
                        end_locus_index = i
                        end_locus_found = True
                        break

                if not end_locus_found:
                    end_locus_index = 0
                    end_locus = self.locus_list[end_locus_index]
            
            # flat.print_log_msg('self.locus_list control output: '+repr(self.locus_list))

            flat.print_log_msg('Running precompute for partition: '+str(p))

            flat.print_log_msg('start_locus: '+repr(start_locus)+' end_locus: '+repr(end_locus)+' end_locus_index '+repr(end_locus_index))
            # This will not include the very last SNP of the complete range, but that shouldn't be too important since the end of the range shouldn't be a defining location for LD
            while curr_locus <= end_locus:                     
                self.add_locus_to_precomputed(curr_locus) # We want snp_bottom to be added here always (for later use). Same thing for snp_top
                
                # flat.print_log_msg('curr_locus: '+repr(curr_locus)+' end_locus: '+repr(end_locus))

                if (curr_locus > self.snp_first or self.initial_breakpoint_index == 0) and (curr_locus <= self.snp_last): # Do not include snp_first in the calculation unless the very first block is being taken into account. Do not calculate anything above snp_last, just insert dummies
                    for key, el in self.matrix[curr_locus].items():
                        # don't take into account anything over snp_top
                        if key <= self.snp_top:                        
                            corr_coeff = self.matrix[curr_locus][key] / math.sqrt( self.matrix[curr_locus][curr_locus] * self.matrix[key][key] )
                            
    #                         if curr_locus != key: # Don't include diagonal! ...although not that important.
                            self.add_val_to_precomputed(decimal.Decimal(corr_coeff**2), curr_locus, key) # If the diagonal is included, it doesn't matter because later we add and subtract is exactly once when adding and subra
    #                         else:
    #                             self.add_val_to_precomputed(decimal.Decimal(0), curr_locus, key)
                else:
                    self.add_val_to_precomputed(decimal.Decimal(0), curr_locus, curr_locus) # Dummy value for snp_first! ...in order to be consistent for some other future use of these data structures
                    
                if curr_locus_index+1 < len(self.locus_list):
                    curr_locus_index+=1
                    curr_locus = self.locus_list[curr_locus_index]
                else:
                    flat.print_log_msg('curr_locus_index out of bounds') # The possibility of this happening is only at the end of the range [usually chromosome] (end of last partition)
                    break
    
            # flat.delete_loci_smaller_than_and_output_matrix_to_file(end_locus, self.matrix, locus_list, locus_list_deleted, cnst.const['out_matrix_filename'])
            if self.dynamic_delete:
                flat.print_log_msg('Deleting loci not required any more')
                flat.delete_loci_smaller_than_leanest(end_locus, self.matrix, self.locus_list)
    
        self.start_locus = start_locus
        self.start_locus_index = start_locus_index
        self.end_locus = end_locus
        self.end_locus_index = end_locus_index
    
        self.init_complete = True

    def add_val_to_precomputed(self, val, curr_locus, key):
        if curr_locus not in self.precomputed['data']:
            self.precomputed['data'][curr_locus] = {}
            self.precomputed['data'][curr_locus]['sum_vert'] = decimal.Decimal(0)
            self.precomputed['data'][curr_locus]['sum_horiz'] = decimal.Decimal(0)
#             self.precomputed['data'][curr_locus]['N_vert'] = decimal.Decimal(0)
#             self.precomputed['data'][curr_locus]['N_horiz'] = decimal.Decimal(0)
            
        if key not in self.precomputed['data']:
            self.precomputed['data'][key] = {}
            self.precomputed['data'][key]['sum_vert'] = decimal.Decimal(0)
            self.precomputed['data'][key]['sum_horiz'] = decimal.Decimal(0)
#             self.precomputed['data'][key]['N_vert'] = decimal.Decimal(0)
#             self.precomputed['data'][key]['N_horiz'] = decimal.Decimal(0)
            

        self.precomputed['data'][curr_locus]['sum_vert'] += val
#         self.precomputed['data'][curr_locus]['N_vert'] += 1
        self.precomputed['data'][key]['sum_horiz'] += val
#         self.precomputed['data'][key]['N_horiz'] += 1
        
    def add_locus_to_precomputed(self, curr_locus):
        self.precomputed['locus_list'].append(curr_locus)
        
    def search(self):
        if not self.init_complete:
            flat.print_log_msg('init_search() must be run before search(). Starting automatically...')
            self.init_search()
            
        flat.print_log_msg('Starting local search...')
        
        # In case the value itself is not in the list:
        try:
            snp_bottom_ind = binsrch.find_ge_ind(self.precomputed['locus_list'], self.snp_bottom)
            snp_top_ind = binsrch.find_le_ind(self.precomputed['locus_list'], self.snp_top)
        except Exception as e:
            flat.print_log_msg('Error2!')
            flat.print_log_msg(repr(e))
            flat.print_log_msg('self.precomputed[\'locus_list\']: '+repr(self.precomputed['locus_list']))
            flat.print_log_msg('self.snp_bottom: '+repr(self.snp_bottom))
            flat.print_log_msg('self.snp_first: '+repr(self.snp_first))
            flat.print_log_msg('self.snp_last: '+repr(self.snp_last))
            flat.print_log_msg('self.snp_top: '+repr(self.snp_top))
            flat.print_log_msg('self.__dict__: '+repr(self.__dict__))
            flat.print_log_msg('Continuing...')
            return self.breakpoints[self.initial_breakpoint_index], None


        # Old:
        # snp_first_ind = self.precomputed['locus_list'].index(self.snp_first) # This should be snp_bottom
        # snp_top_ind = self.precomputed['locus_list'].index(self.snp_top) 
        
        # Start from init breakpoint and search left. Then start from init_breakpoint again and search right.
        # We start from init_breakpoint because that's the initial sum and N that we have -> so we can use the precomputed data to incrementally check for 
        # Find the closest locus to the breakpoint value, because a breakpoint doesn't necessarily have to be in the locus_list
        breakpoint_index_in_locus_list = binsrch.find_le_ind(self.precomputed['locus_list'], self.breakpoints[self.initial_breakpoint_index])
        init_breakpoint_locus = self.precomputed['locus_list'][breakpoint_index_in_locus_list]
        # Old:
        # breakpoint_index_in_locus_list = self.precomputed['locus_list'].index(self.breakpoints[self.initial_breakpoint_index])

        curr_sum = self.total_sum
        curr_N = self.total_N
        
        min_metric = decimal.Decimal(self.total_sum) / decimal.Decimal(self.total_N)
        min_breakpoint = None
        
        min_metric_details = {}
        min_metric_details['sum'] = self.total_sum
        min_metric_details['N_zero'] = self.total_N
        min_distance_right = 0 # because the initial distance of the minimum actually is 0! (until we find a new minima to the RIGHT, or we don't in which case it doesn't matter)

        # Go RIGHT!
        flat.print_log_msg('Searching right...')
        if breakpoint_index_in_locus_list+1 < len(self.precomputed['locus_list']):
            curr_loc_ind = breakpoint_index_in_locus_list+1
            curr_loc = self.precomputed['locus_list'][curr_loc_ind]
            
            while curr_loc <= self.snp_last:
                curr_sum = curr_sum - self.precomputed['data'][curr_loc]['sum_horiz'] + self.precomputed['data'][curr_loc]['sum_vert']
                
                horiz_N = curr_loc_ind-snp_bottom_ind-1
                vert_N = snp_top_ind-curr_loc_ind
                curr_N = curr_N - horiz_N + vert_N
                
                curr_metric = decimal.Decimal(curr_sum) / decimal.Decimal(curr_N)
                
                if curr_metric < min_metric:
                    min_metric = curr_metric
                    min_breakpoint = curr_loc
                    min_metric_details['sum'] = curr_sum
                    min_metric_details['N_zero'] = curr_N
                    min_distance_right = curr_loc - init_breakpoint_locus

                
                if curr_loc_ind+1 < len(self.precomputed['locus_list']):
                    curr_loc_ind += 1
                    curr_loc = self.precomputed['locus_list'][curr_loc_ind]
                else:
                    flat.print_log_msg('curr_locus_index out of bounds') # The possibility of this happening is only at the end of the chromosome (end of last partition)
                    break
        else:
            flat.print_log_msg('Warning: breakpoint_index_in_locus_list+1 < len(self.precomputed["locus_list"]) not satisfied!')
            flat.print_log_msg('Breakpoints: '+repr(self.breakpoints))
            flat.print_log_msg('Locus_list: '+repr(self.precomputed['locus_list']))
            flat.print_log_msg('breakpoint_index_in_locus_list: '+ repr(breakpoint_index_in_locus_list))
        
        # Reset search for left
        curr_sum = self.total_sum
        curr_N = self.total_N

        # Go LEFT!    
        flat.print_log_msg('Searching left...')
        if breakpoint_index_in_locus_list-1 >= 0:
            curr_loc_ind = breakpoint_index_in_locus_list-1
            curr_loc = self.precomputed['locus_list'][curr_loc_ind]
            
            curr_sum = self.total_sum
            curr_N = self.total_N
            
            while curr_loc > self.snp_first: # Don't include previous breakpoint!
                curr_sum = curr_sum + self.precomputed['data'][curr_loc]['sum_horiz'] - self.precomputed['data'][curr_loc]['sum_vert']
                
                horiz_N = curr_loc_ind-snp_bottom_ind-1
                vert_N = snp_top_ind-curr_loc_ind
                curr_N = curr_N + horiz_N - vert_N
                
                curr_metric = decimal.Decimal(curr_sum) / decimal.Decimal(curr_N)
                
                if (curr_metric < min_metric) or (curr_metric == min_metric and (init_breakpoint_locus - curr_loc)<min_distance_right): # min_distance_right is used to compare to RIGHT metric, not within LEFT metric!
                    min_metric = curr_metric
                    min_breakpoint = curr_loc
                    min_metric_details['sum'] = curr_sum
                    min_metric_details['N_zero'] = curr_N

                if curr_loc_ind-1 >= 0:
                    curr_loc_ind -= 1
                    curr_loc = self.precomputed['locus_list'][curr_loc_ind]
                else:
                    flat.print_log_msg('curr_locus_index out of bounds') # The possibility of this happening is only at the beginning of the chromosome (start of first partition)
                    break
        else:
            flat.print_log_msg('Warning: breakpoint_index_in_locus_list-1 >=0 not satisfied!')
            flat.print_log_msg('Breakpoints: '+repr(self.breakpoints))
            flat.print_log_msg('Locus_list: '+repr(self.precomputed['locus_list']))
            flat.print_log_msg('breakpoint_index_in_locus_list: '+ repr(breakpoint_index_in_locus_list))
        
        self.search_complete = True
        
        flat.print_log_msg('Search done')
        
        return min_breakpoint, min_metric_details

def main():
    breakpoints = [20056346,
   23172864,
   26207725,
   27249779,
   29266559,
   29978822,
   31322564,
   33063813,
   33715859,
   35472318,
   36913379,
   39281968,
   40255964,
   42098394,
   43453056,
   44942909,
   46458807,
   48205362,
   50455659,
   51373940,
   53447645,
   54329768,
   55790298,
   58012944,
   58660548,
   59583650,
   61326949,
   63053905,
   65293922,
   67115837,
   68982238,
   70224115,
   71640609,
   74148469,
   75436123,
   78606362,
   81047897,
   84207279,
   87017863,
   88725515,
   90302399,
   92152779,
   93740145,
   94947530,
   96546735,
   97961068,
   99269331,
   100716423,
   102054465]

    breakpoint_index = 33 
    total_sum = decimal.Decimal('41049.603797938148512195858044319004257218538046177')
    total_N = decimal.Decimal('116785159748')

    tmp_begin = int((breakpoints[breakpoint_index-1]+breakpoints[breakpoint_index])/2)
    tmp_end = int((breakpoints[breakpoint_index]+breakpoints[breakpoint_index+1])/2)

    print('tmp_begin', tmp_begin, 'tmp_end', tmp_end)

    local_search_run = LocalSearch('chr15', tmp_begin, tmp_end, breakpoint_index, breakpoints, total_sum, total_N, cnst.const['orig_data_EUR'])       
    
    new_breakpoint, new_metric = local_search_run.search()
    
    print(new_breakpoint, new_metric['sum']/new_metric['N_zero'])
    print(breakpoints[breakpoint_index], total_sum/total_N)

if __name__ == '__main__':
    main()

