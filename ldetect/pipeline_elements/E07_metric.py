'''
Created on Aug 4, 2014

@author: tberisa
'''

import ldetect.baselib.flat_file_consts as cnst 
import ldetect.baselib.flat_file as flat 

import math
import sys
import decimal

class Metric:
    # def __init__(self, name, snp_first, snp_last, input_config, breakpoints):
    def __init__(self, name, input_config, breakpoints, snp_first=-1, snp_last=-1):
        self.breakpoints = breakpoints 

        self.matrix = {}
        self.locus_list = []
        self.locus_list_deleted = []
        
        decimal.getcontext().prec=50
        self.metric = {}
        self.metric['sum'] = decimal.Decimal('0')
        self.metric['N_nonzero'] = decimal.Decimal('0')
        self.metric['N_zero'] = decimal.Decimal('0')
        
        self.name = name
        self.snp_first, self.snp_last = flat.first_last(name, input_config, snp_first, snp_last)
        # self.snp_first = snp_first
        # self.snp_last = snp_last
        self.input_config = input_config
        
        self.partitions = flat.get_final_partitions(self.input_config, self.name, self.snp_first, self.snp_last)
        
        self.dynamic_delete = True
        self.calculation_complete = False
        
        self.start_locus = -1
        self.start_locus_index = -1
        self.end_locus = -1
        self.end_locus_index = -1
        

    def calc_metric(self, lean=True):
        if lean:
            return self.calc_metric_lean()
        else:
            return self.calc_metric_full()

    # Based on calc_vert: 
    def calc_metric_full(self):
        # flat.print_log_msg('Removing existing matrix output file')
        # try:
        #     os.remove(cnst.const['out_matrix_delim'])
        # except OSError:
        #     pass
        
        if not self.dynamic_delete:
            raise Exception('Error: dynamic delete must be True for metric calculation!')

        flat.print_log_msg('Start metric')
        
        curr_breakpoint_index = 0
        block_height = 0
        block_width = 0
        
        total_N_SNPs = decimal.Decimal('0')
        block_width_sum = decimal.Decimal('0')

        # pre-read all relevant partitions at beginning!
        last_p_num = -1
        for p_num_init in range(0, len(self.partitions)-1):
            if self.snp_first >= self.partitions[p_num_init+1][0]:
                flat.print_log_msg('Pre-reading partition: '+str(self.partitions[p_num_init])) 
                flat.read_partition_into_matrix(self.partitions, p_num_init, self.matrix, self.locus_list, self.name, self.input_config, self.snp_first, self.snp_last)
                last_p_num = p_num_init
            else:
                break

        curr_locus = -1
        # for p_num, p in enumerate(self.partitions):
        for p_num in range(last_p_num+1, len(self.partitions)):
            p = self.partitions[p_num]

            flat.print_log_msg('Reading partition: '+str(p))
            flat.read_partition_into_matrix(self.partitions, p_num, self.matrix, self.locus_list, self.name, self.input_config, self.snp_first, self.snp_last)

            # Determine first locus
            if curr_locus<0: # Either first partition or not found in first partition
                # curr_locus = -1 # <- this should have been set to -1 before entering the main for loop
                if len(self.locus_list)>0:
                    # Find first locus >= snp_first
                    for i, locus in enumerate(self.locus_list):
                        if locus >= self.snp_first:
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
                flat.print_log_msg('Comment: This is possibly due to snp_first being very close to end of partition.')
                flat.print_log_msg('Details: ')
                flat.print_log_msg('Partition: '+repr(p))
                flat.print_log_msg('snp_first: '+repr(self.snp_first))
                flat.print_log_msg('curr_locus: '+repr(curr_locus)) 
                continue #continue to next partition 
                # raise Exception('Error: curr_locus not found!')   
            
            # Determine last locus
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
            
            flat.print_log_msg('Running metric for partition: '+str(p))
            # This will not include the very last SNP of the complete range, but that shouldn't be too important since the end of the range shouldn't be a defining location for LD
            while curr_locus <= end_locus:
                if  curr_breakpoint_index<len(self.breakpoints): 
                    if curr_locus > self.breakpoints[curr_breakpoint_index]: # Breakpoint is the last element of the block!
#                         block_height =  len(self.locus_list) - curr_locus_index
                        block_height =  0 - total_N_SNPs # - 1 # ? # this is in accordance with the formula for deferred sum calculation 
                        self.metric['N_zero'] += block_height * block_width
                        block_width_sum += block_width
                        
                        curr_breakpoint_index += 1
                        block_width = 0
                
                if  curr_breakpoint_index>=len(self.breakpoints):
                    break
                
#                 found = False
                try:
                    for key, el in self.matrix[curr_locus]['data'].items():
                        if key > self.breakpoints[curr_breakpoint_index]: # Only add those above the breakpoint!
                            corr_coeff = self.matrix[curr_locus]['data'][key]['shrink'] / math.sqrt( self.matrix[curr_locus]['data'][curr_locus]['shrink'] * self.matrix[key]['data'][key]['shrink'] )
                            self.metric['sum'] += decimal.Decimal(corr_coeff**2)
                            self.metric['N_nonzero'] += 1
#                             found = True
                except IndexError as e:
                    print('Error!')
                    print(e)
                    print(key, el)
                    print(curr_locus)
                    print(self.matrix)
                    print(self.breakpoints)
                    print(curr_breakpoint_index)
                    
#                 if found:
                block_width += 1 # block_width needs to be increased even if it doesn't have values in the outer part of the matrix! 
                    
                if curr_locus_index+1 < len(self.locus_list):
                    curr_locus_index+=1
                    curr_locus = self.locus_list[curr_locus_index]
                    total_N_SNPs += 1
                else:
                    flat.print_log_msg('curr_locus_index out of bounds')
                    break

#             if block_width > 0: # If an LD block hasn't finished, but a new partition must be read into memory
# #                 index_of_breakpoint_in_locus_list = -1
#                 for ind in range(curr_locus_index, len(self.locus_list)):
#                     if self.locus_list[ind] >= self.breakpoints[curr_breakpoint_index]:
# #                         index_of_breakpoint_in_locus_list = ind
#                         break
#                 
#                 num_of_SNPs_to_add = ind - curr_locus_index
#                 
# #                 if index_of_breakpoint_in_locus_list < 0:
# #                     raise Exception('Error: index_of_breakpoint_in_locus_list not found!')
#                 
# #                 block_height =  len(self.locus_list) - index_of_breakpoint_in_locus_list
#                 block_height =  0 - (total_N_SNPs+num_of_SNPs_to_add)
#                 self.metric['N_zero'] += block_height * block_width
#                 
#                 block_width_sum += block_width
#                 block_width = 0
                
            # flat.delete_loci_smaller_than_and_output_matrix_to_file(end_locus, self.matrix, locus_list, locus_list_deleted, cnst.const['out_matrix_filename'])
            if self.dynamic_delete:
                flat.print_log_msg('Deleting loci not required any more')
                flat.delete_loci_smaller_than(end_locus, self.matrix, self.locus_list, self.locus_list_deleted)

        self.start_locus = start_locus
        self.start_locus_index = start_locus_index
        self.end_locus = end_locus
        self.end_locus_index = end_locus_index
        
        self.metric['N_zero'] += total_N_SNPs * block_width_sum # this is in accordance with the formula for deferred sum calculation
        
        print('total_N_SNPs, block_width', total_N_SNPs, block_width)
        print('total_N_SNPs-block_width', total_N_SNPs-block_width)
        print('block_width_sum', block_width_sum)
        
        self.calculation_complete = True
        
        return self.metric

    # Based on calc_vert: 
    def calc_metric_lean(self):
        # flat.print_log_msg('Removing existing matrix output file')
        # try:
        #     os.remove(cnst.const['out_matrix_delim'])
        # except OSError:
        #     pass
        
        if not self.dynamic_delete:
            raise Exception('Error: dynamic delete must be True for metric calculation!')

        flat.print_log_msg('Start metric')
        
        curr_breakpoint_index = 0
        block_height = 0
        block_width = 0
        
        total_N_SNPs = decimal.Decimal('0')
        block_width_sum = decimal.Decimal('0')

        # pre-read all relevant partitions at beginning!
        last_p_num = -1
        for p_num_init in range(0, len(self.partitions)-1):
            if self.snp_first >= self.partitions[p_num_init+1][0]:
                flat.print_log_msg('Pre-reading partition: '+str(self.partitions[p_num_init])) 
                flat.read_partition_into_matrix_lean(self.partitions, p_num_init, self.matrix, self.locus_list, self.name, self.input_config, self.snp_first, self.snp_last)
                last_p_num = p_num_init
            else:
                break

        curr_locus = -1
        # for p_num, p in enumerate(self.partitions):
        for p_num in range(last_p_num+1, len(self.partitions)):
            p = self.partitions[p_num]

            flat.print_log_msg('Reading partition: '+str(p))
            flat.read_partition_into_matrix_lean(self.partitions, p_num, self.matrix, self.locus_list, self.name, self.input_config, self.snp_first, self.snp_last)

            # Determine first locus
            if curr_locus<0: # Either first partition or not found in first partition
                # curr_locus = -1 # <- this should have been set to -1 before entering the main for loop
                if len(self.locus_list)>0:
                    # Find first locus >= snp_first
                    for i, locus in enumerate(self.locus_list):
                        if locus >= self.snp_first:
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
                flat.print_log_msg('Comment: This is possibly due to snp_first being very close to end of partition.')
                flat.print_log_msg('Details: ')
                flat.print_log_msg('Partition: '+repr(p))
                flat.print_log_msg('snp_first: '+repr(self.snp_first))
                flat.print_log_msg('curr_locus: '+repr(curr_locus)) 
                continue #continue to next partition 
                # raise Exception('Error: curr_locus not found!')   
            
            # Determine last locus
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
            
            flat.print_log_msg('Running metric for partition: '+str(p))
            # This will not include the very last SNP of the complete range, but that shouldn't be too important since the end of the range shouldn't be a defining location for LD
            while curr_locus <= end_locus:
                if  curr_breakpoint_index<len(self.breakpoints): 
                    if curr_locus > self.breakpoints[curr_breakpoint_index]: # Breakpoint is the last element of the block!
#                         block_height =  len(self.locus_list) - curr_locus_index
                        block_height =  0 - total_N_SNPs # - 1 # ? # this is in accordance with the formula for deferred sum calculation 
                        self.metric['N_zero'] += block_height * block_width
                        block_width_sum += block_width
                        
                        curr_breakpoint_index += 1
                        block_width = 0
                
                if  curr_breakpoint_index>=len(self.breakpoints):
                    break
                
#                 found = False
                try:
                    for key, el in self.matrix[curr_locus].items():
                        if key > self.breakpoints[curr_breakpoint_index]: # Only add those above the breakpoint!
                            corr_coeff = self.matrix[curr_locus][key] / math.sqrt( self.matrix[curr_locus][curr_locus] * self.matrix[key][key] )
                            self.metric['sum'] += decimal.Decimal(corr_coeff**2)
                            self.metric['N_nonzero'] += 1
#                             found = True
                except IndexError as e:
                    print('Error!')
                    print(e)
                    print(key, el)
                    print(curr_locus)
                    print(self.matrix)
                    print(self.breakpoints)
                    print(curr_breakpoint_index)
                    
#                 if found:
                block_width += 1 # block_width needs to be increased even if it doesn't have values in the outer part of the matrix! 
                    
                if curr_locus_index+1 < len(self.locus_list):
                    curr_locus_index+=1
                    curr_locus = self.locus_list[curr_locus_index]
                    total_N_SNPs += 1
                else:
                    flat.print_log_msg('curr_locus_index out of bounds')
                    break

#             if block_width > 0: # If an LD block hasn't finished, but a new partition must be read into memory
# #                 index_of_breakpoint_in_locus_list = -1
#                 for ind in range(curr_locus_index, len(self.locus_list)):
#                     if self.locus_list[ind] >= self.breakpoints[curr_breakpoint_index]:
# #                         index_of_breakpoint_in_locus_list = ind
#                         break
#                 
#                 num_of_SNPs_to_add = ind - curr_locus_index
#                 
# #                 if index_of_breakpoint_in_locus_list < 0:
# #                     raise Exception('Error: index_of_breakpoint_in_locus_list not found!')
#                 
# #                 block_height =  len(self.locus_list) - index_of_breakpoint_in_locus_list
#                 block_height =  0 - (total_N_SNPs+num_of_SNPs_to_add)
#                 self.metric['N_zero'] += block_height * block_width
#                 
#                 block_width_sum += block_width
#                 block_width = 0
                
            # flat.delete_loci_smaller_than_and_output_matrix_to_file(end_locus, self.matrix, locus_list, locus_list_deleted, cnst.const['out_matrix_filename'])
            if self.dynamic_delete:
                flat.print_log_msg('Deleting loci not required any more')
                # delete_loci_smaller_than does not need to be lean!!!
                flat.delete_loci_smaller_than_leanest(end_locus, self.matrix, self.locus_list)

        self.start_locus = start_locus
        self.start_locus_index = start_locus_index
        self.end_locus = end_locus
        self.end_locus_index = end_locus_index
        
        self.metric['N_zero'] += total_N_SNPs * block_width_sum # this is in accordance with the formula for deferred sum calculation
        
        print('total_N_SNPs, block_width', total_N_SNPs, block_width)
        print('total_N_SNPs-block_width', total_N_SNPs-block_width)
        print('block_width_sum', block_width_sum)
        
        self.calculation_complete = True
        
        return self.metric

def main():
#     def __init__(self, name, snp_first, snp_last, input_config, breakpoints):
#      
    begin = 9411243
    end = 48119216
    
#     begin = 46287140
#     end = 48119216
     
    breakpoints1 = [10148322, 15250019, 15864313, 16491839, 17748811, 18252127, 18912106, 19637870, 20332293, 20929869, 21190923, 21649595, 22318833, 23231365, 24271200, 24774771, 25035980, 26088085, 27431612, 27666047, 28290149, 28485200, 28761470, 29335757, 29790442, 30972911, 32778127, 33370496, 34413058, 35253882, 35614394, 36328018, 37283402, 38078491, 39227880, 39908770, 40259482, 40965403, 41448115, 41676786, 42689700, 43100808, 43345207, 43799567, 44748107, 45265729, 45789905, 46336509, 46883153, 47465743]   
    
    metric = Metric('chr21', cnst.const['orig_data'], breakpoints1, begin, end)
    out = metric.calc_metric()
    print(out)
    print(out['sum']/out['N_zero'])

    breakpoints2 = [i for i in range(begin, end+1, int((end-begin)/(len(breakpoints1)-1)))]
    
    metric = Metric('chr21', cnst.const['orig_data'], breakpoints2, begin, end)
    out = metric.calc_metric()
    print(out)
    print(out['sum']/out['N_zero'])
    
    flat.print_log_msg('Done')

if __name__ == '__main__':
    main()
